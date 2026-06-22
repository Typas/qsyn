/****************************************************************************
  PackageName  [ tableau/optimize ]
  Synopsis     [ Partial-Hadamard gadgetization (gadgetize step) ]
  Copyright    [ Copyright(c) 2024 DVLab, GIEE, NTU, Taiwan ]
****************************************************************************/

#include "tableau/optimize/adaptive_gadget.hpp"

#include <spdlog/spdlog.h>

#include <array>
#include <cstddef>
#include <optional>
#include <queue>
#include <unordered_set>
#include <utility>
#include <vector>

#include "convert/qcir_to_tableau.hpp"
#include "convert/tableau_to_qcir.hpp"
#include "qcir/basic_gate_type.hpp"
#include "qcir/operation.hpp"
#include "qcir/qcir.hpp"
#include "qcir/qcir_gate.hpp"

namespace qsyn::tableau {

namespace {

// Grow `cone` along the circuit DAG from `start`, following predecessors
// (forward == false) or successors (forward == true). Uses only QCir's public
// neighbor accessors, so no qcir-module change is needed.
void add_cone(qcir::QCir const& qcir, qcir::QCirGate* start,
              std::unordered_set<qcir::QCirGate*>& cone, bool forward) {
    std::queue<size_t> queue;
    queue.push(start->get_id());
    cone.insert(start);
    while (!queue.empty()) {
        auto const id = queue.front();
        queue.pop();
        auto const neighbors = forward ? qcir.get_successors(id) : qcir.get_predecessors(id);
        for (auto const& nb : neighbors) {
            if (!nb.has_value()) continue;
            auto* g = qcir.get_gate(nb);
            if (cone.insert(g).second) queue.push(*nb);
        }
    }
}

// The set of *internal* Hadamards: an H with a non-Clifford (T-like) gate
// reachable both backward and forward through the circuit DAG. This is the same
// causal-cone rule qsyn uses for its internal-H statistic
// (qcir.cpp get_gate_statistics) and matches the structural count exactly. A
// naive global first/last-rotation span over get_gates() is WRONG here because
// get_gates() is topologically reordered, so boundary Hadamards on quiet wires
// land between the global first/last rotation and get miscounted as internal.
std::unordered_set<qcir::QCirGate*> internal_hadamards(qcir::QCir const& qcir) {
    std::unordered_set<qcir::QCirGate*> has_rotation_after;   // a non-Clifford lies ahead
    std::unordered_set<qcir::QCirGate*> has_rotation_before;  // a non-Clifford lies behind
    for (auto* gate : qcir.get_gates()) {
        if (is_clifford(gate->get_operation())) continue;
        add_cone(qcir, gate, has_rotation_after, /*forward=*/false);
        add_cone(qcir, gate, has_rotation_before, /*forward=*/true);
    }

    std::unordered_set<qcir::QCirGate*> internal;
    for (auto* gate : qcir.get_gates()) {
        if (gate->get_operation().get_type() == "h" &&
            has_rotation_after.contains(gate) && has_rotation_before.contains(gate)) {
            internal.insert(gate);
        }
    }
    return internal;
}

// Emit the fixed Vandaele Hadamard gadget body onto data wire q and ancilla a
// (mirrors quantum-circuit-optimization/src/circuit.rs:190-196):
//     S(a) . S(q) . CX(q,a) . S(a) . Z(a) . CX(a,q) . CX(q,a)
// The gadget contains no Hadamard; the data wire stays on q (no remap). The
// ancilla's Hadamard prep/uncompute are emitted separately, around the body.
void emit_hadamard_gadget(qcir::QCir& out, size_t q, size_t a) {
    out.append(qcir::SGate(), {a});
    out.append(qcir::SGate(), {q});
    out.append(qcir::CXGate(), {q, a});
    out.append(qcir::SGate(), {a});
    out.append(qcir::ZGate(), {a});
    out.append(qcir::CXGate(), {a, q});
    out.append(qcir::CXGate(), {q, a});
}

}  // namespace

qcir::QCir gadgetize_all_internal_hadamards(qcir::QCir const& qcir) {
    auto const n        = qcir.get_num_qubits();
    auto const internal = internal_hadamards(qcir);

    qcir::QCir out{n + internal.size()};

    // Ancilla Hadamard prep, before the body (c_out = anc + body + anc).
    for (size_t a = 0; a < internal.size(); ++a) {
        out.append(qcir::HGate(), {n + a});
    }

    // Body: copy gates in order, replacing each internal Hadamard with the gadget
    // onto a fresh ancilla. The data wire is never remapped.
    size_t next_ancilla = n;
    for (auto* gate : qcir.get_gates()) {
        if (internal.contains(gate)) {
            emit_hadamard_gadget(out, gate->get_qubits()[0], next_ancilla);
            ++next_ancilla;
        } else {
            out.append(gate->get_operation(), gate->get_qubits());
        }
    }

    // Ancilla Hadamard uncompute, after the body.
    for (size_t a = 0; a < internal.size(); ++a) {
        out.append(qcir::HGate(), {n + a});
    }

    return out;
}

namespace {

// Splice the fixed Vandaele gadget body in place of a Hadamard on data wire q,
// using fresh ancilla a (mirrors circuit.rs:190-196, minus the prep/uncompute H
// which are handled at the boundaries):
//     S(a) . S(q) . CX(q,a) . S(a) . Z(a) . CX(a,q) . CX(q,a)
void append_gadget_body(CliffordOperatorString& ops, size_t q, size_t a) {
    ops.emplace_back(CliffordOperatorType::s, std::array<size_t, 2>{a, 0});
    ops.emplace_back(CliffordOperatorType::s, std::array<size_t, 2>{q, 0});
    ops.emplace_back(CliffordOperatorType::cx, std::array<size_t, 2>{q, a});
    ops.emplace_back(CliffordOperatorType::s, std::array<size_t, 2>{a, 0});
    ops.emplace_back(CliffordOperatorType::z, std::array<size_t, 2>{a, 0});
    ops.emplace_back(CliffordOperatorType::cx, std::array<size_t, 2>{a, q});
    ops.emplace_back(CliffordOperatorType::cx, std::array<size_t, 2>{q, a});
}

size_t count_hadamards(StabilizerTableau const& st) {
    size_t count = 0;
    for (auto const& op : extract_clifford_operators(st)) {
        if (op.first == CliffordOperatorType::h) ++count;
    }
    return count;
}

}  // namespace

bool gadgetize_internal_hadamards(Tableau& tableau) {
    if (tableau.is_empty()) return true;
    auto const n = tableau.n_qubits();

    // Structural internal Hadamards live in Clifford blocks strictly between the
    // first and last phase-polynomial region. Tableau-native: no circuit
    // round-trip (which is lossy for FastTODD), rotations stay symbolic.
    std::optional<size_t> first_rot;
    std::optional<size_t> last_rot;
    for (size_t i = 0; i < tableau.size(); ++i) {
        if (std::holds_alternative<std::vector<PauliRotation>>(tableau[i])) {
            if (!first_rot) first_rot = i;
            last_rot = i;
        }
    }
    if (!first_rot) return true;  // no phase polynomial, nothing to gadgetize

    size_t k = 0;
    for (size_t i = *first_rot + 1; i < *last_rot; ++i) {
        if (auto const* st = std::get_if<StabilizerTableau>(&tableau[i])) {
            k += count_hadamards(*st);
        }
    }
    if (k == 0) return true;  // no internal Hadamards

    // Widen everything to n + k idle ancillae, then gadgetize each internal block:
    // rebuild it from its gate string with every Hadamard replaced by the gadget
    // body onto a fresh ancilla.
    Tableau ext        = tableau.extended_with_ancillae(k);
    size_t next_ancilla = n;
    for (size_t i = *first_rot + 1; i < *last_rot; ++i) {
        auto* st = std::get_if<StabilizerTableau>(&ext[i]);
        if (!st) continue;
        CliffordOperatorString spliced;
        for (auto const& op : extract_clifford_operators(*st)) {
            if (op.first == CliffordOperatorType::h) {
                append_gadget_body(spliced, op.second[0], next_ancilla++);
            } else {
                spliced.push_back(op);
            }
        }
        *st = StabilizerTableau{n + k}.apply(spliced);
    }

    // Bounded interior-merge: merge ONLY the phase-poly regions and the (now
    // Hadamard-free) gadget bodies between the first and last region into one
    // diagonal region, keeping the leading and trailing Clifford blocks as
    // boundaries. Running collapse on the whole tableau would instead thread the
    // trailing block's Hadamards leftward into the region (H Z = X H), making it
    // non-diagonal so FastTODD skips it.
    Tableau interior{n + k};  // seeded with an identity Clifford front
    for (size_t i = *first_rot; i <= *last_rot; ++i) interior.push_back(ext[i]);
    collapse(interior);  // -> [C_int, P_merged]
    auto const& c_int    = std::get<StabilizerTableau>(interior.front());
    auto& p_merged       = std::get<std::vector<PauliRotation>>(interior.back());

    // Front boundary: leading Clifford blocks, then the interior residual Clifford,
    // with the ancilla Hadamard prep prepended (each ancilla starts |+>).
    StabilizerTableau front{n + k};
    for (size_t i = 0; i < *first_rot; ++i) {
        front.apply(extract_clifford_operators(std::get<StabilizerTableau>(ext[i])));
    }
    front.apply(extract_clifford_operators(c_int));
    for (size_t a = n; a < n + k; ++a) front.prepend_h(a);

    // Back boundary: trailing Clifford blocks, then the ancilla Hadamard uncompute.
    StabilizerTableau back{n + k};
    for (size_t i = *last_rot + 1; i < ext.size(); ++i) {
        back.apply(extract_clifford_operators(std::get<StabilizerTableau>(ext[i])));
    }
    for (size_t a = n; a < n + k; ++a) back.h(a);

    tableau = Tableau{front, p_merged, back};
    return true;
}

}  // namespace qsyn::tableau
