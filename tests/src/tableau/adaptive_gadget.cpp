#include "tableau/optimize/adaptive_gadget.hpp"

#include <catch2/catch_test_macros.hpp>
#include <cstddef>
#include <limits>
#include <optional>
#include <variant>
#include <vector>

#include "convert/qcir_to_tableau.hpp"
#include "convert/tableau_to_qcir.hpp"
#include "qcir/basic_gate_type.hpp"
#include "qcir/operation.hpp"
#include "qcir/qcir.hpp"
#include "qcir/qcir_equiv.hpp"
#include "qcir/qcir_io.hpp"
#include "tableau/pauli_rotation.hpp"
#include "tableau/stabilizer_tableau.hpp"
#include "tableau/tableau.hpp"
#include "util/phase.hpp"

using namespace qsyn::tableau;
using qsyn::qcir::CXGate;
using qsyn::qcir::HGate;
using qsyn::qcir::PZGate;
using qsyn::qcir::QCir;
using dvlab::Phase;

namespace {

size_t count_hadamards(QCir const& qcir) {
    size_t count = 0;
    for (auto const& gate : qcir.get_gates()) {
        if (gate->get_operation().get_type() == "h") ++count;
    }
    return count;
}

// Hadamards left in the interior of a data wire: an H on wire w that has a
// rotation both before and after it on that same wire. This is robust to the
// topological tie-breaking of get_gates() (it only depends on per-wire order)
// and captures the property the gadgetizer must establish: no data wire carries
// an interior Hadamard. Ancilla prep/uncompute Hadamards live on rotation-free
// wires, so they are never counted.
size_t count_interior_wire_hadamards(QCir const& qcir) {
    constexpr size_t k_none = std::numeric_limits<size_t>::max();
    std::vector<size_t> first_rot(qcir.get_num_qubits(), k_none);
    std::vector<size_t> last_rot(qcir.get_num_qubits(), k_none);

    size_t idx = 0;
    for (auto const& gate : qcir.get_gates()) {
        if (!is_clifford(gate->get_operation())) {
            for (auto const q : gate->get_qubits()) {
                if (first_rot[q] == k_none) first_rot[q] = idx;
                last_rot[q] = idx;
            }
        }
        ++idx;
    }

    size_t count = 0;
    idx          = 0;
    for (auto const& gate : qcir.get_gates()) {
        if (gate->get_operation().get_type() == "h") {
            auto const w = gate->get_qubits()[0];
            if (first_rot[w] != k_none && first_rot[w] < idx && idx < last_rot[w]) ++count;
        }
        ++idx;
    }
    return count;
}

}  // namespace

TEST_CASE("gadgetize_all_internal_hadamards removes interior Hadamards", "[tableau][adaptive-gadget]") {
    SECTION("an internal H (T . H . T) is gadgetized onto an ancilla") {
        QCir qcir{1};
        qcir.append(PZGate(Phase(1, 4)), {0});  // T
        qcir.append(HGate(), {0});              // internal H
        qcir.append(PZGate(Phase(1, 4)), {0});  // T

        auto const out = gadgetize_all_internal_hadamards(qcir);

        REQUIRE(out.get_num_qubits() == 2);          // one ancilla added
        REQUIRE(count_interior_wire_hadamards(out) == 0);  // interior H absorbed by the gadget
        REQUIRE(count_hadamards(out) == 2);           // only the ancilla prep + uncompute H remain
    }

    SECTION("a leading H (H . T) is left in place") {
        QCir qcir{1};
        qcir.append(HGate(), {0});              // leading H, not internal
        qcir.append(PZGate(Phase(1, 4)), {0});  // T

        auto const out = gadgetize_all_internal_hadamards(qcir);

        REQUIRE(out.get_num_qubits() == 1);   // no ancilla
        REQUIRE(count_hadamards(out) == 1);   // H untouched
    }

    SECTION("a trailing H (T . H) is left in place") {
        QCir qcir{1};
        qcir.append(PZGate(Phase(1, 4)), {0});  // T
        qcir.append(HGate(), {0});              // trailing H, not internal

        auto const out = gadgetize_all_internal_hadamards(qcir);

        REQUIRE(out.get_num_qubits() == 1);
        REQUIRE(count_hadamards(out) == 1);
    }

    SECTION("two internal Hadamards allocate two ancillae") {
        QCir qcir{1};
        qcir.append(PZGate(Phase(1, 4)), {0});  // T
        qcir.append(HGate(), {0});              // internal H #1
        qcir.append(PZGate(Phase(1, 4)), {0});  // T
        qcir.append(HGate(), {0});              // internal H #2
        qcir.append(PZGate(Phase(1, 4)), {0});  // T

        auto const out = gadgetize_all_internal_hadamards(qcir);

        REQUIRE(out.get_num_qubits() == 3);           // two ancillae
        REQUIRE(count_interior_wire_hadamards(out) == 0);
        REQUIRE(count_hadamards(out) == 4);           // 2 prep + 2 uncompute
    }
}

// NOTE: a unitary qcir equivalence check is impossible for a gadgetized circuit
// (the Vandaele gadget is measurement-based, so the result is equivalent to the
// original only under measurement + feedforward, never as a bare unitary). The
// rigorous-feasible correctness checks are the structural invariant below (the
// merged region MUST be diagonal -- a non-diagonal region means the gadget left a
// Hadamard in the T-region, which FastTODD then silently skips) plus reference
// T-count matching against Vandaele in the golden dofiles.

namespace {

bool all_rotations_diagonal(Tableau const& tableau) {
    for (auto const& subtableau : tableau) {
        if (auto const* pr = std::get_if<std::vector<PauliRotation>>(&subtableau)) {
            for (auto const& rotation : *pr) {
                if (!rotation.is_diagonal()) return false;
            }
        }
    }
    return true;
}

}  // namespace

// Count internal Hadamards in a tableau: a Hadamard inside a Clifford block that
// sits strictly between the first and last phase-poly region (the same notion the
// gadgetizer eliminates). After gadgetization this must be 0.
size_t count_internal_block_hadamards(Tableau const& tableau) {
    std::optional<size_t> first_rot;
    std::optional<size_t> last_rot;
    for (size_t i = 0; i < tableau.size(); ++i) {
        if (std::holds_alternative<std::vector<PauliRotation>>(tableau[i])) {
            if (!first_rot) first_rot = i;
            last_rot = i;
        }
    }
    if (!first_rot) return 0;
    size_t count = 0;
    for (size_t i = *first_rot + 1; i < *last_rot; ++i) {
        if (auto const* st = std::get_if<StabilizerTableau>(&tableau[i])) {
            for (auto const& op : extract_clifford_operators(*st)) {
                if (op.first == CliffordOperatorType::h) ++count;
            }
        }
    }
    return count;
}

TEST_CASE("gadgetize_internal_hadamards eliminates internal Hadamards", "[tableau][adaptive-gadget]") {
    SECTION("an internal-block Hadamard is absorbed onto an ancilla; region stays diagonal") {
        // Directly build C0 . T(q0) . [H(q0)] . T(q0) over 2 qubits: the Hadamard
        // sits in an internal Clifford block (a rotation region before and after).
        auto c0 = StabilizerTableau{2};
        auto p0 = std::vector<PauliRotation>{PauliRotation(PauliProduct("ZI"), Phase(1, 4))};
        auto c1 = StabilizerTableau{2};
        c1.h(0);  // the internal Hadamard
        auto p1 = std::vector<PauliRotation>{PauliRotation(PauliProduct("ZI"), Phase(1, 4))};
        auto tableau = Tableau{c0, p0, c1, p1};

        REQUIRE(count_internal_block_hadamards(tableau) == 1);  // precondition

        REQUIRE(gadgetize_internal_hadamards(tableau));

        REQUIRE(tableau.n_qubits() == 3);                       // one ancilla allocated
        REQUIRE(count_internal_block_hadamards(tableau) == 0);  // internal H eliminated
        REQUIRE(all_rotations_diagonal(tableau));               // merged region diagonal
    }

    SECTION("a Hadamard-free tableau is left untouched (no spurious ancillae)") {
        auto c0 = StabilizerTableau{2};
        auto p0 = std::vector<PauliRotation>{
            PauliRotation(PauliProduct("ZI"), Phase(1, 4)),
            PauliRotation(PauliProduct("ZZ"), Phase(1, 4))};
        auto tableau = Tableau{c0, p0};

        REQUIRE(count_internal_block_hadamards(tableau) == 0);  // nothing to gadgetize

        REQUIRE(gadgetize_internal_hadamards(tableau));

        REQUIRE(tableau.n_qubits() == 2);          // no-op
        REQUIRE(all_rotations_diagonal(tableau));
    }
}

// End-to-end gadgetization quality against committed circuits (run from the repo
// root, as `make unit-test` does). Reads tests/tableau/circuits/*.qc directly, so
// it is reproducible from the commit -- unlike dof+ref goldens, whose .log
// references are gitignored. (FastTODD's own correctness/quality is covered in
// fasttodd.cpp.)

namespace {

// T-count of (tmerge -> hopt -> [gadgetize] -> fasttodd) on the synthesized
// circuit.
size_t pipeline_t_count(qsyn::qcir::QCir const& qcir, bool gadgetize) {
    auto tableau = qsyn::tableau::to_tableau(qcir);
    REQUIRE(tableau.has_value());
    merge_rotations(*tableau);
    minimize_internal_hadamards(*tableau);
    if (gadgetize) REQUIRE(gadgetize_internal_hadamards(*tableau));
    optimize_phase_polynomial(*tableau, FastToddPhasePolynomialOptimizationStrategy{});
    auto synthesized = qsyn::tableau::to_qcir(
        *tableau, HOptSynthesisStrategy{}, NaivePauliRotationsSynthesisStrategy{});
    REQUIRE(synthesized.has_value());
    return qsyn::qcir::get_gate_statistics(*synthesized).at("t-family");
}

}  // namespace

TEST_CASE("gadgetization never increases the T-count over the no-gadget baseline", "[tableau][adaptive-gadget]") {
    auto const paths = {
        "tests/tableau/circuits/tof_3.qc",
        "tests/tableau/circuits/gf2^5_mult.qc",  // no internal H -> gadgetize is a no-op
    };
    for (auto const* path : paths) {
        auto const qcir = qsyn::qcir::from_file(path);
        REQUIRE(qcir.has_value());

        REQUIRE(pipeline_t_count(*qcir, /*gadgetize=*/true) <=
                pipeline_t_count(*qcir, /*gadgetize=*/false));
    }
}
