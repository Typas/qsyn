/****************************************************************************
  PackageName  [ tableau/optimize ]
  Synopsis     [ Partial-Hadamard gadgetization (gadgetize step) ]
  Copyright    [ Copyright(c) 2024 DVLab, GIEE, NTU, Taiwan ]
****************************************************************************/

#include "tableau/optimize/adaptive_gadget.hpp"

#include <spdlog/spdlog.h>

#include <unistd.h>

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdio>
#include <cstdlib>
#include <limits>
#include <numeric>
#include <optional>
#include <queue>
#include <span>
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

// Emit the fixed Hadamard gadget body onto data wire q and ancilla a:
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

// Splice the fixed Hadamard gadget body in place of a Hadamard on data wire q,
// using fresh ancilla a (minus the prep/uncompute H, which are handled at the
// boundaries):
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

// One gadgetized merge-group. `clifford` is the residual front Clifford from the bounded merge
// (c_int); `region` is the merged diagonal phase polynomial.
struct GadgetizedSpan {
    StabilizerTableau clifford;
    std::vector<PauliRotation> region;
};

// Gadgetize the internal Clifford blocks strictly between region subtableau indices `lo` and `hi`
// of the already-extended tableau `ext`, splicing each Hadamard onto a fresh ancilla drawn from
// `anc_next` (advanced as ancillae are consumed), then collapse [lo..hi] into one diagonal region.
// Factored out so the gadgetize-everything path and the selective path share identical merge logic.
GadgetizedSpan gadgetize_span(Tableau& ext, size_t lo, size_t hi, size_t& anc_next) {
    auto const n_total = ext.n_qubits();
    for (size_t i = lo + 1; i < hi; ++i) {
        auto* st = std::get_if<StabilizerTableau>(&ext[i]);
        if (st == nullptr) continue;
        CliffordOperatorString spliced;
        for (auto const& op : extract_clifford_operators(*st)) {
            if (op.first == CliffordOperatorType::h) {
                append_gadget_body(spliced, op.second[0], anc_next++);
            } else {
                spliced.push_back(op);
            }
        }
        *st = StabilizerTableau{n_total}.apply(spliced);
    }
    Tableau interior{n_total};
    // Move (not copy) this group's subtableaux out of `ext` -- they are not read again elsewhere
    // (front/back/separators only touch ext indices outside [lo..hi]), so `ext` shrinks as groups are
    // absorbed instead of staying fully resident beside `interior`. Halves the gadgetize peak.
    for (size_t i = lo; i <= hi; ++i) interior.push_back(std::move(ext[i]));
    collapse(interior);  // -> [C_int, P_merged]
    return GadgetizedSpan{std::get<StabilizerTableau>(interior.front()),
                          std::get<std::vector<PauliRotation>>(interior.back())};
}

}  // namespace

bool gadgetize_internal_hadamards(Tableau& tableau) {
    if (tableau.is_empty()) return true;
    auto const n = tableau.n_qubits();

    // Structural internal Hadamards live in Clifford blocks strictly between the
    // first and last phase-polynomial region. Tableau-native: no circuit
    // round-trip (which is lossy for the optimizer), rotations stay symbolic.
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

    // Widen to n + k idle ancillae and gadgetize+merge the whole interior into one region.
    Tableau ext = tableau.extended_with_ancillae(k);
    size_t next_ancilla = n;
    auto const span     = gadgetize_span(ext, *first_rot, *last_rot, next_ancilla);

    // Front boundary: leading Clifford blocks, then the interior residual Clifford,
    // with the ancilla Hadamard prep prepended (each ancilla starts |+>).
    StabilizerTableau front{n + k};
    for (size_t i = 0; i < *first_rot; ++i) {
        front.apply(extract_clifford_operators(std::get<StabilizerTableau>(ext[i])));
    }
    front.apply(extract_clifford_operators(span.clifford));
    for (size_t a = n; a < n + k; ++a) front.prepend_h(a);

    // Back boundary: trailing Clifford blocks, then the ancilla Hadamard uncompute.
    StabilizerTableau back{n + k};
    for (size_t i = *last_rot + 1; i < ext.size(); ++i) {
        back.apply(extract_clifford_operators(std::get<StabilizerTableau>(ext[i])));
    }
    for (size_t a = n; a < n + k; ++a) back.h(a);

    tableau = Tableau{front, span.region, back};
    return true;
}

bool gadgetize_internal_hadamards(Tableau& tableau, std::span<size_t const> selected_boundaries) {
    if (tableau.is_empty() || selected_boundaries.empty()) return false;
    auto const n = tableau.n_qubits();

    // Region subtableau indices, in pipeline order. Boundary b sits between regions[b], regions[b+1].
    std::vector<size_t> regions;
    for (size_t i = 0; i < tableau.size(); ++i) {
        if (std::holds_alternative<std::vector<PauliRotation>>(tableau[i])) regions.push_back(i);
    }
    if (regions.size() < 2) return false;  // no boundary to merge across

    std::unordered_set<size_t> selected;
    for (auto const b : selected_boundaries) {
        if (b + 1 < regions.size()) selected.insert(b);  // ignore out-of-range
    }
    if (selected.empty()) return false;

    // Ancilla count = #internal Hadamards in the Clifford blocks of the selected boundaries.
    size_t k = 0;
    for (auto const b : selected) {
        for (size_t i = regions[b] + 1; i < regions[b + 1]; ++i) {
            if (auto const* st = std::get_if<StabilizerTableau>(&tableau[i])) k += count_hadamards(*st);
        }
    }
    if (k == 0) return false;  // selected boundaries carry no internal Hadamards

    Tableau ext         = tableau.extended_with_ancillae(k);
    size_t next_ancilla = n;

    // Merge-groups: maximal runs of regions joined by selected boundaries. Group spans
    // region-indices [a..b]; unselected boundaries separate groups.
    std::vector<std::pair<size_t, size_t>> groups;
    for (size_t a = 0; a < regions.size();) {
        size_t b = a;
        while (b + 1 < regions.size() && selected.contains(b)) ++b;
        groups.emplace_back(a, b);
        a = b + 1;
    }

    // Front boundary: leading Clifford blocks before the first region, then the first group's
    // residual Clifford, with all ancilla Hadamard preps prepended (each ancilla starts |+>).
    StabilizerTableau front{n + k};
    for (size_t i = 0; i < regions.front(); ++i) {
        front.apply(extract_clifford_operators(std::get<StabilizerTableau>(ext[i])));
    }

    std::vector<SubTableau> middle;  // separators + merged regions, between front and back
    for (size_t g = 0; g < groups.size(); ++g) {
        auto const [a, b] = groups[g];
        auto span         = gadgetize_span(ext, regions[a], regions[b], next_ancilla);
        if (g == 0) {
            front.apply(extract_clifford_operators(span.clifford));
        } else {
            // Separator = unselected boundary Clifford(s) between the previous group's last region
            // and this group's first region, then this group's residual Clifford.
            StabilizerTableau sep{n + k};
            for (size_t i = regions[groups[g - 1].second] + 1; i < regions[a]; ++i) {
                if (auto const* st = std::get_if<StabilizerTableau>(&ext[i])) {
                    sep.apply(extract_clifford_operators(*st));
                }
            }
            sep.apply(extract_clifford_operators(span.clifford));
            middle.emplace_back(std::move(sep));
        }
        middle.emplace_back(std::move(span.region));
    }
    for (size_t a = n; a < n + k; ++a) front.prepend_h(a);

    // Back boundary: trailing Clifford blocks after the last region, then ancilla uncompute.
    StabilizerTableau back{n + k};
    for (size_t i = regions.back() + 1; i < ext.size(); ++i) {
        back.apply(extract_clifford_operators(std::get<StabilizerTableau>(ext[i])));
    }
    for (size_t a = n; a < n + k; ++a) back.h(a);

    // Assemble: [front] [separator, region]... [back]. Move (not copy) each `middle` element into
    // `out`: the result is identical, but it avoids holding a second full-size copy of the widened
    // subtableaux at peak (copying here would add a third full resident copy of the widened set).
    Tableau out{n + k};
    out.erase(out.begin(), out.end());
    out.push_back(SubTableau{std::move(front)});
    for (auto& st : middle) out.push_back(std::move(st));
    out.push_back(SubTableau{std::move(back)});
    tableau = std::move(out);
    return true;
}

// Live resident-set size of this process, read from /proc/self/statm (Linux). The planner uses this as
// the baseline floor instead of a constant: the resident input (QCir DAG + un-widened tableau) varies
// by orders of magnitude with circuit size, so a fixed floor catastrophically under-counts big inputs.
// This is a MEASUREMENT, not a fitted number. Returns 0 off Linux (then the floor is just 0 -> the
// structural terms still bound the growth; macOS is predict-only without an OS backstop anyway).
size_t current_rss_bytes() {
    std::FILE* f = std::fopen("/proc/self/statm", "r");
    if (f == nullptr) return 0;
    unsigned long sz = 0;
    unsigned long rss = 0;
    auto const matched = std::fscanf(f, "%lu %lu", &sz, &rss);
    std::fclose(f);
    return matched == 2 ? static_cast<size_t>(rss) * static_cast<size_t>(sysconf(_SC_PAGESIZE)) : 0;
}

// Structural peak-RSS GROWTH (bytes) over the plan-time baseline, for the whole gadgetize -> phasepoly
// pipeline of a plan that widens the tableau to `n` qubits, where the widened tableau holds `s_clifford`
// StabilizerTableau blocks + `m_total` resident PauliRotations, and the largest single region to be
// phasepoly-optimized has `m_region` terms. The planner adds current_rss_bytes() as the baseline.
//
// EVERY constant below is a sizeof or a paper dimension (no fitted "engineering coefficient"):
//   - StabilizerTableau = 2n PauliProducts, each a sul::dynamic_bitset of (2n+1) bits  -> 2n*(2n+1)/8 B
//     (PauliProduct::_bitset is bit-packed, 1 bit/bit; stabilizer_tableau.hpp / pauli_rotation.hpp).
//   - PauliRotation     = one (2n+1)-bit bitset                                         -> (2n+1)/8 B.
//   - CliffordOperator  = pair<u8, array<size_t,2>>                                     -> 24 B.
//   - FastTODD L-matrix = BooleanMatrix (Row = std::vector<unsigned char>, 1 byte/bit), rows = terms,
//     cols = n + C(n,2)  (fastTodd.cpp build_l_transpose_row_from_term reserves n + n(n-1)/2); the
//     l_matrix + l_matrix_transpose pair (x2) + augmented=identity(m) (m^2).
//
// sizeof a PauliProduct (== its sole member sul::dynamic_bitset, measured 32 B) and a PauliRotation
// (PauliProduct + dvlab::Phase, ~48 B). These struct headers dominate the bit payload at small n and
// were the source of an earlier small-circuit under-prediction, so they are kept as explicit terms.
constexpr size_t kPauliProductBytes  = 32;
constexpr size_t kPauliRotationBytes = 48;

// FIXME: not well optimized for large circuits. This fixed envelope over-predicts more as n grows,
// so utilization (actual/budget) falls at large budgets and big budgets gadgetize fewer Hadamards than
// would fit. Safe (only over-predicts, never OOMs), but a tight fix needs an n-dependent multiplier.
// Collapse runtime also grows ~ s_clifford*n^2.
//
// STRUCTURAL: the selective gadgetize transform keeps every Clifford block as a full-width
// StabilizerTableau, and at peak holds ~3 coexisting full-width copies of that set -- `ext`
// (extended_with_ancillae, being consumed group-by-group), `middle` (the accumulated output
// separators+regions), and `out` (assembled from middle).
constexpr size_t kWidenCopies = 3;
// The one allocator-dependent knob: real RSS = structural bytes x allocator overhead, which no sizeof
// captures and which varies by allocator. User-overridable via env QSYN_GADGET_HEADROOM (a percent).
// Default 175: under the threaded phasepoly executor, glibc retains freed gadgetize arenas, so the
// gadgetize and phasepoly stages sum in RSS instead of max(), and a lower value can leave predicted <
// actual (unsafe) at high thread counts. 175 restores predicted >= actual on glibc; jemalloc/mimalloc
// are safe lower. The glibc margin compresses as the budget grows, so very large budgets should pin a
// higher QSYN_GADGET_HEADROOM; a malloc_trim() at the gadgetize->phasepoly boundary (return arenas,
// break the stage-sum) is the structural alternative to padding -- deferred.
size_t gadget_alloc_headroom_percent() {
    static size_t const pct = []() -> size_t {
        if (char const* e = std::getenv("QSYN_GADGET_HEADROOM")) {
            if (auto const v = std::strtoul(e, nullptr, 10); v > 0) return static_cast<size_t>(v);
        }
        return 175;
    }();
    return pct;
}

size_t predict_peak_bytes(size_t n, size_t s_clifford, size_t m_total, size_t m_region) {
    // One widened StabilizerTableau: 2n PauliProducts, each (2n+1) bits packed PLUS the 32 B struct.
    auto const stab_bytes = (2 * n) * (2 * n + 1) / 8 + (2 * n) * kPauliProductBytes;
    // One widened PauliRotation: (2n+1) bits packed PLUS the struct.
    auto const rot_bytes = (2 * n + 1) / 8 + kPauliRotationBytes;

    // GADGETIZE stage: kWidenCopies resident copies of the width-n tableau (s_clifford StabilizerTableaux
    // + m_total rotations) + synthesize() temporaries (CliffordOperatorString O((2n)^2) x 24 B), all
    // inflated by allocator fragmentation.
    auto const live_gadgetize = kWidenCopies * (s_clifford * stab_bytes + m_total * rot_bytes)
                              + 24 * (2 * n) * (2 * n);
    auto const gadgetize = live_gadgetize / 100 * gadget_alloc_headroom_percent();

    // Structural growth only; the planner adds the live plan-time RSS as the baseline floor.
    return std::max(gadgetize, phasepoly_region_bytes(n, m_region));
}

size_t phasepoly_region_bytes(size_t n, size_t m_region) {
    // PHASEPOLY stage (FastTODD, BooleanMatrix::Row = 1 byte/bit): a region runs at the global width n,
    // so its peak is the region's L-matrices + augmented + phase-poly matrix.
    auto const c_n_2 = n * (n - 1) / 2;
    return 2 * m_region * (n + c_n_2) + m_region * m_region + n * m_region;
}

namespace {

// Greedy boundary selection: gadgetize the highest-benefit region boundaries whose merged region
// still fits the budget, filling toward the ceiling. Benefit proxy = shared z-active qubits between
// adjacent regions (combinable parity structure). Returns the chosen boundary indices (possibly
// empty if even one merge does not fit). Pure: no optimization runs.
std::vector<size_t> plan_boundaries(Tableau const& tableau, size_t budget_bytes) {
    // The budget itself is the ceiling (cgroup memory.max is the real hard enforcer; no fudge margin).
    // Headroom = budget minus the live resident baseline (QCir DAG + un-widened tableau already loaded);
    // predict_peak_bytes returns the structural GROWTH on top, so a plan fits iff growth <= headroom.
    auto const ceiling   = budget_bytes;
    auto const baseline  = current_rss_bytes();
    auto const base_n    = tableau.n_qubits();

    // s_clifford = number of StabilizerTableau blocks. extended_with_ancillae widens ALL of them to the
    // global width, so this (fixed) count is the driver of the gadgetize-stage WIDEN cost.
    size_t s_clifford = 0;
    for (size_t i = 0; i < tableau.size(); ++i) {
        if (std::holds_alternative<StabilizerTableau>(tableau[i])) ++s_clifford;
    }

    // Per-region term count and z-active qubit mask, in pipeline order.
    std::vector<size_t> region_m;
    std::vector<std::vector<bool>> region_z;
    std::vector<size_t> region_idx;
    for (size_t i = 0; i < tableau.size(); ++i) {
        auto const* pr = std::get_if<std::vector<PauliRotation>>(&tableau[i]);
        if (pr == nullptr) continue;
        region_idx.push_back(i);
        region_m.push_back(pr->size());
        std::vector<bool> z(base_n, false);
        for (auto const& rot : *pr) {
            for (size_t q = 0; q < base_n; ++q) {
                if (rot.pauli_product().is_z_set(q)) z[q] = true;
            }
        }
        region_z.push_back(std::move(z));
    }
    auto const R = region_idx.size();
    if (R < 2) return {};

    // boundary b is between region b and b+1. Hadamard cost + benefit proxy per boundary.
    std::vector<size_t> bh(R - 1, 0);
    std::vector<size_t> benefit(R - 1, 0);
    for (size_t b = 0; b + 1 < R; ++b) {
        for (size_t i = region_idx[b] + 1; i < region_idx[b + 1]; ++i) {
            if (auto const* st = std::get_if<StabilizerTableau>(&tableau[i])) bh[b] += count_hadamards(*st);
        }
        for (size_t q = 0; q < base_n; ++q) {
            if (region_z[b][q] && region_z[b + 1][q]) ++benefit[b];
        }
    }

    // prefix sums of region_m for O(1) group term counts.
    std::vector<size_t> prefix_m(R + 1, 0);
    for (size_t r = 0; r < R; ++r) prefix_m[r + 1] = prefix_m[r] + region_m[r];
    auto const m_total = prefix_m[R];  // all resident rotations (widen-resident block)

    // The widen is GLOBAL: extended_with_ancillae(sum of selected boundary Hadamards) widens the whole
    // tableau, so n_total grows with every selected boundary. Each region (merged or not) is later
    // phasepoly-optimized at n_total, so the phasepoly term uses the largest region's term count.
    std::vector<bool> selected(R - 1, false);
    size_t total_anc    = 0;  // accumulated ancillae (= internal Hadamards) of selected boundaries
    size_t max_region_m = *std::max_element(region_m.begin(), region_m.end());
    while (true) {
        bool found          = false;
        size_t best_b       = 0;
        size_t best_benefit = 0;
        size_t best_predict = std::numeric_limits<size_t>::max();
        for (size_t b = 0; b + 1 < R; ++b) {
            if (selected[b]) continue;
            // group that results from also selecting b: extend left/right over already-selected boundaries.
            size_t lo = b;
            while (lo > 0 && selected[lo - 1]) --lo;
            size_t hi = b + 1;
            while (hi + 1 < R && selected[hi]) ++hi;
            auto const cand_group_m = prefix_m[hi + 1] - prefix_m[lo];
            auto const cand_n       = base_n + total_anc + bh[b];  // selecting b adds only bh[b] new ancillae
            auto const cand_pred    = baseline + predict_peak_bytes(cand_n, s_clifford, m_total,
                                                                    std::max(max_region_m, cand_group_m));
            if (cand_pred > ceiling) continue;  // would not fit the budget
            if (!found || benefit[b] > best_benefit ||
                (benefit[b] == best_benefit && cand_pred < best_predict)) {
                found        = true;
                best_b       = b;
                best_benefit = benefit[b];
                best_predict = cand_pred;
            }
        }
        if (!found) break;
        selected[best_b] = true;
        total_anc += bh[best_b];
        // update largest-region term count with the group now containing best_b
        size_t lo = best_b;
        while (lo > 0 && selected[lo - 1]) --lo;
        size_t hi = best_b + 1;
        while (hi + 1 < R && selected[hi]) ++hi;
        max_region_m = std::max(max_region_m, prefix_m[hi + 1] - prefix_m[lo]);
    }

    std::vector<size_t> result;
    for (size_t b = 0; b + 1 < R; ++b) {
        if (selected[b]) result.push_back(b);
    }

    // Developer diagnostic (debug level only -- never shown at the default user-facing level, so it
    // does not surface the predicted-peak / baseline internals). Lets a verification run inspect the
    // planner's decision with `--verbose`.
    auto const total_h   = std::accumulate(bh.begin(), bh.end(), size_t{0});
    auto const predicted = baseline + predict_peak_bytes(base_n + total_anc, s_clifford, m_total, max_region_m);
    constexpr double mib = 1024.0 * 1024.0;
    spdlog::debug(
        "gadgetize budget {:.0f} MiB: selected {} of {} internal Hadamards (n {} -> {}, s_clifford {}, "
        "m_region {}); predicted peak {:.0f} MiB (baseline {:.0f} MiB)",
        static_cast<double>(budget_bytes) / mib, total_anc, total_h, base_n, base_n + total_anc,
        s_clifford, max_region_m,
        static_cast<double>(predicted) / mib, static_cast<double>(baseline) / mib);
    return result;
}

}  // namespace

bool gadgetize_within_budget(Tableau& tableau, size_t budget_bytes) {
    auto const selected = plan_boundaries(tableau, budget_bytes);
    if (selected.empty()) return false;  // budget too small / nothing to gadgetize
    return gadgetize_internal_hadamards(tableau, selected);
}

}  // namespace qsyn::tableau
