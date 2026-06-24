/****************************************************************************
  PackageName  [ tableau/optimize ]
  Synopsis     [ Partial-Hadamard gadgetization.

                 Additive sibling of optimize_phase_polynomial. The existing
                 per-region strategy.optimize(clifford, polynomial) contract and
                 optimize_phase_polynomial are left untouched; this is a new,
                 separately-named path. ]
  Copyright    [ Copyright(c) 2024 DVLab, GIEE, NTU, Taiwan ]
****************************************************************************/

#pragma once

#include <cstddef>
#include <span>

#include "qcir/qcir.hpp"
#include "tableau/tableau.hpp"
#include "tableau/tableau_optimization.hpp"

namespace qsyn::tableau {

// gadget core
// ---------------------------------------------------------------------------

// Gadgetize *every* internal Hadamard in the circuit, mirroring the Vandaele
// oracle (quantum-circuit-optimization/src/circuit.rs:174 hadamard_gadgetization).
// An internal Hadamard is one with a non-Clifford (T-like) rotation both before
// and after it in the circuit (the global `i < last && flag` guard); leading and
// trailing Hadamards are left in place.
//
// Each gadgetized Hadamard on qubit q allocates a fresh ancilla a and emits the
// fixed Clifford gadget body
//     S(a) . S(q) . CX(q,a) . S(a) . Z(a) . CX(a,q) . CX(q,a)
// in place of the Hadamard. The data wire stays on q (no remap). Every ancilla
// gets a Hadamard prep before the whole body and a Hadamard uncompute after it
// (c_out = anc + body + anc). The gadget body contains no Hadamard, so the only
// Hadamards left on data wires are leading/trailing and the T-region is diagonal
// over the extended qubit set.
//
// The result keeps the ancillae (n + #internal Hadamards qubits); it is
// equivalent to the input with the ancillae initialized |0>. This matches the
// reference, which also keeps ancillae rather than uncomputing back to n qubits.
qcir::QCir gadgetize_all_internal_hadamards(qcir::QCir const& qcir);

// Tableau-level gadgetize step: synthesize the tableau to a circuit, gadgetize
// every internal Hadamard, and read it back as a tableau (now with ancillae).
// This is the standalone `gadgetize` pipeline step; the phase-polynomial
// optimizer (phasepoly fasttodd/tohpe/todd) runs as a separate, following step
// so the two levers stay independently selectable for the ablation study.
// Returns false (leaving the tableau unchanged) if the round-trip synthesis
// fails.
bool gadgetize_internal_hadamards(Tableau& tableau);

// Selective gadgetize: gadgetize only the internal Hadamards at the chosen region boundaries,
// leaving the rest as separate non-gadget regions. `selected_boundaries` indexes region *pairs*:
// boundary `i` sits between phase-poly region `i` and region `i+1` (0-based over the regions in
// pipeline order). Adjacent regions joined by selected boundaries are merged into one diagonal
// region (the same bounded-interior-merge the no-arg overload uses, scoped to that group);
// unselected boundaries stay as Clifford separators. Passing every boundary reproduces the no-arg
// (gadgetize-everything) result exactly.
//
// `selected_boundaries` must be NON-EMPTY (there is no "gadgetize nothing" mode -- for the
// non-gadget baseline, simply omit this step). Out-of-range indices are ignored. Returns false
// (leaving the tableau unchanged) if there is nothing to gadgetize.
bool gadgetize_internal_hadamards(Tableau& tableau, std::span<size_t const> selected_boundaries);

// Budget-aware gadgetize: plan which region boundaries to gadgetize so the predicted later
// phase-poly optimization of every resulting region stays under `budget_bytes` (with an internal
// safety margin), then apply the selective transform. Gadget-only -- the phase-poly optimizer runs
// as a separate following step, exactly like the no-arg path. Returns false (tableau unchanged) if
// the budget is too small to gadgetize even one boundary, or there is nothing to gadgetize -- the
// caller should report a "budget too small" error rather than silently producing a non-gadget
// result (for the non-gadget baseline, omit this step entirely).
bool gadgetize_within_budget(Tableau& tableau, size_t budget_bytes);

// memory model
// ---------------------------------------------------------------------------

// Predicted peak RSS (bytes) of the whole gadgetize -> phasepoly pipeline for a plan that widens the
// tableau to `n` qubits, where the widened tableau holds `s_clifford` StabilizerTableau blocks and
// `m_total` resident rotations, and the largest single region to be optimized has `m_region` terms.
//
// Profile-derived 2026-06-23 (RSS probes; see docs/adaptive-gadget-worklog.md). Every constant is a
// sizeof (StabilizerTableau bit-packed 2n*(2n+1)/8; PauliRotation (2n+1)/8; CliffordOperator 24 B;
// BooleanMatrix::Row 1 byte/bit) or a paper dimension (FastTODD L-matrix cols n+C(n,2),
// Heyfron2018/Vandaele2024). The one empirical input is the measured RSS floor; the only flagged
// estimate is the widen-copies factor (=2, upper bound). Used predict-only: the planner keeps every
// plan under the budget because exceeding it is an unrecoverable OOM kill (cgroup memory.max enforces).
size_t predict_peak_bytes(size_t n, size_t s_clifford, size_t m_total, size_t m_region);

// Predicted peak RSS (bytes) of optimizing ONE phase-poly region of `m_region` terms at width `n` --
// the FastTODD L-matrix + augmented + phase-poly matrix term of predict_peak_bytes, broken out for
// sizing a single region. Upper bound (predict >= actual).
size_t phasepoly_region_bytes(size_t n, size_t m_region);

}  // namespace qsyn::tableau
