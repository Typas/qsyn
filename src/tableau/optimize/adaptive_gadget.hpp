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

/**
 * @brief Gadgetize every internal Hadamard, returning a wider circuit with one ancilla per gadget.
 *
 * An internal Hadamard has a non-Clifford rotation both before and after it; leading/trailing
 * Hadamards are left in place. The fixed Clifford gadget body leaves the non-Clifford region
 * diagonal over the extended qubit set. Ancillae are kept (initialized |0>), not uncomputed.
 */
qcir::QCir gadgetize_all_internal_hadamards(qcir::QCir const& qcir);

/**
 * @brief Synthesize, gadgetize every internal Hadamard, and read back as a (wider) tableau.
 * @return false (tableau unchanged) if the round-trip synthesis fails.
 */
bool gadgetize_internal_hadamards(Tableau& tableau);

/**
 * @brief Gadgetize only the internal Hadamards at the selected region boundaries, merging joined regions.
 * @param selected_boundaries non-empty region-pair indices (boundary i joins regions i and i+1);
 *        out-of-range indices are ignored. Passing every boundary equals the no-arg overload.
 * @return false (tableau unchanged) if there is nothing to gadgetize.
 */
bool gadgetize_internal_hadamards(Tableau& tableau, std::span<size_t const> selected_boundaries);

/**
 * @brief Plan and apply the gadgetization whose predicted phase-poly peak stays under budget_bytes.
 * @return false (tableau unchanged) if the budget is too small to gadgetize even one boundary, or
 *         there is nothing to gadgetize; the caller should report a "budget too small" error.
 */
bool gadgetize_within_budget(Tableau& tableau, size_t budget_bytes);

// memory model
// ---------------------------------------------------------------------------

/**
 * @brief Predicted peak RSS (bytes) of the gadgetize -> phasepoly pipeline at width n. Upper bound.
 *
 * The widened tableau holds s_clifford StabilizerTableau blocks and m_total rotations; the largest
 * region to optimize has m_region terms. Used predict-only: the planner keeps every plan under budget.
 */
size_t predict_peak_bytes(size_t n, size_t s_clifford, size_t m_total, size_t m_region);

}  // namespace qsyn::tableau
