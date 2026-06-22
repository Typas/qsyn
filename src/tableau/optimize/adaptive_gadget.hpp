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

}  // namespace qsyn::tableau
