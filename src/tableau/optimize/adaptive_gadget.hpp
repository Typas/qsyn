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

qcir::QCir gadgetize_all_internal_hadamards(qcir::QCir const& qcir);

bool gadgetize_internal_hadamards(Tableau& tableau);

bool gadgetize_internal_hadamards(Tableau& tableau, std::span<size_t const> selected_boundaries);

bool gadgetize_within_budget(Tableau& tableau, size_t budget_bytes);

// memory model
// ---------------------------------------------------------------------------

size_t predict_peak_bytes(size_t n, size_t s_clifford, size_t m_total, size_t m_region);

size_t phasepoly_region_bytes(size_t n, size_t m_region);

}  // namespace qsyn::tableau
