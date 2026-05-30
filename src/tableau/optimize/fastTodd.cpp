/**
 * @file
 * @brief implementation of the fastTodd phase polynomial optimization
 * @copyright Copyright(c) 2024 DVLab, GIEE, NTU, Taiwan
 */

#include "../tableau_optimization.hpp"

namespace qsyn::tableau {
namespace {
using Polynomial = PhasePolynomialOptimizationStrategy::Polynomial;
}  // namespace

std::pair<StabilizerTableau, Polynomial> TohpePhasePolynomialOptimizationStrategy::optimize(StabilizerTableau const& clifford, Polynomial const& polynomial) const {
    return {clifford, polynomial};
}

std::pair<StabilizerTableau, Polynomial> FastToddPhasePolynomialOptimizationStrategy::optimize(StabilizerTableau const& clifford, Polynomial const& polynomial) const {
    return {clifford, polynomial};
}

}  // namespace qsyn::tableau
