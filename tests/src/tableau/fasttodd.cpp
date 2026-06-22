#include <catch2/catch_test_macros.hpp>
#include <cstddef>
#include <optional>

#include "convert/qcir_to_tableau.hpp"
#include "convert/tableau_to_qcir.hpp"
#include "qcir/qcir.hpp"
#include "qcir/qcir_equiv.hpp"
#include "qcir/qcir_io.hpp"
#include "tableau/stabilizer_tableau.hpp"
#include "tableau/tableau.hpp"
#include "tableau/tableau_optimization.hpp"

using namespace qsyn::tableau;

namespace {

// Run the phase-polynomial optimization pipeline (tmerge -> hopt -> phasepoly
// <strategy>) on a circuit and synthesize it back. No gadgetization, so the
// result is a unitary circuit equivalent to the input.
std::optional<qsyn::qcir::QCir> phasepoly_optimize(
    qsyn::qcir::QCir const& qcir,
    PhasePolynomialOptimizationStrategy const& strategy) {
    auto tableau = qsyn::tableau::to_tableau(qcir);
    if (!tableau) return std::nullopt;
    merge_rotations(*tableau);
    minimize_internal_hadamards(*tableau);
    optimize_phase_polynomial(*tableau, strategy);
    return qsyn::tableau::to_qcir(
        *tableau, HOptSynthesisStrategy{}, NaivePauliRotationsSynthesisStrategy{});
}

std::optional<qsyn::qcir::QCir> fasttodd_optimize(qsyn::qcir::QCir const& qcir) {
    return phasepoly_optimize(qcir, FastToddPhasePolynomialOptimizationStrategy{});
}

size_t t_count(qsyn::qcir::QCir const& qcir) {
    return qsyn::qcir::get_gate_statistics(qcir).at("t-family");
}

}  // namespace

TEST_CASE("FastTODD preserves circuit equivalence", "[tableau][fasttodd]") {
    // The optimized circuit must be the same unitary as the input.
    auto const qcir = qsyn::qcir::from_file("tests/tableau/circuits/tof_3.qc");
    REQUIRE(qcir.has_value());

    auto const optimized = fasttodd_optimize(*qcir);
    REQUIRE(optimized.has_value());

    REQUIRE(qsyn::qcir::is_equivalent(*qcir, *optimized));
}

TEST_CASE("TOHPE preserves circuit equivalence", "[tableau][fasttodd]") {
    auto const qcir = qsyn::qcir::from_file("tests/tableau/circuits/tof_3.qc");
    REQUIRE(qcir.has_value());

    auto const optimized =
        phasepoly_optimize(*qcir, TohpePhasePolynomialOptimizationStrategy{});
    REQUIRE(optimized.has_value());

    REQUIRE(qsyn::qcir::is_equivalent(*qcir, *optimized));
}

TEST_CASE("FastTODD reaches at most the paper's T-count", "[tableau][fasttodd]") {
    struct Case {
        char const* path;
        size_t paper_t;
    };
    // T-counts reported in Vandaele et al. 2024 for these circuits; FastTODD must
    // be at least as good (<=). Using <= (not ==) so an improvement never breaks
    // the test, while a regression above the paper does.
    auto const cases = {
        Case{"tests/tableau/circuits/tof_3.qc", 15},
        Case{"tests/tableau/circuits/tof_4.qc", 23},
        Case{"tests/tableau/circuits/tof_5.qc", 31},
        Case{"tests/tableau/circuits/barenco_tof_3.qc", 16},
    };
    for (auto const& c : cases) {
        auto const qcir = qsyn::qcir::from_file(c.path);
        REQUIRE(qcir.has_value());

        auto const optimized = fasttodd_optimize(*qcir);
        REQUIRE(optimized.has_value());

        REQUIRE(t_count(*optimized) <= c.paper_t);
    }
}
