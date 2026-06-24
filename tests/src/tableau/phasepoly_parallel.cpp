#include <catch2/catch_test_macros.hpp>
#include <cstddef>
#include <optional>
#include <variant>
#include <vector>

#include "convert/qcir_to_tableau.hpp"
#include "convert/tableau_to_qcir.hpp"
#include "qcir/basic_gate_type.hpp"
#include "qcir/operation.hpp"
#include "qcir/qcir.hpp"
#include "qcir/qcir_equiv.hpp"
#include "tableau/pauli_rotation.hpp"
#include "tableau/stabilizer_tableau.hpp"
#include "tableau/tableau.hpp"
#include "tableau/tableau_optimization.hpp"
#include "util/phase.hpp"

using namespace qsyn::tableau;
using dvlab::Phase;
using qsyn::qcir::CXGate;
using qsyn::qcir::PZGate;
using qsyn::qcir::QCir;

namespace {

// A diagonal (phase-polynomial) circuit, so its tableau is a single rotation region.
QCir diagonal_circuit() {
    QCir qcir{2};
    for (int rep = 0; rep < 4; ++rep) {
        qcir.append(PZGate(Phase(1, 4)), {0});
        qcir.append(CXGate(), {0, 1});
        qcir.append(PZGate(Phase(1, 4)), {1});
    }
    return qcir;
}

// Tableau conversion collapses every Clifford to one side, so a circuit yields a
// single region. To get the multiple independent regions the parallel executor
// targets, replicate one region behind its own identity-Clifford separator.
Tableau multi_region_tableau(size_t copies) {
    auto base = qsyn::tableau::to_tableau(diagonal_circuit());
    REQUIRE(base.has_value());
    std::vector<PauliRotation> region;
    for (auto const& sub : *base) {
        if (auto const* r = std::get_if<std::vector<PauliRotation>>(&sub)) region = *r;
    }
    REQUIRE_FALSE(region.empty());
    auto const n = region.front().n_qubits();

    Tableau tableau{n};  // starts with one identity StabilizerTableau
    for (size_t i = 0; i < copies; ++i) {
        tableau.push_back(SubTableau{region});
        tableau.push_back(SubTableau{StabilizerTableau{n}});
    }
    return tableau;
}

size_t count_regions(Tableau const& tableau) {
    size_t n = 0;
    for (auto const& sub : tableau) {
        if (std::holds_alternative<std::vector<PauliRotation>>(sub)) ++n;
    }
    return n;
}

bool tableaux_equal(Tableau const& a, Tableau const& b) {
    if (a.size() != b.size()) return false;
    for (size_t i = 0; i < a.size(); ++i) {
        if (auto const* sa = std::get_if<StabilizerTableau>(&a[i])) {
            auto const* sb = std::get_if<StabilizerTableau>(&b[i]);
            if (sb == nullptr || *sa != *sb) return false;
        } else {
            auto const& ra = std::get<std::vector<PauliRotation>>(a[i]);
            auto const* rb = std::get_if<std::vector<PauliRotation>>(&b[i]);
            if (rb == nullptr || ra != *rb) return false;
        }
    }
    return true;
}

constexpr size_t k_ample_budget = size_t{1} << 30;  // 1 GiB: ample for the tiny test regions

}  // namespace

TEST_CASE("parallel phasepoly preserves the unitary", "[tableau][phasepoly-parallel]") {
    // No gadgetize, so the tableau stays a plain unitary -- the optimized circuit
    // must be the same unitary as before optimization.
    auto const original = multi_region_tableau(4);
    auto const before   = qsyn::tableau::to_qcir(
        original, HOptSynthesisStrategy{}, NaivePauliRotationsSynthesisStrategy{});
    REQUIRE(before.has_value());

    auto threaded = original;
    REQUIRE(optimize_phase_polynomial(
        threaded, FastToddPhasePolynomialOptimizationStrategy{}, 4, k_ample_budget));
    auto const after = qsyn::tableau::to_qcir(
        threaded, HOptSynthesisStrategy{}, NaivePauliRotationsSynthesisStrategy{});
    REQUIRE(after.has_value());

    REQUIRE(qsyn::qcir::is_equivalent(*before, *after));
}

TEST_CASE("parallel phasepoly result is byte-identical to sequential", "[tableau][phasepoly-parallel]") {
    auto const original = multi_region_tableau(4);
    REQUIRE(count_regions(original) > 1);  // confirms >1 independent unit runs in parallel

    auto sequential = original;
    optimize_phase_polynomial(sequential, FastToddPhasePolynomialOptimizationStrategy{});

    auto threaded            = original;
    auto const threaded_done = optimize_phase_polynomial(
        threaded, FastToddPhasePolynomialOptimizationStrategy{}, 4, k_ample_budget);

    REQUIRE(threaded_done);
    REQUIRE(tableaux_equal(sequential, threaded));
}

TEST_CASE("parallel phasepoly is thread-count-invariant", "[tableau][phasepoly-parallel]") {
    auto const original = multi_region_tableau(4);

    auto two = original;
    REQUIRE(optimize_phase_polynomial(two, FastToddPhasePolynomialOptimizationStrategy{}, 2, k_ample_budget));
    auto eight = original;
    REQUIRE(optimize_phase_polynomial(eight, FastToddPhasePolynomialOptimizationStrategy{}, 8, k_ample_budget));

    REQUIRE(tableaux_equal(two, eight));
}

TEST_CASE("parallel phasepoly refuses a budget too small for one region", "[tableau][phasepoly-parallel]") {
    auto tableau      = multi_region_tableau(4);
    auto const before = tableau;

    // A 1-byte budget cannot fit any region: refuse up front, leave the tableau untouched.
    auto const ok = optimize_phase_polynomial(
        tableau, FastToddPhasePolynomialOptimizationStrategy{}, 4, std::optional<size_t>{1});

    REQUIRE_FALSE(ok);
    REQUIRE(tableaux_equal(tableau, before));
}
