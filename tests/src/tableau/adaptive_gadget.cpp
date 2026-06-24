#include "tableau/optimize/adaptive_gadget.hpp"

#include <fmt/format.h>

#include <catch2/catch_test_macros.hpp>
#include <cstddef>
#include <limits>
#include <optional>
#include <span>
#include <string>
#include <variant>
#include <array>
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

// predict_peak_bytes is a sum of exact sizeof / paper-dimension terms (no fitted magic number); the
// planner adds the live RSS baseline on top. It must (a) equal the explicit structural formula -- this
// is the "no magic numbers" guard: if anyone slips in an unexplained constant, the recomputation here
// diverges -- and (b) be monotonic, so more gadgetization (larger n / s_clifford / m) predicts a larger
// peak (the basis of the full > adaptive > no-gadget memory ordering).
TEST_CASE("predict_peak_bytes is the exact structural sum and is monotonic", "[tableau][adaptive-gadget]") {
    auto const expected = [](size_t n, size_t s, size_t mt, size_t mr) {
        auto const stab = (2 * n) * (2 * n + 1) / 8 + (2 * n) * size_t{32};  // bits + PauliProduct(32 B)
        auto const rot  = (2 * n + 1) / 8 + size_t{48};                      // bits + PauliRotation(48 B)
        auto const gad  = (3 * (s * stab + mt * rot) + 24 * (2 * n) * (2 * n)) / 100 * 175;  // 3 copies x1.75 frag
        auto const cn2  = n * (n - 1) / 2;
        auto const pp   = 2 * mr * (n + cn2) + mr * mr + n * mr;             // FastTODD L + augmented + A
        return std::max(gad, pp);
    };
    SECTION("equals the explicit sizeof/paper formula (no hidden constant)") {
        for (auto const [n, s, mt, mr] : std::vector<std::array<size_t, 4>>{
                 {16, 50, 200, 47}, {351, 200, 1019, 211}, {1115, 1145, 3517, 3517}}) {
            CHECK(predict_peak_bytes(n, s, mt, mr) == expected(n, s, mt, mr));
        }
    }
    SECTION("monotonic non-decreasing in n, s_clifford, and m_region") {
        CHECK(predict_peak_bytes(100, 50, 200, 100) >= predict_peak_bytes(50, 50, 200, 100));   // n
        CHECK(predict_peak_bytes(100, 80, 200, 100) >= predict_peak_bytes(100, 40, 200, 100));  // s_clifford
        CHECK(predict_peak_bytes(100, 50, 200, 400) >= predict_peak_bytes(100, 50, 200, 100));  // m_region
        CHECK(predict_peak_bytes(1115, 1145, 3517, 3517) > predict_peak_bytes(84, 40, 237, 100));
    }
}

namespace {
// Serialize a tableau to a canonical string for byte-identical comparison.
std::string serialize(Tableau const& t) {
    std::string s;
    for (auto const& sub : t) s += fmt::format("{}\n----\n", sub);
    return s;
}

// Build C0 . T(q0) . [H] . T(q0) . [H] . T(q0) over 2 qubits: 3 regions, 2 internal-H boundaries.
Tableau three_region_two_boundary() {
    auto c0 = StabilizerTableau{2};
    auto p0 = std::vector<PauliRotation>{PauliRotation(PauliProduct("ZI"), Phase(1, 4))};
    auto c1 = StabilizerTableau{2};
    c1.h(0);
    auto p1 = std::vector<PauliRotation>{PauliRotation(PauliProduct("ZI"), Phase(1, 4))};
    auto c2 = StabilizerTableau{2};
    c2.h(0);
    auto p2 = std::vector<PauliRotation>{PauliRotation(PauliProduct("ZI"), Phase(1, 4))};
    return Tableau{c0, p0, c1, p1, c2, p2};
}
}  // namespace

// GOLDEN INVARIANT: selecting EVERY boundary must reproduce the no-arg (gadgetize-everything)
// result byte-for-byte. This pins the selective overload against the proven no-arg transform.
TEST_CASE("selective gadgetize with all boundaries equals the no-arg overload", "[tableau][adaptive-gadget]") {
    auto full = three_region_two_boundary();
    auto sel  = full;  // copy

    REQUIRE(gadgetize_internal_hadamards(full));  // no-arg: gadgetize everything

    auto const boundaries = std::vector<size_t>{0, 1};  // both boundaries
    REQUIRE(gadgetize_internal_hadamards(sel, std::span<size_t const>{boundaries}));

    CHECK(full.n_qubits() == sel.n_qubits());
    CHECK(full.n_pauli_rotations() == sel.n_pauli_rotations());
    CHECK(serialize(full) == serialize(sel));
}

// Selecting a SUBSET gadgetizes only those boundaries: merged where selected, separate elsewhere.
TEST_CASE("selective gadgetize with a subset leaves unselected boundaries separate", "[tableau][adaptive-gadget]") {
    SECTION("select only boundary 0 -> merge regions 0,1; region 2 stays separate") {
        auto t = three_region_two_boundary();
        auto const boundaries = std::vector<size_t>{0};
        REQUIRE(gadgetize_internal_hadamards(t, std::span<size_t const>{boundaries}));

        CHECK(t.n_qubits() == 3);                 // one ancilla (boundary 0 had one H)
        CHECK(all_rotations_diagonal(t));         // every region diagonal
        // boundary 1's Hadamard is NOT gadgetized -> it remains an internal-block Hadamard
        CHECK(count_internal_block_hadamards(t) == 1);
        // still has >1 region (regions not all merged into one)
        CHECK(t.n_pauli_rotations() >= 2);
    }
    SECTION("empty / out-of-range selections do nothing and report false") {
        auto t = three_region_two_boundary();
        CHECK_FALSE(gadgetize_internal_hadamards(t, std::span<size_t const>{}));
        auto const oob = std::vector<size_t>{99};
        CHECK_FALSE(gadgetize_internal_hadamards(t, std::span<size_t const>{oob}));
    }
}

// gadgetize_within_budget: planner + selective transform end to end.
TEST_CASE("gadgetize_within_budget gadgetizes within the budget or reports false", "[tableau][adaptive-gadget]") {
    SECTION("a generous budget reproduces gadgetize-everything") {
        auto t    = three_region_two_boundary();
        auto full = t;
        REQUIRE(gadgetize_within_budget(t, static_cast<size_t>(8) * 1024 * 1024 * 1024));  // 8 GiB
        REQUIRE(gadgetize_internal_hadamards(full));
        CHECK(t.n_qubits() == full.n_qubits());
        CHECK(serialize(t) == serialize(full));
    }
    SECTION("a tiny budget gadgetizes nothing and returns false (tableau unchanged)") {
        auto t = three_region_two_boundary();
        CHECK_FALSE(gadgetize_within_budget(t, static_cast<size_t>(1)));
        CHECK(t.n_qubits() == 2);
        CHECK(count_internal_block_hadamards(t) == 2);  // both internal Hadamards still present
    }
}
