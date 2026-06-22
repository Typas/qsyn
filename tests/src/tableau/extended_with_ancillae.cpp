#include <catch2/catch_test_macros.hpp>

#include "tableau/pauli_rotation.hpp"
#include "tableau/stabilizer_tableau.hpp"
#include "tableau/tableau.hpp"
#include "util/phase.hpp"

using namespace qsyn::tableau;
using dvlab::Phase;

TEST_CASE("PauliProduct::extended_with_ancillae appends identity qubits", "[tableau][substrate]") {
    SECTION("preserves paulis and sign, ancillae are identity") {
        auto const p   = PauliProduct("-XYZ");
        auto const ext = p.extended_with_ancillae(2);

        REQUIRE(ext.n_qubits() == 5);
        REQUIRE(ext == PauliProduct("-XYZII"));
    }

    SECTION("k = 0 is a no-op") {
        auto const p = PauliProduct("XZ");
        REQUIRE(p.extended_with_ancillae(0) == p);
    }

    SECTION("non-negative product stays non-negative") {
        auto const ext = PauliProduct("Y").extended_with_ancillae(1);
        REQUIRE(ext == PauliProduct("YI"));
        REQUIRE_FALSE(ext.is_neg());
    }
}

TEST_CASE("PauliRotation::extended_with_ancillae widens product and keeps phase", "[tableau][substrate]") {
    auto const rot = PauliRotation(PauliProduct("XZ"), Phase(1, 4));
    auto const ext = rot.extended_with_ancillae(2);

    REQUIRE(ext.n_qubits() == 4);
    REQUIRE(ext == PauliRotation(PauliProduct("XZII"), Phase(1, 4)));
    REQUIRE(ext.phase() == Phase(1, 4));
}

TEST_CASE("StabilizerTableau::extended_with_ancillae rebuilds with correct stab/destab split", "[tableau][substrate]") {
    SECTION("identity tableau extends to larger identity tableau") {
        auto const t = StabilizerTableau(3);
        REQUIRE(t.extended_with_ancillae(2) == StabilizerTableau(5));
    }

    SECTION("non-trivial rows keep their qubit, ancillae are Z/X identity") {
        // Build a non-identity 2-qubit tableau via a CX, then widen by 1.
        auto t = StabilizerTableau(2);
        t.cx(0, 1);
        auto const ext = t.extended_with_ancillae(1);

        REQUIRE(ext.n_qubits() == 3);

        // Original qubit rows are the originals widened by one identity qubit.
        REQUIRE(ext.stabilizer(0) == t.stabilizer(0).extended_with_ancillae(1));
        REQUIRE(ext.stabilizer(1) == t.stabilizer(1).extended_with_ancillae(1));
        REQUIRE(ext.destabilizer(0) == t.destabilizer(0).extended_with_ancillae(1));
        REQUIRE(ext.destabilizer(1) == t.destabilizer(1).extended_with_ancillae(1));

        // Fresh ancilla qubit is Z/X identity, matching a fresh tableau's row.
        REQUIRE(ext.stabilizer(2) == StabilizerTableau(3).stabilizer(2));
        REQUIRE(ext.destabilizer(2) == StabilizerTableau(3).destabilizer(2));
    }
}

TEST_CASE("Tableau::extended_with_ancillae widens all subtableaux and n_qubits", "[tableau][substrate]") {
    auto const clifford = StabilizerTableau(2);
    auto const rotations = std::vector<PauliRotation>{
        PauliRotation(PauliProduct("ZZ"), Phase(1, 4)),
        PauliRotation(PauliProduct("XI"), Phase(1, 2)),
    };
    auto const t = Tableau{clifford, rotations};

    auto const ext = t.extended_with_ancillae(1);

    REQUIRE(ext.n_qubits() == 3);
    REQUIRE(ext.size() == t.size());

    // Clifford block widened.
    REQUIRE(std::get<StabilizerTableau>(ext[0]) == clifford.extended_with_ancillae(1));

    // Rotation block widened term-by-term, phases preserved.
    auto const& ext_rots = std::get<std::vector<PauliRotation>>(ext[1]);
    REQUIRE(ext_rots.size() == rotations.size());
    REQUIRE(ext_rots[0] == rotations[0].extended_with_ancillae(1));
    REQUIRE(ext_rots[1] == rotations[1].extended_with_ancillae(1));
}
