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

}  // namespace qsyn::tableau
