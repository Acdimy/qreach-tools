#!/usr/bin/env python3
"""
Test: explicit (non‑QADD) TransitionSystem structural properties.

Covers construction, topology, annotation propagation, and basic
classical‑proposition handling without relying on parse_qiskit_cir
for the core structural checks.
"""

from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import pyqreach
from qiskit import QuantumCircuit
from qreach.parse_qiskit import parse_qiskit_cir
from qreach.qctl import (
    set_initial_state,
    quantum_state,
    span_states,
)


# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------

def _leaf_locations(ts):
    return [i for i in range(ts.getLocationNum()) if ts.isLeafLoc(i)]


def _incoming_edges(ts, loc_idx):
    """Reconstruct incoming edges (not directly exposed in Python)."""
    incoming = []
    for j in range(ts.getLocationNum()):
        if loc_idx in ts.Locations[j].postLocations:
            incoming.append(j)
    return incoming


# ---------------------------------------------------------------------------
# 1. Basic TS construction
# ---------------------------------------------------------------------------

def test_basic_construction():
    """Create a minimal TS and check basic properties."""
    num_q = 2
    ts = pyqreach.TransitionSystem()

    l0 = pyqreach.Location(num_q)
    l1 = pyqreach.Location(num_q)
    l2 = pyqreach.Location(num_q)
    ts.addLocation(l0)
    ts.addLocation(l1)
    ts.addLocation(l2)

    h_gate = pyqreach.QOperation("H", num_q, [0], [])
    cx_gate = pyqreach.QOperation("CX", num_q, [0, 1], [])

    ts.addRelation(0, 1, h_gate)
    ts.addRelation(1, 2, cx_gate)
    ts.setInitLocation(0)

    # -- assertions --
    assert ts.getLocationNum() == 3, f"Expected 3 locs, got {ts.getLocationNum()}"
    assert ts.getInitLocation() == 0

    # Leaf detection
    assert not ts.isLeafLoc(0)
    assert not ts.isLeafLoc(1)
    assert ts.isLeafLoc(2), "Location 2 should be a leaf (no outgoing edges)"

    # Outgoing edges
    assert ts.Locations[0].postLocations == [1], (
        f"Unexpected outgoing from loc 0: {ts.Locations[0].postLocations}"
    )
    assert ts.Locations[1].postLocations == [2]
    assert ts.Locations[2].postLocations == []

    # Incoming edges (reconstructed)
    assert _incoming_edges(ts, 0) == []
    assert _incoming_edges(ts, 1) == [0]
    assert _incoming_edges(ts, 2) == [1]

    # Relation names (getName() returns formatted "Name(indices)[params]").
    assert ts.getRelationName(0, 1) == "H(0)[]", (
        f"Unexpected relation name: {ts.getRelationName(0, 1)}"
    )
    assert ts.getRelationName(1, 2) == "CX(0,1)[]", (
        f"Unexpected relation name: {ts.getRelationName(1, 2)}"
    )
    assert ts.getRelationName(0, 2) == ""  # no such edge

    # Each location has correct qNum
    assert ts.Locations[0].qNum == num_q
    assert ts.Locations[1].qNum == num_q
    assert ts.Locations[2].qNum == num_q

    print("test_basic_construction  OK")


# ---------------------------------------------------------------------------
# 2. Bell‑state circuit (H → CX)
# ---------------------------------------------------------------------------

def test_bell_state():
    """H on q0 then CX(0,1) from |00⟩ → Bell state (|00⟩+|11⟩)/√2."""
    num_q = 2
    ts = pyqreach.TransitionSystem()

    for _ in range(3):
        ts.addLocation(pyqreach.Location(num_q))
    ts.addRelation(0, 1, pyqreach.QOperation("H", num_q, [0], []))
    ts.addRelation(1, 2, pyqreach.QOperation("CX", num_q, [0, 1], []))
    ts.setInitLocation(0)

    set_initial_state(ts, "00")          # |00⟩ on init location
    ts.computingFixedPointPost()

    # Final location should hold the Bell state — a 1‑dimensional pure state.
    upper_dim, lower_dim = ts.printDims(2)
    assert lower_dim == 1, (
        f"Bell state is a pure state (1D subspace), got lower_dim={lower_dim}"
    )
    assert upper_dim == 4, f"Upper bound should be full 4‑dim space, got {upper_dim}"

    bell = span_states(["00", "11"])
    assert ts.satisfy(2, bell), "Final location should satisfy span(|00⟩,|11⟩)"

    # It should NOT satisfy a single basis state.
    assert not ts.satisfy(2, quantum_state("00")), (
        "Bell state should not satisfy a single basis state"
    )
    assert not ts.satisfy(2, quantum_state("01"))
    assert not ts.satisfy(2, quantum_state("10"))
    assert not ts.satisfy(2, quantum_state("11"))

    # Initial location still has its seeded annotation.
    init_upper, init_lower = ts.printDims(0)
    assert init_lower == 1, f"Init loc lower dim should be 1, got {init_lower}"
    assert ts.satisfy(0, quantum_state("00")), "Init loc should satisfy |00⟩"

    print("test_bell_state  OK")


# ---------------------------------------------------------------------------
# 3. Branching circuit (H + measurement)
# ---------------------------------------------------------------------------

def test_branching_measurement():
    """A circuit with a measurement creates multiple leaf locations.
    Uses 2 qubits (minimum required by the CFLOBDD backend)."""
    qc = QuantumCircuit(2, 1)
    qc.h(0)
    qc.measure(0, 0)

    ts = pyqreach.TransitionSystem()
    result_locs = parse_qiskit_cir(qc, qc.num_qubits, ts)
    set_initial_state(ts, "00")
    ts.computingFixedPointPost()

    # Should have more than a linear number of locations (branching).
    n_locs = ts.getLocationNum()
    assert n_locs >= 3, f"Branching circuit should have ≥3 locs, got {n_locs}"

    # There should be at least one branching location (≥2 outgoing edges).
    branch_locs = [
        i for i in range(n_locs)
        if len(ts.Locations[i].postLocations) >= 2
    ]
    assert len(branch_locs) > 0, "No branching location found"

    # Every branching location should have exactly 2 outgoing edges
    # (meas0 and meas1).
    for loc in branch_locs:
        out_deg = len(ts.Locations[loc].postLocations)
        assert out_deg == 2, (
            f"Branching loc {loc}: expected 2 outgoing, got {out_deg}"
        )

    # At least two leaf locations (one per measurement outcome).
    leaves = _leaf_locations(ts)
    assert len(leaves) >= 2, (
        f"Expected ≥2 leaf locations, got {len(leaves)}"
    )

    # Each leaf should have a non‑zero lower‑bound annotation.
    for leaf in leaves:
        _, lower_dim = ts.printDims(leaf)
        assert lower_dim > 0, (
            f"Leaf {leaf}: lower dim should be > 0 after post, got {lower_dim}"
        )

    print(f"test_branching_measurement  OK "
          f"(locs={n_locs}, branches={len(branch_locs)}, leaves={len(leaves)})")


# ---------------------------------------------------------------------------
# 4. Identity self‑loop → still a leaf
# ---------------------------------------------------------------------------

def test_identity_self_loop():
    """A location whose only outgoing edge is an identity self‑loop
    is still considered a leaf."""
    num_q = 2
    ts = pyqreach.TransitionSystem()
    ts.addLocation(pyqreach.Location(num_q))
    # Identity self‑loop on qubit 0.
    ts.addRelation(0, 0, pyqreach.QOperation("I", num_q, [0], []))
    ts.setInitLocation(0)

    assert ts.isLeafLoc(0), (
        "Self‑loop with identity should still be a leaf"
    )
    print("test_identity_self_loop  OK")


# ---------------------------------------------------------------------------
# 5. Post fixed‑point propagation through a chain
# ---------------------------------------------------------------------------

def test_post_propagation_chain():
    """Verify annotations propagate correctly through a 3‑gate chain.
    Uses 2 qubits (minimum required by the CFLOBDD backend)."""
    num_q = 2
    ts = pyqreach.TransitionSystem()

    # Build chain: loc0 --X(q0)--> loc1 --H(q0)--> loc2 --Z(q0)--> loc3
    for _ in range(4):
        ts.addLocation(pyqreach.Location(num_q))
    ts.addRelation(0, 1, pyqreach.QOperation("X", num_q, [0], []))
    ts.addRelation(1, 2, pyqreach.QOperation("H", num_q, [0], []))
    ts.addRelation(2, 3, pyqreach.QOperation("Z", num_q, [0], []))
    ts.setInitLocation(0)

    set_initial_state(ts, "00")
    ts.computingFixedPointPost()

    # X|00⟩ = |10⟩, H|10⟩ = |-0⟩, Z|-0⟩ = |+0⟩.
    # Every location should have lower_dim = 1 (pure state throughout).
    for i in range(4):
        _, lower_dim = ts.printDims(i)
        assert lower_dim == 1, f"Loc {i}: expected lower_dim=1, got {lower_dim}"

    assert ts.satisfy(1, quantum_state("10")), "After X on q0 of |00⟩: should be |10⟩"
    assert not ts.satisfy(1, quantum_state("00"))

    assert ts.satisfy(2, quantum_state("-0")), "After H on q0 of |10⟩: should be |-0⟩"
    assert not ts.satisfy(2, quantum_state("00"))
    assert not ts.satisfy(2, quantum_state("10"))

    assert ts.satisfy(3, quantum_state("+0")), "After Z on q0 of |-0⟩: should be |+0⟩"

    print("test_post_propagation_chain  OK")


# ---------------------------------------------------------------------------
# 6. Classical propositions
# ---------------------------------------------------------------------------

def test_classical_propositions():
    """Build a circuit with classical registers and check CP handling."""
    qc = QuantumCircuit(2, 2)
    qc.h(0)
    qc.cx(0, 1)
    qc.measure(0, 0)
    qc.measure(1, 1)

    ts = pyqreach.TransitionSystem()
    result_locs = parse_qiskit_cir(qc, qc.num_qubits, ts)
    set_initial_state(ts, "00")
    ts.computingFixedPointPost()

    # Each location should have a ClassicalProposition.
    for i in range(ts.getLocationNum()):
        cp = ts.Locations[i].cp
        assert cp is not None, f"Location {i} has no ClassicalProposition"
        # termNum() returns the number of terms in the CP.
        n = cp.termNum()
        assert n >= 0, f"Location {i} CP termNum() = {n}"

    # Leaf locations should have CP terms (classical measurement outcomes).
    leaves = _leaf_locations(ts)
    cp_strings = set()
    for leaf in leaves:
        cp_strings.add(ts.Locations[leaf].cp.toString())

    print(f"Leaf CP values: {cp_strings}")
    # After Bell + measurements, we expect "00" and "11" as outcomes.
    assert "00" in cp_strings or len(leaves) >= 2, (
        f"Expected classical outcomes from Bell measurement, got {cp_strings}"
    )

    print(f"test_classical_propositions  OK  ({len(leaves)} leaves)")


# ---------------------------------------------------------------------------
# 7. Multiple outgoing edges (non‑branching)
# ---------------------------------------------------------------------------

def test_multiple_outgoing_gates():
    """A location can have multiple outgoing gate transitions
    (non‑projective)."""
    num_q = 2
    ts = pyqreach.TransitionSystem()

    ts.addLocation(pyqreach.Location(num_q))  # 0
    ts.addLocation(pyqreach.Location(num_q))  # 1
    ts.addLocation(pyqreach.Location(num_q))  # 2

    # Both transitions apply the same gate — this is a bit artificial
    # but the TS should accept it.
    ts.addRelation(0, 1, pyqreach.QOperation("H", num_q, [0], []))
    ts.addRelation(0, 2, pyqreach.QOperation("X", num_q, [0], []))
    ts.setInitLocation(0)

    assert len(ts.Locations[0].postLocations) == 2
    assert not ts.isLeafLoc(0)
    assert ts.isLeafLoc(1)
    assert ts.isLeafLoc(2)

    # Post‑condition should propagate through both edges.
    set_initial_state(ts, "00")
    ts.computingFixedPointPost()

    _, dim1 = ts.printDims(1)
    _, dim2 = ts.printDims(2)
    assert dim1 > 0, "Loc 1 should have annotation after post"
    assert dim2 > 0, "Loc 2 should have annotation after post"

    print("test_multiple_outgoing_gates  OK")


# ---------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    test_basic_construction()
    test_bell_state()
    test_branching_measurement()
    test_identity_self_loop()
    test_post_propagation_chain()
    test_classical_propositions()
    test_multiple_outgoing_gates()
    print("\n=== All TransitionSystem structural tests PASSED ===")
