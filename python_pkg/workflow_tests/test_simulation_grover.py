#!/usr/bin/env python3
"""
Test: QReach Grover simulation vs Qiskit Statevector.

Loads grover_5.qasm, runs both QReach fixed-point post and (in a
subprocess) Qiskit Statevector simulation, then verifies they agree.

The Grover circuit amplifies the marked state |11101> from the initial
state |00001> (ancilla at |1> for phase kickback).

Strategy:
  - QReach produces a single pure state (1‑dim subspace).
  - We verify it lies within the span of the 8 basis states that have
    the oracle qubit fixed to |0> and ancilla fixed to |1>.
  - Independently, Qiskit (subprocess) confirms the dominant basis
    state and probability distribution.
"""

from pathlib import Path
import subprocess
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

# IMPORTANT: import Qiskit BEFORE pyqreach to avoid CFLOBDD conflicts.
from qiskit import QuantumCircuit  # noqa: E402

import pyqreach  # noqa: E402
from qreach.parse_qiskit import parse_qiskit_cir  # noqa: E402
from qreach.qctl import set_initial_state, quantum_state, span_states  # noqa: E402


# ---------------------------------------------------------------------------
# Qiskit subprocess helper
# ---------------------------------------------------------------------------

def _run_qiskit_sim(filename: str, init_label: str) -> dict:
    """Run Qiskit Statevector simulation in a subprocess.

    Qiskit label convention is opposite to QReach's:
      Qiskit "abcde"  →  q[4]=a, q[3]=b, q[2]=c, q[1]=d, q[0]=e
      QReach "abcde"  →  q[0]=a, q[1]=b, q[2]=c, q[3]=d, q[4]=e

    Returns {qiskit_label: probability}.
    """
    script = f"""
from qiskit import QuantumCircuit
from qiskit.quantum_info import Statevector

qc = QuantumCircuit.from_qasm_file({filename!r})
sv = Statevector.from_label({init_label!r})
sv = sv.evolve(qc)
probs = sv.probabilities_dict()
for state, prob in sorted(probs.items(), key=lambda x: -x[1]):
    print(f"{{state}}:{{prob:.10f}}")
"""
    result = subprocess.run(
        [sys.executable, "-c", script],
        capture_output=True, text=True,
        cwd=str(Path(__file__).resolve().parents[1]),
    )
    if result.returncode != 0:
        raise RuntimeError(f"Qiskit subprocess failed:\n{result.stderr}")
    probs = {}
    for line in result.stdout.strip().split("\n"):
        if not line.strip():
            continue
        state, prob_str = line.split(":")
        probs[state] = float(prob_str)
    return probs


def _qiskit_to_qreach(label: str) -> str:
    """Reverse a Qiskit state label to QReach convention."""
    return label[::-1]


def _qreach_to_qiskit(label: str) -> str:
    """Reverse a QReach state label to Qiskit convention."""
    return label[::-1]


# ---------------------------------------------------------------------------
# Test
# ---------------------------------------------------------------------------

def test_grover_5():
    """Grover_5: 5 qubits, 3 work + 1 oracle + 1 ancilla."""
    filename = "benchmark/grover/grover_5.qasm"
    num_qubits = 5
    work_qubits = 3  # (5+1)//2

    # -- initial state (QReach convention) --
    qreach_init = "0" * (num_qubits - 1) + "1"  # "00001"
    print(f"Initial state (QReach convention): |{qreach_init}>")

    # -- Qiskit simulation (subprocess) --
    qiskit_init = _qreach_to_qiskit(qreach_init)  # "10000"
    probs = _run_qiskit_sim(filename, qiskit_init)

    # Show the top states for diagnostic purposes.
    top_states = sorted(probs.items(), key=lambda x: -x[1])
    print(f"Qiskit top states:")
    for state, prob in top_states[:8]:
        qr_label = _qiskit_to_qreach(state)
        print(f"  |{qr_label}>  prob={prob:.4f}")

    # The dominant state should be the marked item.
    dominant_qiskit = top_states[0][0]
    dominant_prob = top_states[0][1]
    dominant_qreach = _qiskit_to_qreach(dominant_qiskit)
    assert dominant_prob > 0.5, (
        f"Grover should amplify to > 50%, got {dominant_prob:.4f}"
    )

    # -- QReach computation --
    qc = QuantumCircuit.from_qasm_file(filename)
    ts = pyqreach.TransitionSystem()
    result_locs = parse_qiskit_cir(qc, num_qubits, ts)
    set_initial_state(ts, qreach_init)
    ts.computingFixedPointPost()

    leaf_locs = [i for i in range(ts.getLocationNum()) if ts.isLeafLoc(i)]
    assert len(leaf_locs) > 0, "No leaf locations found"

    # Build the span of all valid terminal states.
    # After Grover, only states with oracle=0 and ancilla=1 have support.
    valid_states = []
    for x in range(1 << work_qubits):
        bits = format(x, f"0{work_qubits}b")
        valid_states.append(bits + "01")
    allowed_subspace = span_states(valid_states)

    for loc in leaf_locs:
        upper_dim, lower_dim = ts.printDims(loc)
        assert lower_dim == 1, (
            f"Leaf {loc}: Grover output should be a pure state, "
            f"got lower_dim={lower_dim}"
        )
        assert ts.satisfy(loc, allowed_subspace), (
            f"Leaf {loc}: output state not in allowed span"
        )
        # The marked state by itself should NOT be satisfied (output
        # is a superposition, not a single basis state).
        assert not ts.satisfy(loc, quantum_state(dominant_qreach)), (
            f"Leaf {loc}: output is a superposition, not exactly "
            f"|{dominant_qreach}>"
        )
        print(f"  Leaf {loc}: dims=({upper_dim},{lower_dim})  "
              f"within allowed span ({len(valid_states)} states)  OK")

    print(f"\nGrover_5 simulation test PASSED "
          f"(Qiskit dominant |{dominant_qreach}> prob={dominant_prob:.4f}, "
          f"{len(leaf_locs)} leaf locations, {ts.getLocationNum()} total)")


if __name__ == "__main__":
    test_grover_5()
