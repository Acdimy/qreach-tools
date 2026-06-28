#!/usr/bin/env python3
"""
Test: QReach QFT simulation.

Loads qft_5.qasm, runs both QReach fixed-point post and (in a subprocess)
Qiskit Statevector simulation, then verifies they agree.

QFT|0…0⟩ = H^⊗n |0…0⟩ = |+…+⟩ (equal superposition over all basis states).
Since the circuit is unitary with no measurements, the output is a single
pure state, which we can check directly with satisfy().
"""

from pathlib import Path
import subprocess
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

# IMPORTANT: import Qiskit BEFORE pyqreach to avoid CFLOBDD conflicts.
from qiskit import QuantumCircuit  # noqa: E402

import pyqreach  # noqa: E402
from parse_qiskit import parse_qiskit_cir  # noqa: E402
from qctl import set_initial_state, quantum_state  # noqa: E402


def _run_qiskit_sim(filename: str, init_label: str) -> dict:
    """Run Qiskit Statevector simulation in a subprocess.

    Returns {qiskit_label: probability}.
    """
    script = f"""
from qiskit import QuantumCircuit
from qiskit.quantum_info import Statevector

qc = QuantumCircuit.from_qasm_file({filename!r})
sv = Statevector.from_label({init_label!r})
sv = sv.evolve(qc)
probs = sv.probabilities_dict()
for state, prob in sorted(probs.items()):
    if prob > 1e-15:
        print(f"{{state}}:{{prob:.10f}}")
print(f"TOTAL_STATES:{{len(probs)}}")
"""
    result = subprocess.run(
        [sys.executable, "-c", script],
        capture_output=True, text=True,
        cwd=str(Path(__file__).resolve().parents[1]),
    )
    if result.returncode != 0:
        raise RuntimeError(f"Qiskit subprocess failed:\n{result.stderr}")
    probs = {}
    total = 0
    for line in result.stdout.strip().split("\n"):
        if not line.strip():
            continue
        if line.startswith("TOTAL_STATES:"):
            total = int(line.split(":")[1])
            continue
        state, prob_str = line.split(":")
        probs[state] = float(prob_str)
    return probs, total


def test_qft_5():
    """QFT_5: 5‑qubit QFT from |00000⟩ → equal superposition |+++++⟩."""
    filename = "benchmark/qft/qft_5.qasm"
    num_qubits = 5

    init_state = "0" * num_qubits  # |00000⟩
    print(f"Initial state: |{init_state}>")

    # -- Qiskit simulation (subprocess) --
    probs, total = _run_qiskit_sim(filename, init_state)
    expected_uniform = 1.0 / (2**num_qubits)
    eps = 1e-12

    assert total == 2**num_qubits, (
        f"Expected {2**num_qubits} non-zero states, got {total}"
    )
    for state, prob in probs.items():
        assert abs(prob - expected_uniform) < eps, (
            f"QFT|00000⟩ should be uniform, |{state}⟩ prob={prob:.6f}"
        )
    print(f"Qiskit: all {total} basis states uniform at {expected_uniform:.6f}  OK")

    # -- QReach computation --
    qc = QuantumCircuit.from_qasm_file(filename)
    ts = pyqreach.TransitionSystem()
    result_locs = parse_qiskit_cir(qc, num_qubits, ts)
    set_initial_state(ts, init_state)
    ts.computingFixedPointPost()

    leaf_locs = [i for i in range(ts.getLocationNum()) if ts.isLeafLoc(i)]
    assert len(leaf_locs) > 0, "No leaf locations found"

    # QFT|0…0⟩ = |+…+⟩ — the equal superposition pure state.
    expected_op = quantum_state("+" * num_qubits)

    for loc in leaf_locs:
        upper_dim, lower_dim = ts.printDims(loc)
        assert lower_dim == 1, (
            f"Leaf {loc}: QFT is unitary, expected lower_dim=1, got {lower_dim}"
        )
        assert ts.satisfy(loc, expected_op), (
            f"Leaf {loc}: QFT output must be |+++++>"
        )
        print(f"  Leaf {loc}: dims=({upper_dim},{lower_dim})  "
              f"satisfy(|+++++>) = True  OK")

    print(f"\nQFT_5 simulation test PASSED "
          f"(Qiskit confirms uniform {total} states, "
          f"{len(leaf_locs)} leaf locations, {ts.getLocationNum()} total)")


if __name__ == "__main__":
    test_qft_5()
