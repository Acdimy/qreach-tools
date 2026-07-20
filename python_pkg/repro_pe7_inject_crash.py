#!/usr/bin/env python3
"""Minimal reproduction of pe_7 injected crash (SIGBUS/SIGSEGV, code -10/-11).

Usage:
    cd python_pkg
    ../.venv/bin/python repro_pe7_inject_crash.py

The crash occurs in the verification step (_run_debug_check) which does a
SECOND parse of the original circuit via _simulate_final_operation.
CFLOBDD global state interaction between the two parses is the suspected trigger.
"""

import copy
import random
import sys
from pathlib import Path
from time import perf_counter

import pyqreach
from qiskit import QuantumCircuit
from qiskit.circuit import CircuitInstruction
from qiskit.circuit.library import XGate, YGate, ZGate

from parse_qiskit import parse_qiskit_cir_lazy
from qctl import modelChecking, quantum_state, set_initial_state, span_qops, tsLabelling

PYTHON_PKG = Path(__file__).resolve().parent
_NUSMV_PATH = str(
    (PYTHON_PKG.parent.parent / "NuSMV-2.7.0-macos-universal" / "bin" / "NuSMV").resolve()
)

# ---------------------------------------------------------------------------
# Config — matches experiment run (seed=1, pe_7 at file_index=5)
# ---------------------------------------------------------------------------
QASM_FILE = PYTHON_PKG / "benchmark" / "pe" / "pe_7.qasm"
SEED = 1
FILE_INDEX = 5  # pe_2=0, pe_3=1, pe_4=2, pe_5=3, pe_6=4, pe_7=5


def insert_random_pauli(qc, rng):
    new_qc = copy.deepcopy(qc)
    num_qubits = new_qc.num_qubits
    num_ops = len(new_qc.data)
    pos = rng.randint(0, num_ops)
    q = rng.randint(0, num_qubits - 1)
    gate_name = rng.choice(["x", "y", "z"])
    gate_map = {"x": XGate(), "y": YGate(), "z": ZGate()}
    instr = CircuitInstruction(gate_map[gate_name], [new_qc.qubits[q]], [])
    new_qc.data.insert(pos, instr)
    return new_qc, (pos, q, gate_name)


def simulate_final(qc, initial_state):
    """Replicate _simulate_final_operation — second CFLOBDD parse."""
    ts = pyqreach.TransitionSystem()
    result = parse_qiskit_cir_lazy(
        qc, qc.num_qubits, ts,
        initial_state=initial_state,
        return_metadata=True,
    )
    ts.computingFixedPointPost()
    if len(result.result_locations) == 1:
        return ts.Locations[result.result_locations[0]].lowerBound
    return span_qops([ts.Locations[loc].lowerBound for loc in result.result_locations])


def main():
    pyqreach.initializeTransitionSystem()

    # Inference of initial state
    n = 7
    rng_seed = SEED * 10000 + n + 3_000_000
    initial_state = "".join(random.Random(rng_seed).choice(["0", "1"]) for _ in range(8))
    print(f"[1] Initial state: {initial_state}")

    # Load circuit
    qc_orig = QuantumCircuit.from_qasm_file(str(QASM_FILE))
    print(f"[2] Loaded {QASM_FILE.name}:  {qc_orig.num_qubits} qubits, {len(qc_orig.data)} gates")

    # Insert Pauli error
    rng = random.Random(SEED + FILE_INDEX)
    qc_injected, (pos, q, gate) = insert_random_pauli(qc_orig, rng)
    print(f"[3] Injected {gate.upper()} at pos {pos}/{len(qc_injected.data)} on qubit {q}")

    # ---- FIRST parse (injected circuit) ----
    print("[4] FIRST parse (injected circuit) ...", flush=True)
    t0 = perf_counter()
    ts1 = pyqreach.TransitionSystem()
    r1 = parse_qiskit_cir_lazy(qc_injected, qc_injected.num_qubits, ts1,
                               initial_state=initial_state, return_metadata=True)
    print(f"    {ts1.getLocationNum()} locs, {len(r1.result_locations)} result  "
          f"({perf_counter()-t0:.3f}s)")

    print("[5] FIRST fixed-point post ...", flush=True)
    t0 = perf_counter()
    ts1.computingFixedPointPost()
    print(f"    Done ({perf_counter()-t0:.3f}s)")

    # ---- SECOND parse (original circuit for verification) ----
    # This is _simulate_final_operation — the suspected crash site
    print("[6] SECOND parse (original circuit for verification) ...", flush=True)
    t0 = perf_counter()
    expected = simulate_final(qc_orig, initial_state)
    print(f"    Done ({perf_counter()-t0:.3f}s)")

    # ---- Model checking ----
    print("[7] Model checking ...", flush=True)
    all_leaf = list(r1.result_locations)
    for loc in all_leaf:
        ts1.setLabel(loc, "final")
    tsLabelling(ts1, expected, "debug_op", locList=all_leaf)
    t0 = perf_counter()
    result = modelChecking(ts1, "AG (debug_op <-> final)", nusmv_path=_NUSMV_PATH)
    print(f"    satisfied={result.get('satisfied')}  ({perf_counter()-t0:.3f}s)")

    print("\nSUCCESS: no crash.")


if __name__ == "__main__":
    main()
