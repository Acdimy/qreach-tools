#!/usr/bin/env python3
"""Test dqc_qft_2 under Pauli injection — does it crash?"""
import sys, copy, random
from pathlib import Path

PYTHON_PKG = Path(__file__).resolve().parent
sys.path.insert(0, str(PYTHON_PKG))

import pyqreach
from qiskit import QuantumCircuit
from qiskit.circuit import CircuitInstruction
from qiskit.circuit.library import XGate, YGate, ZGate
from parse_qiskit import parse_qiskit_cir_lazy
from qctl import span_qops

def insert_random_pauli(qc, rng):
    pauli = [XGate(), YGate(), ZGate()]
    gate = rng.choice(pauli)
    pos = rng.randint(0, len(qc.data))
    qb = rng.randint(0, qc.num_qubits - 1)
    qc.data.insert(pos, CircuitInstruction(gate, [qc.qubits[qb]], []))
    return qc, {"position": pos, "qubit": qb, "gate": gate.name}

def run_one(qasm_path, initial_state, label, seed=1):
    pyqreach.initializeTransitionSystem()
    qc = QuantumCircuit.from_qasm_file(str(qasm_path))
    print(f"\n=== {label} === nq={qc.num_qubits} ngates={len(qc.data)} init={initial_state}")

    rng = random.Random(seed)
    original_qc = copy.deepcopy(qc)
    qc_injected, info = insert_random_pauli(qc, rng)
    print(f"  Injected {info['gate']} at pos={info['position']} qb={info['qubit']}")

    # Injected TS
    ts = pyqreach.TransitionSystem()
    pres = parse_qiskit_cir_lazy(qc_injected, qc_injected.num_qubits, ts,
                                 initial_state=initial_state, return_metadata=True)
    ts.computingFixedPointPost()
    print(f"  Injected TS: {ts.getLocationNum()} locs, {len(pres.result_locations)} result_locs")

    # Reference TS
    ref_ts = pyqreach.TransitionSystem()
    ref_pres = parse_qiskit_cir_lazy(original_qc, original_qc.num_qubits, ref_ts,
                                     initial_state=initial_state, return_metadata=True)
    ref_ts.computingFixedPointPost()
    print(f"  Ref TS: {ref_ts.getLocationNum()} locs, {len(ref_pres.result_locations)} result_locs")

    for loc in ref_pres.result_locations:
        print(f"  ref loc {loc}: dims={ref_ts.printDims(loc)}")
        ref_ts.Locations[loc].lowerBound.printFormal()

    print(f"  Building expected via span_qops...")
    try:
        expected = span_qops([ref_ts.Locations[loc].lowerBound
                              for loc in ref_pres.result_locations])
        print(f"  expected: qNum={expected.qNum}, dim={expected.dim}")
        expected.printFormal()
    except Exception as e:
        print(f"  span_qops CRASHED: {e}")
        import traceback
        traceback.print_exc()
        return

    all_leaf = list(pres.result_locations) + list(getattr(pres, "lazy_pruned_locations", []) or [])
    print(f"  all_leaf = {all_leaf}")

    for loc in all_leaf:
        d = ts.printDims(loc)
        print(f"  leaf loc {loc}: dims={d}")

    print(f"  Calling satisfy...")
    try:
        result = ts.Locations[all_leaf[0]].satisfy(expected)
        print(f"  satisfy OK: {result}")
    except Exception as e:
        print(f"  satisfy CRASHED: {e}")

if __name__ == "__main__":
    # dqc_qft_2 with seed=1 (init="01")
    run_one(PYTHON_PKG / "benchmark" / "dqc_qft" / "dqc_qft_2.qasm",
            "01", "dqc_qft_2 init=01 (seed=1)", seed=1)

    # dqc_qft_2 with seed=2 (init="11")
    run_one(PYTHON_PKG / "benchmark" / "dqc_qft" / "dqc_qft_2.qasm",
            "11", "dqc_qft_2 init=11 (seed=2)", seed=2)

    # dqc_qft_2 with seed=3 (init="00")
    run_one(PYTHON_PKG / "benchmark" / "dqc_qft" / "dqc_qft_2.qasm",
            "00", "dqc_qft_2 init=00 (seed=3)", seed=3)
