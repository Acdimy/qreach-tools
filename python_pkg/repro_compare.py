#!/usr/bin/env python3
"""Compare dqc_pe_2 vs dqc_qft_2 under identical Pauli injection.

Goal: understand why dqc_qft sometimes doesn't crash.
"""
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

    # Print dims for ref result locations
    for loc in ref_pres.result_locations:
        d = ref_ts.printDims(loc)
        print(f"  ref loc {loc}: dims={d}")

    # Build expected from reference
    print(f"  Building expected via span_qops...")
    try:
        ref_ops = [ref_ts.Locations[loc].lowerBound for loc in ref_pres.result_locations]
        print(f"  ref_ops: {len(ref_ops)} ops")
        for i, op in enumerate(ref_ops):
            print(f"    op[{i}]: qNum={op.qNum}, dim={op.dim}")
            op.printFormal()
        expected = span_qops(ref_ops)
        print(f"  expected: qNum={expected.qNum}, dim={expected.dim}")
        expected.printFormal()
    except Exception as e:
        print(f"  span_qops CRASHED: {e}")
        import traceback
        traceback.print_exc()
        return

    # CRASH POINT
    all_leaf = list(pres.result_locations) + list(getattr(pres, "lazy_pruned_locations", []) or [])
    print(f"  all_leaf = {all_leaf}")
    print(f"  Injected leaf[0] lowerBound:")
    ts.Locations[all_leaf[0]].lowerBound.printFormal()

    print(f"  Calling satisfy...")
    try:
        result = ts.Locations[all_leaf[0]].satisfy(expected)
        print(f"  satisfy OK: {result}")
    except Exception as e:
        print(f"  satisfy CRASHED: {e}")
        import traceback
        traceback.print_exc()

if __name__ == "__main__":
    # dqc_pe_2 — known to crash
    run_one(PYTHON_PKG / "benchmark" / "dqc_pe" / "dqc_pe_2.qasm",
            "001", "dqc_pe_2", seed=1)
