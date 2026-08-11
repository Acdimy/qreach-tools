#!/usr/bin/env python3
"""Run benchpress grover benchmarks."""
import sys, time, os, re
from pathlib import Path

PYTHON_PKG = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PYTHON_PKG))

from parse_qiskit import parse_qiskit_cir_lazy
from qiskit import QuantumCircuit
import pyqreach
from qctl import quantum_state, span_qops

CONVERTED = PYTHON_PKG / "benchmark" / "converted_qasm"
circuits = sorted([f for f in os.listdir(CONVERTED) if "grover" in f and f.endswith(".qasm")])

header = f"{'Circuit':<48} {'Qubits':>6} {'Status':>10} {'Time':>8}"
print(header)
print("-" * 76)

for cfn in circuits:
    path = CONVERTED / cfn
    qc = QuantumCircuit.from_qasm_file(str(path))
    n = qc.num_qubits
    m = re.search(r"grover(\d+)", cfn)
    work = int(m.group(1)) if m else n

    t0 = time.perf_counter()
    try:
        pyqreach.initializeTransitionSystem()
        ts = pyqreach.TransitionSystem()
        init = "0" * (n - 1) + "1"
        r = parse_qiskit_cir_lazy(qc, n, ts, initial_state=init, return_metadata=False)
        elapsed = time.perf_counter() - t0

        final_loc = r[-1] if isinstance(r, list) else r
        if isinstance(final_loc, list):
            final_loc = final_loc[0]
        grover_init = ts.Locations[work].lowerBound
        grover_good = quantum_state("1" * work + "0" * (n - work - 1) + "1")
        grover_final = span_qops([grover_init, grover_good])
        ok = ts.Locations[final_loc].satisfy(grover_final)
        status = "PASS" if ok else "WRONG"
    except Exception as e:
        elapsed = time.perf_counter() - t0
        status = f"ERROR({str(e)[:20]})"

    print(f"{cfn:<48} {n:>6} {status:>10} {elapsed:>7.1f}s")
