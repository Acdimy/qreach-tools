#!/usr/bin/env python3
"""Test linear grover circuits with post-image verification."""
import sys, os
from pathlib import Path
PYTHON_PKG = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PYTHON_PKG))

from parse_qiskit import parse_qiskit_cir_lazy
from qiskit import QuantumCircuit
import pyqreach
from qctl import quantum_state, span_qops

CONV = PYTHON_PKG / "benchmark" / "converted_qasm"

for cfn in [
    "single-it-grover32-plus-linear.qasm",
    "single-it-grover32-zero-linear.qasm",
    "single-it-grover64-plus-linear.qasm",
    "single-it-grover64-zero-linear.qasm",
]:
    path = CONV / cfn
    if not path.exists():
        continue
    qc = QuantumCircuit.from_qasm_file(str(path))
    n = qc.num_qubits
    ndata = (n + 1) // 2
    nanc = n // 2

    pyqreach.initializeTransitionSystem()
    ts = pyqreach.TransitionSystem()
    r = parse_qiskit_cir_lazy(qc, n, ts, initial_state="0" * n, return_metadata=True)
    fl = r.result_locations[-1]
    if isinstance(fl, list):
        fl = fl[0]

    bz = quantum_state("0" * n)
    if n < 31:
        bp = quantum_state("+" * ndata + "0" * nanc)
    else:
        bp = pyqreach.QOperation(["0" * n])
        for i in range(ndata):
            bp = bp.post_image(pyqreach.QOperation("H", n, [i], []))
    target = span_qops([bz, bp])
    ok = ts.Locations[fl].satisfy(target)
    print(f"{cfn:<45} n={n:>3} satisfy={ok}")
