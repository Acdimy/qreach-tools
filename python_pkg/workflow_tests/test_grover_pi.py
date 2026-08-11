#!/usr/bin/env python3
"""Quick verification: post_image chain match and test on grover circuits."""
import sys, os
from pathlib import Path
PYTHON_PKG = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PYTHON_PKG))

from parse_qiskit import parse_qiskit_cir_lazy
from qiskit import QuantumCircuit
import pyqreach
from qctl import quantum_state, span_qops

# Step 1: verify post_image chain produces correct state
print("=== Step 1: Verify post_image ===")
for n in [3, 7, 15]:
    ndata = (n + 1) // 2
    pyqreach.initializeTransitionSystem()
    ref = quantum_state("+" * ndata + "0" * (n - ndata))
    state = pyqreach.QOperation(["0" * n])
    for i in range(ndata):
        state = state.post_image(pyqreach.QOperation("H", n, [i], []))
    sp = span_qops([ref, state])
    match = "OK" if sp.dim() == 1 else f"FAIL(dim={sp.dim()})"
    print(f"  n={n} ndata={ndata}: {match}")

# Step 2: test small grover circuits using post_image
print("\n=== Step 2: Small grover with post_image ===")
CONV = PYTHON_PKG / "benchmark" / "converted_qasm"
for cfn in sorted(os.listdir(CONV)):
    if not cfn.startswith("single-it-grover") or not cfn.endswith(".qasm"):
        continue
    path = CONV / cfn
    qc = QuantumCircuit.from_qasm_file(str(path))
    n = qc.num_qubits
    if n > 63:
        continue
    ndata = (n + 1) // 2

    pyqreach.initializeTransitionSystem()
    ts = pyqreach.TransitionSystem()
    r = parse_qiskit_cir_lazy(qc, n, ts, initial_state="0" * n, return_metadata=True)
    fl = r.result_locations[-1]
    if isinstance(fl, list):
        fl = fl[0]

    bz = quantum_state("0" * n)
    bp = pyqreach.QOperation(["0" * n])
    for i in range(ndata):
        bp = bp.post_image(pyqreach.QOperation("H", n, [i], []))
    target = span_qops([bz, bp])
    ok = ts.Locations[fl].satisfy(target)
    print(f"  {cfn:<45} n={n:>3} satisfy={ok}")
