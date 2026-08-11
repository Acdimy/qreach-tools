#!/usr/bin/env python3
"""Quick diagnosis of grover failures on small circuits."""
import sys, os, time
from pathlib import Path

PYTHON_PKG = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PYTHON_PKG))

from parse_qiskit import parse_qiskit_cir_lazy
from qiskit import QuantumCircuit
import pyqreach
from qctl import quantum_state, span_qops

CONV = PYTHON_PKG / "benchmark" / "converted_qasm"
TESTS = [
    "single-it-grover8-plus.qasm",
    "single-it-grover8-zero.qasm",
    "single-it-grover16-plus.qasm",
    "single-it-grover32-plus-linear.qasm",
]

for cfn in TESTS:
    import re
    path = CONV / cfn
    qc = QuantumCircuit.from_qasm_file(str(path))
    n = qc.num_qubits
    m = re.search(r"grover(\d+)", cfn)
    work = int(m.group(1))

    pyqreach.initializeTransitionSystem()
    ts = pyqreach.TransitionSystem()
    init = "0" * n if "zero" in cfn else "0" * (n - 1) + "1"
    r = parse_qiskit_cir_lazy(qc, n, ts, initial_state=init, return_metadata=True)

    fl = r.result_locations[-1]
    if isinstance(fl, list):
        fl = fl[0]

    lb = ts.Locations[fl].lowerBound
    ub = ts.Locations[fl].upperBound

    print(f"--- {cfn} n={n} work={work} locs={len(ts.Locations)} ---")
    print(f"  lowerBound dim={lb.dim()}")
    print(f"  upperBound dim={ub.dim()}")
    print(f"  satisfy(self)={ts.Locations[fl].satisfy(lb)}")

    # Check the debug condition
    gi = ts.Locations[work].lowerBound
    gg = quantum_state("1" * work + "0" * (n - work - 1) + "1")
    print(f"  init dim={gi.dim()}, good dim={gg.dim()}")
    gf = span_qops([gi, gg])
    print(f"  final_span dim={gf.dim()}, satisfy(final_span)={ts.Locations[fl].satisfy(gf)}")

    # Check some intermediate locations
    import random
    for i in random.sample(range(min(10, len(ts.Locations))), min(3, len(ts.Locations))):
        lb_i = ts.Locations[i].lowerBound
        print(f"  loc[{i}]: lb.dim={lb_i.dim()} self_ok={ts.Locations[i].satisfy(lb_i)}")
    print()
