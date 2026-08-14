from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import pyqreach
### Verifiable quantum secret sharing
from qiskit import QuantumCircuit, QuantumRegister, ClassicalRegister, transpile
from qiskit.quantum_info import Statevector
from qiskit_aer import Aer
import numpy as np
from math import pi
import random
from time import time
from qreach.parse_qiskit import *
from qreach.qctl import *
from qreach.circ_utils import *

# Correct
# qc = QuantumCircuit(3, 1)
# qc.x(2)
# qc.measure(2, 0)
# qc.reset(1)
# with qc.while_loop((0, 0b1)):
#     qc.reset(0)
#     qc.h(0)
#     qc.t(0)
#     qc.cx(0, 1)
#     qc.h(0)
#     qc.cx(0, 1)
#     qc.t(0)
#     qc.h(0)
#     qc.measure(0, 0)

# Buggy
qc = QuantumCircuit(3, 1)
qc.x(2)
qc.measure(2, 0)
qc.reset(1)
with qc.while_loop((0, 0b1)):
    # qc.reset(0)
    qc.h(0)
    qc.t(0)
    qc.cx(0, 1)
    qc.h(0)
    qc.cx(0, 1)
    qc.t(0)
    qc.h(0)
    qc.measure(0, 0)

ts = pyqreach.TransitionSystem()
parse_result = parse_qiskit_cir(qc, qc.num_qubits, ts, return_metadata=True)
set_initial_state(ts, "000")
ts.computingFixedPointPost()

# Specification 1: EF (Zero), where Zero is a zero-dimensional subspace.
# EF (Zero) is true iff there exists a reachable location whose quantum state
# became the zero subspace after fixed-point computation (i.e. an unreachable
# location that is graph-reachable from the initial state).
zero_op = zero_subspace(qc.num_qubits)
tsLabelling(ts, zero_op, "Zero")
result1 = modelChecking(ts, "EF (Zero)")
print("Spec 1: EF (Zero)")
print("  satisfied:", result1['satisfied'])
if result1.get('analysis'):
    from qreach.qctl import _format_counterexample_analysis
    print(_format_counterexample_analysis(result1['analysis']))

# Specification 2: AF (outloop -> s), where s is the subspace of
# (|001> + i\sqrt(2)|011>), and outloop is the first location exiting the loop.

# Find the outloop location (EW = exit-while, with non-zero quantum state).
outloop_locs = [
    loc for loc in range(ts.getLocationNum())
    if "EW" in (ts.Locations[loc].getIdentifier() or "")
    and ts.printDims(loc)[1] > 0
]
if not outloop_locs:
    raise RuntimeError("No outloop (EW) location with non-zero quantum state found")
# Label only the non-zero outloop location(s)
for loc in outloop_locs:
    ts.setLabel(loc, "outloop")
print(f"outloop locations: {outloop_locs}")

# Target subspace: |001> + i*sqrt(2)|011>  (unnormalized)
target_s = amplitude_state({"001": 1.0 + 0j, "011": 1j * IRR_SQRT2})
tsLabelling(ts, target_s, "s")

# AG (outloop -> s): at the outloop location, the quantum state must satisfy s.
result2 = modelChecking(ts, "AG (outloop -> s)")
print("Spec 2: AF (outloop -> s)")
print("  satisfied:", result2['satisfied'])
if result2.get('analysis'):
    from qreach.qctl import _format_counterexample_analysis
    print(_format_counterexample_analysis(result2['analysis']))