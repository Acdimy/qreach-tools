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
from parse_qiskit import *
from qctl import *
from circ_utils import *

# Correct
# qc = QuantumCircuit(3, 1)
# qc.x(2)
# qc.measure(2, 0)
# qc.initialize(0, 1)
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
    from qctl import _format_counterexample_analysis
    print(_format_counterexample_analysis(result1['analysis']))

# Specification 2: AF (outloop -> s), where s is the subspace of
# (|001> + i\sqrt(2)|011>), and outloop is the first location exiting the loop.