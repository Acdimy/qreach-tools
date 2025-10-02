import pyqreach
### Verifiable quantum secret sharing
from qiskit import QuantumCircuit, QuantumRegister, ClassicalRegister, transpile
from qiskit.quantum_info import Statevector
import numpy as np
from time import time
from qiskit_aer import AerSimulator

from parse_qiskit import *
from qctl import *
from circ_utils import *

qc = QuantumCircuit(3, 2)
qc.h(2)
qc.measure(2, 0)
qc.h(2)
qc.measure(2, 1)

qc.h(0)
qc.cx(0, 1)
with qc.if_test((0, 1)): # b1 == 1
    qc.z(0)
with qc.if_test((1, 1)): # b2 == 1
    qc.x(0)
qc.cx(0, 1)
qc.h(0)

ts1 = pyqreach.TransitionSystem()
ts2 = pyqreach.TransitionSystem(False)
resultList = parse_qiskit_cir(qc, 3, ts2)
