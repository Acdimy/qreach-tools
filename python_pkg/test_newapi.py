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
resultList = parse_qiskit_cir(qc, qc.num_qubits, ts)
opinit = pyqreach.QOperation(["000"])
ts.setAnnotation([[0, opinit]])
ts.computingFixedPointPost()
annotated_loc = annotate(ts, ["leaf", "reached"])
print(annotated_loc)
