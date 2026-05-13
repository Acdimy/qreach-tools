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

qc = QuantumCircuit(16, 16)
qc.h(0)
for i in range(12):
    qc.measure(i, i)

ts1 = pyqreach.initializeTransitionSystem()
ts2 = pyqreach.TransitionSystem(False)
resultList = parse_qiskit_cir(qc, 16, ts2)
opinit = pyqreach.QOperation(["0000000000000000"])
ts2.setAnnotation([[0, opinit]])
ts2.computingFixedPointPost()
ts2.printSupp(20)


nx2Graph_hierarchical(dict2NX(ts2Dict(ts2)), "test_graph")
