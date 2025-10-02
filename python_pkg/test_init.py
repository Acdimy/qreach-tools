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

# Unit Test Qiskit initialize
q = QuantumRegister(3, "q")
c = ClassicalRegister(3, "c")
circ = QuantumCircuit(q, c)
# circ.h(0)
# circ.cx(0, 1)
circ.initialize([1/np.sqrt(2), 1/np.sqrt(2)], [0])
circ.h(0)

ts = pyqreach.TransitionSystem()
resultList = parse_qiskit_cir(circ, circ.num_qubits, ts)
# visualize_transition_system(ts, "unit_test_initialize")
op00 = pyqreach.QOperation(["000"])
ts.setAnnotation([[0, op00]])
start_time = time()
ts.computingFixedPointPost()
end_time = time()
print("Fixed point computation time: ", end_time - start_time)
ts.printSupp(1)
ts.printSupp(2)
ts.printSupp(3)
ts.printSupp(4)
ts.printSupp(5)
