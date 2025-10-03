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

circ = QuantumCircuit(3, 3)
circ.h(0)
circ.x(1)
circ.z(2)
circ.reset([0, 1, 2])
ts = pyqreach.TransitionSystem()
resultList = parse_qiskit_cir(circ, ts)
print("Transition System Locations:", ts.getLocationNum())
print("Result List size:", len(resultList))

op00 = pyqreach.QOperation(["000"])
ts.setAnnotation([[0, op00]])
start_time = time()
ts.computingFixedPointPost()
end_time = time()
print(f"Time taken for computing fixed point post: {end_time - start_time:.2f} seconds")
ts.printSupp(ts.getLocationNum()-1)

visualize_transition_system(ts, 'unit_test_reset')