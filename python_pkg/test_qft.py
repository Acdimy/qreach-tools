import pyqreach
### Verifiable quantum secret sharing
from qiskit import QuantumCircuit, QuantumRegister, ClassicalRegister, transpile
from qiskit.quantum_info import Statevector
from qiskit_aer import Aer
import numpy as np
import math
import random
from time import time
from parse_qiskit import *

# Unit Test Measurements
qubits = QuantumRegister(5, 'q')
c = ClassicalRegister(2, 'c')
qc = QuantumCircuit(qubits, c)
qc.h(0)
qc.h(1)
qc.cx(0, 2)
qc.cx(1, 3)
qc.measure(0, 0)
qc.measure(1, 1)

ts = pyqreach.TransitionSystem()
resultList = parse_qiskit_cir(qc, ts)

visualize_transition_system(ts, 'unit_test')
