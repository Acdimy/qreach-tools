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
from qctl import *
from circ_utils import *

qubits = QuantumRegister(2)
clbits = ClassicalRegister(2)
circuit = QuantumCircuit(qubits, clbits)

q0, q1 = qubits
c0, c1 = clbits

circuit.h([q0, q1])
circuit.measure(q0, c0)
circuit.measure(q1, c1)
# with circuit.while_loop((clbits, 0b11)):
#     circuit.h([q0, q1])
#     circuit.measure(q0, c0)
#     circuit.measure(q1, c1)
with circuit.if_test((clbits, 0b01)):
    circuit.x(q1)

ts = pyqreach.TransitionSystem()
resultList = parse_qiskit_cir(circuit, 2, ts)
for i, loc in enumerate(ts.Locations):
    print(i, loc.getIdentifier())
visualize_transition_system(ts, "while_identifier_test")
