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

# Unit Test Measurements
# qubits = QuantumRegister(5, 'q')
# c = ClassicalRegister(2, 'c')
# qc = QuantumCircuit(qubits, c)
# qc.h(0)
# qc.h(1)
# qc.cx(0, 2)
# qc.cx(1, 3)
# qc.measure(0, 0)
# qc.measure(1, 1)

# ts = pyqreach.TransitionSystem()
# resultList = parse_qiskit_cir(qc, ts)

# Unit Test if-else condition
# q = QuantumRegister(3, "q")
# c = ClassicalRegister(3, "c")
# circ = QuantumCircuit(q, c)
# circ.h(q[0])
# circ.cx(q[0], q[1])
# circ.measure(q[0], c[0])
 
# # Qiskit 2.0 uses if_test instead of c_if
# with circ.if_test((c, 4)):
#     circ.x(q[1])

# ts = pyqreach.TransitionSystem()
# resultList = parse_qiskit_cir(circ, ts)

# Unit Test while loop
# qubits = QuantumRegister(2)
# clbits = ClassicalRegister(2)
# circuit = QuantumCircuit(qubits, clbits)

# q0, q1 = qubits
# c0, c1 = clbits

# circuit.h(0)
# circuit.h(1)
# circuit.measure(0, 0)
# circuit.measure(1, 1)
# with circuit.while_loop((clbits, 0b11)):
#     circuit.h(0)
#     circuit.h(1)
#     circuit.measure(0, 0)
#     circuit.measure(1, 1)

# ts = pyqreach.TransitionSystem()
# resultList = parse_qiskit_cir(circuit, ts)

# Integration Test: Quantum Bernoulli Factory
qubits = QuantumRegister(3)
midBits = ClassicalRegister(2, name='mid')
ancBits = ClassicalRegister(1, name='anc')
circ = QuantumCircuit(qubits, midBits, ancBits)
theta = 2 * np.arccos(np.sqrt(0.2))
circ.reset([0, 1, 2])
for i in range(2):
    with circ.if_test((midBits, 0b00)):
        # circ.ry(theta, 0); circ.ry(theta, 1)
        circ.u(theta, 0, 0, 0)
        circ.u(theta, 0, 0, 1)
        circ.cx(0, 1)
        circ.h(0)
        circ.measure(0,0)
        circ.measure(1,1)
    with circ.if_test((midBits, 0b00)):
        circ.x(2)
        circ.measure(2,2)
    with circ.if_test((midBits, 0b11)):
        circ.x(2)
        circ.measure(2,2)
    with circ.while_loop((ancBits, 0b1)):
        circ.reset([0, 1, 2])
        
        # Necessary to re-measure qubit 2 to reset the ancBits
        circ.measure(2,2)
        # circ.ry(theta, 0); circ.ry(theta, 1)
        circ.u(theta, 0, 0, 0)
        circ.u(theta, 0, 0, 1)

        circ.cx(0, 1)
        circ.h(0)
        circ.measure(0,0)
        circ.measure(1,1)
        with circ.if_test((midBits, 0b00)):
            circ.x(2)
        with circ.if_test((midBits, 0b11)):
            circ.x(2)
        circ.measure(2,2)
    with circ.if_test((midBits, 0b01)):
        circ.reset([0, 1, 2])
        circ.measure(0,0)
        circ.measure(1,1)
ts = pyqreach.TransitionSystem()
resultList = parse_qiskit_cir(circ, circ.num_qubits, ts)
print("Transition System Locations:", ts.getLocationNum())
print("Result List size:", len(resultList))

op00 = pyqreach.QOperation(["000"])
ts.setAnnotation([[0, op00]])
start_time = time()
ts.computingFixedPointPost()
end_time = time()
print(f"Time taken for computing fixed point post: {end_time - start_time:.2f} seconds")

# visualize_transition_system(ts, 'unit_test')
