import pyqreach
### Verifiable quantum secret sharing
from qiskit import QuantumCircuit, QuantumRegister, ClassicalRegister, transpile
from qiskit.quantum_info import Statevector
# from qiskit_aer import AerSimulator
import numpy as np

from parse_qiskit import *

def prepare_5perfect(circ, idx):
    circ.z(idx[0])

    circ.h(idx[1])
    circ.cz(idx[1],idx[2])
    circ.cz(idx[1],idx[4])
    circ.cx(idx[1],idx[0])

    circ.h(idx[4])
    circ.cz(idx[4],idx[3])
    circ.cz(idx[4],idx[1])
    circ.cx(idx[4],idx[0])

    circ.h(idx[3])
    circ.cz(idx[3],idx[2])
    circ.cz(idx[3],idx[0])
    circ.cx(idx[3],idx[4])

    circ.h(idx[2])
    circ.cz(idx[2],idx[1])
    circ.cz(idx[2],idx[4])
    circ.cx(idx[2],idx[3])

circ = QuantumCircuit(6, 6)
idx = [1,2,3,4,5]
for i in idx:
    circ.h(i)
prepare_5perfect(circ, idx)
circ.x(0)

# print([circ.data[i][0].name for i in range(len(circ.data))])

ts = parse_qiskit(circ)

# print the relations in the transition system


op00 = pyqreach.QOperation(["000000"])
ts.setAnnotation([[0, op00]])
ts.computingFixedPointPost()
ts.printDims(0)
