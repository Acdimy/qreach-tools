import pyqreach
### Verifiable quantum secret sharing
from qiskit import QuantumCircuit, QuantumRegister, ClassicalRegister, transpile
from qiskit.quantum_info import Statevector
# from qiskit_aer import AerSimulator
import numpy as np

from parse_qiskit import *

def prepare_5perfect_code(circ, idx):
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

circ = QuantumCircuit(10, 10)
idx0 = list(range(5))
idx1 = list(range(5,10))


circ.h(0)
circ.t(0)
circ.x(0)

circ.h(5)
prepare_5perfect_code(circ, idx0)
prepare_5perfect_code(circ, idx1)

for i in range(5):
    circ.cx(i, 5+i)

ts = parse_qiskit(circ)

# print the relations in the transition system

measList = applyFinalMeasurement(ts, 'ZZZZZ', idx1, 10)

print(ts.getLocationNum())

op00 = pyqreach.QOperation(["0000000000"])
ts.setAnnotation([[0, op00]])
ts.computingFixedPointPost()

measOutput = '11001'


for l in measList:
    ts.printDims(l)
    # ts.printSupp(l)
ts.printSupp(measList[int(measOutput, 2)])
