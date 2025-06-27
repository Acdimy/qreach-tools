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
gateLength = circ.depth() - 1

prepare_5perfect_code(circ, idx1)

for i in range(5):
    circ.cx(i, 5+i)

ts = parse_qiskit(circ)

locBeforeMeas = ts.getLocationNum() - 6

# print the relations in the transition system

measList = applyMeasureAndReset(ts, 'ZZZZZ', idx1, 10)
# measList = applyFinalMeasurement(ts, 'ZZZZZ', idx1, 10)

# ts.addLocation(pyqreach.Location(10, 0))

# locMerge = ts.getLocationNum() - 1

# id = pyqreach.QOperation("I", 10, [0], [])
# for locMeas in measList:
#     ts.addRelation(locMeas, locMerge, id)

# ts.addRelation(gateLength, locMerge, id)

op00 = pyqreach.QOperation(["0000000000"])
ts.setAnnotation([[0, op00]])
ts.computingFixedPointPost()
outputMeas = '11100'
ts.printSupp(measList[int(outputMeas, 2)])
# for measLoc in measList[:1]:
#     ts.printDims(measLoc)
#     ts.printSupp(measLoc)

# ts.printDims(locMerge)

visualize_transition_system(ts)
