import pyqreach
### Verifiable quantum secret sharing
from qiskit import QuantumCircuit, QuantumRegister, ClassicalRegister, transpile
from qiskit.quantum_info import Statevector
# from qiskit_aer import AerSimulator
import numpy as np

from parse_qiskit import *

circ = QuantumCircuit(10, 10)

# circ.x(0)

prepare_5perfect_code(circ, [0,1,2,3,4])
gateLength = len(circ.data)

circ.h(5)
prepare_5perfect_code(circ, [5,6,7,8,9])

for i in range(5):
    circ.cx(i, 5+i)

ts = parse_qiskit(circ)

# print the relations in the transition system

measList = applyMeasureAndReset(ts, 'ZZZZZ', [5,6,7,8,9], 10)
# measList = applyFinalMeasurement(ts, 'ZZZZZ', [5,6,7,8,9], 10)

ts.addLocation(pyqreach.Location(10, 0))
locMerge = ts.getLocationNum() - 1
id = pyqreach.QOperation("I", 10, [0], [])
for locMeas in measList[7:8]:
    ts.addRelation(locMeas, locMerge, id)

ts.addRelation(gateLength, locMerge, id)

op00 = pyqreach.QOperation(["0000000000"])
ts.setAnnotation([[0, op00]])
ts.computingFixedPointPost()
ts.printDims(locMerge)
ts.printSupp(locMerge)
ts.printSupp(10)
# outputMeas = '11100'
# ts.printSupp(measList[int(outputMeas, 2)])
# for measLoc in measList[:1]:
#     ts.printDims(measLoc)
#     ts.printSupp(measLoc)

# ts.printDims(locMerge)

# visualize_transition_system(ts)
