import pyqreach
### Verifiable quantum secret sharing
from qiskit import QuantumCircuit, QuantumRegister, ClassicalRegister, transpile
from qiskit.quantum_info import Statevector
import numpy as np

from parse_qiskit import *

circ = QuantumCircuit(14, 14)

prepare_steane_code(circ, list(range(7)))
gateLength = len(circ.data)
circ.h(7)
prepare_steane_code(circ, list(range(7, 14)))
for i in range(7):
    circ.cx(i, 7+i)
ts = parse_qiskit(circ)
measList = applyMeasureAndReset(ts, 'ZZZZZZZ', list(range(7, 14)), 14)


op00 = pyqreach.QOperation(["00000000000000"])
ts.setAnnotation([[0, op00]])
ts.computingFixedPointPost()
print(measList)
# for measLoc in measList:
#     dims = ts.printDims(measLoc)
#     if dims[1] > 0:
#         print(f"Location {measLoc}: upperBound dimension = {dims[0]}, lowerBound dimension = {dims[1]}")
#         ts.printSupp(measLoc)

visualize_transition_system(ts, 'test_steane')
