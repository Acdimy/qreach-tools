import pyqreach
### Verifiable quantum secret sharing
from qiskit import QuantumCircuit, QuantumRegister, ClassicalRegister, transpile
from qiskit.quantum_info import Statevector
import numpy as np
from time import time

from parse_qiskit import *

"""
Target: Verify the feasibility of using the Steane code to share a quantum secret.
The setting in this test: Verifiable quantum secret sharing protocol using the Steane code.

Result: The measurement outcomes are only correct codewords of the [7,4] Hamming code. And the data qubits remain unchanged.
"""

circ = QuantumCircuit(14, 14)

prepare_steane_code(circ, list(range(7)))
gateLength = len(circ.data)
circ.h(13)
prepare_steane_code(circ, list(range(7, 14)))
for i in range(7):
    circ.cx(i, 7+i)
ts = parse_qiskit(circ)
measList = applyMeasureAndReset(ts, 'ZZZZZZZ', list(range(7, 14)), 14)
# measList = applySelectiveFinalMeasurement(ts, 'ZZZZZZZ', list(range(7, 14)), 14)

op00 = pyqreach.QOperation(["00000000000000"])
ts.setAnnotation([[0, op00]])
start_time = time()
ts.computingFixedPointPost()
end_time = time()
print(f"Time taken for computing fixed point post: {end_time - start_time:.2f} seconds")

for i,measLoc in enumerate(measList):
    dims = ts.printDims(measLoc)
    if dims[1] > 0:
        print(f"Location {getBinary(i, 7)}: upperBound dimension = {dims[0]}, lowerBound dimension = {dims[1]}")
        ts.printSupp(measLoc)

# visualize_transition_system(ts, 'test_steane')
