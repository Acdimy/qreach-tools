import pyqreach
### Verifiable quantum secret sharing
from qiskit import QuantumCircuit, QuantumRegister, ClassicalRegister, transpile
from qiskit.quantum_info import Statevector
import numpy as np
from time import time

from parse_qiskit import *

"""
Target: Verify the infeasibility of using the 5-perfect code to share a quantum secret.
The setting in this test: Trying to use the error-detection procedure of 5-perfect code. Measure 'ZXXZI' on the distributed ancilla qubits.

Result: The data qubits also collapsed.
"""

circ = QuantumCircuit(10, 10)

# circ.x(0)

prepare_5perfect_code(circ, [0,1,2,3,4])
prepare_5perfect_code(circ, [5,6,7,8,9])

# for i in range(5):
#     circ.cx(i, 5+i)

circ.cx(0,5)
circ.h(6)
circ.cx(6,1)
circ.h(6)
circ.h(7)
circ.cx(7,2)
circ.h(7)
circ.cx(3,8)
circ.cx(4,9)

ts = parse_qiskit(circ)
measList = applyFinalMeasurement(ts, 'ZZZZI', [5,6,7,8,9], 10)

ts.addLocation(pyqreach.Location(10, 0))
locMerge = ts.getLocationNum() - 1
id = pyqreach.QOperation("I", 10, [0], [])
recoverLoc = [0,3,5,6,9,10,12,15]
for i in recoverLoc:
    locMeas = measList[i]
    ts.addRelation(locMeas, locMerge, id)

op00 = pyqreach.QOperation(["0000000000"])
ts.setAnnotation([[0, op00]])
start_time = time()
ts.computingFixedPointPost()
end_time = time()
print(f"Time taken for computing fixed point post: {end_time - start_time:.2f} seconds")
ts.printDims(locMerge)
ts.printSupp(58)

# visualize_transition_system(ts, 'test_513_bug2')
