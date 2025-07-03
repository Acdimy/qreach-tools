import pyqreach
### Verifiable quantum secret sharing
from qiskit import QuantumCircuit, QuantumRegister, ClassicalRegister, transpile
from qiskit.quantum_info import Statevector
# from qiskit_aer import AerSimulator
import numpy as np
from time import time

from parse_qiskit import *

"""
Target: Entanglement distillation test.
The setting in this test: 1. Forward. Over the right input, the protocol Ensure the correctness of the first pair of entangled states 
at the expense of the second pair of entangled states
2. Backward.

"""

circ = QuantumCircuit(4,4)
circ.h(0)
circ.cx(0,1)
circ.h(2)
circ.cx(2,3)

circ.cx(0,2)
circ.cx(1,3)
ts = parse_qiskit(circ)
measList = applyFinalMeasurement(ts, 'ZZ', [2,3], 4)

# Set annotations
# op00 = pyqreach.QOperation(["0000"])
# ts.setAnnotation([[0, op00]])
opI = pyqreach.CreateIdentityQO(4)
opO = pyqreach.CreateZeroQO(4)
for measLoc in measList:
    ts.setAnnotation([[measLoc, opO]])
ts.setAnnotation([[measList[1], opI]])

start_time = time()
ts.computingFixedPointPre()
end_time = time()
print(f"Time taken for computing fixed point pre: {end_time - start_time:.2f} seconds")
for locMeas in measList:
    ts.printSupp(locMeas)
ts.printSupp(4)
ts.printDims(4)

visualize_transition_system(ts, 'test_ent_distillation')
