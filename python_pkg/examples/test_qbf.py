from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import pyqreach
### Verifiable quantum secret sharing
from qiskit import QuantumCircuit, QuantumRegister, ClassicalRegister
import numpy as np
from time import time
from qreach.parse_qiskit import *
from qreach.qctl import *
from qreach.circ_utils import *

# Integration Test: Quantum Bernoulli Factory
qubits = QuantumRegister(3)
midBits = ClassicalRegister(2, name='mid')
ancBits = ClassicalRegister(1, name='anc')
circ = QuantumCircuit(qubits, midBits, ancBits)
theta = 2 * np.arccos(np.sqrt(0.2))
circ.reset([0, 1, 2])
k = 2
for i in range(k):
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
        with circ.if_test((midBits, 0b00)): # Change here to simulate a bug
            circ.x(2)
        with circ.if_test((midBits, 0b11)):
            circ.x(2)
        circ.measure(2,2)
    if i < k-1:
        with circ.if_test((midBits, 0b10)): # The first register is 0, the second is 1
            circ.reset([0, 1, 2])
            circ.measure(0,0)
            circ.measure(1,1)
            circ.measure(2,2)

ts = pyqreach.TransitionSystem()
parse_qiskit_cir(circ, circ.num_qubits, ts, return_metadata=True)
set_initial_state(ts, "000")
ts.computingFixedPointPost()
annotate(ts, ["leaf", "loop", "reached"])
annotate_classical(ts, "t", "10")
annotate_classical(ts, "f", ["00", "11"])
annotate_classical(ts, "h", "01")
result = modelChecking(ts, 'AG ((t & valid) -> ! E [valid U (valid & leaf & ! t)])')
print("Output: ", result["output"])
print("Model checking result:", result['satisfied'])
