import pyqreach
### Verifiable quantum secret sharing
from qiskit import QuantumCircuit, QuantumRegister, ClassicalRegister, transpile
from qiskit.quantum_info import Statevector
import numpy as np
from time import time

from parse_qiskit import *
from qctl import *
from circ_utils import *


circ = QuantumCircuit(10, 10)

prepare_5perfect_code(circ, [0,1,2,3,4])
circ.h(5)
prepare_5perfect_code(circ, [5,6,7,8,9])

for i in range(5):
    circ.cx(i, 5+i)
for i in range(5, 10):
    circ.measure(i, i)
for i in range(5, 10):
    with circ.if_test((i, 1)):
        circ.x(i)

# Parse the Qiskit circuit into a transition system
cons_start_time = time()
ts = pyqreach.TransitionSystem()
resultList = parse_qiskit_cir(circ, circ.num_qubits, ts)
cons_end_time = time()
print(f"Time taken for constructing transition system: {cons_end_time - cons_start_time:.2f} seconds")
print("Transition System Locations:", ts.getLocationNum())
print("Result List size:", len(resultList))

# Set the initial state (annotation) of the transition system
op00 = pyqreach.QOperation(["0000000000"])
ts.setAnnotation([[0, op00]])
start_time = time()
ts.computingFixedPointPost()
end_time = time()
print(f"Time taken for computing fixed point post: {end_time - start_time:.2f} seconds")

# Specify atomic propositions and the CTL formula
# Atomic proposition p
ap_start_time = time()
prop0 = ts.Locations[17].lowerBound
tsLabelling(ts, prop0, "p")
# Atomic proposition q: leaf nodes
for loc in range(ts.getLocationNum()):
    if ts.isLeafLoc(loc):
        ts.setLabel(loc, "q")
ap_end_time = time()
print(f"Time taken for atomic proposition labelling: {ap_end_time - ap_start_time:.2f} seconds")

# Model checking
model_start_time = time()
result = modelChecking(ts, 'AG (q -> p)')
model_end_time = time()
print(f"Time taken for model checking: {model_end_time - model_start_time:.2f} seconds")
print("Model checking result:", result['satisfied'])
print("Model checking counterexample:", result['counterexample'])

# visualize_transition_system(ts, "ts_513_qiskit")

