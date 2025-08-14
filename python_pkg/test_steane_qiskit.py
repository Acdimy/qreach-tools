import pyqreach
### Verifiable quantum secret sharing
from qiskit import QuantumCircuit, QuantumRegister, ClassicalRegister, transpile
from qiskit.quantum_info import Statevector
import numpy as np
from time import time

from parse_qiskit import *
from qctl import tsLabelling


circ = QuantumCircuit(14, 14)
prepare_steane_code(circ, list(range(7)))
circ.h(13)
prepare_steane_code(circ, list(range(7, 14)))
for i in range(7):
    circ.cx(i, 7+i)
for i in range(7, 14):
    circ.measure(i, i)
for i in range(7, 14):
    with circ.if_test((i, 1)):
        circ.x(i)

ts = pyqreach.TransitionSystem()
resultList = parse_qiskit_cir(circ, circ.num_qubits, ts)
print("Transition System Locations:", ts.getLocationNum())
print("Result List size:", len(resultList))

op00 = pyqreach.QOperation(["00000000000000"])
ts.setAnnotation([[0, op00]])
start_time = time()
ts.computingFixedPointPost()
end_time = time()
print(f"Time taken for computing fixed point post: {end_time - start_time:.2f} seconds")

# visualize_transition_system(ts, 'steane_dynamic_qiskit')

# Atomic proposition p
prop0 = ts.Locations[14].lowerBound
tsLabelling(ts, prop0, "p")
# Atomic proposition q: leaf nodes
for loc in range(ts.getLocationNum()):
    if ts.isLeafLoc(loc):
        ts.setLabel(loc, "q")
# Atomic proposition r: is a valid [7,4,3] Hamming code
graph_dic = ts2Dict(ts)
# graph_nx = dict2NX(graph_dic)
# nx2Graph_hierarchical(graph_nx, 'steane_dynamic_qiskit')
smv_content = dict2SMV(graph_dic, 'AF p')
