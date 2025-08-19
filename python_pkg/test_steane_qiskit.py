import pyqreach
### Verifiable quantum secret sharing
from qiskit import QuantumCircuit, QuantumRegister, ClassicalRegister, transpile
from qiskit.quantum_info import Statevector
import numpy as np
from time import time

from parse_qiskit import *
from qctl import tsLabelling, tsLabellingClRegList

# Qiskit program
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

# Parse the Qiskit circuit into a transition system
cons_start_time = time()
ts = pyqreach.TransitionSystem()
resultList = parse_qiskit_cir(circ, circ.num_qubits, ts)
cons_end_time = time()
print(f"Time taken for constructing transition system: {cons_end_time - cons_start_time:.2f} seconds")
print("Transition System Locations:", ts.getLocationNum())
print("Result List size:", len(resultList))

# Set the initial state (annotation) of the transition system
op00 = pyqreach.QOperation(["00000000000000"])
ts.setAnnotation([[0, op00]])
start_time = time()
ts.computingFixedPointPost()
end_time = time()
print(f"Time taken for computing fixed point post: {end_time - start_time:.2f} seconds")


# Specify atomic propositions and the CTL formula
# Atomic proposition p
prop0 = ts.Locations[14].lowerBound
tsLabelling(ts, prop0, "p")
# Atomic proposition q: leaf nodes
for loc in range(ts.getLocationNum()):
    if ts.isLeafLoc(loc):
        ts.setLabel(loc, "q")
# Atomic proposition r: is a valid [7,4,3] Hamming code
pat = ["1111111", "0000111", "1001011", "0110011", "0101101", "1010101", "0011001", "1100001",
       "0011110", "1100110", "0101010", "1010010", "1001100", "0110100", "1111000", "0000000"]
# padding pat with 0s to match the number of clbits
pat = [p.ljust(14, '0') for p in pat]
tsLabellingClRegList(ts, pat, "r")

# Convert the transition system to a dictionary and then to a SMV file
graph_dic = ts2Dict(ts)
# graph_nx = dict2NX(graph_dic)
# nx2Graph_hierarchical(graph_nx, 'steane_dynamic_qiskit')
smv_content = dict2SMV(graph_dic, 'AG ((q & r) -> p)')
with open('unit_test.smv', 'w') as f:
    f.write(smv_content)
