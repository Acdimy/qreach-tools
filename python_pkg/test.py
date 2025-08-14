import pyqreach
### Verifiable quantum secret sharing
from qiskit import QuantumCircuit, QuantumRegister, ClassicalRegister, transpile
from qiskit.quantum_info import Statevector
import numpy as np
from time import time

from parse_qiskit import *
from qctl import tsLabelling


circ = QuantumCircuit(4, 4)

circ.h(0)
circ.h(1)
circ.cx(0, 2)
circ.cx(1, 3)
circ.measure(0, 0)
with circ.if_test((0, 1)):
    circ.x(2)

ts = pyqreach.TransitionSystem()
resultList = parse_qiskit_cir(circ, circ.num_qubits, ts)
print("Transition System Locations:", ts.getLocationNum())
print("Result List size:", len(resultList))
op00 = pyqreach.QOperation(["0000"])
ts.setAnnotation([[0, op00]])
start_time = time()
ts.computingFixedPointPost()
end_time = time()
print(f"Time taken for computing fixed point post: {end_time - start_time:.2f} seconds")
# for resloc in resultList:
#     ts.printSupp(resloc)
# visualize_transition_system(ts, 'unit_test')
prop0 = ts.Locations[2].lowerBound
# print(prop0.normalized)
tsLabelling(ts, prop0, "p")

graph_dic = ts2Dict(ts)
graph_nx = dict2NX(graph_dic)
# nx2Graph_hierarchical(graph_nx, 'unit_test_nx')
smv_content = dict2SMV(graph_dic, 'AF p')
# Save the SMV content to a file
with open('unit_test.smv', 'w') as f:
    f.write(smv_content)
