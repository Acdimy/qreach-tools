import pyqreach
### Verifiable quantum secret sharing
from qiskit import QuantumCircuit, QuantumRegister, ClassicalRegister, transpile
from qiskit.quantum_info import Statevector
import numpy as np
from time import time
from qiskit_aer import AerSimulator

from parse_qiskit import *
from qctl import *
from circ_utils import *

qubits = QuantumRegister(3)
midBits = ClassicalRegister(2, name='mid')
ancBits = ClassicalRegister(1, name='anc')
circ = QuantumCircuit(qubits, midBits, ancBits)

theta = 2 * np.arccos(np.sqrt(0.4))
circ.u(theta, 0, 0, 0)
circ.u(theta, 0, 0, 1)
circ.cx(0, 1)
circ.h(0)
circ.measure(0,0)
circ.measure(1,1)


ts = pyqreach.TransitionSystem()
resultList = parse_qiskit_cir(circ, circ.num_qubits, ts)
print("Transition System Locations:", ts.getLocationNum())
print("Result List size:", len(resultList))

op00 = pyqreach.QOperation(["000"])
ts.setAnnotation([[0, op00]])
start_time = time()
ts.computingFixedPointPost()
end_time = time()
print(f"Time taken for computing fixed point post: {end_time - start_time:.2f} seconds")

print(ts.printDims(9))

graph_dic = ts2Dict(ts)
graph_nx = dict2NX(graph_dic)
nx2Graph_hierarchical(graph_nx, 'unit_test_qbf')
