import pyqreach
### Verifiable quantum secret sharing
from qiskit import QuantumCircuit, QuantumRegister, ClassicalRegister, transpile
from qiskit.quantum_info import Statevector
from qiskit_aer import Aer
import numpy as np
from math import pi
import random
from time import time
from parse_qiskit import *
from qctl import *
from circ_utils import *

# qnum = 8

# # Initialize the transition system
# ts = pyqreach.TransitionSystem()
# ts.setInitLocation(0)

# # Add locations to the transition system
# # Here we create 5 locations with indices from 0 to 4
# loc_list = []
# for i in range(4):
#     loc = pyqreach.Location(qnum,i)
#     ts.addLocation(loc)
#     loc_list.append(loc)

# # Define the operations
# op0 = pyqreach.QOperation(["00000000"])
# op1 = pyqreach.QOperation(["10000000"])
# oph = pyqreach.QOperation("H", qnum, [0], [])
# opx = pyqreach.QOperation("X", qnum, [0], [])
# opy = pyqreach.QOperation("Y", qnum, [0], [])
# opz = pyqreach.QOperation("Z", qnum, [0], [])
# opm0 = pyqreach.QOperation("meas0", qnum, [0], [])
# opm1 = pyqreach.QOperation("meas1", qnum, [0], [])
# opi = pyqreach.QOperation("I", qnum, [0], [])

# # Add the operations to the transition system
# # The topological order of the transition system: 
# # 0 -> 1, 1 -> 2, 1 -> 3, 2 -> 4, 3 -> 4, 4 -> 4
# ts.addRelation(0, 1, oph)
# ts.addRelation(1, 2, opm0)
# ts.addRelation(1, 3, opm1)
# ts.addRelation(2, 1, oph)
# ts.addRelation(3, 3, opi)

# ts.setAnnotation([[3, op1]])

# ts.computingFixedPointPre()
# ts.printDims(1)
# ts.printDims(2)
# ts.printDims(3)
# ts.printSupp(0)

qc = QuantumCircuit(3, 1)
qc.x(2)
qc.measure(2, 0)
qc.initialize([1/np.sqrt(2), 1/np.sqrt(2)], 1)
with qc.while_loop((0, 0b1)):
    qc.reset(0)
    qc.h(0)
    qc.t(0)
    qc.cx(0, 1)
    qc.h(0)
    qc.cx(0, 1)
    qc.t(0)
    qc.h(0)
    qc.measure(0, 0)

ts = pyqreach.TransitionSystem()
resultList = parse_qiskit_cir(qc, qc.num_qubits, ts)
opinit = pyqreach.QOperation(["000"])
ts.setAnnotation([[0, opinit]])
ts.computingFixedPointPost()
