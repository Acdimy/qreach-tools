import pyqreach
### Verifiable quantum secret sharing
from qiskit import QuantumCircuit, QuantumRegister, ClassicalRegister, transpile
from qiskit.circuit import CircuitInstruction
from qiskit.circuit.library import XGate, YGate, ZGate
from qiskit.quantum_info import Statevector
from qiskit_aer import Aer
import numpy as np
from math import pi
import random
from time import time
from parse_qiskit import *
from qctl import *
from circ_utils import *
import pandas as pd

filename = "benchmark/grover/grover_5.qasm"
qc = QuantumCircuit.from_qasm_file(filename)
num_qubits = qc.num_qubits

pyqreach.initializeTransitionSystem()

ts = pyqreach.TransitionSystem(False)
resultList = parse_qiskit_cir(qc, num_qubits, ts)
op00 = pyqreach.QOperation(["0"*num_qubits])
ts.setAnnotation([[0, op00]])
time_post_start = time()
ts.computingFixedPointPost()
time_post_end = time()

work_qubits = int((num_qubits+1)/2)
grover_init = ts.Locations[work_qubits].lowerBound # superposition |++...+00...01>
grover_good = pyqreach.QOperation(["1"*work_qubits + "0"*(num_qubits - work_qubits - 1) + "1"]) # good state |11...100...01>
ts_temp = pyqreach.TransitionSystem(False)
loc0, loc1, loc2 = pyqreach.Location(num_qubits,0), pyqreach.Location(num_qubits,1), pyqreach.Location(num_qubits,2)
ts_temp.addLocation(loc0)
ts_temp.addLocation(loc1)
ts_temp.addLocation(loc2)
ts_temp.addRelation(0, 2, pyqreach.QOperation("I", num_qubits, [0], []))
ts_temp.addRelation(1, 2, pyqreach.QOperation("I", num_qubits, [0], []))
ts_temp.setAnnotation([[0, grover_init], [1, grover_good]])
ts_temp.computingFixedPointPost()
grover_final = ts_temp.Locations[2].lowerBound

ts_temp.printSupp(2)

# ts.printSupp(0)
# ts.printSupp(resultList[-1])

ts.resetLocationBounds()
ts.setAnnotation([[resultList[-1], grover_final]])
# ts.setAnnotation([[0, pyqreach.CreateIdentityQO(num_qubits)]])
time_pre_start = time()
ts.computingFixedPointPre()
# ts.computingFixedPointPost()
time_pre_end = time()
print(f"Time for post computation: {time_post_end - time_post_start} seconds")
print(f"Time for all computation: {time_pre_end - time_pre_start} seconds")
ts.printSupp(0)
ts.printSupp(resultList[-1])
