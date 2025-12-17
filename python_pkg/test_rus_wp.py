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

qc = QuantumCircuit(3,1)
qc.x(2)
qc.measure(2, 0)
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
    
pyqreach.initializeTransitionSystem()
ts = pyqreach.TransitionSystem(False)
resultList = parse_qiskit_cir(qc, qc.num_qubits, ts)
opinit = pyqreach.QOperation(["000"])
ts.setAnnotation([[0, opinit]])
ts.computingFixedPointPost()
op_final = ts.Locations[resultList[-1]].lowerBound
ts.printSupp(resultList[-1])

qc_check = QuantumCircuit(3,1)
qc_check.x(2)
qc_check.measure(2, 0)
with qc_check.while_loop((0, 0b1)):
    qc_check.h(0)
    qc_check.t(0)
    qc_check.cx(0, 1)
    qc_check.h(0)
    qc_check.cx(0, 1)
    qc_check.t(0)
    qc_check.h(0)
    qc_check.measure(0, 0)
ts_check = pyqreach.TransitionSystem(False)
resultList_check = parse_qiskit_cir(qc_check, qc_check.num_qubits, ts_check)
ts_check.setAnnotation([[resultList_check[-1], op_final]])
visualize_transition_system(ts_check, "RUS_fix")
ts_check.computingFixedPointPre()
ts_check.printSupp(0)
