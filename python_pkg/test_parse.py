import pyqreach
from qiskit import QuantumCircuit, QuantumRegister, ClassicalRegister, transpile
from qiskit.quantum_info import Statevector
import numpy as np

from parse_qiskit import *


circ = QuantumCircuit(2, 2)

circ.h(0)
circ.t(0)
circ.cx(0, 1)

ts = parse_qiskit(circ)

locNum = ts.getLocationNum()
print("Number of locations in the transition system:", locNum)

resultsPerEntry, errorLoc = applySelectiveFinalMeasurement(ts, 'ZZ', [0,1], 2, [locNum-1], ['00','11'])

visualize_transition_system(ts)

