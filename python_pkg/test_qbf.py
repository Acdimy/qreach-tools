import pyqreach
### Verifiable quantum secret sharing
from qiskit import QuantumCircuit, QuantumRegister, ClassicalRegister, transpile
from qiskit.quantum_info import Statevector
# from qiskit_aer import AerSimulator
import numpy as np

from parse_qiskit import *

circ = QuantumCircuit(6, 6)
idx0 = list(range(5))

prepare_5perfect_code(circ, idx0)
gateLength = len(circ.data)
circ.cx(0,5)
circ.h(5)
circ.cx(5,1)
circ.cx(5,2)
circ.h(5)
circ.cx(3,5)

ts = parse_qiskit(circ)
measList = applyFinalMeasurement(ts, 'Z', [5], 5)
# measList = applyFinalMeasurement(ts, 'ZXXZI', [0,1,2,3,4], 5)

op00 = pyqreach.QOperation(["000000"])
ts.setAnnotation([[0, op00]])
ts.computingFixedPointPost()
ts.printSupp(gateLength)
ts.printSupp(measList[0])
