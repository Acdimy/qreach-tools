import pyqreach
from qiskit import QuantumCircuit, QuantumRegister, ClassicalRegister, transpile
from qiskit.quantum_info import Statevector
import numpy as np

from parse_qiskit import *


circ = QuantumCircuit(2, 2)

circ.h(0)
circ.cx(0, 1)

ts = parse_qiskit(circ)

# measOp = pyqreach.QOperation('meas0', 2, [0], [])
# ts.addLocation(pyqreach.Location(2, 0))
# ts.addRelation(2, 3, measOp)

results = applyFinalMeasurement(ts, 'ZZ', [0, 1], 2)


# locNum = ts.getLocationNum()
# print("Number of locations in the transition system:", locNum)
# resultsPerEntry, errorLoc = applySelectiveFinalMeasurement(ts, 'ZZ', [0,1], 2, [locNum-1], ['00','11'])
op00 = pyqreach.QOperation(["00"])
ts.setAnnotation([[0, op00]])
ts.computingFixedPointPost()
for measLoc in results:
    ts.printSupp(measLoc)

visualize_transition_system(ts)
print("Results of final measurement:", results)
