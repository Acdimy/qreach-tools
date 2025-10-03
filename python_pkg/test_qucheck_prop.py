import pyqreach
from qiskit import QuantumCircuit, QuantumRegister, ClassicalRegister, transpile
from qiskit.quantum_info import Statevector
from qiskit_aer import Aer
from qiskit.quantum_info import Pauli
import numpy as np
from math import pi
import random
from time import time
from parse_qiskit import *
from qctl import *
from circ_utils import *

# Pack the qucheck assertions: assert_equal, assert_entangled, assert_separable, etc.

# assert_equal: For all resultList of the transition system of the program, the state should be equal to the expected state
def assert_equal(circuit1, circuit2):
    qubits1 = circuit1.num_qubits
    qubits2 = circuit2.num_qubits
    if qubits1 != qubits2:
        return False
    ts1 = pyqreach.TransitionSystem()
    resultList1 = parse_qiskit_cir(circuit1, qubits1, ts1)
    # For all locations in resultList1, add a new location and connect with identity transition
    if len(resultList1) > 1:
        ts1.addLocation(pyqreach.Location(qubits1, 0))
        for resultLoc in resultList1:
            ts1.addRelation(resultLoc, ts1.getLocationNum()-1, pyqreach.QOperation("I", qubits1, [0], []))
    opinit1 = pyqreach.QOperation(["0"*qubits1])
    ts1.setAnnotation([[0, opinit1]])

    circ1EndNode = ts1.getLocationNum()-1
    ts1.addLocation(pyqreach.Location(qubits1, 0))
    newStartLoc = ts1.getLocationNum()-1


    resultList2 = parse_qiskit_cir(circuit2, qubits2, ts1, startNodes=[newStartLoc])
    # For all locations in resultList2, add a new location and connect with identity transition, and extract the lowerBound of the location as expectedState
    if len(resultList2) > 1:
        ts1.addLocation(pyqreach.Location(qubits2, 0))
        for resultLoc in resultList2:
            ts1.addRelation(resultLoc, ts1.getLocationNum()-1, pyqreach.QOperation("I", qubits2, [0], []))
    opinit2 = pyqreach.QOperation(["0"*qubits2])
    ts1.setAnnotation([[newStartLoc, opinit2]])

    
    ts1.computingFixedPointPost()

    prop0 = ts1.Locations[ts1.getLocationNum()-1].lowerBound
    # ts1.printSupp(3)
    # ts1.printSupp(5)
    # visualize_transition_system(ts1, "unit_test_assert_equal")
    if ts1.Locations[circ1EndNode].satisfy(prop0):
        return True
    return False
    

def assert_entangled(qubit_indices, circuit, basis='Z'):
    if basis not in ['X', 'Z']:
        raise ValueError("Basis must be 'X' or 'Z'")
    qubits = circuit.num_qubits
    ts = pyqreach.TransitionSystem()
    resultList = parse_qiskit_cir(circuit, qubits, ts)
    resultList = applyFinalMeasurement(ts, basis*len(qubit_indices), qubit_indices, qubits)

    op0 = pyqreach.QOperation(["0"*qubits])
    ts.setAnnotation([[0, op0]])
    ts.computingFixedPointPost()
    print(ts.printDims(2))
    ent_type = ts.printDims(resultList[0])[1] > 0
    for i, loc in enumerate(resultList):
        if i % 4 == 0 or i % 4 == 3:
            if ts.printDims(loc)[1] != ent_type:
                return False
        else:
            if ts.printDims(loc)[1] == ent_type:
                return False
    return True




# Entanglement test
# qc = QuantumCircuit(2)
# qc.h(0)
# qc.cx(0, 1)
# qc.x(1)
# print(assert_entangled([0, 1], qc, basis='Z'))

# Equality test


circ1 = QuantumCircuit(2)
circ1.z(0)
circ1.y(0)
circ1.append(Pauli('-iX'), [0])

circ2 = QuantumCircuit(2)
circ2.id(0)

print(assert_equal(circ1, circ2))
