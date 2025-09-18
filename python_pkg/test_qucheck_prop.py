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

    ts2 = pyqreach.TransitionSystem()
    resultList2 = parse_qiskit_cir(circuit2, qubits2, ts2)
    # For all locations in resultList2, add a new location and connect with identity transition, and extract the lowerBound of the location as expectedState
    if len(resultList2) > 1:
        ts2.addLocation(pyqreach.Location(qubits2, 0))
        for resultLoc in resultList2:
            ts2.addRelation(resultLoc, ts2.getLocationNum()-1, pyqreach.QOperation("I", qubits2, [0], []))
    opinit2 = pyqreach.QOperation(["0"*qubits2])
    ts2.setAnnotation([[0, opinit2]])

    ts1.computingFixedPointPost()
    ts2.computingFixedPointPost()

    prop0 = ts2.Locations[ts2.getLocationNum()-1].lowerBound
    if ts1.Locations[ts1.getLocationNum()-1].satisfy(prop0):
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

from qiskit.quantum_info import Pauli
def quantum_fourier_transform_em(qubits, swap=True):
    # build circuit
    qft = QuantumCircuit(qubits, qubits)

    # modify phase
    for qubit in range(qubits):
        # insert the initial hadamard gate on all qubits in the register

        # semantic preserving changes loc 1, gates 12, index 0
        qft.z(0)
        qft.y(0)
        qft.append(Pauli('-iX'), [0])

        qft.h(qubit)

        # iterate across all indexes to get the appropriate controlled gates
        for offset in range(1, qubits - qubit):
            control_index = qubit + offset
            target_index = qubit
            rotation_amount = (np.pi / 2 ** offset)
            qft.cp(rotation_amount, control_index, target_index)

    # do swaps
    if swap:
        for qubit in range(qubits // 2):
            qft.swap(qubit, qubits - 1 - qubit)
    return qft

def quantum_fourier_transform(qubits, swap=True):
    # build circuit
    qft = QuantumCircuit(qubits, qubits)

    # modify phase
    for qubit in range(qubits):
        # insert the initial hadamard gate on all qubits in the register
        qft.h(qubit)

        # iterate across all indexes to get the appropriate controlled gates
        for offset in range(1, qubits - qubit):
            control_index = qubit + offset
            target_index = qubit
            rotation_amount = (np.pi / 2 ** offset)
            qft.cp(rotation_amount, control_index, target_index)

    # do swaps
    if swap:
        for qubit in range(qubits // 2):
            qft.swap(qubit, qubits - 1 - qubit)
    return qft


# Entanglement test
# qc = QuantumCircuit(2)
# qc.h(0)
# qc.cx(0, 1)
# qc.x(1)
# print(assert_entangled([0, 1], qc, basis='Z'))

# Equality test
circ1 = quantum_fourier_transform(8)
circ2 = quantum_fourier_transform_em(8)
print(assert_equal(circ1, circ2))
