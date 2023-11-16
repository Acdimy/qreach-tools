import numpy as np
from qiskit import QuantumCircuit
from qiskit import Aer, execute
from qiskit.quantum_info import Operator
from qiskit.quantum_info import Statevector
import copy

def getCirEntry(cir, i, j):
    # bin_j = bin(j).replace('0b', '')
    qnum = cir.num_qubits
    circ = copy.deepcopy(cir)
    initVector = [0 for _ in range(2**qnum)]
    initVector[j] = 1
    # circ.initialize(initVector, list(range(qnum)))
    circ.initialize([1,0], 0)
    simulator = Aer.get_backend('statevector_simulator')
    result = execute(circ, simulator).result()
    statevector = result.get_statevector(circ)
    return statevector[i]

def getUnitary(cir):
    U = Operator(cir)
    return U.data

def simulate(stateInt:int, cir):
    qnum = cir.num_qubits
    state = Statevector.from_int(stateInt, 2**qnum)
    state = state.evolve(cir)
    return state
