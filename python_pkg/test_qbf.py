import pyqreach
### Verifiable quantum secret sharing
from qiskit import QuantumCircuit, QuantumRegister, ClassicalRegister, transpile
from qiskit.quantum_info import Statevector
# from qiskit_aer import AerSimulator
import numpy as np

from parse_qiskit import *

def p_coin(circ: QuantumCircuit, p: float, k: int):
    par = 2 * np.arccos(np.sqrt(p))
    for i in range(k):
        circ.ry(par, 0)
        circ.ry(par, 1)
        circ.cx(0, 1)
        circ.h(0)
        circ.measure(0, 0)
        circ.measure(1, 1)
        circ.reset(0)
        circ.reset(1)


circ = QuantumCircuit(2, 2)

