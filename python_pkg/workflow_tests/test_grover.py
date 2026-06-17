from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

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

qc = QuantumCircuit.from_qasm_file("benchmark/grover/grover_5.qasm")
ts = pyqreach.TransitionSystem()
parse_result = parse_qiskit_cir(qc, qc.num_qubits, ts, return_metadata=True)
set_initial_state(ts, "00001")
ts.computingFixedPointPost()
work_qubits = int((qc.num_qubits+1)/2)
# target_subspace = span_states(["1"*work_qubits + "0"*(qc.num_qubits - work_qubits - 1) + "1", "+"*work_qubits + "0"*(qc.num_qubits - work_qubits - 1) + "1"])
target_subspace = span_states(["+++01", "11101"])
annotate(ts, ["leaf"])
# res = ts.Locations[parse_result[-1]].satisfy(target_subspace)
# print("Grover benchmark result:", res)

# the formula: leaf -> target_subspace
