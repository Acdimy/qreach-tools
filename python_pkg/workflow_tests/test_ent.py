from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import pyqreach
### Verifiable quantum secret sharing
from inline_annotations import QReachCircuit
from qiskit.quantum_info import Statevector
# from qiskit_aer import AerSimulator
import numpy as np
from time import time

from parse_qiskit import *
from qctl import *
from circ_utils import *

"""
Target: Entanglement distillation test.
The setting in this test: 1. Forward. Over the right input, the protocol Ensure the correctness of the first pair of entangled states 
at the expense of the second pair of entangled states
2. Backward.

"""

circ = QReachCircuit(4,4)
circ.h(0)
circ.cx(0,1)
circ.h(2)
circ.cx(2,3)

circ.cx(0,2)
circ.cx(1,3)
circ.measure(2,2)
circ.measure(3,3)

ts = pyqreach.TransitionSystem()
parse_result = parse_qiskit_cir(circ, circ.num_qubits, ts, return_metadata=True)
opI = pyqreach.CreateIdentityQO(4)
opO = pyqreach.CreateZeroQO(4)
#TODO: immitate the codes in python_pkg/test_ent.py, but use the new api.