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
from qasm_workflow_runner import insert_random_pauli

_BENCH_DIR = Path(__file__).resolve().parents[1] / "benchmark"
qc = QuantumCircuit.from_qasm_file(str(_BENCH_DIR / "grover" / "grover_5.qasm"))
# set random seed for reproducibility
rng = random.Random(42)
qc, _ = insert_random_pauli(qc, rng)
ts = pyqreach.TransitionSystem()
parse_result = parse_qiskit_cir(qc, qc.num_qubits, ts, return_metadata=True)
set_initial_state(ts, "00001")
ts.computingFixedPointPost()
work_qubits = int((qc.num_qubits+1)/2)
# target_subspace = span_states(["1"*work_qubits + "0"*(qc.num_qubits - work_qubits - 1) + "1", "+"*work_qubits + "0"*(qc.num_qubits - work_qubits - 1) + "1"])
target_subspace = span_states(["+++01", "11101"])
annotate(ts, ["leaf"])
tsLabelling(ts, target_subspace, "target_subspace", locList=leaf_locations(ts))
result = modelChecking(ts, "AG (leaf -> target_subspace)")
print("Output: ", result["output"])
print("Model checking result:", result['satisfied'])
if result.get('analysis'):
    from qctl import _format_counterexample_analysis
    print(_format_counterexample_analysis(result['analysis']))
# res = ts.Locations[parse_result[-1]].satisfy(target_subspace)
# print("Grover benchmark result:", res)

# the formula: leaf -> target_subspace
