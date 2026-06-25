from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

### Verifiable quantum secret sharing
from inline_annotations import QReachCircuit
from time import time

from parse_qiskit import *
from qctl import *
from circ_utils import *

FILE_PATH = Path(__file__).resolve().parents[1].joinpath("benchmark/benchpress-medium/supported/bv_n14.qasm")
qc = QReachCircuit.from_qasm_file(FILE_PATH)
start_time = time()
ts = pyqreach.TransitionSystem()
parse_result = parse_qiskit_cir(qc, qc.num_qubits, ts, return_metadata=True)
set_initial_state(ts, "0"*qc.num_qubits)
build_time = time() - start_time
ts.computingFixedPointPost()
check_time = time() - start_time - build_time
print(f"Transition System Locations: {ts.getLocationNum()}")
print(f"Time taken for building transition system: {build_time:.2f} seconds")
print(f"Time taken for model checking: {check_time:.2f} seconds")
