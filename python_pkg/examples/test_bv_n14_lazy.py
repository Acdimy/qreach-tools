from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from qreach.inline_annotations import QReachCircuit
from time import time

from qreach.parse_qiskit import *
from qreach.qctl import *
from qreach.circ_utils import *

FILE_PATH = Path(__file__).resolve().parents[1].joinpath("benchmark/benchpress-medium/supported/bv_n14.qasm")
qc = QReachCircuit.from_qasm_file(FILE_PATH)
start_time = time()
ts = pyqreach.TransitionSystem()
parse_result = parse_qiskit_cir_lazy(qc, qc.num_qubits, ts, initial_state="0"*qc.num_qubits, return_metadata=True)
build_time = time() - start_time
ts.computingFixedPointPost()
check_time = time() - start_time - build_time
print(f"Transition System Locations: {ts.getLocationNum()}")
print(f"Lazy pruned locations: {len(parse_result.lazy_pruned_locations)}")
print(f"Time taken for building transition system: {build_time:.2f} seconds")
print(f"Time taken for model checking: {check_time:.2f} seconds")
