from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from inline_annotations import QReachCircuit
from time import time

from parse_qiskit import *
from qctl import *
from circ_utils import *

# Qiskit program
circ = QReachCircuit(14, 14)
prepare_steane_code(circ, list(range(7)))
circ.mark('enc')
circ.h(13)
prepare_steane_code(circ, list(range(7, 14)))
for i in range(7):
    circ.cx(i, 7+i)
for i in range(7, 14):
    circ.measure(i, i)
for i in range(7, 14):
    with circ.if_test((i, 1)):
        circ.x(i)

start_time = time()
ts = pyqreach.TransitionSystem()
parse_result = parse_qiskit_cir_lazy(circ, circ.num_qubits, ts, initial_state="00000000000000", return_metadata=True)
build_time = time() - start_time
ts.computingFixedPointPost()
label_snapshot(ts, parse_result, "enc", "target")
annotate(ts, ["leaf"])
result = modelChecking(ts, 'AG (leaf -> target)')
check_time = time() - start_time - build_time
print(f"Transition System Locations: {ts.getLocationNum()}")
print(f"Lazy pruned locations: {len(parse_result.lazy_pruned_locations)}")
print(f"Time taken for building transition system: {build_time:.2f} seconds")
print(f"Time taken for model checking: {check_time:.2f} seconds")
print("Output: ", result["output"])
print("Model checking result:", result['satisfied'])
