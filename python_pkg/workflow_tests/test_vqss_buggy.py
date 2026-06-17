from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

### Verifiable quantum secret sharing
from inline_annotations import QReachCircuit
from time import time

from parse_qiskit import *
from qctl import *
from circ_utils import *

circ = QReachCircuit(10, 10)

prepare_5perfect_code(circ, [0,1,2,3,4])
circ.mark('enc')
circ.h(5)
prepare_5perfect_code(circ, [5,6,7,8,9])

for i in range(5):
    circ.cx(i, 5+i)
for i in range(5, 10):
    circ.measure(i, i)
for i in range(5, 10):
    with circ.if_test((i, 1)):
        circ.x(i)

ts = pyqreach.TransitionSystem()
parse_result = parse_qiskit_cir(circ, circ.num_qubits, ts, return_metadata=True)
set_initial_state(ts, "0000000000")
ts.computingFixedPointPost()
label_snapshot(ts, parse_result, "enc", "target")
annotate(ts, ["leaf"])
result = modelChecking(ts, 'AG (leaf -> target)')
print("Output: ", result["output"])
print("Model checking result:", result['satisfied'])
