from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

### Verifiable quantum secret sharing
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

ts = pyqreach.TransitionSystem()
parse_result = parse_qiskit_cir(circ, circ.num_qubits, ts, return_metadata=True)
set_initial_state(ts, "00000000000000")
ts.computingFixedPointPost()
label_snapshot(ts, parse_result, "enc", "target")
annotate(ts, ["leaf"])
result = modelChecking(ts, 'AG (leaf -> target)')
print("Output: ", result["output"])
print("Model checking result:", result['satisfied'])
