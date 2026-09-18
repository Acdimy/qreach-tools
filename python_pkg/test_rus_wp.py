"""Legacy RUS weakest-precondition example.

This script now follows the same single-transition-system workflow as the
maintained workflow test, so it exercises the repaired reset semantics instead
of the older two-TS pattern.
"""

import pyqreach
from time import time

from qiskit import QuantumCircuit

from parse_qiskit import parse_qiskit_cir
from qctl import set_initial_state


qc = QuantumCircuit(3, 1)
qc.x(2)
qc.measure(2, 0)
with qc.while_loop((0, 0b1)):
    qc.reset(0)
    qc.h(0)
    qc.t(0)
    qc.cx(0, 1)
    qc.h(0)
    qc.cx(0, 1)
    qc.t(0)
    qc.h(0)
    qc.measure(0, 0)

pyqreach.initializeTransitionSystem()
ts = pyqreach.TransitionSystem(False)
result_list = parse_qiskit_cir(qc, qc.num_qubits, ts)
set_initial_state(ts, "000")
ts.computingFixedPointPost()

final_loc = result_list[-1]
op_final = ts.Locations[final_loc].lowerBound
print("=== Post-condition at Final Location ===")
ts.printSupp(final_loc)

ts.resetLocationBounds()
ts.setAnnotation([[final_loc, op_final]])

time_pre_start = time()
ts.computingFixedPointPre()
time_pre_end = time()

print(f"Backward reachability (pre): {time_pre_end - time_pre_start:.3f}s")
print("\n=== Weakest Pre-condition at Initial Location ===")
ts.printSupp(0)
