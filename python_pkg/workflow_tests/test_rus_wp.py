"""Weakest-precondition verification of an RUS (Repeat-Until-Success) circuit.

This is a modernized version of ``test_rus_wp.py`` using the current
``qctl.py`` APIs.

The workflow:
  1. Parse the RUS circuit and compute forward reachability.
  2. Snapshot the quantum post-condition at the final location.
  3. Reset the TS bounds and attach that post-condition back to the leaf.
  4. Run ``computingFixedPointPre`` on the same explicit TS.

This single-TS workflow intentionally avoids the older two-parse pattern,
which is still susceptible to unrelated CFLOBDD global-state crashes when a
``QOperation`` computed from one transition system is transferred into
another in the same process.
"""

import pyqreach
from time import time

from qiskit import QuantumCircuit

from parse_qiskit import parse_qiskit_cir
from qctl import (
    set_initial_state,
)


# ============================================================================
# Step 1: Forward reachability on the RUS circuit
# ============================================================================
print("=" * 60)
print("Step 1: Forward reachability on original RUS circuit")
print("=" * 60)

qc = QuantumCircuit(3, 1)
qc.x(2)
qc.measure(2, 0)
with qc.while_loop((0, 0b1)):
    qc.reset(0)          # ← present in the original circuit
    qc.h(0)
    qc.t(0)
    qc.cx(0, 1)
    qc.h(0)
    qc.cx(0, 1)
    qc.t(0)
    qc.h(0)
    qc.measure(0, 0)

ts = pyqreach.TransitionSystem()
parse_result = parse_qiskit_cir(qc, qc.num_qubits, ts, return_metadata=True)

set_initial_state(ts, "000")
ts.computingFixedPointPost()

# The post-condition is the lower-bound subspace at the final location.
final_locs = parse_result.result_locations
op_final = ts.Locations[final_locs[-1]].lowerBound
print("Post-condition at final location(s):")
for loc in final_locs:
    ts.printSupp(loc)

# ============================================================================
# Step 2: Backward reachability on the same TS
# ============================================================================
print("\n" + "=" * 60)
print("Step 2: Backward reachability on the same RUS circuit")
print("=" * 60)

ts.resetLocationBounds()
ts.setAnnotation([[loc, op_final] for loc in final_locs])

time_pre_start = time()
ts.computingFixedPointPre()
time_pre_end = time()
print(f"Backward reachability (pre): {time_pre_end - time_pre_start:.3f}s")

# ============================================================================
# Step 3: Inspect weakest pre-condition
# ============================================================================
print("\n=== Weakest Pre-condition at Initial Location ===")
ts.printSupp(0)

print("\n=== Post-condition at Final Location(s) ===")
for loc in final_locs:
    ts.printSupp(loc)
