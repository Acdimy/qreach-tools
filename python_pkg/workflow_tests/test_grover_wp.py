"""Weakest-precondition verification of Grover's search algorithm using QReach.

This is a modernized version of ``test_grover_wp.py`` that uses the current
``qctl.py`` APIs instead of the older direct-annotation / temporary-TS
workaround patterns.

The workflow:
  1. Parse a Grover circuit into a transition system.
  2. Run forward reachability (post-image fixed point).
  3. Construct the desired quantum post-condition as the span of the initial
     superposition and the target (good) state.
  4. Set the post-condition on the leaf location(s).
  5. Run backward reachability (pre-image fixed point) via
     ``computingFixedPointPre``.
  6. Inspect the weakest pre-condition at the initial location.
"""

import pyqreach
from time import time

from qiskit import QuantumCircuit

from parse_qiskit import parse_qiskit_cir
from qctl import (
    quantum_state,
    set_zero_initial_state,
    span_qops,
    annotate_leaf_operation,
)

# ---------------------------------------------------------------------------
# 1. Load the Grover-5 benchmark circuit
# ---------------------------------------------------------------------------
filename = "benchmark/grover/grover_5.qasm"
qc = QuantumCircuit.from_qasm_file(filename)
num_qubits = qc.num_qubits
print(f"Loaded {filename}: {num_qubits} qubits")

# ---------------------------------------------------------------------------
# 2. Parse into explicit transition system and compute forward reachability
# ---------------------------------------------------------------------------
ts = pyqreach.TransitionSystem()
parse_result = parse_qiskit_cir(qc, num_qubits, ts, return_metadata=True)

set_zero_initial_state(ts)

time_post_start = time()
ts.computingFixedPointPost()
time_post_end = time()
print(f"Forward reachability (post): {time_post_end - time_post_start:.3f}s")

# ---------------------------------------------------------------------------
# 3. Construct the desired post-condition (Grover target subspace)
# ---------------------------------------------------------------------------
# After post-image fixed point, the lowerBound at intermediate locations
# contains the reachable quantum states at that program point.
# For Grover-5 (5 qubits total, 3 work qubits):
#   - grover_init: the uniform superposition after initial H gates
#   - grover_good: the target marked state |11100>
work_qubits = int((num_qubits + 1) / 2)

# Extract the superposition state from the TS after the H-gate layer
grover_init = ts.Locations[work_qubits].lowerBound

# Build the target "good" state using the product-state string syntax
grover_good = quantum_state(
    "1" * work_qubits + "0" * (num_qubits - work_qubits - 1) + "1"
)

# Span the two states to create the post-condition subspace
grover_final = span_qops([grover_init, grover_good])
print(f"Post-condition subspace dimension: {grover_final.dim()}")

# ---------------------------------------------------------------------------
# 4. Set post-condition on leaf locations and compute weakest pre-condition
# ---------------------------------------------------------------------------
ts.resetLocationBounds()

# Attach the post-condition to the leaf location(s)
leaf_locs = parse_result.result_locations
annotate_leaf_operation(ts, grover_final, loc_list=leaf_locs)

time_pre_start = time()
ts.computingFixedPointPre()
time_pre_end = time()
print(f"Backward reachability (pre): {time_pre_end - time_pre_start:.3f}s")

# ---------------------------------------------------------------------------
# 5. Inspect results
# ---------------------------------------------------------------------------
print("\n=== Weakest Pre-condition at Initial Location ===")
ts.printSupp(0)

print("\n=== Post-condition at Final Location ===")
for leaf_loc in leaf_locs:
    ts.printSupp(leaf_loc)

print(f"\nTotal wall time: {time_pre_end - time_post_start:.3f}s")
