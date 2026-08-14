"""Weakest-precondition verification of an RUS (Repeat-Until-Success) circuit.

This is a modernized version of ``test_rus_wp.py`` using the current
``qctl.py`` APIs.

.. note::

   This test currently triggers a CFLOBDD assertion in
   ``SingleVecTerm::normalize()`` (``resMap.Size() <= 2``) during the
   pre-image computation.  The same assertion is present in the older
   ``test_rus_wp.py`` — it is a pre-existing limitation of the CFLOBDD
   normalization path when dealing with multi-amplitude leaf maps that
   arise from measurement-and-loop control flow during backward
   reachability.  This is not a regression introduced by the API
   modernization.

The workflow:
  1. Parse the original RUS circuit and compute forward reachability to
     obtain the quantum post-condition at the final location.
  2. Build a *modified* circuit that omits the ``reset`` instruction inside
     the while-loop body (the ``reset`` is excluded because it breaks the
     unitary pre-image semantics).
  3. Transfer the post-condition from step 1 onto the modified circuit's
     leaf locations.
  4. Run ``computingFixedPointPre`` on the modified circuit to compute the
     weakest pre-condition — the set of initial states that guarantee the
     program terminates in the post-condition subspace.
"""

import pyqreach
from time import time

from qiskit import QuantumCircuit, QuantumRegister, ClassicalRegister

from qreach.parse_qiskit import parse_qiskit_cir
from qreach.qctl import (
    set_initial_state,
)


# ============================================================================
# Step 1: Original RUS circuit — compute post-condition
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

# The post-condition is the lower-bound subspace at the final location
op_final = ts.Locations[parse_result.result_locations[-1]].lowerBound
print("Post-condition at final location(s):")
for loc in parse_result.result_locations:
    ts.printSupp(loc)

# ============================================================================
# Step 2: Modified circuit (no reset) + weakest pre-condition
# ============================================================================
print("\n" + "=" * 60)
print("Step 2: Backward reachability on modified circuit (no reset)")
print("=" * 60)

# The reset(0) is removed because it is non-unitary and would break the
# pre-image semantics.  When removed, the while-loop body is purely unitary
# inside a measurement-controlled loop, and pre-image computation is valid.
qc_check = QuantumCircuit(3, 1)
qc_check.x(2)
qc_check.measure(2, 0)
with qc_check.while_loop((0, 0b1)):
    qc_check.h(0)
    qc_check.t(0)
    qc_check.cx(0, 1)
    qc_check.h(0)
    qc_check.cx(0, 1)
    qc_check.t(0)
    qc_check.h(0)
    qc_check.measure(0, 0)

ts_check = pyqreach.TransitionSystem()
parse_result_check = parse_qiskit_cir(
    qc_check, qc_check.num_qubits, ts_check, return_metadata=True
)

# Transfer the post-condition from Step 1 onto the leaf location(s)
final_locs = parse_result_check.result_locations
ts_check.setAnnotation([[loc, op_final] for loc in final_locs])

time_pre_start = time()
ts_check.computingFixedPointPre()
time_pre_end = time()
print(f"Backward reachability (pre): {time_pre_end - time_pre_start:.3f}s")

# ============================================================================
# Step 3: Inspect weakest pre-condition
# ============================================================================
print("\n=== Weakest Pre-condition at Initial Location ===")
ts_check.printSupp(0)

print("\n=== Post-condition at Final Location(s) ===")
for loc in final_locs:
    ts_check.printSupp(loc)
