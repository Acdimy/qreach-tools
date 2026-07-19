#!/usr/bin/env python3
"""Isolate: bit-span (rows 2 vs 10) vs opposite-sign.

Test A: same sign, rows 2+10 → |001> + |101>  (H(2) on |001>)
Test B: opposite sign, rows 2+6  → |001> - |011>  (Z(1) after creating |001>+|011>)
"""

import sys
from pathlib import Path

PYTHON_PKG = Path(__file__).resolve().parent
sys.path.insert(0, str(PYTHON_PKG))

import pyqreach
from qiskit import QuantumCircuit
from parse_qiskit import parse_qiskit_cir_lazy
from qctl import span_qops

pyqreach.initializeTransitionSystem()

# --- Test A: rows 2 and 10, SAME sign ---
# X(0): |001>, H(2): (|001> + |101>)/sqrt(2)
print("=== Test A: rows 2 and 10, SAME sign ===")
qc_a = QuantumCircuit(3)
qc_a.x(0)
qc_a.h(2)
ts_a = pyqreach.TransitionSystem()
pres_a = parse_qiskit_cir_lazy(qc_a, 3, ts_a, initial_state="000", return_metadata=True)
ts_a.computingFixedPointPost()
ref_vecs_a = [ts_a.Locations[loc].lowerBound for loc in pres_a.result_locations]
print(f"  {len(ref_vecs_a)} location(s)")

# Collect all vectors into span_qops to trigger GramSchmidt/normalize
if len(ref_vecs_a) >= 2:
    try:
        span_qops(ref_vecs_a)
        print("  span_qops OK")
    except Exception as e:
        print(f"  span_qops crashed: {e}")
elif len(ref_vecs_a) == 1:
    print("  only 1 vector, trying disjunction with self to trigger GS...")
    try:
        ref_vecs_a[0].disjunction(ref_vecs_a[0])
        print("  disjunction OK")
    except Exception as e:
        print(f"  disjunction crashed: {e}")

# --- Test B: rows 2 and 6, OPPOSITE sign ---
# Need: (|001> - |011>)/sqrt(2) = entries at rows 2 and 6 with opposite signs
# X(0): |001>, H(1): (|001> + |011>)/sqrt(2), Z(1): (|001> - |011>)/sqrt(2)
print("\n=== Test B: rows 2 and 6, OPPOSITE sign ===")
qc_b = QuantumCircuit(3)
qc_b.x(0)
qc_b.h(1)
qc_b.z(1)
ts_b = pyqreach.TransitionSystem()
pres_b = parse_qiskit_cir_lazy(qc_b, 3, ts_b, initial_state="000", return_metadata=True)
ts_b.computingFixedPointPost()
ref_vecs_b = [ts_b.Locations[loc].lowerBound for loc in pres_b.result_locations]
print(f"  {len(ref_vecs_b)} location(s)")

if len(ref_vecs_b) >= 2:
    try:
        span_qops(ref_vecs_b)
        print("  span_qops OK")
    except Exception as e:
        print(f"  span_qops crashed: {e}")
else:
    print("  only 1 vector")
