#!/usr/bin/env python3
"""Minimal reproduction of CFLOBDD MatrixTranspose bug at level 8+ (qNum >= 128).

At CFLOBDD level >= 8, MatrixTranspose destroys row 0 of the transposed vector,
causing dot() to return 0 for inner products that should be non-zero.

The bug manifests in Gram-Schmidt via span_qops: even when |0>^n is one of the
spanning vectors, the satisfy check fails because dot() returns 0 for the
projection coefficients.

Usage:
    cd python_pkg
    ../.venv/bin/python workflow_tests/repro_transpose_bug.py
"""

from __future__ import annotations

import sys
from pathlib import Path

PYTHON_PKG = Path(__file__).resolve().parents[1]
if str(PYTHON_PKG) not in sys.path:
    sys.path.insert(0, str(PYTHON_PKG))

import pyqreach
from qiskit import QuantumCircuit
from parse_qiskit import parse_qiskit_cir_lazy


def make_state(n_qubits: int, plus_indices: list[int] | None = None):
    """Build |0>^n with H applied to specified qubits."""
    s = pyqreach.QOperation(["0" * n_qubits])
    if plus_indices:
        for i in plus_indices:
            s = s.post_image(pyqreach.QOperation("H", n_qubits, [i], []))
    return s


def _get_zero_location(n: int):
    """Get a Location whose lowerBound = upperBound = |0>^n.

    Uses a simple X-gate circuit: loc 0 = |0>^n, loc 1 = |1>|0>^(n-1).
    Returns (TransitionSystem, location_id).
    """
    qc = QuantumCircuit(n)
    qc.x(0)
    ts = pyqreach.TransitionSystem()
    parse_qiskit_cir_lazy(qc, n, ts, initial_state="0" * n, return_metadata=True)
    return ts, 0  # location 0 = |0>^n (before X gate)


def test_self_consistency(n: int) -> bool:
    """Check: does |0>^n satisfy span(|0>^n, |0>^n)?"""
    pyqreach.initializeTransitionSystem()
    z1 = make_state(n)
    z2 = make_state(n)
    target = pyqreach.span_qops([z1, z2])
    ts, loc = _get_zero_location(n)
    return ts.Locations[loc].satisfy(target)


def test_partial_plus(n: int, k: int) -> bool:
    """Check: does |0>^n satisfy span(|0>^n, |+>^k|0>^(n-k))?"""
    pyqreach.initializeTransitionSystem()
    z = make_state(n)
    p = make_state(n, list(range(k)))
    target = pyqreach.span_qops([z, p])
    ts, loc = _get_zero_location(n)
    return ts.Locations[loc].satisfy(target)


def main():
    print("=" * 70)
    print("CFLOBDD MatrixTranspose Bug — Minimal Reproduction")
    print("=" * 70)
    print()
    print("The bug: at CFLOBDD level >= 8, MatrixTranspose destroys row 0.")
    print("dot(|0>, |0>) returns 0 instead of 1.0.")
    print("Gram-Schmidt projection coefficients become 0 instead of correct values.")
    print()

    # --- Test 1: self-consistency across levels ---
    print("--- Test 1: |0>^n ∈ span(|0>^n, |0>^n) ---")
    print("  (Two identical vectors → Gram-Schmidt → 1D span)")
    print("  (Should ALWAYS be True — it's one of the spanning vectors)")
    print()
    all_ok = True
    for n in [32, 64, 128, 256]:
        ok = test_self_consistency(n)
        status = "✓" if ok else "✗ TRANSPOSE BUG"
        if not ok:
            all_ok = False
        print(f"  n={n:>4} (qNum={n}, level={n.bit_length()}): {status}")
    print()

    # --- Test 2: partial |+⟩ states ---
    print("--- Test 2: |0>^n ∈ span(|0>^n, |+>^k|0>^(n-k)) ---")
    print("  (|0>^n IS the first spanning vector — must be True)")
    print()
    for n in [128, 256]:
        print(f"  n={n}:")
        for k in [1, 16, 32, 64, 128]:
            if k > n:
                break
            ok = test_partial_plus(n, k)
            status = "✓" if ok else "✗"
            if not ok:
                all_ok = False
            print(f"    k={k:>4}: {status}")
    print()

    # --- Summary ---
    print("=" * 70)
    if all_ok:
        print("All tests PASSED — transpose bug not triggered at these sizes.")
    else:
        print("Some tests FAILED — transpose bug is active.")
        print()
        print("Root cause: MatrixTranspose at CFLOBDD level >= 8 destroys row 0.")
        print("The [0,0] entry of the transposed matrix is 0 even via")
        print("EvaluateIteratively.  No threshold adjustment can fix this — it is")
        print("a structural DAG corruption in the transpose operation.")
        print()
        print("Workaround: use the Identity-multiply trick (I × content) in dot(),")
        print("the same way normalize() already does.")
    print("=" * 70)


if __name__ == "__main__":
    main()
