#!/usr/bin/env python3
"""Diagnose MatrixMultiplyV4WithInfo DAG structure at levels 7-9.

Builds |+>^n via sequential post_image(H) step by step,
checking self-consistency at each step to find where the
DAG first becomes corrupted.

The C++ diagnostic output (from matrix1234_complex_float_boost_top_node.cpp)
prints [WithInfo L=X] lines to stderr for each gate application at level >= 7.

Usage:
    cd python_pkg
    ../.venv/bin/python workflow_tests/diag_withinfo_dag.py [n]
"""

from __future__ import annotations

import sys
import os
from pathlib import Path

PYTHON_PKG = Path(__file__).resolve().parents[1]
if str(PYTHON_PKG) not in sys.path:
    sys.path.insert(0, str(PYTHON_PKG))

import pyqreach


def _check_self(op: pyqreach.QOperation, n: int) -> bool:
    """Check whether a QOperation satisfies itself (self-consistency)."""
    ts = pyqreach.TransitionSystem()
    loc = pyqreach.Location(n)
    loc.lowerBound = op
    ts.addLocation(loc)
    return ts.Locations[0].satisfy(op)


def diagnose(n: int, check_every: int = 1):
    """Build |+>^n via post_image, checking self-consistency at each step."""
    print(f"=== Diagnosing n={n} (check every {check_every} H gates) ===")
    print()

    pyqreach.initializeTransitionSystem()

    # Start with |0>^n
    states = [pyqreach.QOperation(["0" * n])]
    ok_start = _check_self(states[0], n)
    print(f"  step  0: |0>^{n}  satisfy(self)={ok_start}")

    # Apply H gates one by one
    first_bad = -1
    for i in range(n):
        gate = pyqreach.QOperation("H", n, [i], [])
        states.append(states[-1].post_image(gate))

        if (i + 1) % check_every == 0 or i == n - 1:
            ok = _check_self(states[-1], n)
            marker = "  <-- CORRUPTED" if not ok else ""
            print(f"  step {i+1:>3}: H(q{i}) satisfy(self)={ok}{marker}")
            sys.stdout.flush()
            if not ok and first_bad < 0:
                first_bad = i + 1

    print()
    if first_bad >= 0:
        print(f"*** First corruption at step {first_bad} (H on qubit {first_bad-1}) ***")
        print(f"    State: |+>^{first_bad}|0>^{n-first_bad}")
    else:
        print("All steps passed — no corruption detected.")
    print()
    return first_bad


def main():
    if len(sys.argv) > 1:
        ns = [int(sys.argv[1])]
    else:
        ns = [64, 128, 256]

    for n in ns:
        diagnose(n, check_every=max(1, n // 16))


if __name__ == "__main__":
    main()
