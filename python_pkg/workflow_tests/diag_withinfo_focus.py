#!/usr/bin/env python3
"""Pinpoint the exact WithInfo call where DAG corruption first occurs.

Focuses on step ~208 of n=256 |+>^n construction where the first
corruption was detected.

Usage:
    cd python_pkg
    ../.venv/bin/python workflow_tests/diag_withinfo_dag.py --focus 256 205 212
"""

from __future__ import annotations

import sys
from pathlib import Path

PYTHON_PKG = Path(__file__).resolve().parents[1]
if str(PYTHON_PKG) not in sys.path:
    sys.path.insert(0, str(PYTHON_PKG))

import pyqreach


def _check_self(op: pyqreach.QOperation, n: int) -> bool:
    ts = pyqreach.TransitionSystem()
    loc = pyqreach.Location(n)
    loc.lowerBound = op
    ts.addLocation(loc)
    return ts.Locations[0].satisfy(op)


def main():
    if len(sys.argv) >= 4:
        n = int(sys.argv[1])
        start = int(sys.argv[2])
        end = int(sys.argv[3])
    else:
        n = 256
        start = 200
        end = 212

    print(f"=== Focusing n={n}, steps {start}-{end} ===")
    print()

    pyqreach.initializeTransitionSystem()
    state = pyqreach.QOperation(["0" * n])

    # Fast-forward to just before start
    for i in range(start):
        gate = pyqreach.QOperation("H", n, [i], [])
        state = state.post_image(gate)

    print(f"  After step {start}: |+>^{start}|0>^{n-start}")
    ok = _check_self(state, n)
    print(f"    satisfy(self)={ok}")
    if not ok:
        print("    *** Already corrupted before focus window!")
        return

    # Check each remaining step
    for i in range(start, end):
        print(f"\n  >>> Step {i+1}: applying H to qubit {i} <<<", flush=True)
        gate = pyqreach.QOperation("H", n, [i], [])
        state = state.post_image(gate)
        ok = _check_self(state, n)
        marker = "  <-- CORRUPTED!" if not ok else ""
        print(f"    satisfy(self)={ok}{marker}", flush=True)
        if not ok:
            print(f"\n*** First corruption at step {i+1} ***")
            print(f"    State: |+>^{i+1}|0>^{n-i-1}")
            break

    print()
    print("Done.")


if __name__ == "__main__":
    main()
