#!/usr/bin/env python3
"""Minimal reproduction: post_image DAG corruption at CFLOBDD level >= 8.

Building |+>^n via sequential post_image (as the lazy parser does for
gate applications) produces a CFLOBDD state that fails self-consistency
checks (satisfy(self) == False) for n >= 64.

At small n (<= 16), the same state built via post_image or directly via
QOperation(['+'*n]) yields identical, self-consistent results — confirming
that post_image itself is the source of corruption at high levels.

Usage:
    cd python_pkg
    ../.venv/bin/python workflow_tests/repro_postimage_corruption.py
"""

from __future__ import annotations

import sys
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


def _build_via_post_image(n: int) -> list[pyqreach.QOperation]:
    """Build |+>^n via sequential post_image(H), keeping all intermediates."""
    states = [pyqreach.QOperation(["0" * n])]
    for i in range(n):
        gate = pyqreach.QOperation("H", n, [i], [])
        states.append(states[-1].post_image(gate))
    return states


def main():
    print("=" * 70)
    print("post_image DAG Corruption — Minimal Reproduction")
    print("=" * 70)

    # --- Small n: verify post_image produces correct |+>^n ---
    print()
    print("--- Small n: verify post_image chain correctness ---")
    print("  (direct QOperation(['+'*n]) is known-correct for n <= 16)")
    print()
    for n in [2, 4, 8, 16]:
        pyqreach.initializeTransitionSystem()
        chained = _build_via_post_image(n)
        direct = pyqreach.QOperation(["+" * n])

        # Both should be the same state: span(chained, direct) should be 1D
        span = pyqreach.span_qops([direct, chained[-1]])
        chain_ok = _check_self(chained[-1], n)
        direct_ok = _check_self(direct, n)
        print(f"  n={n:>4}: dim(span)={span.dim()}  chain_self={chain_ok}  direct_self={direct_ok}")

    # --- Large n: post_image chain breaks self-consistency ---
    print()
    print("--- Large n: post_image chain self-consistency ---")
    print("  A valid state must satisfy itself. False means DAG corrupted.")
    print()
    all_ok = True
    for n in [32, 64, 96, 128, 256]:
        pyqreach.initializeTransitionSystem()
        chained = _build_via_post_image(n)

        # Check self-consistency at final step
        chain_ok = _check_self(chained[-1], n)

        status = "OK" if chain_ok else "CORRUPTED"
        if not chain_ok:
            all_ok = False
        print(f"  n={n:>4}: satisfy(self) = {status}")

    print()
    print("=" * 70)
    if all_ok:
        print("All sizes pass — no corruption detected.")
    else:
        print("post_image chain produces corrupted DAG at level >= 8.")
        print()
        print("Root cause: MatrixMultiplyV4 (gate application / post_image)")
        print("corrupts CFLOBDD internal return-map routing for certain DAG")
        print("topologies at high CFLOBDD levels.  The H×content trick in")
        print("dot()/normalize() mitigates Gram-Schmidt but does NOT protect")
        print("the post_image / gate-application path.")
    print("=" * 70)


if __name__ == "__main__":
    main()
