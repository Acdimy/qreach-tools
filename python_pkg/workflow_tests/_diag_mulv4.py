#!/usr/bin/env python3
"""Trigger MatrixMultiplyV4 computation with corrupted |+>^64 and print diag."""
import sys
from pathlib import Path
PYTHON_PKG = Path('/Users/ftdac/thu/qreach-tools/python_pkg')
sys.path.insert(0, str(PYTHON_PKG))
import pyqreach

def check_self(op, n):
    ts = pyqreach.TransitionSystem()
    loc = pyqreach.Location(n)
    loc.lowerBound = op
    ts.addLocation(loc)
    return ts.Locations[0].satisfy(op)

n = 64
pyqreach.initializeTransitionSystem()

# Build corrupted |+>^64
state = pyqreach.QOperation(['0' * n])
for i in range(n):
    state = state.post_image(pyqreach.QOperation('H', n, [i], []))

print("=== check_self(corrupted |+>^64) ===", flush=True)
result = check_self(state, n)
print(f"Result: {result}", flush=True)
