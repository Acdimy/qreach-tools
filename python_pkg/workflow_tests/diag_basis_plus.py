#!/usr/bin/env python3
"""Test basis_plus construction and single-it-grover debug check."""
import sys, os
from pathlib import Path
PYTHON_PKG = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PYTHON_PKG))

import pyqreach
from qctl import quantum_state, span_qops

# Test basis_plus construction
for n in [15, 31, 64, 128]:
    half = n // 2
    pyqreach.initializeTransitionSystem()
    try:
        z = quantum_state("0" * n)
        p = quantum_state("0" * half + "+" * (n - half))
        sp = span_qops([z, p])
        print(f"n={n} half={half}: dim={sp.dim()}")
    except Exception as e:
        print(f"n={n} half={half}: ERROR: {e}")
