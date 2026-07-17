#!/usr/bin/env python3
"""Minimal repro: pinpoint exactly where the crash happens.

Run:  cd python_pkg && ../.venv/bin/python repro_minimal.py
"""
import sys, copy, random
from pathlib import Path

PYTHON_PKG = Path(__file__).resolve().parent
sys.path.insert(0, str(PYTHON_PKG))

import pyqreach
from qiskit import QuantumCircuit
from qiskit.circuit import CircuitInstruction
from qiskit.circuit.library import XGate, YGate, ZGate
from parse_qiskit import parse_qiskit_cir_lazy
from qctl import span_qops

def insert_random_pauli(qc, rng):
    pauli = [XGate(), YGate(), ZGate()]
    gate = rng.choice(pauli)
    pos = rng.randint(0, len(qc.data))
    qb = rng.randint(0, qc.num_qubits - 1)
    qc.data.insert(pos, CircuitInstruction(gate, [qc.qubits[qb]], []))
    return qc

qasm_path = PYTHON_PKG / "benchmark" / "dqc_pe" / "dqc_pe_2.qasm"
initial_state = "001"

pyqreach.initializeTransitionSystem()

qc = QuantumCircuit.from_qasm_file(str(qasm_path))

# reference TS
ref_ts = pyqreach.TransitionSystem()
ref_pres = parse_qiskit_cir_lazy(qc, qc.num_qubits, ref_ts,
                                 initial_state=initial_state, return_metadata=True)
ref_ts.computingFixedPointPost()
ref_vecs = [ref_ts.Locations[loc].lowerBound for loc in ref_pres.result_locations]

# injected TS
rng = random.Random(1)
qc_injected = insert_random_pauli(copy.deepcopy(qc), rng)
ts = pyqreach.TransitionSystem()
pres = parse_qiskit_cir_lazy(qc_injected, qc_injected.num_qubits, ts,
                             initial_state=initial_state, return_metadata=True)
ts.computingFixedPointPost()

all_leaf = list(pres.result_locations) + list(getattr(pres, "lazy_pruned_locations", []) or [])

# Step 1: span_qops
print("1. span_qops(ref_vecs)...", flush=True)
expected = span_qops(ref_vecs)
print("1. OK", flush=True)

# Step 2: get leaf location
print("2. get leaf location...", flush=True)
loc = ts.Locations[all_leaf[0]]
print("2. OK", flush=True)

# Step 3: satisfy
print("3. loc.satisfy(expected)...", flush=True)
loc.satisfy(expected)
print("3. OK", flush=True)
