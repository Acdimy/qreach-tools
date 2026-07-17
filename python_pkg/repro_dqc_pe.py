#!/usr/bin/env python3
"""Minimal repro: dqc_pe_2 with Pauli injection + comparison debug.
Crashes in span_qops inside _simulate_final_operation.
"""
import sys, copy, random
from pathlib import Path

PYTHON_PKG = Path(__file__).resolve().parent
sys.path.insert(0, str(PYTHON_PKG))

import pyqreach
from qiskit import QuantumCircuit
from qiskit.circuit.library import XGate, YGate, ZGate
from parse_qiskit import parse_qiskit_cir_lazy
from qctl import span_qops, modelChecking, tsLabelling

# ---------------------------------------------------------------------------
def insert_random_pauli(qc, rng):
    """Insert one random Pauli gate at a random position."""
    pauli = [XGate(), YGate(), ZGate()]
    gate = rng.choice(pauli)
    pos = rng.randint(0, len(qc.data))
    qb = rng.randint(0, qc.num_qubits - 1)
    from qiskit.circuit import CircuitInstruction
    qc.data.insert(pos, CircuitInstruction(gate, [qc.qubits[qb]], []))
    return qc, {"position": pos, "qubit": qb, "gate": gate.name}

# ---------------------------------------------------------------------------
qasm_path = PYTHON_PKG / "benchmark" / "dqc_pe" / "dqc_pe_2.qasm"
initial_state = "001"  # correct: "0"*2 + "1"

pyqreach.initializeTransitionSystem()

qc = QuantumCircuit.from_qasm_file(str(qasm_path))
print(f"[1] Loaded: {qc.num_qubits} qubits, {len(qc.data)} gates")

# Inject Pauli
rng = random.Random(1)
original_qc = copy.deepcopy(qc)
qc_injected, info = insert_random_pauli(qc, rng)
print(f"[2] Injected {info['gate']} at pos={info['position']} qubit={info['qubit']}")
print(f"    injected circuit: {len(qc_injected.data)} gates")

# Build injected TS (this works)
ts = pyqreach.TransitionSystem()
pres = parse_qiskit_cir_lazy(qc_injected, qc_injected.num_qubits, ts,
                             initial_state=initial_state, return_metadata=True)
ts.computingFixedPointPost()
print(f"[3] Injected TS: {ts.getLocationNum()} locs, "
      f"{len(pres.result_locations)} result_locs")

# Build reference TS from original circuit → span_qops CRASHES HERE
print(f"[4] Building reference TS (original circuit)...")
ref_ts = pyqreach.TransitionSystem()
ref_pres = parse_qiskit_cir_lazy(original_qc, original_qc.num_qubits, ref_ts,
                                 initial_state=initial_state, return_metadata=True)
ref_ts.computingFixedPointPost()
print(f"    Ref TS: {ref_ts.getLocationNum()} locs, "
      f"{len(ref_pres.result_locations)} result_locs")

# Show each result location
for loc in ref_pres.result_locations:
    print(f"    loc {loc}")
    # ref_ts.printSupp(loc)

print(f"[5] span_qops(...)")
expected = span_qops([ref_ts.Locations[loc].lowerBound
                      for loc in ref_pres.result_locations])
print(f"    expected type: {type(expected).__name__}")
expected.printFormal()

# --- apply label & model-check (assertion triggers here) ---
print(f"[6] tsLabelling + modelChecking — THIS MAY CRASH:")
_rl = list(pres.result_locations)
_pl = list(getattr(pres, "lazy_pruned_locations", []) or [])
all_leaf = _rl + _pl

ts.printSupp(all_leaf[0])
ts.Locations[all_leaf[0]].satisfy(expected)

# for loc in all_leaf:
#     ts.setLabel(loc, "final")
# tsLabelling(ts, expected, "debug_op", locList=all_leaf)

# _NUSMV = str((PYTHON_PKG.parent.parent / "NuSMV-2.7.0-macos-universal" / "bin" / "NuSMV").resolve())
# result = modelChecking(ts, "AG (debug_op <-> final)", nusmv_path=_NUSMV)
# print(f"    satisfied={result.get('satisfied')}, status={result.get('status', 'ok')}")
