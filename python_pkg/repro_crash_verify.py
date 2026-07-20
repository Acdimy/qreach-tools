#!/usr/bin/env python3
"""Rerun all 9 crash cases WITH verify, single-process, 300s timeout.

Usage:  cd python_pkg && ../.venv/bin/python repro_crash_verify.py
"""
import copy, csv, random, sys, time, traceback
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import pyqreach
from qiskit import QuantumCircuit
from qiskit.circuit import CircuitInstruction
from qiskit.circuit.library import XGate, YGate, ZGate
from parse_qiskit import parse_qiskit_cir_lazy
from qctl import modelChecking, span_qops, tsLabelling

PY = Path(__file__).resolve().parent
NUSMV = str((PY.parent.parent / "NuSMV-2.7.0-macos-universal" / "bin" / "NuSMV").resolve())
SEED = 1
TIMEOUT = 300.0

CASES = [
    ("pe", "pe_7", 7, 5, 8),
    ("pe", "pe_8", 8, 6, 9),
    ("pe", "pe_9", 9, 7, 10),
    ("pe", "pe_10", 10, 8, 11),
    ("pe", "pe_11", 11, 9, 12),
    ("pe", "pe_12", 12, 10, 13),
    ("qft", "qft_3", 3, 1, 3),
    ("qft", "qft_7", 7, 5, 7),
    ("qft", "qft_8", 8, 6, 8),
]

print("=" * 60)
print("Crash-case verify rerun  (timeout=300s)")
print("=" * 60)

pyqreach.initializeTransitionSystem()

results = []
for idx, (family, stem, n, file_index, num_qubits) in enumerate(CASES):
    label = f"[{idx+1}/{len(CASES)}] {stem}"
    print(f"\n{label}  n={n}  idx={file_index} ...", flush=True)
    t_start = time.perf_counter()
    qasm = PY / "benchmark" / family / f"{stem}.qasm"

    try:
        # Init state
        import re
        if family == "pe":
            rng_seed = SEED * 10000 + n + 3_000_000
        else:
            rng_seed = SEED * 10000 + n + 2_000_000
        init = "".join(random.Random(rng_seed).choice(["0", "1"]) for _ in range(num_qubits))

        # Load
        qc_orig = QuantumCircuit.from_qasm_file(str(qasm))

        # Inject
        rng = random.Random(SEED + file_index)
        qc_inj = copy.deepcopy(qc_orig)
        pos = rng.randint(0, len(qc_inj.data))
        q = rng.randint(0, qc_inj.num_qubits - 1)
        gate_name = rng.choice(["x", "y", "z"])
        gate_map = {"x": XGate(), "y": YGate(), "z": ZGate()}
        qc_inj.data.insert(pos, CircuitInstruction(gate_map[gate_name], [qc_inj.qubits[q]], []))
        print(f"  init={init}  inject={gate_name}@pos{pos}/q{q}", flush=True)

        # Parse injected
        print(f"  [1] parse+fp injected ...", flush=True)
        ts = pyqreach.TransitionSystem()
        r = parse_qiskit_cir_lazy(qc_inj, qc_inj.num_qubits, ts, initial_state=init, return_metadata=True)
        ts.computingFixedPointPost()
        print(f"  [1] done ({ts.getLocationNum()} locs, {time.perf_counter()-t_start:.1f}s)", flush=True)

        # Parse clean
        print(f"  [2] parse+fp clean ...", flush=True)
        ts2 = pyqreach.TransitionSystem()
        r2 = parse_qiskit_cir_lazy(qc_orig, qc_orig.num_qubits, ts2, initial_state=init, return_metadata=True)
        ts2.computingFixedPointPost()
        expected = ts2.Locations[r2.result_locations[0]].lowerBound
        print(f"  [2] done ({time.perf_counter()-t_start:.1f}s)", flush=True)

        # tsLabelling
        print(f"  [3] tsLabelling ...", flush=True)
        for loc in r.result_locations:
            ts.setLabel(loc, "final")
        tsLabelling(ts, expected, "debug_op", locList=list(r.result_locations))
        print(f"  [3] done", flush=True)

        # NuSMV
        print(f"  [4] NuSMV ...", flush=True)
        mc = modelChecking(ts, "AG (debug_op <-> final)", nusmv_path=NUSMV)
        sat = mc.get("satisfied")
        print(f"  [4] done  sat={sat}  ({time.perf_counter()-t_start:.1f}s)", flush=True)

        results.append((stem, "ok", sat, round(time.perf_counter() - t_start, 2)))
    except Exception as e:
        results.append((stem, f"error: {e}", None, round(time.perf_counter() - t_start, 2)))
        print(f"  EXCEPTION: {e}", flush=True)

# Summary
print(f"\n{'=' * 60}")
print(f"{'case':10s} {'status':20s} {'sat':8s} {'time':>8s}")
print("-" * 50)
ok = err = 0
for stem, status, sat, t in results:
    print(f"{stem:10s} {str(status):20s} {str(sat):8s} {str(t):>7s}s")
    if status == "ok":
        ok += 1
    else:
        err += 1
print(f"\n{ok} ok, {err} error")
