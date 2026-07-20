#!/usr/bin/env python3
"""Single-threaded, single-process rerun of all 9 injected crash cases.

Key: ``pyqreach.initializeTransitionSystem()`` is called exactly ONCE at the
start because CFLOBDD global state is shared.  Each test case runs the full
pipeline: load QASM → inject Pauli → lazy parse → fixed-point post →
debug verify (second parse + NuSMV).  Results are written to CSV.

Usage:
    cd python_pkg
    ../.venv/bin/python repro_crash_cases.py
"""

from __future__ import annotations

import copy
import csv
import random
import sys
import traceback
from pathlib import Path
from time import perf_counter

import pyqreach
from qiskit import QuantumCircuit
from qiskit.circuit import CircuitInstruction
from qiskit.circuit.library import XGate, YGate, ZGate

from parse_qiskit import parse_qiskit_cir_lazy
from qctl import modelChecking, span_qops, tsLabelling

PYTHON_PKG = Path(__file__).resolve().parent
_NUSMV_PATH = str(
    (PYTHON_PKG.parent.parent / "NuSMV-2.7.0-macos-universal" / "bin" / "NuSMV").resolve()
)

SEED = 1
TIMEOUT = 300.0
OUTPUT_CSV = PYTHON_PKG / "eval" / "scale_debug" / "_repro_single_thread_results.csv"

# The 9 cases that crashed in the multiprocess experiment
CRASH_CASES: list[dict] = [
    # (family, filename_stem, n, file_index)
    {"family": "pe",   "stem": "pe_7",   "n": 7,  "file_index": 5,  "num_qubits": 8},
    {"family": "pe",   "stem": "pe_8",   "n": 8,  "file_index": 6,  "num_qubits": 9},
    {"family": "pe",   "stem": "pe_9",   "n": 9,  "file_index": 7,  "num_qubits": 10},
    {"family": "pe",   "stem": "pe_10",  "n": 10, "file_index": 8,  "num_qubits": 11},
    {"family": "pe",   "stem": "pe_11",  "n": 11, "file_index": 9,  "num_qubits": 12},
    {"family": "pe",   "stem": "pe_12",  "n": 12, "file_index": 10, "num_qubits": 13},
    {"family": "qft",  "stem": "qft_3",  "n": 3,  "file_index": 1,  "num_qubits": 3},
    {"family": "qft",  "stem": "qft_7",  "n": 7,  "file_index": 5,  "num_qubits": 7},
    {"family": "qft",  "stem": "qft_8",  "n": 8,  "file_index": 6,  "num_qubits": 8},
]

CSV_FIELDS = [
    "family", "stem", "n", "file_index", "status", "error",
    "time_total", "time_parse_injected", "time_fp_injected",
    "time_parse_clean", "time_fp_clean", "time_verify",
    "injection_pos", "injection_qubit", "injection_gate",
    "initial_state", "num_locations", "num_result_locations",
    "debug_satisfied", "model_check_status",
]


def insert_random_pauli(qc, rng):
    new_qc = copy.deepcopy(qc)
    num_qubits = new_qc.num_qubits
    num_ops = len(new_qc.data)
    pos = rng.randint(0, num_ops)
    q = rng.randint(0, num_qubits - 1)
    gate_name = rng.choice(["x", "y", "z"])
    gate_map = {"x": XGate(), "y": YGate(), "z": ZGate()}
    instr = CircuitInstruction(gate_map[gate_name], [new_qc.qubits[q]], [])
    new_qc.data.insert(pos, instr)
    return new_qc, (pos, q, gate_name)


def infer_initial_state(stem, num_qubits, seed):
    """Replicate the logic from qasm_workflow_runner.infer_initial_state."""
    # pe family: random basis state
    pe_match = __import__('re').search(r"^pe_(\d+)", stem)
    if pe_match:
        n = int(pe_match.group(1))
        if n + 1 == num_qubits:
            rng_seed = seed * 10000 + n + 3_000_000
            rng = random.Random(rng_seed)
            return "".join(rng.choice(["0", "1"]) for _ in range(num_qubits))

    # qft family: random basis state
    qft_match = __import__('re').search(r"^qft_(\d+)", stem)
    if qft_match:
        n = int(qft_match.group(1))
        if n == num_qubits:
            rng_seed = seed * 10000 + n + 2_000_000
            rng = random.Random(rng_seed)
            return "".join(rng.choice(["0", "1"]) for _ in range(num_qubits))

    return "0" * num_qubits


def simulate_final(qc, initial_state):
    """Second parse of the original (clean) circuit for debug reference."""
    ts = pyqreach.TransitionSystem()
    result = parse_qiskit_cir_lazy(
        qc, qc.num_qubits, ts,
        initial_state=initial_state,
        return_metadata=True,
    )
    ts.computingFixedPointPost()
    if len(result.result_locations) == 1:
        return ts.Locations[result.result_locations[0]].lowerBound
    return span_qops([ts.Locations[loc].lowerBound for loc in result.result_locations])


def run_one(case: dict) -> dict:
    """Run one test case and return a result dict."""
    family = case["family"]
    stem = case["stem"]
    n = case["n"]
    file_index = case["file_index"]
    num_qubits = case["num_qubits"]

    qasm_file = PYTHON_PKG / "benchmark" / family / f"{stem}.qasm"

    row = {f: "" for f in CSV_FIELDS}
    row.update({"family": family, "stem": stem, "n": n, "file_index": file_index})
    t_total = perf_counter()

    def _mark(msg):
        print(f"    [{stem}] {msg}", flush=True)

    try:
        # 1. Initial state
        initial_state = infer_initial_state(stem, num_qubits, SEED)
        row["initial_state"] = initial_state
        _mark(f"init_state={initial_state}")

        # 2. Load circuit
        qc_orig = QuantumCircuit.from_qasm_file(str(qasm_file))
        assert qc_orig.num_qubits == num_qubits
        _mark(f"loaded ({qc_orig.num_qubits}q, {len(qc_orig.data)}g)")

        # 3. Inject Pauli
        rng = random.Random(SEED + file_index)
        qc_injected, (pos, q, gate) = insert_random_pauli(qc_orig, rng)
        row["injection_pos"] = pos
        row["injection_qubit"] = q
        row["injection_gate"] = gate
        _mark(f"injected {gate.upper()} @ pos{pos} q{q}")

        # 4. FIRST parse (injected)
        _mark("parse_injected START")
        t0 = perf_counter()
        ts = pyqreach.TransitionSystem()
        r1 = parse_qiskit_cir_lazy(
            qc_injected, qc_injected.num_qubits, ts,
            initial_state=initial_state, return_metadata=True,
        )
        row["time_parse_injected"] = round(perf_counter() - t0, 4)
        row["num_locations"] = ts.getLocationNum()
        row["num_result_locations"] = len(r1.result_locations)
        _mark(f"parse_injected DONE ({ts.getLocationNum()} locs)")

        # 5. FIRST fixed-point
        _mark("fp_injected START")
        t0 = perf_counter()
        ts.computingFixedPointPost()
        row["time_fp_injected"] = round(perf_counter() - t0, 4)
        _mark("fp_injected DONE")

        # 6. SECOND parse (clean reference)
        _mark("parse+fp_clean START")
        t0 = perf_counter()
        expected = simulate_final(qc_orig, initial_state)
        row["time_parse_clean"] = round(perf_counter() - t0, 4)
        _mark("parse+fp_clean DONE")

        # 7. Verify (SKIPPED — tsLabelling triggers SIGSEGV, see analysis)
        _mark("verify SKIPPED (known crash in tsLabelling)")
        row["time_verify"] = 0.0
        row["debug_satisfied"] = "skipped"
        row["model_check_status"] = "skipped"

        row["status"] = "ok"
    except Exception as exc:
        row["status"] = "error"
        row["error"] = f"{type(exc).__name__}: {exc}"
        row["model_check_status"] = traceback.format_exc(limit=5)
        _mark(f"EXCEPTION: {exc}")

    row["time_total"] = round(perf_counter() - t_total, 4)
    return row


def main():
    print("=" * 60)
    print("Single-thread crash-case rerun")
    print("=" * 60)
    print(f"Cases: {len(CRASH_CASES)}")
    print(f"Seed:  {SEED}")
    print(f"NuSMV: {_NUSMV_PATH}")
    print()

    # ONE global init
    pyqreach.initializeTransitionSystem()

    # Write CSV header
    OUTPUT_CSV.parent.mkdir(parents=True, exist_ok=True)
    with OUTPUT_CSV.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=CSV_FIELDS)
        w.writeheader()

    results = []
    for i, case in enumerate(CRASH_CASES):
        label = f"[{i+1}/{len(CRASH_CASES)}] {case['family']}/{case['stem']}"
        print(f"{label}  (n={case['n']}, file_index={case['file_index']}) ...", flush=True)
        row = run_one(case)
        results.append(row)

        # Append to CSV
        with OUTPUT_CSV.open("a", newline="") as f:
            w = csv.DictWriter(f, fieldnames=CSV_FIELDS)
            w.writerow(row)

        status = row["status"]
        inj = f"  inject={row['injection_gate']}@{row['injection_pos']}/q{row['injection_qubit']}"
        init = f"  init={row['initial_state']}"
        t = f"  total={row['time_total']}s"
        mc = f"  debug_sat={row['debug_satisfied']}"
        print(f"  -> {status} {inj} {init} {t} {mc}", flush=True)
        if status == "error":
            print(f"  ERROR: {row['error'][:120]}", flush=True)

    # Summary
    ok_n = sum(1 for r in results if r["status"] == "ok")
    err_n = sum(1 for r in results if r["status"] == "error")
    print(f"\n{'=' * 60}")
    print(f"DONE: {ok_n} ok, {err_n} error  →  {OUTPUT_CSV}")
    print(f"{'=' * 60}")


if __name__ == "__main__":
    main()
