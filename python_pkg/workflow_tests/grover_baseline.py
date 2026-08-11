#!/usr/bin/env python3
"""Quick grover benchmark baseline: which circuits pass, timeout, or produce wrong results.

Usage:
    cd python_pkg
    ../.venv/bin/python workflow_tests/grover_baseline.py
"""

from __future__ import annotations

import sys
import time
import signal
from pathlib import Path

PYTHON_PKG = Path(__file__).resolve().parents[1]
if str(PYTHON_PKG) not in sys.path:
    sys.path.insert(0, str(PYTHON_PKG))

from parse_qiskit import parse_qiskit_cir_lazy
from qiskit import QuantumCircuit
import pyqreach


class TimeoutError(Exception):
    pass


def _handler(signum, frame):
    raise TimeoutError()


GROVER_DIR = PYTHON_PKG / "benchmark" / "grover"
TIMEOUT_SEC = 300  # 5 minutes per circuit


def check_one(qasm_path: Path) -> dict:
    """Run one grover circuit and return results."""
    result = {
        "name": qasm_path.stem,
        "n_qubits": 0,
        "status": "unknown",
        "time_s": 0.0,
        "satisfied": None,
        "error": None,
    }

    qc = QuantumCircuit.from_qasm_file(str(qasm_path))
    n = qc.num_qubits
    result["n_qubits"] = n

    # Determine search width from filename
    import re
    m = re.search(r"grover_(\d+)", qasm_path.stem)
    work_qubits = int(m.group(1)) if m else n

    try:
        pyqreach.initializeTransitionSystem()
        ts = pyqreach.TransitionSystem()

        t0 = time.perf_counter()
        signal.signal(signal.SIGALRM, _handler)
        signal.alarm(TIMEOUT_SEC)

        # Initial state: |0>^(n-1)|1>
        init = "0" * (n - 1) + "1"
        result_locs = parse_qiskit_cir_lazy(qc, n, ts, initial_state=init, return_metadata=False)

        signal.alarm(0)
        elapsed = time.perf_counter() - t0
        result["time_s"] = elapsed

        # Check: final state should be in span(|init>, |good>)
        # where |good> = |1>^work |0>^(n-work-1)|1>
        grover_init = ts.Locations[work_qubits].lowerBound
        grover_good = pyqreach.QOperation(["1" * work_qubits + "0" * (n - work_qubits - 1) + "1"])
        grover_final = pyqreach.span_qops([grover_init, grover_good])

        final_loc = result_locs[-1] if isinstance(result_locs, list) else result_locs
        if isinstance(final_loc, list):
            final_loc = final_loc[0]

        satisfied = ts.Locations[final_loc].satisfy(grover_final)
        result["satisfied"] = satisfied
        result["status"] = "PASS" if satisfied else "WRONG"

    except TimeoutError:
        signal.alarm(0)
        elapsed = TIMEOUT_SEC
        result["time_s"] = elapsed
        result["status"] = "TIMEOUT"
    except Exception as e:
        signal.alarm(0)
        result["status"] = "ERROR"
        result["error"] = str(e)[:200]

    return result


def main():
    files = sorted(GROVER_DIR.glob("grover_*.qasm"), key=lambda p: int(p.stem.split("_")[1]))
    print(f"{'Circuit':<24} {'Qubits':>6} {'Status':>10} {'Time(s)':>8} {'Satisfied':>10}")
    print("-" * 62)

    results = []
    for f in files:
        print(f"  {f.stem:<22} ...", end=" ", flush=True)
        r = check_one(f)
        results.append(r)
        print(f"{r['n_qubits']:>6}  {r['status']:>10}  {r['time_s']:>7.1f}s"
              + (f"  {str(r['satisfied']):>10}" if r['satisfied'] is not None else ""))

    # Summary
    passed = [r for r in results if r["status"] == "PASS"]
    wrong = [r for r in results if r["status"] == "WRONG"]
    timeout = [r for r in results if r["status"] == "TIMEOUT"]
    errors = [r for r in results if r["status"] == "ERROR"]

    print()
    print("=" * 62)
    print(f"Summary: {len(passed)} PASS, {len(wrong)} WRONG, {len(timeout)} TIMEOUT, {len(errors)} ERROR")
    if wrong:
        print(f"  Wrong: {[r['name'] for r in wrong]}")
    if timeout:
        print(f"  Timeout: {[r['name'] for r in timeout]}")
    if errors:
        print(f"  Errors: {[r['name'] for r in errors]}")


if __name__ == "__main__":
    main()
