"""Test all -linear QASM variants for performance.

Linear variants have symmetric CCX qubit ordering and are expected to be fast.
However, non-power-of-2 qubit counts (grover150: 300→512) may still hit issues.
"""

from __future__ import annotations

import sys
from pathlib import Path
from time import perf_counter

PYTHON_PKG = Path(__file__).resolve().parents[1]
if str(PYTHON_PKG) not in sys.path:
    sys.path.insert(0, str(PYTHON_PKG))

import pyqreach
from qiskit import QuantumCircuit
from parse_qiskit import parse_qiskit_cir_lazy

BENCHMARK_DIR = PYTHON_PKG / "benchmark" / "converted_qasm"

TIMEOUT_PER_TEST = 60  # seconds — linear variants should all be fast


def run_one(qasm_file: str, timeout: int = TIMEOUT_PER_TEST) -> dict:
    qasm_path = BENCHMARK_DIR / qasm_file
    name = qasm_path.stem
    t0 = perf_counter()

    qc = QuantumCircuit.from_qasm_file(str(qasm_path))
    t_load = perf_counter() - t0

    initial_state = "0" * qc.num_qubits
    qnum = qc.num_qubits
    cflobdd_qnum = 1 << (qnum.bit_length() - 1)
    if cflobdd_qnum < qnum:
        cflobdd_qnum <<= 1
    is_power_of_2 = cflobdd_qnum == qnum

    pyqreach.initializeTransitionSystem()
    ts = pyqreach.TransitionSystem()

    t1 = perf_counter()
    result = parse_qiskit_cir_lazy(
        qc, qnum, ts,
        initial_state=initial_state,
        return_metadata=True,
    )
    t_parse = perf_counter() - t1
    total = perf_counter() - t0

    if t_parse > timeout:
        status = "TIMEOUT"
    elif t_parse > 1.0:
        status = f"SLOW ({t_parse:.1f}s)"
    else:
        status = f"OK ({t_parse:.3f}s)"

    return {
        "name": name,
        "qubits": qnum,
        "cflobdd_qnum": cflobdd_qnum,
        "pow2": is_power_of_2,
        "gates": len(qc.data),
        "locations": ts.getLocationNum(),
        "result_locs": len(result.result_locations),
        "parse_time": t_parse,
        "load_time": t_load,
        "total_time": total,
        "status": status,
    }


def main():
    import glob

    files = sorted(glob.glob(str(BENCHMARK_DIR / "*-linear*")))
    files = [Path(f).name for f in files]

    print(f"Testing {len(files)} linear variants...")
    print()

    results = []
    for f in files:
        name = Path(f).stem
        qbits_str = name.split("-")[0].replace("grover", "").replace("singleit", "")
        print(f"  {name:<45} ", end="", flush=True)
        r = run_one(f)
        results.append(r)
        print(f"{r['status']:<20} ({r['qubits']}q→{r['cflobdd_qnum']}, {r['gates']}g, {r['locations']}locs)")

    # Summary table
    print()
    print(f"{'='*90}")
    print(f"{'Circuit':<38} {'Qubits':>6} {'CFLOBDD':>8} {'2^n':>4} {'Gates':>6} {'Parse':>10} {'Status'}")
    print(f"{'-'*90}")
    all_fast = True
    for r in results:
        pow2 = "✓" if r["pow2"] else "✗"
        print(f"{r['name']:<38} {r['qubits']:>6} {r['cflobdd_qnum']:>8} {pow2:>4} {r['gates']:>6} {r['parse_time']:>9.3f}s {r['status']}")
        if "TIMEOUT" in r["status"] or "SLOW" in r["status"]:
            all_fast = False

    print()
    if all_fast:
        print("✓ All linear variants completed within budget.")
    else:
        print("✗ Some linear variants did NOT complete within budget.")


if __name__ == "__main__":
    main()
