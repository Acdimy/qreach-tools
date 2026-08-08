"""Quick reproduction of Grover32 timeout behavior after BV fix.

Tests the corrected behavior:
- Grover32-plus: should TIMEOUT (Reduce fragmentation at level=7)
- Grover32-plus-linear: should be FAST (symmetric CCX, power-of-2)
- Grover32-zero: should TIMEOUT (H gates in diffusion trigger same path)

The initial state is always |0>^n; H^⊗n superposition is created by the circuit's
first H-gate layer (for "plus" variants) or the diffusion operator (for "zero").
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
TIMEOUT = 30  # seconds per test


def run_one(name: str, qasm_file: str, timeout: int = TIMEOUT) -> dict:
    qasm_path = BENCHMARK_DIR / qasm_file
    print(f"\n{'='*60}")
    print(f"Testing: {name} ({qasm_file})")
    print(f"{'='*60}")

    t0 = perf_counter()
    try:
        qc = QuantumCircuit.from_qasm_file(str(qasm_path))
        t_load = perf_counter() - t0
        print(f"  Load QASM: {t_load:.3f}s")
        print(f"  Qubits: {qc.num_qubits}, Gates: {len(qc.data)}")

        initial_state = "0" * qc.num_qubits
        print(f"  Initial state: |0>^{qc.num_qubits}")

        pyqreach.initializeTransitionSystem()
        ts = pyqreach.TransitionSystem()

        t1 = perf_counter()
        result = parse_qiskit_cir_lazy(
            qc, qc.num_qubits, ts,
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
            status = f"FAST ({t_parse:.3f}s)"

        print(f"  Parse time: {t_parse:.3f}s")
        print(f"  Locations: {ts.getLocationNum()}")
        print(f"  Result locations: {len(result.result_locations)}")
        print(f"  Total time: {total:.3f}s")
        print(f"  Status: {status}")

        return {"name": name, "status": status, "parse_time": t_parse, "total_time": total,
                "locations": ts.getLocationNum()}

    except Exception as e:
        total = perf_counter() - t0
        if total > timeout:
            status = "TIMEOUT"
        else:
            status = f"ERROR: {e}"
        print(f"  Total time: {total:.3f}s")
        print(f"  Status: {status}")
        return {"name": name, "status": status, "parse_time": total, "total_time": total, "locations": -1}


def main():
    tests = [
        # The only FAST case: 64 qubits, power-of-2, symmetric CCX
        ("Grover32-plus-linear", "single-it-grover32-plus-linear.qasm", 10),
        # Should TIMEOUT with BV fix: nonlinear CCX at level=7
        ("Grover32-plus", "single-it-grover32-plus.qasm", 30),
        # Should TIMEOUT with BV fix: H gates in diffusion trigger path
        ("Grover32-zero", "single-it-grover32-zero.qasm", 30),
    ]

    # Run fast case first to confirm the build works
    results = []
    for name, qasm_file, timeout in tests:
        r = run_one(name, qasm_file, timeout)
        results.append(r)

    print(f"\n{'='*60}")
    print("SUMMARY")
    print(f"{'='*60}")
    print(f"{'Circuit':<25} {'Parse time':>12} {'Locations':>10} {'Status'}")
    print("-" * 60)
    for r in results:
        print(f"{r['name']:<25} {r['parse_time']:>10.3f}s {r['locations']:>10} {r['status']}")

    # Expected outcomes after BV fix:
    expected = {
        "Grover32-plus-linear": "FAST",
        "Grover32-plus": "TIMEOUT",
        "Grover32-zero": "TIMEOUT",
    }
    print(f"\n{'='*60}")
    print("VERIFICATION (against expected outcomes after BV fix)")
    print(f"{'='*60}")
    all_ok = True
    for r in results:
        actual = "FAST" if "FAST" in r["status"] else ("TIMEOUT" if "TIMEOUT" in r["status"] else "OTHER")
        exp = expected.get(r["name"], "?")
        match = "✓" if actual == exp else "✗"
        if actual != exp:
            all_ok = False
        print(f"  {match} {r['name']}: expected={exp}, actual={actual}")

    if all_ok:
        print("\n✓ All results match expected outcomes.")
    else:
        print("\n✗ Some results DO NOT match expected outcomes.")


if __name__ == "__main__":
    main()
