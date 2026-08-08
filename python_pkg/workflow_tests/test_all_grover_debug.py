"""Complete debug-mode verification of ALL Grover benchmarks with CORRECTED target.

BUG FOUND in existing `_run_debug_check` (qasm_workflow_runner.py:281):
The old target `|0>^half |+>^(n-half)` puts `+` on the SECOND half.
But the circuit applies H to the FIRST {search} qubits (the search register).

CORRECTED target: span{ |0>^n , |+>^{search} |0>^{helper} }
  where search = diffusion-H qubit count (from circuit structure).

The diffusion H covers only the search qubits, not all qubits — this is a
simplified Grover where the diffusion operates only on the search register.
"""

from __future__ import annotations

import re, sys
from pathlib import Path
from time import perf_counter

PYTHON_PKG = Path(__file__).resolve().parents[1]
if str(PYTHON_PKG) not in sys.path:
    sys.path.insert(0, str(PYTHON_PKG))

import pyqreach
from qiskit import QuantumCircuit
from parse_qiskit import parse_qiskit_cir_lazy

BENCHMARK_DIR = PYTHON_PKG / "benchmark" / "converted_qasm"
TIMEOUT_PER_TEST = 120


def _make_state(n_qubits: int, plus_indices: list[int] | None = None) -> pyqreach.QOperation:
    """Create |0>^n with H applied to each qubit in `plus_indices`."""
    state = pyqreach.QOperation(["0" * n_qubits])
    if plus_indices:
        for i in plus_indices:
            h_gate = pyqreach.QOperation("H", n_qubits, [i], [])
            state = state.post_image(h_gate)
    return state


def _get_search_qubits(qc: QuantumCircuit) -> int:
    """Return the number of search qubits from the diffusion-H layer."""
    # Count last consecutive H gates (diffusion): they cover search qubits
    count = 0
    for inst in reversed(qc.data):
        if inst.operation.name == 'h':
            count += 1
        else:
            break
    return count if count > 0 else 1  # fallback


def run_one(qasm_file: str) -> dict | None:
    qasm_path = BENCHMARK_DIR / qasm_file
    name = qasm_path.stem
    t0 = perf_counter()

    qc = QuantumCircuit.from_qasm_file(str(qasm_path))
    n = qc.num_qubits
    search = _get_search_qubits(qc)
    half = n // 2

    initial_state = "0" * n
    pyqreach.initializeTransitionSystem()
    ts = pyqreach.TransitionSystem()

    t1 = perf_counter()
    parse_result = parse_qiskit_cir_lazy(
        qc, n, ts, initial_state=initial_state, return_metadata=True,
    )
    t_parse = perf_counter() - t1

    if t_parse > TIMEOUT_PER_TEST:
        print(f"TIMEOUT ({t_parse:.0f}s)")
        return None

    # Build CORRECTED target: span{ |0>^n , |+>^search |0>^helper }
    t2 = perf_counter()
    zero_n = _make_state(n)
    plus_search = _make_state(n, list(range(search)))
    target = pyqreach.span_qops([zero_n, plus_search])
    t_subspace = perf_counter() - t2

    # Check
    result_locs = list(parse_result.result_locations or [])
    all_sat = all(ts.Locations[loc].satisfy(target) for loc in result_locs)

    # Diagnostic
    final_loc = result_locs[-1] if result_locs else -1
    final_in_zero = ts.Locations[final_loc].satisfy(pyqreach.span_qops([zero_n])) if final_loc >= 0 else None

    return {
        "name": name, "qubits": n, "half": half, "search": search,
        "gates": len(qc.data), "locations": ts.getLocationNum(),
        "result_locs": len(result_locs),
        "parse_time": t_parse, "subspace_time": t_subspace,
        "satisfied": all_sat, "final_in_zero": final_in_zero,
    }


def main():
    import glob

    all_files = sorted(glob.glob(str(BENCHMARK_DIR / "single-it-grover*.qasm")))
    # Separate into linear and non-linear
    linear = sorted([f for f in all_files if '-linear' in f])
    nonlinear = sorted([f for f in all_files if '-linear' not in f])

    print(f"Debug-verifying Grover benchmarks with CORRECTED target:")
    print(f"  span{{ |0>^n , |+>^{{search}} |0>^{{helper}} }}")
    print(f"  ({len(linear)} linear + {len(nonlinear)} non-linear)")
    print()

    results = []
    for lbl, group in [("LINEAR", linear), ("NON-LINEAR", nonlinear)]:
        print(f"--- {lbl} ---")
        for f in group:
            fname = Path(f).name
            name = Path(f).stem

            # Skip very slow non-linear circuits
            m = re.search(r'grover(\d+)', name)
            n_val = int(m.group(1)) if m else 999
            if lbl == "NON-LINEAR" and n_val > 16:
                print(f"  {name:<45} SKIP (too slow)")
                continue

            print(f"  {name:<45} ", end="", flush=True)
            r = run_one(fname)
            if r is None:
                continue
            results.append(r)

            status = "✓" if r["satisfied"] else "✗ FAILED"
            extra = f"final∈|0>:{r['final_in_zero']}"
            print(f"{status:<10} ({r['qubits']}q, search={r['search']}, {r['gates']}g, {r['parse_time']:.3f}s parse)  [{extra}]")
        print()

    # Summary
    print(f"{'='*100}")
    print(f"{'Circuit':<38} {'Q':>3} {'S':>3} {'satisfied':>10} {'∈|0>':>6}  {'Parse':>8}")
    print(f"{'-'*100}")
    all_ok = True
    for r in results:
        chk = "✓" if r["satisfied"] else "✗"
        print(f"{r['name']:<38} {r['qubits']:>3} {r['search']:>3} {chk:>10} {str(r['final_in_zero']):>6}  {r['parse_time']:>7.3f}s")
        if not r["satisfied"]:
            all_ok = False

    print()
    if all_ok:
        print("✓ All tested benchmarks pass the corrected debug subspace check.")
    else:
        print("✗ Some benchmarks FAILED.")


if __name__ == "__main__":
    main()
