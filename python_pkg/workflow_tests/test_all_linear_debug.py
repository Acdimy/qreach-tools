"""Debug-mode verification of all -linear Grover benchmarks.

Verifies that every result location's final quantum state lies in the expected
subspace for single-iteration Grover circuits:

    span{ |0...0>,  |0...0, +...+> }

where the first half qubits are all zeros and the second half are uniform
superpositions.

Matches the existing `_run_debug_check` logic in qasm_workflow_runner.py
(lines 278-284), but avoids the broken `+`-string path for >= 32 qubits by
constructing H^⊗k states via `post_image` with single-qubit H gates.
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
TIMEOUT_PER_TEST = 120


def _make_basis_state(bitstring: str) -> pyqreach.QOperation:
    """Create a basis state.  Must contain only 0/1 (no +/-)."""
    return pyqreach.QOperation([bitstring])


def _make_plus_state(n_qubits: int, qubit_indices: list[int]) -> pyqreach.QOperation:
    """Create |0>^n with H applied to each qubit in `qubit_indices`.

    Produces a state where the designated qubits are in |+> and the rest in |0>.
    Avoids the `+`-string path which overflows for n >= 32.
    """
    state = pyqreach.QOperation(["0" * n_qubits])
    for i in qubit_indices:
        h_gate = pyqreach.QOperation("H", n_qubits, [i], [])
        state = state.post_image(h_gate)
    return state


def _construct_target_subspace(n_qubits: int) -> pyqreach.QOperation:
    """Construct span{ |0>^n , |0>^(n/2) |+>^(n/2) }.

    This is the target from _run_debug_check for single-it-grover circuits.
    """
    half = n_qubits // 2
    basis_zero = _make_basis_state("0" * n_qubits)
    basis_plus_half = _make_plus_state(n_qubits, list(range(half, n_qubits)))
    return pyqreach.span_qops([basis_zero, basis_plus_half])


def run_debug_one(qasm_file: str) -> dict:
    qasm_path = BENCHMARK_DIR / qasm_file
    name = qasm_path.stem
    t0 = perf_counter()

    qc = QuantumCircuit.from_qasm_file(str(qasm_path))
    n = qc.num_qubits
    half = n // 2
    initial_state = "0" * n

    pyqreach.initializeTransitionSystem()
    ts = pyqreach.TransitionSystem()

    # Parse (lazy)
    t1 = perf_counter()
    parse_result = parse_qiskit_cir_lazy(
        qc, n, ts,
        initial_state=initial_state,
        return_metadata=True,
    )
    t_parse = perf_counter() - t1

    # Build target subspace
    t2 = perf_counter()
    target = _construct_target_subspace(n)
    t_subspace = perf_counter() - t2

    # Check: every result location satisfies the target subspace
    result_locs = list(parse_result.result_locations or [])
    all_satisfied = True
    failed_locs = []
    for loc in result_locs:
        if not ts.Locations[loc].satisfy(target):
            all_satisfied = False
            failed_locs.append(loc)

    # Also check what the final state actually is (diagnostic)
    final_loc = result_locs[-1] if result_locs else -1
    final_in_zero = False
    final_in_plus = False
    if final_loc >= 0:
        zero_n = _make_basis_state("0" * n)
        plus_n = _make_plus_state(n, list(range(n)))
        final_in_zero = ts.Locations[final_loc].satisfy(pyqreach.span_qops([zero_n]))
        final_in_plus = ts.Locations[final_loc].satisfy(pyqreach.span_qops([plus_n]))

    t_total = perf_counter() - t0

    return {
        "name": name,
        "qubits": n,
        "half": half,
        "gates": len(qc.data),
        "locations": ts.getLocationNum(),
        "result_locs": len(result_locs),
        "parse_time": t_parse,
        "subspace_time": t_subspace,
        "total_time": t_total,
        "satisfied": all_satisfied,
        "failed_locs": failed_locs,
        "final_in_zero": final_in_zero,
        "final_in_plus": final_in_plus,
    }


def main():
    import glob

    files = sorted(glob.glob(str(BENCHMARK_DIR / "*-linear*")))
    files = [Path(f).name for f in files]

    print(f"Debug-verifying {len(files)} linear variants...")
    print(f"  Target subspace: span{{ |0...0>, |0...0, +...+> }}")
    print()

    results = []
    for f in files:
        name = Path(f).stem
        print(f"  {name:<45} ", end="", flush=True)
        r = run_debug_one(f)
        results.append(r)
        status = "✓ SATISFIED" if r["satisfied"] else f"✗ FAILED"
        extra = f"final∈|0>: {r['final_in_zero']}, final∈|+>: {r['final_in_plus']}"
        print(f"{status:<20} ({r['qubits']}q, {r['parse_time']:.3f}s parse, {r['subspace_time']:.3f}s subspace)  [{extra}]")

    # Summary
    print()
    print(f"{'='*105}")
    print(f"{'Circuit':<38} {'Qubits':>6} {'Gates':>6} {'Parse':>8} {'Subsp':>8} {'Sat':>5} {'∈|0>':>6} {'∈|+>':>6}")
    print(f"{'-'*105}")
    all_ok = True
    for r in results:
        chk = "✓" if r["satisfied"] else "✗"
        print(f"{r['name']:<38} {r['qubits']:>6} {r['gates']:>6} "
              f"{r['parse_time']:>7.3f}s {r['subspace_time']:>7.3f}s "
              f"{chk:>5} {str(r['final_in_zero']):>6} {str(r['final_in_plus']):>6}")
        if not r["satisfied"]:
            all_ok = False

    print()
    if all_ok:
        print("✓ All linear variants pass the debug subspace check.")
    else:
        print("✗ Some linear variants FAILED the debug subspace check.")
        print("  Note: 'plus' variants have initial H-gate layers that change the final-state subspace.")
        print("  This check is designed for 'zero' variant circuits.")


if __name__ == "__main__":
    main()
