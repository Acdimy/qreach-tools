"""Debug-mode verification of smaller Grover benchmarks (non-linear).

Tests grover4/8/16 variants where inner products (2^{-n/2}) are large enough
to avoid the checkifzero 1e-8 threshold issue.

Uses the H-gate workaround for constructing |+> states to avoid the
simpleProductStateAmplitudes overflow for qubits >= 32 (grover16 has 31→32 qubits).
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
TIMEOUT_PER_TEST = 60


def _make_basis_state(bitstring: str) -> pyqreach.QOperation:
    return pyqreach.QOperation([bitstring])


def _make_plus_state(n_qubits: int, qubit_indices: list[int]) -> pyqreach.QOperation:
    """Create |0>^n with H applied to each qubit in `qubit_indices`."""
    state = pyqreach.QOperation(["0" * n_qubits])
    for i in qubit_indices:
        h_gate = pyqreach.QOperation("H", n_qubits, [i], [])
        state = state.post_image(h_gate)
    return state


def _construct_target_subspace(n_qubits: int):
    """Construct span{ |0>^n , |0>^(n/2) |+>^(n/2) } using span_qops (Gram-Schmidt)."""
    half = n_qubits // 2
    basis_zero = _make_basis_state("0" * n_qubits)
    basis_plus_half = _make_plus_state(n_qubits, list(range(half, n_qubits)))
    return pyqreach.span_qops([basis_zero, basis_plus_half])


def _construct_target_disjunction(n_qubits: int):
    """Same target but using disjunction (semantically should be equivalent)."""
    half = n_qubits // 2
    basis_zero = _make_basis_state("0" * n_qubits)
    basis_plus_half = _make_plus_state(n_qubits, list(range(half, n_qubits)))
    return basis_zero.disjunction(basis_plus_half)


def run_one(qasm_file: str) -> dict:
    qasm_path = BENCHMARK_DIR / qasm_file
    name = qasm_path.stem
    t0 = perf_counter()

    qc = QuantumCircuit.from_qasm_file(str(qasm_path))
    n = qc.num_qubits
    half = n // 2

    # Always use H-gate workaround to avoid simpleProductStateAmplitudes issues
    use_string_path = False

    initial_state = "0" * n
    pyqreach.initializeTransitionSystem()
    ts = pyqreach.TransitionSystem()

    t1 = perf_counter()
    parse_result = parse_qiskit_cir_lazy(
        qc, n, ts, initial_state=initial_state, return_metadata=True,
    )
    t_parse = perf_counter() - t1

    # Build target via H-gate workaround (always works)
    t2 = perf_counter()
    target_span = _construct_target_subspace(n)
    target_disj = _construct_target_disjunction(n)
    t_subspace = perf_counter() - t2

    result_locs = list(parse_result.result_locations or [])
    final_loc = result_locs[-1] if result_locs else -1

    # Check using span_qops
    sat_span = all(ts.Locations[loc].satisfy(target_span) for loc in result_locs) if result_locs else False
    # Check using disjunction
    sat_disj = all(ts.Locations[loc].satisfy(target_disj) for loc in result_locs) if result_locs else False

    # Compare with string path if available
    sat_span_str = None
    sat_disj_str = None
    if use_string_path:
        sat_span_str = all(ts.Locations[loc].satisfy(target_span_str) for loc in result_locs)
        sat_disj_str = all(ts.Locations[loc].satisfy(target_disj_str) for loc in result_locs)

    # Diagnostic: what is the final state?
    zero_n = _make_basis_state("0" * n)
    plus_n = _make_plus_state(n, list(range(n)))
    final_in_zero = ts.Locations[final_loc].satisfy(pyqreach.span_qops([zero_n])) if final_loc >= 0 else None
    final_in_plus = ts.Locations[final_loc].satisfy(pyqreach.span_qops([plus_n])) if final_loc >= 0 else None

    return {
        "name": name,
        "qubits": n,
        "half": half,
        "gates": len(qc.data),
        "locations": ts.getLocationNum(),
        "parse_time": t_parse,
        "subspace_time": t_subspace,
        "sat_span": sat_span,
        "sat_disj": sat_disj,
        "sat_span_str": sat_span_str,
        "sat_disj_str": sat_disj_str,
        "final_in_zero": final_in_zero,
        "final_in_plus": final_in_plus,
        "use_string_path": use_string_path,
    }


def main():
    files = [
        "single-it-grover4-plus.qasm", "single-it-grover4-zero.qasm",
        "single-it-grover8-plus.qasm", "single-it-grover8-zero.qasm",
        "single-it-grover16-plus.qasm", "single-it-grover16-zero.qasm",
    ]

    print(f"Debug-verifying {len(files)} smaller Grover benchmarks...")
    print(f"  Target subspace: span{{ |0...0>, |0...0, +...+> }}")
    print()

    results = []
    for f in files:
        name = Path(f).stem
        print(f"  {name:<35} ", end="", flush=True)
        r = run_one(f)
        results.append(r)

        parts = []
        if r["sat_span"]:
            parts.append("span_qops:✓")
        else:
            parts.append("span_qops:✗")
        if r["sat_disj"]:
            parts.append("disj:✓")
        else:
            parts.append("disj:✗")

        if r["use_string_path"]:
            if r["sat_span_str"]:
                parts.append("str_span:✓")
            else:
                parts.append("str_span:✗")

        extra = f"final∈|0>:{r['final_in_zero']}, ∈|+>:{r['final_in_plus']}"
        print(f"{' '.join(parts):<45} ({r['qubits']}q→{1<<((r['qubits']-1).bit_length())}, {r['parse_time']:.3f}s)  [{extra}]")

    # Summary
    print()
    print(f"{'='*100}")
    header = f"{'Circuit':<30} {'Q':>3} {'span_qops':>10} {'disj':>6} {'∈|0>':>6} {'∈|+>':>6}  {'Parse':>8}"
    print(header)
    print(f"{'-'*100}")
    for r in results:
        print(f"{r['name']:<30} {r['qubits']:>3} {str(r['sat_span']):>10} {str(r['sat_disj']):>6} "
              f"{str(r['final_in_zero']):>6} {str(r['final_in_plus']):>6}  {r['parse_time']:>7.3f}s")

    # Check consistency
    print()
    all_span_ok = all(r["sat_span"] for r in results)
    all_disj_ok = all(r["sat_disj"] for r in results)
    print(f"span_qops all-pass: {all_span_ok}")
    print(f"disjunction all-pass: {all_disj_ok}")


if __name__ == "__main__":
    main()
