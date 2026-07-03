from __future__ import annotations

import argparse
import multiprocessing as mp
import queue
import sys
import traceback
from dataclasses import dataclass
from pathlib import Path
from time import perf_counter
from typing import Any

PYTHON_PKG = Path(__file__).resolve().parents[1]
REPO_ROOT = PYTHON_PKG.parent
if str(PYTHON_PKG) not in sys.path:
    sys.path.insert(0, str(PYTHON_PKG))

import pyqreach  # noqa: E402
from qiskit import QuantumCircuit  # noqa: E402
from qiskit.circuit import CircuitInstruction  # noqa: E402

from parse_qiskit import parse_qiskit_cir_lazy  # noqa: E402


@dataclass
class PrefixResult:
    status: str
    prefix_len: int
    num_qubits: int | str = ""
    num_gates: int | str = ""
    num_locations: int | str = ""
    num_result_locations: int | str = ""
    time_load_qasm: float | str = ""
    time_build_prefix: float | str = ""
    time_parse: float | str = ""
    time_fixed_post: float | str = ""
    time_total: float | str = ""
    error: str = ""


def _resolve_qasm(path: str) -> Path:
    qasm_path = Path(path).expanduser()
    if not qasm_path.is_absolute():
        candidate_from_cwd = Path.cwd() / qasm_path
        candidate_from_python_pkg = PYTHON_PKG / qasm_path
        if candidate_from_cwd.exists():
            qasm_path = candidate_from_cwd
        else:
            qasm_path = candidate_from_python_pkg
    qasm_path = qasm_path.resolve()
    if not qasm_path.exists():
        raise FileNotFoundError(f"QASM file does not exist: {qasm_path}")
    return qasm_path


def make_prefix_circuit(qc: QuantumCircuit, max_instructions: int | None) -> QuantumCircuit:
    prefix_len = len(qc.data) if max_instructions is None else min(max_instructions, len(qc.data))
    prefix = QuantumCircuit(*qc.qregs, *qc.cregs, name=f"{qc.name}_prefix_{prefix_len}")
    for instruction in qc.data[:prefix_len]:
        prefix.append(
            CircuitInstruction(
                instruction.operation,
                instruction.qubits,
                instruction.clbits,
            )
        )
    return prefix


def run_prefix_once(
    qasm_path: Path,
    *,
    max_instructions: int | None,
    initial_state: str | None,
    run_fixed_post: bool,
) -> PrefixResult:
    total_start = perf_counter()
    load_start = perf_counter()
    qc = QuantumCircuit.from_qasm_file(str(qasm_path))
    time_load_qasm = perf_counter() - load_start

    build_start = perf_counter()
    prefix = make_prefix_circuit(qc, max_instructions)
    time_build_prefix = perf_counter() - build_start

    state = initial_state or ("0" * prefix.num_qubits)
    if len(state) != prefix.num_qubits:
        raise ValueError(f"Initial state length {len(state)} does not match {prefix.num_qubits} qubits")

    pyqreach.initializeTransitionSystem()
    ts = pyqreach.TransitionSystem()

    parse_start = perf_counter()
    parse_result = parse_qiskit_cir_lazy(
        prefix,
        prefix.num_qubits,
        ts,
        initial_state=state,
        return_metadata=True,
    )
    time_parse = perf_counter() - parse_start

    time_fixed_post: float | str = ""
    if run_fixed_post:
        fixed_start = perf_counter()
        ts.computingFixedPointPost()
        time_fixed_post = perf_counter() - fixed_start

    return PrefixResult(
        status="ok",
        prefix_len=len(prefix.data),
        num_qubits=prefix.num_qubits,
        num_gates=len(prefix.data),
        num_locations=ts.getLocationNum(),
        num_result_locations=len(parse_result.result_locations),
        time_load_qasm=time_load_qasm,
        time_build_prefix=time_build_prefix,
        time_parse=time_parse,
        time_fixed_post=time_fixed_post,
        time_total=perf_counter() - total_start,
    )


def _worker(payload: dict[str, Any], out_queue: mp.Queue) -> None:
    try:
        result = run_prefix_once(
            Path(payload["qasm_path"]),
            max_instructions=payload["max_instructions"],
            initial_state=payload["initial_state"],
            run_fixed_post=payload["run_fixed_post"],
        )
    except Exception as exc:  # noqa: BLE001 - diagnostic script should report failures
        result = PrefixResult(
            status="error",
            prefix_len=payload["max_instructions"] if payload["max_instructions"] is not None else -1,
            error=f"{type(exc).__name__}: {exc}\n{traceback.format_exc(limit=5)}",
        )
    out_queue.put(result)


def run_prefix_with_timeout(
    qasm_path: Path,
    *,
    max_instructions: int | None,
    initial_state: str | None,
    run_fixed_post: bool,
    timeout_seconds: float,
) -> PrefixResult:
    ctx = mp.get_context("spawn")
    out_queue = ctx.Queue()
    payload = {
        "qasm_path": str(qasm_path),
        "max_instructions": max_instructions,
        "initial_state": initial_state,
        "run_fixed_post": run_fixed_post,
    }
    proc = ctx.Process(target=_worker, args=(payload, out_queue))
    proc.start()
    proc.join(timeout_seconds)
    prefix_len = max_instructions if max_instructions is not None else -1
    if proc.is_alive():
        proc.terminate()
        proc.join()
        return PrefixResult(status="timeout", prefix_len=prefix_len, time_total=timeout_seconds)
    try:
        return out_queue.get_nowait()
    except queue.Empty:
        return PrefixResult(
            status="error",
            prefix_len=prefix_len,
            error=f"worker exited with code {proc.exitcode} without returning a result",
        )


def format_result(result: PrefixResult) -> str:
    fields = [
        f"status={result.status}",
        f"prefix_len={result.prefix_len}",
        f"num_qubits={result.num_qubits}",
        f"num_gates={result.num_gates}",
        f"num_locations={result.num_locations}",
        f"num_result_locations={result.num_result_locations}",
        f"time_load_qasm={result.time_load_qasm}",
        f"time_build_prefix={result.time_build_prefix}",
        f"time_parse={result.time_parse}",
        f"time_fixed_post={result.time_fixed_post}",
        f"time_total={result.time_total}",
    ]
    if result.error:
        fields.append(f"error={result.error}")
    return " ".join(fields)


def bisect_prefix(
    qasm_path: Path,
    *,
    high: int,
    initial_state: str | None,
    run_fixed_post: bool,
    timeout_seconds: float,
) -> tuple[int | None, PrefixResult | None, PrefixResult | None]:
    low = 0
    failing: PrefixResult | None = None
    passing: PrefixResult | None = None

    while low < high:
        mid = (low + high) // 2
        if mid == low:
            mid = high
        result = run_prefix_with_timeout(
            qasm_path,
            max_instructions=mid,
            initial_state=initial_state,
            run_fixed_post=run_fixed_post,
            timeout_seconds=timeout_seconds,
        )
        print(f"[bisect] {format_result(result)}", flush=True)
        if result.status == "timeout":
            failing = result
            high = mid
        elif result.status == "ok":
            passing = result
            low = mid
        else:
            failing = result
            high = mid
        if high - low <= 1:
            break

    smallest_failing = high if failing is not None else None
    if smallest_failing is not None:
        failing = run_prefix_with_timeout(
            qasm_path,
            max_instructions=smallest_failing,
            initial_state=initial_state,
            run_fixed_post=run_fixed_post,
            timeout_seconds=timeout_seconds,
        )
        if smallest_failing > 0:
            passing = run_prefix_with_timeout(
                qasm_path,
                max_instructions=smallest_failing - 1,
                initial_state=initial_state,
                run_fixed_post=run_fixed_post,
                timeout_seconds=timeout_seconds,
            )
    return smallest_failing, passing, failing


def main() -> None:
    parser = argparse.ArgumentParser(description="Debug lazy QASM prefix timeouts for converted Grover circuits")
    parser.add_argument("--qasm", required=True, help="Path to a QASM file")
    parser.add_argument("--max-instructions", type=int, default=None, help="Maximum prefix instruction count")
    parser.add_argument("--timeout-seconds", type=float, default=30.0, help="Subprocess timeout for prefix runs")
    parser.add_argument("--bisect", action="store_true", help="Find the smallest prefix that exceeds the timeout")
    parser.add_argument("--skip-fixed-post", action="store_true", help="Only time lazy parse, not fixed-point post")
    parser.add_argument("--initial-state", default=None, help="Initial product state; defaults to all zeroes")
    args = parser.parse_args()

    qasm_path = _resolve_qasm(args.qasm)
    run_fixed_post = not args.skip_fixed_post

    if args.bisect:
        full_qc = QuantumCircuit.from_qasm_file(str(qasm_path))
        high = args.max_instructions if args.max_instructions is not None else len(full_qc.data)
        smallest, passing, failing = bisect_prefix(
            qasm_path,
            high=high,
            initial_state=args.initial_state,
            run_fixed_post=run_fixed_post,
            timeout_seconds=args.timeout_seconds,
        )
        print(f"smallest_failing_prefix={smallest}")
        print(f"largest_adjacent_passing={format_result(passing) if passing else ''}")
        print(f"smallest_failing_result={format_result(failing) if failing else ''}")
        return

    result = run_prefix_with_timeout(
        qasm_path,
        max_instructions=args.max_instructions,
        initial_state=args.initial_state,
        run_fixed_post=run_fixed_post,
        timeout_seconds=args.timeout_seconds,
    )
    print(format_result(result))


if __name__ == "__main__":
    main()
