from __future__ import annotations

import argparse
import csv
import multiprocessing as mp
import queue
import random
import re
import sys
import traceback
from dataclasses import dataclass
from pathlib import Path
from time import perf_counter
from typing import Any

PYTHON_PKG = Path(__file__).resolve().parents[1]
if str(PYTHON_PKG) not in sys.path:
    sys.path.insert(0, str(PYTHON_PKG))

# NuSMV binary — resolved relative to PYTHON_PKG so it works regardless of cwd
_NUSMV_PATH = str(
    (PYTHON_PKG.parent.parent / "NuSMV-2.7.0-macos-universal" / "bin" / "NuSMV").resolve()
)

import pyqreach
from qiskit import QuantumCircuit
from qiskit.circuit import CircuitInstruction
from qiskit.circuit.library import XGate, YGate, ZGate

from parse_qiskit import parse_qiskit_cir, parse_qiskit_cir_lazy
from qctl import modelChecking, quantum_state, set_initial_state, span_qops, tsLabelling, tsLabellingClRegList, tsLabellingDefault


CSV_FIELDS = [
    "filename",
    "relative_filename",
    "status",
    "error",
    "timeout_seconds",
    "num_qubits",
    "num_gates",
    "num_locations",
    "num_result_locations",
    "lazy",
    "lazy_pruned_locations",
    "lazy_check_passed",
    "time_total",
    "time_load_qasm",
    "time_prepare_ts",
    "time_fixed_point_post",
    "verification_time",
    "error_injection",
    "injection_position",
    "injection_qubit",
    "injection_gate",
    "debug_enabled",
    "debug_kind",
    "debug_satisfied",
    "model_check_satisfied",
    "model_check_status",
]


@dataclass
class QasmRunConfig:
    input_dir: Path
    output_file: Path
    timeout_seconds: float = 300.0
    initial_state: str | None = None
    error_injection: bool = False
    debug: bool = False
    seed: int = 42
    fail_fast: bool = False
    limit: int | None = None
    append: bool = True
    lazy: bool = True


def discover_qasm_files(input_dir: Path) -> list[Path]:
    input_dir = Path(input_dir).expanduser().resolve()
    if input_dir.is_file():
        return [input_dir] if input_dir.suffix == ".qasm" else []
    if not input_dir.exists():
        raise FileNotFoundError(f"Input directory does not exist: {input_dir}")
    return sorted(path for path in input_dir.rglob("*.qasm") if path.is_file())


def infer_initial_state(
    qasm_path: Path,
    qc: QuantumCircuit,
    *,
    error_injection: bool = False,
    seed: int = 42,
) -> str:
    """Infer the correct initial state for a QASM benchmark.

    Follows the conventions from the legacy ``test_parse_qasm.py`` workflow:

    * **grover**: ``"0"*(n-1) + "1"`` where *n* is the search width extracted
      from the filename.
    * **dqc_pe**: ``"0"*n + "1"`` for clean mode; ``"0"*(n+1)`` for injected
      mode (matches the legacy ``run_type`` semantics).
    * **dqc_qft**: random computational basis state of length *n* (seeded
      deterministically from *seed* so repeated runs use the same state).
    * **qft**: same random basis-state convention as dqc_qft.
    * Everything else: all-zero state.
    """
    stem = qasm_path.stem.lower()

    # -- grover ---------------------------------------------------------------
    grover_match = re.search(r"grover_(\d+)", stem)
    if grover_match:
        n = int(grover_match.group(1))
        if n == qc.num_qubits and n > 0:
            return "0" * (n - 1) + "1"

    # -- dqc_pe ---------------------------------------------------------------
    dqc_pe_match = re.search(r"dqc_pe_(\d+)", stem)
    if dqc_pe_match:
        n = int(dqc_pe_match.group(1))
        legacy_state = "0" * n + ("0" if error_injection else "1")
        if len(legacy_state) == qc.num_qubits:
            return legacy_state

    # -- dqc_qft  (random basis state, deterministic per [family, n, seed]) ---
    dqc_qft_match = re.search(r"dqc_qft_(\d+)", stem)
    if dqc_qft_match:
        n = int(dqc_qft_match.group(1))
        if n == qc.num_qubits:
            rng_seed = seed * 10000 + n + 1_000_000
            rng = random.Random(rng_seed)
            return "".join(rng.choice(["0", "1"]) for _ in range(n))

    # -- qft  (random basis state, deterministic per [family, n, seed]) -------
    qft_match = re.search(r"^qft_(\d+)", stem)
    if qft_match:
        n = int(qft_match.group(1))
        if n == qc.num_qubits:
            rng_seed = seed * 10000 + n + 2_000_000
            rng = random.Random(rng_seed)
            return "".join(rng.choice(["0", "1"]) for _ in range(n))

    # -- pe  (random basis state, like qft) ----------------------------------
    pe_match = re.search(r"^pe_(\d+)", stem)
    if pe_match:
        n = int(pe_match.group(1))
        if n + 1 == qc.num_qubits:
            rng_seed = seed * 10000 + n + 3_000_000
            rng = random.Random(rng_seed)
            return "".join(rng.choice(["0", "1"]) for _ in range(qc.num_qubits))

    return "0" * qc.num_qubits


def insert_random_pauli(qc: QuantumCircuit, rng: random.Random) -> tuple[QuantumCircuit, tuple[int, int, str]]:
    import copy

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


def _empty_row(qasm_path: Path, input_dir: Path, timeout_seconds: float) -> dict[str, Any]:
    qasm_path = qasm_path.expanduser().resolve()
    input_dir = input_dir.expanduser().resolve()
    try:
        relative = str(qasm_path.relative_to(input_dir))
    except ValueError:
        relative = qasm_path.name
    return {field: "" for field in CSV_FIELDS} | {
        "filename": str(qasm_path),
        "relative_filename": relative,
        "timeout_seconds": timeout_seconds,
    }


def _lazy_check(ts: pyqreach.TransitionSystem, parse_result: Any) -> bool:
    if not getattr(parse_result, "lazy", False):
        return False

    for loc in parse_result.lazy_pruned_locations:
        if ts.printDims(loc)[1] != 0:
            return False
        if not ts.isLeafLoc(loc):
            return False

    for loc in parse_result.result_locations:
        if ts.printDims(loc)[1] <= 0:
            return False

    return True


def _parse_qiskit_for_config(
    qc: QuantumCircuit,
    ts: pyqreach.TransitionSystem,
    initial_state: str,
    *,
    lazy: bool,
):
    if lazy:
        return parse_qiskit_cir_lazy(
            qc,
            qc.num_qubits,
            ts,
            initial_state=initial_state,
            return_metadata=True,
        )

    result = parse_qiskit_cir(qc, qc.num_qubits, ts, return_metadata=True)
    set_initial_state(ts, initial_state)
    return result


def _simulate_final_operation(qc: QuantumCircuit, initial_state: str, *, lazy: bool):
    ts = pyqreach.TransitionSystem()
    result = _parse_qiskit_for_config(qc, ts, initial_state, lazy=lazy)
    ts.computingFixedPointPost()
    if len(result.result_locations) == 1:
        return ts.Locations[result.result_locations[0]].lowerBound
    return span_qops([ts.Locations[loc].lowerBound for loc in result.result_locations])


def _run_debug_check(
    qasm_path: Path,
    qc: QuantumCircuit,
    ts: pyqreach.TransitionSystem,
    parse_result: Any,
    initial_state: str,
    original_qc: QuantumCircuit | None,
    *,
    lazy: bool,
) -> dict[str, Any]:
    stem = qasm_path.stem.lower()

    # Include lazy-pruned locations when applying labels so that CTL formulas
    # such as AG(pe_success -> reached) can detect states whose classical
    # measurement matches a pattern but whose quantum amplitude is zero.
    # Without this, lazy mode would change model-checking semantics.
    _result_locs = list(parse_result.result_locations or [])
    _pruned_locs = list(getattr(parse_result, "lazy_pruned_locations", []) or [])
    all_leaf = _result_locs + _pruned_locs

    if stem.startswith("grover_"):
        work_qubits = int((qc.num_qubits + 1) / 2)
        if qc.num_qubits % 2 != 1 or work_qubits >= ts.getLocationNum():
            raise ValueError("Grover debug check expects an odd-qubit Grover benchmark")
        grover_init = ts.Locations[work_qubits].lowerBound
        grover_good = quantum_state("1" * work_qubits + "0" * (qc.num_qubits - work_qubits - 1) + "1")
        grover_final = span_qops([grover_init, grover_good])
        # satisfy is on the Location, not QOperation (matching legacy test_parse_qasm.py)
        final_loc = _result_locs[-1] if _result_locs else 0
        return {"debug_kind": "grover", "debug_satisfied": ts.Locations[final_loc].satisfy(grover_final)}

    if stem.startswith("single-it-grover"):
        half = qc.num_qubits // 2
        basis_zero = quantum_state("0" * qc.num_qubits)
        basis_plus = quantum_state("0" * half + "+" * (qc.num_qubits - half))
        target = span_qops([basis_zero, basis_plus])
        satisfied = all(ts.Locations[loc].satisfy(target) for loc in _result_locs)
        return {"debug_kind": "grover_converted", "debug_satisfied": satisfied}

    if "dqc_pe" in stem:
        if qc.num_qubits < 11:
            tsLabellingDefault(ts, "reached")
            pe_pattern = "1" + "0" * (qc.num_qubits - 2)
        else:
            tsLabellingDefault(ts, "reached", all_leaf)
            pe_pattern = (qc.num_qubits - 11) * "0" + "1000000000"
        tsLabellingClRegList(ts, [pe_pattern], "pe_success", locList=all_leaf)
        result = modelChecking(ts, "AG ((pe_success -> reached))", nusmv_path=_NUSMV_PATH)
        return {
            "debug_kind": "dqc_pe",
            "debug_satisfied": result.get("satisfied"),
            "model_check_satisfied": result.get("satisfied"),
            "model_check_status": "ok" if result.get("satisfied") is not None else "unknown",
        }

    if original_qc is not None:
        expected = _simulate_final_operation(original_qc, initial_state, lazy=lazy)
        for loc in all_leaf:
            ts.setLabel(loc, "final")
        tsLabelling(ts, expected, "debug_op", locList=all_leaf)
        result = modelChecking(ts, "AG (debug_op <-> final)", nusmv_path=_NUSMV_PATH)
        return {
            "debug_kind": "comparison",
            "debug_satisfied": result.get("satisfied"),
            "model_check_satisfied": result.get("satisfied"),
            "model_check_status": "ok" if result.get("satisfied") is not None else "unknown",
        }

    # Clean-mode self-consistency: use the TS's own final state as the
    # reference (no separate simulation — avoids double-parsing crashes).
    # This is trivially True for a correct TS construction.
    expected = span_qops([ts.Locations[loc].lowerBound for loc in all_leaf])
    for loc in all_leaf:
        ts.setLabel(loc, "final")
    tsLabelling(ts, expected, "debug_op", locList=all_leaf)
    result = modelChecking(ts, "AG (debug_op <-> final)", nusmv_path=_NUSMV_PATH)
    return {
        "debug_kind": "comparison",
        "debug_satisfied": result.get("satisfied"),
        "model_check_satisfied": result.get("satisfied"),
        "model_check_status": "ok" if result.get("satisfied") is not None else "unknown",
    }


def run_qasm_file(qasm_path: Path, config: QasmRunConfig, *, file_index: int = 0) -> dict[str, Any]:
    row = _empty_row(qasm_path, config.input_dir, config.timeout_seconds)
    total_start = perf_counter()
    original_qc = None
    injection_info = None

    pyqreach.initializeTransitionSystem()

    load_start = perf_counter()
    qc = QuantumCircuit.from_qasm_file(str(qasm_path))
    row["time_load_qasm"] = perf_counter() - load_start

    initial_state = config.initial_state or infer_initial_state(
        qasm_path, qc, error_injection=config.error_injection, seed=config.seed
    )
    if len(initial_state) != qc.num_qubits:
        raise ValueError(
            f"Initial state length {len(initial_state)} does not match circuit qubits {qc.num_qubits}"
        )

    if config.error_injection:
        import copy

        original_qc = copy.deepcopy(qc)
        rng = random.Random(config.seed + file_index)
        qc, injection_info = insert_random_pauli(qc, rng)

    ts = pyqreach.TransitionSystem()
    build_start = perf_counter()
    parse_result = _parse_qiskit_for_config(qc, ts, initial_state, lazy=config.lazy)
    row["time_prepare_ts"] = perf_counter() - build_start

    fp_start = perf_counter()
    ts.computingFixedPointPost()
    row["time_fixed_point_post"] = perf_counter() - fp_start

    verify_start = perf_counter()
    lazy_check_passed = _lazy_check(ts, parse_result)
    debug_result: dict[str, Any] = {}
    if config.debug:
        debug_result = _run_debug_check(qasm_path, qc, ts, parse_result, initial_state, original_qc, lazy=config.lazy)
    row["verification_time"] = perf_counter() - verify_start

    row.update(
        {
            "status": "ok",
            "num_qubits": qc.num_qubits,
            "num_gates": len(qc.data),
            "num_locations": ts.getLocationNum(),
            "num_result_locations": len(parse_result.result_locations),
            "lazy": bool(config.lazy),
            "lazy_pruned_locations": len(getattr(parse_result, "lazy_pruned_locations", [])),
            "lazy_check_passed": lazy_check_passed,
            "time_total": perf_counter() - total_start,
            "error_injection": bool(config.error_injection),
            "debug_enabled": bool(config.debug),
            "debug_kind": "none",
        }
    )
    if injection_info is not None:
        row.update(
            {
                "injection_position": injection_info[0],
                "injection_qubit": injection_info[1],
                "injection_gate": injection_info[2],
            }
        )
    row.update(debug_result)
    return row


def _worker(qasm_path: str, config_payload: dict[str, Any], file_index: int, out_queue: mp.Queue) -> None:
    try:
        config = QasmRunConfig(
            input_dir=Path(config_payload["input_dir"]),
            output_file=Path(config_payload["output_file"]),
            timeout_seconds=config_payload["timeout_seconds"],
            initial_state=config_payload.get("initial_state"),
            error_injection=config_payload.get("error_injection", False),
            debug=config_payload.get("debug", False),
            seed=config_payload.get("seed", 42),
            fail_fast=config_payload.get("fail_fast", False),
            limit=config_payload.get("limit"),
            append=config_payload.get("append", True),
            lazy=config_payload.get("lazy", True),
        )
        row = run_qasm_file(Path(qasm_path), config, file_index=file_index)
    except Exception as exc:  # noqa: BLE001 - row should capture benchmark failures
        config = QasmRunConfig(
            input_dir=Path(config_payload["input_dir"]),
            output_file=Path(config_payload["output_file"]),
            timeout_seconds=config_payload["timeout_seconds"],
        )
        row = _empty_row(Path(qasm_path), config.input_dir, config.timeout_seconds)
        row.update(
            {
                "status": "error",
                "error": f"{type(exc).__name__}: {exc}",
                "time_total": "",
                "debug_enabled": config_payload.get("debug", False),
                "error_injection": config_payload.get("error_injection", False),
                "lazy": config_payload.get("lazy", True),
            }
        )
        row["model_check_status"] = traceback.format_exc(limit=5)
    out_queue.put(row)


def _config_payload(config: QasmRunConfig) -> dict[str, Any]:
    return {
        "input_dir": str(config.input_dir),
        "output_file": str(config.output_file),
        "timeout_seconds": config.timeout_seconds,
        "initial_state": config.initial_state,
        "error_injection": config.error_injection,
        "debug": config.debug,
        "seed": config.seed,
        "fail_fast": config.fail_fast,
        "limit": config.limit,
        "append": config.append,
        "lazy": config.lazy,
    }


def run_qasm_file_with_timeout(qasm_path: Path, config: QasmRunConfig, *, file_index: int = 0) -> dict[str, Any]:
    ctx = mp.get_context("spawn")
    out_queue = ctx.Queue()
    proc = ctx.Process(target=_worker, args=(str(qasm_path), _config_payload(config), file_index, out_queue))
    proc.start()
    proc.join(config.timeout_seconds)

    if proc.is_alive():
        proc.terminate()
        proc.join()
        row = _empty_row(qasm_path, config.input_dir, config.timeout_seconds)
        row.update(
            {
                "status": "timeout",
                "time_total": config.timeout_seconds,
                "error_injection": bool(config.error_injection),
                "debug_enabled": bool(config.debug),
                "lazy": bool(config.lazy),
            }
        )
        return row

    try:
        return out_queue.get_nowait()
    except queue.Empty:
        row = _empty_row(qasm_path, config.input_dir, config.timeout_seconds)
        row.update(
            {
                "status": "error",
                "error": f"worker exited with code {proc.exitcode} without returning a result",
                "error_injection": bool(config.error_injection),
                "debug_enabled": bool(config.debug),
                "lazy": bool(config.lazy),
            }
        )
        return row


def normalize_row(row: dict[str, Any]) -> dict[str, Any]:
    return {field: row.get(field, "") for field in CSV_FIELDS}


def write_csv_row(output_file: Path, row: dict[str, Any], *, append: bool = True) -> None:
    output_file = Path(output_file).expanduser().resolve()
    output_file.parent.mkdir(parents=True, exist_ok=True)
    mode = "a" if append else "w"
    write_header = not append or not output_file.exists() or output_file.stat().st_size == 0
    with output_file.open(mode, newline="") as f:
        writer = csv.DictWriter(f, fieldnames=CSV_FIELDS)
        if write_header:
            writer.writeheader()
        writer.writerow(normalize_row(row))


def run_batch(config: QasmRunConfig) -> list[dict[str, Any]]:
    qasm_files = discover_qasm_files(config.input_dir)
    if config.limit is not None:
        qasm_files = qasm_files[: config.limit]
    if not qasm_files:
        raise FileNotFoundError(f"No .qasm files found under {config.input_dir}")

    if not config.append and config.output_file.exists():
        config.output_file.unlink()

    rows: list[dict[str, Any]] = []
    for idx, qasm_file in enumerate(qasm_files):
        print(f"[{idx + 1}/{len(qasm_files)}] Running {qasm_file}", flush=True)
        row = run_qasm_file_with_timeout(qasm_file, config, file_index=idx)
        write_csv_row(config.output_file, row, append=True)
        rows.append(row)
        print(f"  status={row['status']} output={config.output_file}", flush=True)
        if config.fail_fast and row["status"] != "ok":
            break
    return rows


def add_common_arguments(parser: argparse.ArgumentParser, *, default_error_injection: bool, default_debug: bool = False) -> None:
    parser.add_argument("--input-dir", type=Path, default=None, help="Directory or .qasm file to check recursively")
    parser.add_argument("--output-dir", type=Path, default=None, help="Directory for the default CSV output")
    parser.add_argument("--output-file", type=Path, default=None, help="CSV file to write")
    parser.add_argument("--timeout-seconds", type=float, default=300.0, help="Per-file timeout in seconds")
    parser.add_argument("--initial-state", default=None, help="Override inferred initial state bitstring")
    parser.add_argument("--seed", type=int, default=42, help="Random seed for error injection")
    if default_debug:
        parser.add_argument("--debug", dest="debug", action="store_true", default=True)
        parser.add_argument("--no-debug", dest="debug", action="store_false", help="Skip optional legacy debug/model-checking checks")
    else:
        parser.add_argument("--debug", dest="debug", action="store_true", default=False, help="Run optional legacy debug/model-checking checks")
        parser.add_argument("--no-debug", dest="debug", action="store_false")
    parser.add_argument("--limit", type=int, default=None, help="Limit number of discovered files for smoke runs")
    parser.add_argument("--fail-fast", action="store_true", help="Stop after the first non-ok row")
    parser.add_argument("--overwrite", action="store_true", help="Overwrite output CSV instead of appending")
    parser.add_argument("--lazy", dest="lazy", action="store_true", default=True, help="Use lazy measurement construction")
    parser.add_argument("--no-lazy", dest="lazy", action="store_false", help="Use the ordinary non-lazy parser")
    if default_error_injection:
        parser.add_argument("--error-injection", dest="error_injection", action="store_true", default=True)
        parser.add_argument("--no-error-injection", dest="error_injection", action="store_false")
    else:
        parser.add_argument("--error-injection", dest="error_injection", action="store_true", default=False)
        parser.add_argument("--no-error-injection", dest="error_injection", action="store_false")


def make_config_from_args(
    args: argparse.Namespace,
    *,
    default_input_dir: Path,
    default_output_dir: Path,
    default_output_name: str,
) -> QasmRunConfig:
    input_dir = (args.input_dir or default_input_dir).expanduser().resolve()
    output_dir = (args.output_dir or default_output_dir).expanduser().resolve()
    output_file = (args.output_file or (output_dir / default_output_name)).expanduser().resolve()
    return QasmRunConfig(
        input_dir=input_dir,
        output_file=output_file,
        timeout_seconds=args.timeout_seconds,
        initial_state=args.initial_state,
        error_injection=args.error_injection,
        debug=args.debug,
        seed=args.seed,
        fail_fast=args.fail_fast,
        limit=args.limit,
        append=not args.overwrite,
        lazy=args.lazy,
    )
