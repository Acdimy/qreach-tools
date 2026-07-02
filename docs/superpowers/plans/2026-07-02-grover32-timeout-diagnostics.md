# Grover32 Timeout Diagnostics Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add opt-in diagnostics that localize the Grover32 lazy QASM timeout before attempting any root-cause fix.

**Architecture:** Add a small parser profiling helper directly in `python_pkg/parse_qiskit.py`, enabled only by environment variables, and wrap the main instruction loop with timing. Add a focused debug script in `python_pkg/workflow_tests/` that truncates QASM circuits, runs lazy parse/fixed-post in subprocesses when requested, and bisects the smallest slow prefix. After the minimal prefix is known, add a focused regression test using that minimized circuit or prefix.

**Tech Stack:** Python 3 via `.venv`/`uv`, Qiskit 1.4.2, `pyqreach` pybind11 extension, existing explicit lazy parser `parse_qiskit_cir_lazy(...)`, C++/CFLOBDD only after profiling localizes the issue.

## Global Constraints

- Work in repository: `/Users/ftdac/thu/qreach-tools`.
- Use `.venv` / `uv` workflow locally; do not assume conda or Docker.
- Preserve existing workflow APIs in `python_pkg/qctl.py`, `python_pkg/parse_qiskit.py`, and `python_pkg/workflow_tests/qasm_workflow_runner.py`.
- Do not extend product-state `QOperation(["..."])` syntax.
- Default batch timeout is 300 seconds; timed-out rows must remain `status=timeout`.
- New debugging/profiling switches must be opt-in and must not change default CSV output semantics.
- Do not implement a workflow workaround before localizing the root cause.
- Do not assume CCX is the final root cause until profiling evidence points there.
- Do not add a 300-second regression test.

---

## File Structure

- Modify: `python_pkg/parse_qiskit.py`
  - Responsibility: expose opt-in parser instruction profiling without changing default behavior.
  - New internal helpers: `_parse_profile_config()`, `_parse_profile_emit(...)`.
- Create: `python_pkg/workflow_tests/debug_grover32_ccx_pathology.py`
  - Responsibility: load QASM, build truncated circuits, time lazy parse/fixed-post, and bisect prefixes using subprocess timeouts.
- Create after minimization: `python_pkg/workflow_tests/test_grover32_ccx_pathology.py`
  - Responsibility: focused regression for the minimized prefix/circuit once known.
- Do not modify: `python_pkg/workflow_tests/qasm_workflow_runner.py` in the diagnostic phase unless localization proves runner-level metadata is needed.
- Inspect only after localization: `quantum_operation.hpp`, `cflobdd/CFLOBDD/matrix1234_node.cpp`.

---

### Task 1: Add opt-in parser instruction profiling

**Files:**
- Modify: `python_pkg/parse_qiskit.py:1-10`
- Modify: `python_pkg/parse_qiskit.py:530-977`
- Test: manual commands below; no committed test required because feature is env-controlled diagnostic output.

**Interfaces:**
- Consumes: existing `parse_qiskit_cir(...)` loop and `parse_qiskit_cir_lazy(...)` wrapper.
- Produces:
  - Env var `QREACH_PARSE_PROFILE=1` enables output.
  - Env var `QREACH_PARSE_PROFILE_THRESHOLD=<float seconds>` controls slow-line threshold; default `0.25`.
  - Env var `QREACH_PARSE_PROFILE_VERBOSE=1` prints every instruction.
  - Output lines begin with `[parse-profile]`.

- [ ] **Step 1: Add imports**

At the top of `python_pkg/parse_qiskit.py`, change the imports from:

```python
import pyqreach
from qiskit import QuantumCircuit, QuantumRegister, ClassicalRegister
```

to:

```python
import os
from time import perf_counter

import pyqreach
from qiskit import QuantumCircuit, QuantumRegister, ClassicalRegister
```

- [ ] **Step 2: Add parser profiling helper functions**

Insert this block after the `LazyPostContext` dataclass:

```python

def _parse_profile_config() -> tuple[bool, float, bool]:
    enabled = os.environ.get("QREACH_PARSE_PROFILE") == "1"
    verbose = os.environ.get("QREACH_PARSE_PROFILE_VERBOSE") == "1"
    threshold_raw = os.environ.get("QREACH_PARSE_PROFILE_THRESHOLD", "0.25")
    try:
        threshold = float(threshold_raw)
    except ValueError:
        threshold = 0.25
    return enabled, threshold, verbose


def _parse_profile_emit(
    *,
    enabled: bool,
    threshold: float,
    verbose: bool,
    instruction_index: int,
    op_name: str,
    qubits: list[int],
    locations_in: int,
    locations_out: int,
    elapsed: float,
    total_locations: int,
) -> None:
    if not enabled:
        return
    if not verbose and elapsed < threshold:
        return
    print(
        "[parse-profile] "
        f"idx={instruction_index} "
        f"op={op_name} "
        f"qubits={qubits} "
        f"locations_in={locations_in} "
        f"locations_out={locations_out} "
        f"elapsed={elapsed:.6f}s "
        f"total_locations={total_locations}",
        flush=True,
    )
```

- [ ] **Step 3: Read profiling config once per parse call**

In `parse_qiskit_cir(...)`, after:

```python
    instructions = qc.data[pivot:pivotend] if pivotend != 1000000 else qc.data[pivot:]
    currLoc = startNodes
    resultLocs = []
    pruning_resets = False
```

add:

```python
    parse_profile_enabled, parse_profile_threshold, parse_profile_verbose = _parse_profile_config()
```

- [ ] **Step 4: Wrap each instruction in a `try/finally` profile block**

Replace the start of the loop:

```python
    for _,gate in enumerate(instructions):
        # Assume each Locs in currLoc has different classical APs (In the current abstractlevel==1)
        op_name = gate.operation.name
```

with:

```python
    for _, gate in enumerate(instructions):
        profile_instruction_index = pivot + _
        profile_start = perf_counter() if parse_profile_enabled else 0.0
        profile_locations_in = len(currLoc)
        profile_op_name = gate.operation.name
        profile_qubits: list[int] = []
        try:
            # Assume each Locs in currLoc has different classical APs (In the current abstractlevel==1)
            op_name = profile_op_name
```

Then indent the existing loop body one level under the `try:` through the end of the loop body.

At the very end of the loop body, immediately before the loop continues to the next instruction, add this `finally:` block aligned with `try:`:

```python
        finally:
            if parse_profile_enabled:
                _parse_profile_emit(
                    enabled=parse_profile_enabled,
                    threshold=parse_profile_threshold,
                    verbose=parse_profile_verbose,
                    instruction_index=profile_instruction_index,
                    op_name=profile_op_name,
                    qubits=profile_qubits,
                    locations_in=profile_locations_in,
                    locations_out=len(currLoc),
                    elapsed=perf_counter() - profile_start,
                    total_locations=ts.getLocationNum(),
                )
```

- [ ] **Step 5: Store qubit list for profiling**

In the loop, after the existing line:

```python
        qubits = get_global_qb_index(gate.qubits) if gate.qubits else []
```

add:

```python
        profile_qubits = qubits
```

If Step 4's indentation causes this line to be inside the `try:`, keep it inside the `try:`.

- [ ] **Step 6: Run syntax check**

Run from repo root:

```bash
./.venv/bin/python -m py_compile python_pkg/parse_qiskit.py
```

Expected: command exits 0 with no output.

- [ ] **Step 7: Verify default behavior has no profile output**

Run from repo root:

```bash
cd python_pkg
PYTHONPATH=. ../.venv/bin/python workflow_tests/test_lazy_measurement.py
```

Expected output includes:

```text
All lazy measurement tests PASSED
```

Expected output does not include:

```text
[parse-profile]
```

- [ ] **Step 8: Verify profiling output on a small circuit**

Run from repo root:

```bash
QREACH_PARSE_PROFILE=1 QREACH_PARSE_PROFILE_VERBOSE=1 ./.venv/bin/python - <<'PY'
import sys
sys.path.insert(0, 'python_pkg')
import pyqreach
from qiskit import QuantumCircuit
from parse_qiskit import parse_qiskit_cir_lazy
pyqreach.initializeTransitionSystem()
qc = QuantumCircuit(2)
qc.h(0)
qc.cx(0, 1)
ts = pyqreach.TransitionSystem()
parse_qiskit_cir_lazy(qc, qc.num_qubits, ts, initial_state='00', return_metadata=True)
PY
```

Expected output contains two lines similar to:

```text
[parse-profile] idx=0 op=h qubits=[0] locations_in=1 locations_out=1 elapsed=...
[parse-profile] idx=1 op=cx qubits=[0, 1] locations_in=1 locations_out=1 elapsed=...
```

- [ ] **Step 9: Commit Task 1**

Run:

```bash
git add python_pkg/parse_qiskit.py
git commit -m "Add opt-in parser profiling

Co-Authored-By: Claude <noreply@anthropic.com>"
```

---

### Task 2: Add QASM prefix timing and bisect debug script

**Files:**
- Create: `python_pkg/workflow_tests/debug_grover32_ccx_pathology.py`
- Test: manual commands below.

**Interfaces:**
- Consumes: `parse_qiskit_cir_lazy(qc, qc.num_qubits, ts, initial_state=..., return_metadata=True)`.
- Produces CLI:
  - `--qasm PATH` required.
  - `--max-instructions N` optional; default full circuit.
  - `--timeout-seconds S` for subprocess prefix checks; default `30`.
  - `--bisect` finds smallest prefix whose subprocess exceeds timeout.
  - `--skip-fixed-post` skips `ts.computingFixedPointPost()`; default is to run fixed-post after parse.
  - `--initial-state STATE` optional; default `0 * qc.num_qubits`.

- [ ] **Step 1: Create the script**

Create `python_pkg/workflow_tests/debug_grover32_ccx_pathology.py` with this content:

```python
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
    parser.add_argument("--timeout-seconds", type=float, default=30.0, help="Subprocess timeout for --bisect")
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

    result = run_prefix_once(
        qasm_path,
        max_instructions=args.max_instructions,
        initial_state=args.initial_state,
        run_fixed_post=run_fixed_post,
    )
    print(format_result(result))


if __name__ == "__main__":
    main()
```

- [ ] **Step 2: Run syntax check**

Run from repo root:

```bash
./.venv/bin/python -m py_compile python_pkg/workflow_tests/debug_grover32_ccx_pathology.py
```

Expected: command exits 0 with no output.

- [ ] **Step 3: Run a short known-fast prefix**

Run from repo root:

```bash
cd python_pkg
../.venv/bin/python workflow_tests/debug_grover32_ccx_pathology.py \
  --qasm benchmark/converted_qasm/single-it-grover32-plus.qasm \
  --max-instructions 10 \
  --skip-fixed-post
```

Expected output begins with:

```text
status=ok prefix_len=10 num_qubits=63 num_gates=10
```

- [ ] **Step 4: Run bisect with parse-only timeout**

Run from repo root:

```bash
cd python_pkg
QREACH_PARSE_PROFILE=1 QREACH_PARSE_PROFILE_THRESHOLD=0.05 ../.venv/bin/python workflow_tests/debug_grover32_ccx_pathology.py \
  --qasm benchmark/converted_qasm/single-it-grover32-plus.qasm \
  --max-instructions 140 \
  --timeout-seconds 20 \
  --skip-fixed-post \
  --bisect
```

Expected: output contains `[bisect]` lines and ends with:

```text
smallest_failing_prefix=<integer or None>
largest_adjacent_passing=...
smallest_failing_result=...
```

If `smallest_failing_prefix=None`, increase `--max-instructions` or lower `--timeout-seconds` in a later investigation run; do not modify code for this case.

- [ ] **Step 5: Commit Task 2**

Run:

```bash
git add python_pkg/workflow_tests/debug_grover32_ccx_pathology.py
git commit -m "Add Grover32 QASM prefix debugger

Co-Authored-By: Claude <noreply@anthropic.com>"
```

---

### Task 3: Run localization experiments and record evidence

**Files:**
- No source modification expected.
- Optional create: `docs/agent-handoffs/grover32-timeout-localization-notes.md` if evidence should be preserved for the next agent.

**Interfaces:**
- Consumes: parser profiling from Task 1 and prefix debugger from Task 2.
- Produces: exact slow instruction/prefix evidence used to decide whether to inspect Python parser, `QOperation` construction, lazy post-image, or C++ CFLOBDD.

- [ ] **Step 1: Run Grover64 profiling as working reference**

Run from repo root:

```bash
QREACH_PARSE_PROFILE=1 QREACH_PARSE_PROFILE_THRESHOLD=0.05 ./.venv/bin/python - <<'PY'
import sys, time
from pathlib import Path
sys.path.insert(0, 'python_pkg')
import pyqreach
from qiskit import QuantumCircuit
from parse_qiskit import parse_qiskit_cir_lazy
path = Path('python_pkg/benchmark/converted_qasm/single-it-grover64-plus.qasm')
pyqreach.initializeTransitionSystem()
qc = QuantumCircuit.from_qasm_file(str(path))
ts = pyqreach.TransitionSystem()
start = time.perf_counter()
res = parse_qiskit_cir_lazy(qc, qc.num_qubits, ts, initial_state='0' * qc.num_qubits, return_metadata=True)
print('parsed', ts.getLocationNum(), len(res.result_locations), time.perf_counter() - start)
PY
```

Expected: command completes quickly, and any `[parse-profile]` lines identify slow-ish but non-timeout instructions.

- [ ] **Step 2: Run Grover32 under short timeout and capture last profile lines**

Run from repo root:

```bash
./.venv/bin/python - <<'PY'
import os, subprocess
code = r'''
import sys, time
from pathlib import Path
sys.path.insert(0, 'python_pkg')
import pyqreach
from qiskit import QuantumCircuit
from parse_qiskit import parse_qiskit_cir_lazy
path = Path('python_pkg/benchmark/converted_qasm/single-it-grover32-plus.qasm')
pyqreach.initializeTransitionSystem()
qc = QuantumCircuit.from_qasm_file(str(path))
ts = pyqreach.TransitionSystem()
start = time.perf_counter()
parse_qiskit_cir_lazy(qc, qc.num_qubits, ts, initial_state='0' * qc.num_qubits, return_metadata=True)
print('parsed', ts.getLocationNum(), time.perf_counter() - start)
'''
try:
    result = subprocess.run(
        ['./.venv/bin/python', '-c', code],
        cwd='.',
        env={**os.environ, 'QREACH_PARSE_PROFILE': '1', 'QREACH_PARSE_PROFILE_THRESHOLD': '0.05'},
        capture_output=True,
        text=True,
        timeout=20,
    )
    print(result.stdout)
    print(result.stderr)
except subprocess.TimeoutExpired as exc:
    print(exc.stdout or '')
    print(exc.stderr or '')
    print('TIMEOUT')
PY
```

Expected: command prints profile lines before `TIMEOUT`. The last completed line plus the next QASM instruction identify the suspected slow instruction.

- [ ] **Step 3: Run prefix bisect for Grover32 plus**

Run from repo root:

```bash
cd python_pkg
QREACH_PARSE_PROFILE=1 QREACH_PARSE_PROFILE_THRESHOLD=0.05 ../.venv/bin/python workflow_tests/debug_grover32_ccx_pathology.py \
  --qasm benchmark/converted_qasm/single-it-grover32-plus.qasm \
  --max-instructions 180 \
  --timeout-seconds 20 \
  --skip-fixed-post \
  --bisect
```

Expected: records the smallest failing prefix or shows that the failing prefix is beyond 180 instructions.

- [ ] **Step 4: Run prefix bisect for Grover32 zero**

Run from repo root:

```bash
cd python_pkg
QREACH_PARSE_PROFILE=1 QREACH_PARSE_PROFILE_THRESHOLD=0.05 ../.venv/bin/python workflow_tests/debug_grover32_ccx_pathology.py \
  --qasm benchmark/converted_qasm/single-it-grover32-zero.qasm \
  --max-instructions 180 \
  --timeout-seconds 20 \
  --skip-fixed-post \
  --bisect
```

Expected: records whether the zero variant fails at the same or analogous prefix.

- [ ] **Step 5: If useful, save localization notes**

If the evidence is non-trivial, create `docs/agent-handoffs/grover32-timeout-localization-notes.md` with this structure:

```markdown
# Grover32 Timeout Localization Notes

## Commands Run

- `<command>`

## Grover64 Reference

- Result: `<summary>`

## Grover32 Plus

- Last completed profile line: `<line>`
- Smallest failing prefix: `<prefix>`
- Largest adjacent passing prefix: `<prefix>`

## Grover32 Zero

- Last completed profile line: `<line>`
- Smallest failing prefix: `<prefix>`
- Largest adjacent passing prefix: `<prefix>`

## Current Hypothesis

`<hypothesis based on evidence, not pre-existing assumptions>`

## Next Inspection Target

`<Python parser / QOperation construction / lazy post-image / C++ CCX construction>`
```

- [ ] **Step 6: Commit notes if created**

Run only if notes were created:

```bash
git add docs/agent-handoffs/grover32-timeout-localization-notes.md
git commit -m "Record Grover32 timeout localization evidence

Co-Authored-By: Claude <noreply@anthropic.com>"
```

---

### Task 4: Add focused regression after minimization

**Files:**
- Create: `python_pkg/workflow_tests/test_grover32_ccx_pathology.py`
- May import from: `python_pkg/workflow_tests/debug_grover32_ccx_pathology.py`

**Interfaces:**
- Consumes: `make_prefix_circuit(qc, max_instructions)` from the debug script.
- Produces: `test_minimized_grover32_prefix_parses_under_budget()`.

- [ ] **Step 1: Choose prefix from Task 3 evidence**

Set `MINIMIZED_PREFIX_LEN` to the smallest failing prefix discovered in Task 3 after a fix target is known. If Task 3 did not find a failing prefix under 20 seconds, do not create this test yet; instead return to root-cause investigation with finer profiling.

- [ ] **Step 2: Create regression test file**

Create `python_pkg/workflow_tests/test_grover32_ccx_pathology.py` with this content, replacing `MINIMIZED_PREFIX_LEN = 0` with the evidence-backed prefix from Task 3:

```python
from __future__ import annotations

import sys
from pathlib import Path
from time import perf_counter

try:
    import pyqreach
    from qiskit import QuantumCircuit
except ImportError as exc:  # pragma: no cover - script-style regression skip
    print(f"SKIP Grover32 pathology test: {exc}")
    raise SystemExit(0)

PYTHON_PKG = Path(__file__).resolve().parents[1]
if str(PYTHON_PKG) not in sys.path:
    sys.path.insert(0, str(PYTHON_PKG))

from parse_qiskit import parse_qiskit_cir_lazy  # noqa: E402
from workflow_tests.debug_grover32_ccx_pathology import make_prefix_circuit  # noqa: E402

MINIMIZED_PREFIX_LEN = 0
PARSE_BUDGET_SECONDS = 10.0


def test_minimized_grover32_prefix_parses_under_budget() -> None:
    if MINIMIZED_PREFIX_LEN <= 0:
        print("SKIP Grover32 pathology test: minimized prefix length has not been recorded")
        return

    qasm_path = PYTHON_PKG / "benchmark" / "converted_qasm" / "single-it-grover32-plus.qasm"
    qc = QuantumCircuit.from_qasm_file(str(qasm_path))
    prefix = make_prefix_circuit(qc, MINIMIZED_PREFIX_LEN)

    pyqreach.initializeTransitionSystem()
    ts = pyqreach.TransitionSystem()
    start = perf_counter()
    result = parse_qiskit_cir_lazy(
        prefix,
        prefix.num_qubits,
        ts,
        initial_state="0" * prefix.num_qubits,
        return_metadata=True,
    )
    elapsed = perf_counter() - start

    assert result.lazy is True
    assert ts.getLocationNum() == len(prefix.data) + 1
    assert len(result.result_locations) >= 1
    assert elapsed < PARSE_BUDGET_SECONDS, f"parse took {elapsed:.3f}s"


if __name__ == "__main__":
    test_minimized_grover32_prefix_parses_under_budget()
    print("Grover32 pathology test PASSED")
```

- [ ] **Step 3: Run the focused regression**

Run from repo root:

```bash
cd python_pkg
PYTHONPATH=. ../.venv/bin/python workflow_tests/test_grover32_ccx_pathology.py
```

Expected before a fix: fails or exceeds the short budget if `MINIMIZED_PREFIX_LEN` is set to the failing prefix.

Expected after a fix: output includes:

```text
Grover32 pathology test PASSED
```

- [ ] **Step 4: Commit Task 4 after the test is meaningful**

Run only after `MINIMIZED_PREFIX_LEN` has been set to an evidence-backed value:

```bash
git add python_pkg/workflow_tests/test_grover32_ccx_pathology.py
git commit -m "Add Grover32 timeout regression

Co-Authored-By: Claude <noreply@anthropic.com>"
```

---

### Task 5: Decide next root-cause inspection target

**Files:**
- No mandatory modifications.
- Inspect one of:
  - `python_pkg/parse_qiskit.py`
  - `quantum_operation.hpp`
  - `cflobdd/CFLOBDD/matrix1234_node.cpp`

**Interfaces:**
- Consumes: Task 3 evidence.
- Produces: a single stated hypothesis for the next implementation plan or fix task.

- [ ] **Step 1: Classify where time is spent**

Use the evidence from Task 3 to choose one category:

```text
A. Entire instruction is slow, but operation construction vs post-image is unknown.
B. The slow instruction is clearly CCX or a CCX-heavy prefix.
C. The slow instruction is not CCX.
D. The failure depends on accumulated prefix state rather than one instruction.
```

- [ ] **Step 2: If category A, add a second diagnostic plan before fixing**

Do not fix yet. Write a follow-up plan to split profiling around:

```python
op = pyqreach.QOperation(...)
_add_post_and_propagate(ts, cLoc, new_loc, op, lazy_ctx)
```

- [ ] **Step 3: If category B, inspect C++ CCX path**

Start with:

```text
quantum_operation.hpp:772-831
cflobdd/CFLOBDD/matrix1234_node.cpp:907-1238
```

Check whether the slow qubit tuple maps to one of the split cases in `MkCCNOTNode(...)` and whether the ordering path in `quantum_operation.hpp` uses swap-conjugation.

- [ ] **Step 4: If category C, inspect the actual gate path**

Find the matching `op_name` branch in `python_pkg/parse_qiskit.py` and its corresponding C++ `QuantumGateTerm::concretize()` branch before forming a fix hypothesis.

- [ ] **Step 5: If category D, compare neighboring prefixes**

Run the prefix debugger on:

```text
single-it-grover32-plus.qasm at N-1 and N
single-it-grover32-plus-linear.qasm at analogous N
single-it-grover64-plus.qasm at analogous structural point
```

Use only evidence from this comparison to refine the hypothesis.

- [ ] **Step 6: Write the next fix-specific plan**

Do not bundle a C++ fix into this diagnostic plan. Once the root-cause target is known, create a new small implementation plan for that fix with its own failing test and validation commands.
