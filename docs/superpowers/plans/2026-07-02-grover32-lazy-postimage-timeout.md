# Grover32 Lazy Post-Image Timeout Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Make lazy parsing of the pathological Grover32 plus prefix complete quickly by optimizing the accumulated-state post-image path used by repeated X/CCX operations.

**Architecture:** The diagnostic evidence points to accumulated `QOperation`/CFLOBDD state, not Python dispatch or CCX count alone. Add a focused regression harness that times the known failing prefix and preserves fast controls, then optimize the C++ post-image/gate-application path without changing parser semantics.

**Tech Stack:** C++17, CFLOBDD matrix/vector backend, pybind11 `pyqreach`, Python 3 via `uv` venv, Qiskit 1.4.2.

## Global Constraints

- Work in repository `/Users/ftdac/thu/qreach-tools`.
- Do not change user-facing lazy parser semantics; `parse_qiskit_cir_lazy(..., return_metadata=True)` must still return `ParseResult` with the same location/result-location counts for these prefixes.
- Preserve existing workflow APIs listed in `/Users/ftdac/thu/qreach-tools/CLAUDE.md`.
- Prefer explicit/Qiskit/qCTL workflow improvements over symbolic/QADD work.
- Keep the validation focused on the known pathology and nearby controls; do not bundle broad CFLOBDD rewrites.

---

### Task 1: Add a focused lazy-prefix performance regression

**Files:**
- Create: `/Users/ftdac/thu/qreach-tools/python_pkg/workflow_tests/test_grover32_lazy_prefix_performance.py`
- Uses existing helper: `/Users/ftdac/thu/qreach-tools/python_pkg/workflow_tests/debug_grover32_ccx_pathology.py`

**Interfaces:**
- Consumes: `run_prefix_once(qasm_path: pathlib.Path, *, max_instructions: int | None, initial_state: str | None, run_fixed_post: bool) -> PrefixResult`
- Produces: pytest tests that assert the pathological prefix completes under a conservative budget and that control prefixes remain fast.

- [ ] **Step 1: Write the failing regression test**

Create `/Users/ftdac/thu/qreach-tools/python_pkg/workflow_tests/test_grover32_lazy_prefix_performance.py` with this exact content:

```python
from __future__ import annotations

from pathlib import Path

from debug_grover32_ccx_pathology import PYTHON_PKG, run_prefix_once


def _qasm(name: str) -> Path:
    return PYTHON_PKG / "benchmark" / "converted_qasm" / name


def test_grover32_plus_prefix_102_lazy_parse_under_budget() -> None:
    result = run_prefix_once(
        _qasm("single-it-grover32-plus.qasm"),
        max_instructions=102,
        initial_state=None,
        run_fixed_post=False,
    )

    assert result.status == "ok"
    assert result.num_locations == 103
    assert result.num_result_locations == 1
    assert result.time_parse != ""
    assert float(result.time_parse) < 5.0


def test_lazy_parse_controls_stay_fast_near_same_structure() -> None:
    linear = run_prefix_once(
        _qasm("single-it-grover32-plus-linear.qasm"),
        max_instructions=102,
        initial_state=None,
        run_fixed_post=False,
    )
    grover64 = run_prefix_once(
        _qasm("single-it-grover64-plus.qasm"),
        max_instructions=133,
        initial_state=None,
        run_fixed_post=False,
    )

    assert linear.status == "ok"
    assert linear.num_locations == 103
    assert linear.num_result_locations == 1
    assert linear.time_parse != ""
    assert float(linear.time_parse) < 1.0

    assert grover64.status == "ok"
    assert grover64.num_locations == 134
    assert grover64.num_result_locations == 1
    assert grover64.time_parse != ""
    assert float(grover64.time_parse) < 1.0
```

- [ ] **Step 2: Run the test to verify it fails on the pathology**

Run:

```bash
cd /Users/ftdac/thu/qreach-tools/python_pkg
../.venv/bin/python -m pytest workflow_tests/test_grover32_lazy_prefix_performance.py -q
```

Expected before the fix: `test_grover32_plus_prefix_102_lazy_parse_under_budget` fails because `time_parse` is about 22-23 seconds, while `test_lazy_parse_controls_stay_fast_near_same_structure` passes.

- [ ] **Step 3: Commit the failing test**

```bash
cd /Users/ftdac/thu/qreach-tools
git add python_pkg/workflow_tests/test_grover32_lazy_prefix_performance.py
git commit -m "test: add Grover32 lazy prefix performance regression

Co-Authored-By: Claude <noreply@anthropic.com>"
```

---

### Task 2: Split C++ profiling inside post-image gate application

**Files:**
- Modify: `/Users/ftdac/thu/qreach-tools/quantum_operation.hpp:959-978`
- Modify: `/Users/ftdac/thu/qreach-tools/quantum_operation.hpp:1815-1868`
- Inspect if needed: `/Users/ftdac/thu/qreach-tools/python_pkg/parse_qiskit.py:133-141`

**Interfaces:**
- Consumes: `QOperation::postImage(const QOperation& other) const` and `SingleVecTerm::applyGate(const QuantumGateTerm& other, bool direction) const`
- Produces: optional environment-gated timing output that separates gate concretization from matrix-vector multiply for the failing regression.

- [ ] **Step 1: Add temporary environment-gated timing around `SingleVecTerm::applyGate`**

In `/Users/ftdac/thu/qreach-tools/quantum_operation.hpp`, include timing guarded by an environment variable such as `QREACH_GATE_PROFILE`. The profile must report at least gate name, qubit indexes, `qNum`, concretize time, multiply time, and direction. Keep it off unless the environment variable is set.

- [ ] **Step 2: Rebuild bindings**

Run:

```bash
cd /Users/ftdac/thu/qreach-tools/python_pkg
../.venv/bin/python -m invoke build-qreach
../.venv/bin/python -m invoke build-pybind11
```

Expected: both commands complete successfully.

- [ ] **Step 3: Run profiling on only the failing prefix**

Run:

```bash
cd /Users/ftdac/thu/qreach-tools/python_pkg
QREACH_GATE_PROFILE=1 QREACH_PARSE_PROFILE=1 QREACH_PARSE_PROFILE_THRESHOLD=0.5 \
  ../.venv/bin/python workflow_tests/debug_grover32_ccx_pathology.py \
  --qasm benchmark/converted_qasm/single-it-grover32-plus.qasm \
  --max-instructions 102 \
  --skip-fixed-post
```

Expected: output identifies whether the multi-second cost is dominated by `QuantumGateTerm::concretize()` or by `Matrix1234ComplexFloatBoost::MatrixMultiplyV4WithInfo(operand, this->content)`.

- [ ] **Step 4: Do not commit profiling-only code unless it is cleanly gated and useful**

If the instrumentation is temporary, remove it before the optimization commit. If it is retained, keep it environment-gated and document it in the commit message.

---

### Task 3: Optimize the confirmed hot path without changing parser behavior

**Files:**
- Modify one or more of:
  - `/Users/ftdac/thu/qreach-tools/quantum_operation.hpp:560-836`
  - `/Users/ftdac/thu/qreach-tools/quantum_operation.hpp:959-978`
  - `/Users/ftdac/thu/qreach-tools/quantum_operation.hpp:1815-1868`
  - `/Users/ftdac/thu/qreach-tools/cflobdd/CFLOBDD/matrix1234_node.cpp:907-1238`
- Test: `/Users/ftdac/thu/qreach-tools/python_pkg/workflow_tests/test_grover32_lazy_prefix_performance.py`

**Interfaces:**
- Consumes: profiling result from Task 2.
- Produces: unchanged `QOperation::postImage(...)` semantics with faster execution on the Grover32 plus prefix.

- [ ] **Step 1: Choose exactly one optimization target from the profiling result**

Use this decision rule:

```text
If concretize time dominates: optimize or cache gate construction for repeated X/CCX terms.
If matrix multiply time dominates: optimize the matrix-vector application path for singleton normalized subspaces.
If GramSchmidt/normalization dominates: optimize the singleton case in postImage without broad subspace rewrites.
```

- [ ] **Step 2: Implement the smallest semantic-preserving change**

Apply one focused change only. Do not alter `python_pkg/parse_qiskit.py` unless profiling proves Python-level dispatch is material, because current evidence shows the slow instruction has `locations_in=1`, `locations_out=1`, and spends time after the `QOperation` object is handed to `_add_post_and_propagate(...)`.

- [ ] **Step 3: Rebuild bindings**

Run:

```bash
cd /Users/ftdac/thu/qreach-tools/python_pkg
../.venv/bin/python -m invoke build-qreach
../.venv/bin/python -m invoke build-pybind11
```

Expected: both commands complete successfully.

- [ ] **Step 4: Verify the focused regression passes**

Run:

```bash
cd /Users/ftdac/thu/qreach-tools/python_pkg
../.venv/bin/python -m pytest workflow_tests/test_grover32_lazy_prefix_performance.py -q
```

Expected after the fix: both tests pass; the Grover32 plus prefix 102 parse time is under 5 seconds.

- [ ] **Step 5: Verify existing workflow smoke tests still pass**

Run:

```bash
cd /Users/ftdac/thu/qreach-tools/python_pkg
../.venv/bin/python workflow_tests/test_grover.py
../.venv/bin/python workflow_tests/test_newapi.py
```

Expected: both scripts complete without errors.

- [ ] **Step 6: Commit the optimization**

```bash
cd /Users/ftdac/thu/qreach-tools
git add quantum_operation.hpp cflobdd/CFLOBDD/matrix1234_node.cpp python_pkg/workflow_tests/test_grover32_lazy_prefix_performance.py
git commit -m "fix: speed up Grover32 lazy post-image path

Co-Authored-By: Claude <noreply@anthropic.com>"
```

---

## Self-Review

- Spec coverage: The plan starts from the known failing Grover32 plus prefix 102, includes fast controls for the linear and Grover64 variants, and targets only the C++ post-image/gate-application path indicated by Task 5 evidence.
- Placeholder scan: No placeholder implementation steps are included; the test file content and validation commands are explicit.
- Type consistency: The plan uses existing `run_prefix_once(...) -> PrefixResult`, `QOperation::postImage(...)`, and `SingleVecTerm::applyGate(...)` interfaces as currently present in the repository.
