# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project direction

QReach is a quantum model-checking / reachability-analysis tool for quantum Markov chains and Qiskit programs. The current phase focuses on six priorities:

### 1. Lazy construction optimization (primary C++ focus)

SymTS/QADD symbolic work is **shelved**. The active C++ path is `qts_naive::TransitionSystem` with lazy construction (`parse_qiskit_cir_lazy`), which propagates quantum post-images during parsing instead of deferring everything to `computingFixedPointPost()`.

Key optimizations already applied:
- `knownUnitNorm` normalization skip for singleton post-image with unitary gates.
- `unordered_map` for relation containers, `inQueue` set optimization.
- Post-relation map lookup optimization.

Known performance bottlenecks (see `docs/agent-handoffs/`):
- **grover32-timeout-issue.md / grover32-timeout-localization-notes.md**: Root cause — CFLOBDD `Reduce` operation fragments at level=7 with H^⊗n initial state + nonlinear CCX qubit ordering + non-power-of-2 qubit count. `MatrixMultiplyV4WithInfoTopNode` evaluation loop (`c1_sz × c2_sz`) also contributes.
- **benchpress-float-amplitude-explosion.md**: Root cause — parameterized rotation gates (rz/ry with float params) create incommensurate amplitudes that CFLOBDD cannot merge, causing `2^k` leaf-amplitude explosion. Distinct from Grover32: Grover32's bottleneck is Reduce, benchpress's is the evaluation loop.

Both are deep CFLOBDD DAG issues — initial tactical fixes may help, but comprehensive solutions require architectural changes.

### 2. Workflow automation & model-checking tooling (primary Python focus)

Push toward a more automated, user-friendly quantum model checker:
- Streamline Qiskit→TS→labelling→property-checking pipelines.
- Expand `qctl.py` helpers for common model-checking patterns.
- **Counterexample analysis**: Add automated counterexample-trace interpretation from NuSMV results back into QReach locations/states.

### 3. Benchmark experiments

Scale experiments already run on `benchpress-medium` and `converted_qasm`. Extend to remaining benchmarks under `python_pkg/benchmark/`:
- `dqc_pe`, `dqc_qft`, `pe`, `qft` (phase estimation / QFT families).
- `grover`, `quantum_teleportation`, `superdense_coding`.
- `qrw_*`, `rus_*`, `testbigcliff`.

Primary runner scripts: `python_pkg/workflow_tests/run_qasm_benchmarks_lazy.py`, `run_benchpress_qasm_lazy.py`.

### 4. Performance bottleneck investigation

See `docs/agent-handoffs/` — Grover32 CFLOBDD Reduce fragmentation and benchpress float-amplitude explosion. These are hard problems; initial tactical explorations (profiling, pattern documentation, workaround detection) are worthwhile even without full fixes.

### 5. Tool packaging

Package the project as an installable/distributable tool (pip, CLI entry points, or containerized).

### 6. Additional backends

Currently only CFLOBDD. Explore/integrate alternative quantum decision-diagram or simulation backends to broaden applicability (e.g., circuits with parameterized gates that CFLOBDD handles poorly).

### Deprioritized

- SymTS / QADD symbolic optimization (`transition_system_qadd.hpp`, `qadd.hpp`).
- Large symbolic post-image redesigns.
- New QADD internals unless explicitly requested.

## Current architecture

### C++ semantic layer

- `quantum_operation.hpp`
  - Core quantum terminal algebra over the CFLOBDD backend.
  - Defines `QOperation`, `QuantumGateTerm`, `SingleVecTerm`, subspace comparison, conjunction/disjunction/difference, pre/post image, and Gram-Schmidt normalization.
  - Current workflow-facing additions include simple product-state construction from strings containing `0`, `1`, `+`, `-`, and direct span construction via `SpanQOperations(...)`.
  - `knownUnitNorm` flag, `isNormPreservingGate()`, and normalization-skip path for lazy post-image optimization.
- `transition_system.hpp`
  - Explicit baseline namespace `qts_naive`.
  - `qts_naive::TransitionSystem` stores explicit `Location` objects, relation maps, pre/post adjacency, lower/upper quantum bounds, labels, and classical propositions.
  - **Active optimization target** — relation containers now use `unordered_map`; lazy construction propagates post-images during parsing.
- `transition_system_qadd.hpp` — **SHELVED**: symbolic/QADD namespace `qts`. Maintain regressions only.
- `qadd.hpp` — **SHELVED**: QADD node/unique-table/compute-table layer. Maintain regressions only.
- `cl_proposition.hpp`
  - Classical propositions attached to transition-system locations.
- `cflobdd/`
  - Imported/modified CFLOBDD backend. Performance bottlenecks in `MatrixMultiplyV4WithInfoTopNode` (Reduce + evaluation loop) documented in `docs/agent-handoffs/`.

### Python workflow layer

- `python_pkg/qreach_python_wrapper.cpp`
  - Builds the `pyqreach` extension with pybind11.
  - Exposes `QOperation`, `ClassicalProposition`, explicit `TransitionSystem`/`Location`, symbolic `SymTS`, `pyqreach.span_qops(...)`, and selected `QOperation` methods such as `disjunction(...)`.
- `python_pkg/parse_qiskit.py`
  - Lowers Qiskit `QuantumCircuit` objects into explicit `pyqreach.TransitionSystem` via `parse_qiskit_cir(...)` and `parse_qiskit_cir_lazy(...)` (lazy construction — propagates post-images during parsing).
  - Also contains symbolic parser support via `parse_qiskit_cir_sym(...)`.
  - Handles measurements, resets, classical-register bookkeeping, and Qiskit control flow such as `if_else`, `while_loop`, `for_loop`, and `switch_case`.
  - `ParseResult` metadata with inline marks when `return_metadata=True`.
  - `QREACH_PARSE_PROFILE` env-gated per-instruction profiling.
- `python_pkg/qctl.py`
  - User-facing workflow utilities: labelling, initial-state helpers, snapshot labelling, subspace span helpers, SMV/CTL generation, NuSMV invocation, and result parsing.
- `python_pkg/annotations.py`
  - Built-in annotation keywords and `AnnotationRegistry`.
- `python_pkg/inline_annotations.py`
  - Lightweight Qiskit wrapper `QReachCircuit` and `mark(...)` helper.
- `python_pkg/workflow_tests/`
  - Workflow-oriented API tests and benchmark runners:
    - `run_qasm_benchmarks_lazy.py` — batch runner for all benchmark families.
    - `run_benchpress_qasm_lazy.py` — benchpress-medium batch runner.
    - `qasm_workflow_runner.py` — shared utilities for QASM experiments.
    - `debug_grover32_ccx_pathology.py` — prefix bisect/profiling for Grover32.
    - `test_grover32_lazy_prefix_performance.py` — pytest regression for lazy parse timing.
    - `test_grover.py`, `test_newapi.py`, `test_RUS.py` — workflow API smoke tests.
    - `test_lazy_measurement.py` — lazy mode measurement regression.
    - `test_bv_n14.py`, `test_bv_n14_lazy.py` — Bernstein-Vazirani workflows.
    - `test_simulation_grover.py`, `test_simulation_qft.py` — simulation-style comparisons.
- `python_pkg/benchmark/`
  - QASM benchmark circuits organized by family: `grover/`, `dqc_pe/`, `dqc_qft/`, `pe/`, `qft/`, `benchpress-medium/`, `converted_qasm/`, `quantum_teleportation/`, `superdense_coding/`, `qrw_*.qasm`, `rus_*.qasm`, `testbigcliff.qasm`.

## Build and test commands

### C++ build

From the repository root:

```bash
make all          # builds libqreach.so
make test         # builds test_qreach
./test_qreach 8   # run the C++ benchmark/test executable with an example argument
make clean        # remove libqreach.so and object files
```

The current local macOS setup uses Boost at:

```bash
../BOOST/boost_1_81_0
```

`Makefile` defaults `BOOST_PATH` to that location. Override when needed:

```bash
BOOST_PATH=/path/to/boost_1_81_0 make test
```

### Python environment

Do **not** assume conda or Docker are available. The current machine uses `uv` due to compliance constraints.

From the repository root:

```bash
uv venv
uv pip install -r requirements.txt
```

Current required Qiskit versions:

```text
qiskit==1.4.2
qiskit-aer==0.17.1
```

If upgrading from an old Qiskit 0.x environment, make sure stale `qiskit-terra` is removed; mixed Qiskit 0.x/1.x installs are invalid.

### Python extension build

Python dependencies and bindings are manual, not packaged through `setup.py`/`pyproject.toml`.

From `python_pkg`:

```bash
cd python_pkg
../.venv/bin/python -m invoke build-qreach
../.venv/bin/python -m invoke build-pybind11
```

Useful environment overrides:

```bash
export PYTHON=/path/to/python
export BOOST_PATH=/path/to/boost_1_81_0
export PYTHON_INCLUDE=$(python -c "from sysconfig import get_paths as gp; print(gp()['include'])")
```

`python_pkg/tasks.py` now derives Python/pybind11 include paths from the active Python where possible and has Darwin/Linux linker handling. It should remain compatible with Linux; avoid hardcoding macOS-only paths or flags unless guarded by `platform.system()`.

### Common Python regression/debug scripts

Run from `python_pkg` after building `pyqreach`:

```bash
../.venv/bin/python workflow_tests/test_newapi.py       # annotation API smoke test
../.venv/bin/python test_qiskit_grammar.py          # fuller qCTL/Qiskit workflow example
../.venv/bin/python test_513_qiskit.py              # snapshot-style quantum proposition example
../.venv/bin/python test_parse_qasm.py              # QASM/Grover-style debug workflow
../.venv/bin/python workflow_tests/test_grover.py       # Grover workflow/proposition API smoke test
../.venv/bin/python test_RUS.py                     # RUS Qiskit/debugging workflow test
../.venv/bin/python test_symts_minimal.py           # minimal explicit-vs-symbolic post check
../.venv/bin/python test_symts_control_flow.py      # reduced symbolic control-flow checks
../.venv/bin/python test_symts_steane_regression.py # bounded Steane symbolic regression
```

Useful profiling/debug commands:

```bash
TS_PROFILE=1 ./test_qreach 8
cd python_pkg
TS_PROFILE=1 ../.venv/bin/python debug_steane_measure_only.py 1
TS_PROFILE=1 ../.venv/bin/python test_symts_steane_post.py sym 1 full
TS_MAX_POST_ITER=5 ../.venv/bin/python test_symts_minimal.py
QOP_PROFILE=1 ../.venv/bin/python <script.py>
```

NuSMV integration in `python_pkg/qctl.py` defaults to:

```text
../NuSMV-2.7.0-linux64/bin/NuSMV
```

Pass `nusmv_path` to `modelChecking(...)` or ensure `NuSMV` is available there/on `PATH`.

## Current workflow APIs to preserve

The following APIs are part of the current practical workflow surface. When editing them, preserve backward compatibility where possible.

### Initial state helpers in `qctl.py`

```python
quantum_state(bitstring: str) -> pyqreach.QOperation
set_initial_state(ts, bitstring: str, loc: int | None = None) -> pyqreach.QOperation
set_initial_operation(ts, op: pyqreach.QOperation, loc: int | None = None) -> pyqreach.QOperation
set_zero_initial_state(ts, qnum: int | None = None, loc: int | None = None) -> pyqreach.QOperation
```

`quantum_state(...)` delegates to `pyqreach.QOperation([bitstring])`. The string syntax currently supports simple product states over `0`, `1`, `+`, `-`. `set_initial_operation(...)` is the operation-level dual of `set_initial_state(...)` for already-constructed propositions/subspaces.

### Simple product-state `QOperation` strings

`QOperation(std::vector<std::string>)` supports strings such as:

```python
pyqreach.QOperation(["000"])
pyqreach.QOperation(["+++"])
pyqreach.QOperation(["+01-"])
pyqreach.QOperation(["1-0+"])
```

Semantics:

```text
0 = |0>
1 = |1>
+ = H|0> = (|0> + |1>) / sqrt(2)
- = H|1> = (|0> - |1>) / sqrt(2)
```

Do not extend this syntax to arbitrary amplitude expressions without a deliberate design; general amplitudes should use a separate builder/API.

### Multi-dimensional subspace span

Use direct span APIs instead of temporary transition systems:

```python
from qctl import quantum_state, span_states, span_qops

span_states(["00", "11"])
span_qops([quantum_state("+0"), quantum_state("-1")])
```

C++/binding surface:

```python
pyqreach.span_qops([...])
op_a.disjunction(op_b)
```

Prefer `qctl.span_states(...)` / `qctl.span_qops(...)` in user-facing Python code.

This replaces older workaround patterns such as creating a temporary `TransitionSystem`, adding identity transitions into a sink, setting annotations, and calling `computingFixedPointPost()` only to force Gram-Schmidt/span construction.

### Quantum annotations on special locations

`qctl.py` also has helpers for setting quantum annotations on common location sets:

```python
annotate_leaf_state(ts, bitstring: str, *, loc_list=None) -> list[int]
annotate_leaf_operation(ts, op: pyqreach.QOperation, *, loc_list=None) -> list[int]
marker_locations(parse_result, marker: str) -> list[int]
annotate_marker_state(ts, parse_result, marker: str, bitstring: str) -> list[int]
annotate_marker_operation(ts, parse_result, marker: str, op: pyqreach.QOperation) -> list[int]
```

Use these when a workflow needs to attach a quantum state/proposition to leaf locations or to the locations recorded for an inline marker. They set transition-system annotations, not classical AP labels, and return the affected location IDs.

### Built-in annotation keywords

Implemented in `python_pkg/annotations.py` and re-exported by `qctl.py`:

```python
annotate(ts, names=None, *, loc_list=None)
annotate_where(ts, name, predicate, *, loc_list=None)
annotate_identifier(ts, name, pattern, *, loc_list=None)
annotate_classical(ts, name, bit_patterns, *, bit_indices=None, loc_list=None)
default_registry()
```

Current built-ins:

- `reached`: must follow `qctl.tsLabellingDefault`, i.e. `ts.printDims(loc)[1] > 0` for explicit TS.
- `valid`: alias for `reached`.
- `leaf`
- `init`
- `loop`
- `deadend`
- `measured`
- `branch`

Important: do not implement explicit `reached` using `SymTS.locationHasNonZeroAnnotation(...)`; that method belongs to the symbolic transition system and was previously corrected.

### Inline Qiskit marks and parser metadata

`python_pkg/inline_annotations.py` provides:

```python
QReachCircuit
QReachCircuit.from_qasm_file(path)
QReachCircuit.from_qasm_str(qasm_str)
mark(qc, name)
```

`QReachCircuit.from_qasm_file(...)` mirrors Qiskit's `QuantumCircuit.from_qasm_file(...)` for ordinary OpenQASM 2 files without QReach-specific inline annotations. The returned wrapper can still receive `.mark(name)` calls before parsing.

Marks are labelled Qiskit barriers with prefix `qreach:mark:` and are semantic no-ops for transition-system construction.

`parse_qiskit_cir(..., return_metadata=True)` returns a `ParseResult` with:

```python
result_locations
markers
instruction_locations
identifier_index
```

Backward compatibility: ordinary calls without `return_metadata=True` should continue returning the final location list.

### Snapshot labelling

`qctl.label_snapshot(...)` obtains a `QOperation` from a marked location, then follows `tsLabelling` semantics: scan all locations or `locList`, and label every location satisfying that proposition.

```python
label_snapshot(
    ts,
    parse_result,
    marker,
    label=None,
    *,
    bound="lower",
    merge="single",
    locList=None,
)
```

Supported `bound` values:

```text
lower, upper, annotation
```

Supported `merge` values for multi-location marks:

```text
single       # require exactly one location
first        # use first location
join         # span all selected location bounds using span_qops
per_location # create separate labels per marked location
```

## Typical end-to-end flow

1. Build `libqreach.so` and `python_pkg/pyqreach...so`.
2. Create or load a Qiskit circuit.
3. Optionally use `QReachCircuit` / `mark(...)` for inline debug points.
4. Lower with `parse_qiskit_cir(...)` for explicit transition-system work.
5. Set initial state with `set_zero_initial_state(...)` or `set_initial_state(...)`.
6. Run `ts.computingFixedPointPost()`.
7. Add labels using `annotate(...)`, `tsLabelling(...)`, `label_snapshot(...)`, `span_states(...)`, and/or `span_qops(...)`.
8. Generate SMV/CTL and run `modelChecking(...)`, or inspect labels/counterexamples directly.
9. Use symbolic `parse_qiskit_cir_sym(...)` / `SymTS` only when maintaining symbolic regressions or explicitly requested.

## Notes for future changes

- **Primary focus**: lazy construction optimization (C++) and workflow automation + counterexample analysis (Python).
- SymTS / QADD / symbolic post-image code is **shelved** — maintain regressions, do not extend. Do not add new QADD internals or SymTS features without explicit user request.
- The main workflow bottleneck is often obtaining the desired `QOperation`, not merely labelling once it exists. Prioritize APIs that make useful quantum propositions/subspaces easy to construct.
- CFLOBDD performance bottlenecks are documented in `docs/agent-handoffs/`. Tactical fixes (profiling, workaround detection, pattern documentation) are valuable; architectural changes (symbolic leaves, DAG sharing improvements) are high-risk and should be explicitly approved.
- When modifying parser behavior, check explicit parser behavior (both lazy and non-lazy) first. Check symbolic helper branches only if the symbolic regression surface is affected.
- Benchmark runs use `workflow_tests/run_qasm_benchmarks_lazy.py` and `run_benchpress_qasm_lazy.py`. When adding new benchmarks, integrate into these runner scripts.
- Some standalone `pyqreach.QOperation(["..."])` calls may assert if the CFLOBDD/QReach environment has not been initialized. In normal parser/transition-system workflows this is usually avoided. Be cautious when writing tiny standalone smoke tests.
- IDE diagnostics may report missing Boost/pybind11 includes or `BIG_COMPLEX_FLOAT` fallback issues if clangd lacks the Makefile include paths. Trust actual `make test` / pybind11 builds over unconfigured clangd diagnostics.
- Keep README and CLAUDE.md aligned when adding user-facing workflow APIs.
