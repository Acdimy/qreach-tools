# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project direction

QReach is a quantum model-checking / reachability-analysis tool for quantum Markov chains and Qiskit programs. The original project direction emphasized CFLOBDD/QADD symbolic-DD optimization, but the current near-term goal is practical Qiskit program debugging and workflow tooling.

Prioritize:

- Practical Python APIs for constructing useful `QOperation` propositions and subspaces.
- Qiskit parsing, inline marks, location metadata, labelling, qCTL/SMV generation, and regression workflows.
- Explicit transition-system usability and pragmatic optimization around `qts_naive::TransitionSystem`.
- Small, testable workflow improvements that reduce boilerplate in scripts such as `python_pkg/test_qiskit_grammar.py`, `python_pkg/test_513_qiskit.py`, `python_pkg/test_parse_qasm.py`, and the newer `python_pkg/workflow_tests/` scripts.

Deprioritize:

- New QADD symbolic optimization work from `plan.md` unless the user explicitly asks for it.
- Large symbolic post-image redesigns before the explicit/Qiskit debugging workflow is comfortable.

`transition_system_qadd.hpp`, `qadd.hpp`, and symbolic post-image code remain useful references/regression baselines, but do not assume they are the preferred implementation target.

## Current architecture

### C++ semantic layer

- `quantum_operation.hpp`
  - Core quantum terminal algebra over the CFLOBDD backend.
  - Defines `QOperation`, `QuantumGateTerm`, `SingleVecTerm`, subspace comparison, conjunction/disjunction/difference, pre/post image, and Gram-Schmidt normalization.
  - Current workflow-facing additions include simple product-state construction from strings containing `0`, `1`, `+`, `-`, and direct span construction via `SpanQOperations(...)`.
- `transition_system.hpp`
  - Explicit baseline namespace `qts_naive`.
  - `qts_naive::TransitionSystem` stores explicit `Location` objects, relation maps, pre/post adjacency, lower/upper quantum bounds, labels, and classical propositions.
  - Treat this as the semantic reference and the preferred near-term optimization target.
- `transition_system_qadd.hpp`
  - Symbolic/QADD namespace `qts`.
  - Provides `qts::TransitionSystem` / Python `SymTS` and QADD-backed relation/annotation handling.
  - Maintain regressions, but do not extend this route by default.
- `qadd.hpp`
  - QADD node/unique-table/compute-table layer and generic `Apply` operations.
- `cl_proposition.hpp`
  - Classical propositions attached to transition-system locations.
- `cflobdd/`
  - Imported/modified CFLOBDD backend.

### Python workflow layer

- `python_pkg/qreach_python_wrapper.cpp`
  - Builds the `pyqreach` extension with pybind11.
  - Exposes `QOperation`, `ClassicalProposition`, explicit `TransitionSystem`/`Location`, symbolic `SymTS`, `pyqreach.span_qops(...)`, and selected `QOperation` methods such as `disjunction(...)`.
- `python_pkg/parse_qiskit.py`
  - Lowers Qiskit `QuantumCircuit` objects into explicit `pyqreach.TransitionSystem` via `parse_qiskit_cir(...)`.
  - Also contains symbolic parser support via `parse_qiskit_cir_sym(...)`.
  - Handles measurements, resets, classical-register bookkeeping, and Qiskit control flow such as `if_else`, `while_loop`, `for_loop`, and `switch_case`.
  - Current explicit parser can return `ParseResult` metadata with inline marks when `return_metadata=True`.
- `python_pkg/qctl.py`
  - User-facing workflow utilities: labelling, initial-state helpers, snapshot labelling, subspace span helpers, SMV/CTL generation, NuSMV invocation, and result parsing.
- `python_pkg/annotations.py`
  - Built-in annotation keywords and `AnnotationRegistry`.
- `python_pkg/inline_annotations.py`
  - Lightweight Qiskit wrapper `QReachCircuit` and `mark(...)` helper.
- `python_pkg/test_*.py`
  - Mostly executable regression/debug scripts, not pytest-style unit tests. Prefer running the script closest to the touched workflow.

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

- Prefer extending explicit/Qiskit/qCTL workflows before adding new QADD internals.
- The main workflow bottleneck is often obtaining the desired `QOperation`, not merely labelling once it exists. Prioritize APIs that make useful quantum propositions/subspaces easy to construct.
- `python_pkg/test_parse_qasm.py::generate_debug_info` and Grover-style examples are important references for proposition-builder improvements.
- When modifying parser behavior, check explicit parser behavior first. Check symbolic helper branches only if the symbolic regression surface is affected.
- Python-side explicit-vs-symbolic comparisons can be sensitive to repeated global transition-system initialization; use subprocess isolation when comparing both in one workflow.
- Some standalone `pyqreach.QOperation(["..."])` calls may assert if the CFLOBDD/QReach environment has not been initialized. In normal parser/transition-system workflows this is usually avoided. Be cautious when writing tiny standalone smoke tests.
- IDE diagnostics may report missing Boost/pybind11 includes or `BIG_COMPLEX_FLOAT` fallback issues if clangd lacks the Makefile include paths. Trust actual `make test` / pybind11 builds over unconfigured clangd diagnostics.
- Keep README and CLAUDE.md aligned when adding user-facing workflow APIs.
