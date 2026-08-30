# QReach: A Reachability Analysis Tool for Quantum Markov Chains

## Introduction

QReach is a quantum model-checking / reachability-analysis tool for quantum Markov chains and Qiskit programs.

Quantum Markov chains (QMCs) are a fundamental model for quantum information processing systems, including quantum communication protocols and quantum-program semantics. Reachability analysis is a core step in classical model checking, and it also appears in quantum communication, quantum control, and termination analysis of quantum programs.

The original QReach prototype uses Context-Free-Language Ordered Binary Decision Diagrams (CFLOBDDs) as a quantum decision-diagram backend. Current development also emphasizes practical Qiskit program debugging: parsing Qiskit circuits, constructing transition systems, labelling locations with classical/quantum propositions, and checking qCTL/CTL properties through NuSMV.

Some previous works are referred to in this project: [Quasimodo simulator](https://github.com/trishullab/Quasimodo), [Original CFLOBDD](https://github.com/trishullab/cflobdd). QReach builds on these programs with important modifications.

A usable QReach artifact is available at https://doi.org/10.5281/zenodo.10931240.

## Build and installation

This repository is still under active development. Python bindings are built manually rather than through a packaged `setup.py`/`pyproject.toml` workflow.

### C++ build

From the repository root:

```bash
make all
make test
./test_qreach 8
```

If Boost is not in the default path, set `BOOST_PATH`:

```bash
BOOST_PATH=/path/to/boost_1_81_0 make test
```

### Python dependencies

A typical local workflow uses a virtual environment with the dependencies in `requirements.txt`:

```bash
uv venv
uv pip install -r requirements.txt
```

The current tested Qiskit dependency versions are:

```text
qiskit==1.4.2
qiskit-aer==0.17.1
```

### Build the Python extension

From `python_pkg`:

```bash
cd python_pkg
../.venv/bin/python -m invoke build-qreach
../.venv/bin/python -m invoke build-pybind11
```

Useful environment overrides:

```bash
export BOOST_PATH=/path/to/boost_1_81_0
export PYTHON=/path/to/python
export PYTHON_INCLUDE=$(python -c "from sysconfig import get_paths as gp; print(gp()['include'])")
```

## Basic Qiskit workflow

A typical workflow is:

1. Build `libqreach.so` and the `pyqreach` Python extension.
2. Construct or load a Qiskit `QuantumCircuit`.
3. Parse it into a QReach transition system.
4. Set the initial quantum state.
5. Run fixed-point reachability.
6. Add classical/quantum labels.
7. Check CTL/qCTL properties through NuSMV or inspect labels directly.

Minimal example:

```python
import pyqreach
from qiskit import QuantumCircuit
from qreach.parse_qiskit import parse_qiskit_cir
from qreach.qctl import set_zero_initial_state, annotate

qc = QuantumCircuit(2, 0)
qc.h(0)
qc.cx(0, 1)

ts = pyqreach.TransitionSystem()
result_locs = parse_qiskit_cir(qc, qc.num_qubits, ts)

set_zero_initial_state(ts)
ts.computingFixedPointPost()

labels = annotate(ts, ["reached", "leaf"])
print(labels)
```

## Convenient Python APIs

Recent workflow work adds several convenience APIs around common QReach operations. Most are available from `python_pkg/qreach/qctl.py`.

### Initial quantum-state helpers

Instead of manually writing:

```python
opinit = pyqreach.QOperation(["000"])
ts.setAnnotation([[0, opinit]])
```

use:

```python
from qreach.qctl import quantum_state, set_initial_state, set_zero_initial_state

op = quantum_state("000")
set_initial_state(ts, "000")
set_zero_initial_state(ts)
```

APIs:

```python
quantum_state(bitstring: str) -> pyqreach.QOperation
set_initial_state(ts, bitstring: str, loc: int | None = None) -> pyqreach.QOperation
set_initial_operation(ts, op: pyqreach.QOperation, loc: int | None = None) -> pyqreach.QOperation
set_zero_initial_state(ts, qnum: int | None = None, loc: int | None = None) -> pyqreach.QOperation
```

`set_initial_state` defaults to `ts.getInitLocation()`. `set_initial_operation` is the operation-level dual for callers that already constructed a `QOperation`, for example with `span_states`, `span_qops`, or `snapshot_operation`. `set_zero_initial_state` infers the number of qubits from the initial location when possible.

### Simple product-state `QOperation` strings

`pyqreach.QOperation([state])` and `qctl.quantum_state(state)` now support simple product-state strings over:

```text
0, 1, +, -
```

where:

```text
0 = |0>
1 = |1>
+ = H|0> = (|0> + |1>) / sqrt(2)
- = H|1> = (|0> - |1>) / sqrt(2)
```

Examples:

```python
quantum_state("000")
quantum_state("+++")
quantum_state("+01-")
quantum_state("1-0+")
```

For example:

```python
quantum_state("+0")
```

represents:

```text
|+0> = (|00> + |10>) / sqrt(2)
```

This interface is intended for convenient construction of common one-dimensional quantum propositions. General amplitude expressions are not yet part of this string syntax.

### Multi-dimensional subspace span

Use `span_states` or `span_qops` to construct the normalized span of several one-dimensional subspaces.

```python
from qreach.qctl import quantum_state, span_states, span_qops

subspace = span_states(["00", "11"])

op0 = quantum_state("+0")
op1 = quantum_state("-1")
subspace = span_qops([op0, op1])
```

APIs:

```python
span_states(states) -> pyqreach.QOperation
span_qops(ops) -> pyqreach.QOperation
```

These functions call the C++ `QOperation` span/Gram-Schmidt path directly. They replace older workflows that built a temporary transition system and called `computingFixedPointPost()` only to force a union/span of quantum subspaces.

For example, code like:

```python
ts_temp = pyqreach.TransitionSystem(False)
# add locations, identity transitions, annotations ...
ts_temp.computingFixedPointPost()
grover_final = ts_temp.Locations[2].lowerBound
```

can often become:

```python
grover_final = span_qops([grover_init, grover_good])
```

The underlying C++ binding also exposes:

```python
pyqreach.span_qops([...])
op_a.disjunction(op_b)
```

For workflow code, prefer `qctl.span_states` and `qctl.span_qops`.

### Quantum annotations on special locations

In addition to classical proposition labels, qctl provides helpers for setting quantum annotations on common location sets.

Leaf locations:

```python
from qreach.qctl import annotate_leaf_state, annotate_leaf_operation, span_states

annotate_leaf_state(ts, "11")
annotate_leaf_operation(ts, span_states(["00", "11"]))
```

Marker locations:

```python
from qreach.qctl import annotate_marker_state, annotate_marker_operation, quantum_state

annotate_marker_state(ts, parse_result, "after_h", "+0")
annotate_marker_operation(ts, parse_result, "after_h", quantum_state("+0"))
```

APIs:

```python
leaf_locations(ts, loc_list=None) -> list[int]
annotate_leaf_state(ts, bitstring: str, *, loc_list=None) -> list[int]
annotate_leaf_operation(ts, op: pyqreach.QOperation, *, loc_list=None) -> list[int]
marker_locations(parse_result, marker: str) -> list[int]
annotate_marker_state(ts, parse_result, marker: str, bitstring: str) -> list[int]
annotate_marker_operation(ts, parse_result, marker: str, op: pyqreach.QOperation) -> list[int]
```

These functions return the locations whose quantum annotations were set.

### Built-in annotation keywords

Built-in annotation keywords label transition-system locations with common classical propositions. They are implemented through `AnnotationRegistry` in `python_pkg/qreach/annotations.py` and re-exported by `qctl.py`.

```python
from qreach.qctl import annotate, default_registry

labels = annotate(ts, ["reached", "leaf"])
print(labels)
```

Currently supported built-ins include:

| Keyword | Meaning |
| --- | --- |
| `reached` | Location has a non-zero fixed-point annotation. This follows `ts.printDims(loc)[1] > 0`, matching `qctl.tsLabellingDefault`. |
| `valid` | Alias for `reached`, kept for compatibility with older scripts. |
| `leaf` | Terminal/leaf transition-system location. |
| `init` | Initial transition-system location. |
| `loop` | Location whose identifier indicates a while-loop region. Currently heuristic-based. |
| `deadend` | Graph location with no outgoing transition. |
| `measured` | Location inferred from an incoming `meas0`/`meas1` transition. |
| `branch` | Location with multiple outgoing graph edges. |

Useful APIs:

```python
annotate(ts, names=None, *, loc_list=None)
annotate_where(ts, name, predicate, *, loc_list=None)
annotate_identifier(ts, name, pattern, *, loc_list=None)
annotate_classical(ts, name, bit_patterns, *, bit_indices=None, loc_list=None)
default_registry()
```

Example:

```python
from qreach.qctl import annotate, annotate_classical, annotate_identifier

annotate(ts, ["reached", "leaf"])
annotate_classical(ts, "accept", "10")
annotate_identifier(ts, "loop_body", "*W*")
```

### Inline Qiskit marks

`QReachCircuit` is a lightweight wrapper around Qiskit's `QuantumCircuit`. It delegates normal Qiskit syntax to the underlying circuit and adds `.mark(name)`.

```python
from qreach.inline_annotations import QReachCircuit
from qreach.parse_qiskit import parse_qiskit_cir

qc = QReachCircuit(2, 0)
qc.h(0)
qc.mark("after_h")
qc.cx(0, 1)

ts = pyqreach.TransitionSystem()
parse_result = parse_qiskit_cir(
    qc.unwrap(),
    qc.num_qubits,
    ts,
    return_metadata=True,
)

print(parse_result.markers)
```

Marks are represented as labelled Qiskit barriers, so they are semantic no-ops for the transition-system construction. The parser records the current transition-system locations for each mark.

You can also load a normal OpenQASM 2 file through the wrapper:

```python
qc = QReachCircuit.from_qasm_file("benchmark/grover/grover_5.qasm")
qc.mark("after_load")
```

`QReachCircuit.from_qasm_file(...)` mirrors Qiskit's `QuantumCircuit.from_qasm_file(...)`; the file is assumed to be ordinary OpenQASM without QReach-specific inline annotations.

You can also mark a normal `QuantumCircuit`:

```python
from qreach.inline_annotations import mark

mark(qc, "after_oracle")
```

When `return_metadata=True`, `parse_qiskit_cir` returns a `ParseResult` instead of only a list of final locations. `ParseResult` includes:

```python
parse_result.result_locations
parse_result.markers
parse_result.instruction_locations
parse_result.identifier_index
```

For backward compatibility, ordinary calls still return the final location list:

```python
result_locs = parse_qiskit_cir(qc, qc.num_qubits, ts)
```

### Snapshot labelling from marks

A common debugging task is to label every location whose quantum state satisfies the quantum proposition represented by a marked program point.

Manual version:

```python
prop = ts.Locations[17].lowerBound
tsLabelling(ts, prop, "sp17")
```

Convenience version:

```python
from qreach.qctl import label_snapshot

label_snapshot(ts, parse_result, "after_h", "sp_after_h")
```

This follows the existing `tsLabelling` semantics:

1. obtain a `QOperation` from the marked location's `lowerBound`, `upperBound`, or annotation;
2. scan all transition-system locations, or a supplied `locList`;
3. label every location satisfying that quantum proposition.

API:

```python
label_snapshot(
    ts,
    parse_result,
    marker: str,
    label: str | None = None,
    *,
    bound: str = "lower",
    merge: str = "single",
    locList: list | None = None,
) -> dict[str, list[int]]
```

`bound` may be:

```text
lower, upper, annotation
```

`merge` controls marks that correspond to multiple locations:

| Merge mode | Meaning |
| --- | --- |
| `single` | Require exactly one marked location. This is the default and safest mode. |
| `first` | Use the first marked location. Useful for quick debugging. |
| `join` | Construct one quantum proposition spanning all marked location bounds. Uses `span_qops`. |
| `per_location` | Label separately for each marked location, e.g. `label_17`, `label_18`. |

Example:

```python
from qreach.qctl import set_zero_initial_state, label_snapshot

set_zero_initial_state(ts)
ts.computingFixedPointPost()

labelled = label_snapshot(ts, parse_result, "after_h", "sp_after_h")
print(labelled)
```

### Full small example

```python
import pyqreach
from qreach.inline_annotations import QReachCircuit
from qreach.parse_qiskit import parse_qiskit_cir
from qreach.qctl import (
    set_zero_initial_state,
    annotate,
    label_snapshot,
    span_states,
    tsLabelling,
)

qc = QReachCircuit(2, 0)
qc.h(0)
qc.mark("after_h")
qc.cx(0, 1)

ts = pyqreach.TransitionSystem()
parse_result = parse_qiskit_cir(qc.unwrap(), qc.num_qubits, ts, return_metadata=True)

set_zero_initial_state(ts)
ts.computingFixedPointPost()

annotate(ts, ["reached", "leaf"])
label_snapshot(ts, parse_result, "after_h", "sp_after_h")

bell_subspace = span_states(["00", "11"])
tsLabelling(ts, bell_subspace, "bell_span")
```

## Test and regression scripts

Tests and scripts are organized under `python_pkg/`:

- `python_pkg/workflow_tests/` — real unit/regression tests (with assertions) and the
  benchmark runner infrastructure. Symbolic SymTS tests live in
  `python_pkg/workflow_tests/symbolic/` and auto-skip under the LimTDD backend.
  `test_backend_crosscheck.py` is the backend-agnostic dense cross-check oracle
  (small unitary circuits vs Qiskit `Statevector`).
- `python_pkg/examples/` — runnable workflow validation examples (no assertions).
- `python_pkg/eval/` — benchmark evaluation/experiment runners.
- `python_pkg/plots/` — plotting scripts.

Current examples:

```bash
cd python_pkg
../.venv/bin/python workflow_tests/test_ts_structure.py
../.venv/bin/python workflow_tests/test_lazy_measurement.py
../.venv/bin/python workflow_tests/test_backend_crosscheck.py
../.venv/bin/python examples/test_RUS.py
```

## NuSMV integration

`python_pkg/qreach/qctl.py` can emit SMV models and call NuSMV through:

```python
modelChecking(ts, ctl_formula, nusmv_path=None)
```

By default, the code looks for NuSMV at:

```text
../NuSMV-2.7.0-linux64/bin/NuSMV
```

Pass `nusmv_path` explicitly or ensure `NuSMV` is available at the expected path or on `PATH`.

Example:

```python
from qreach.qctl import modelChecking

result = modelChecking(ts, "AG (leaf -> reached)")
print(result)
```

## Development notes

The project currently contains two transition-system directions:

- `qts_naive::TransitionSystem` in `transition_system.hpp`: explicit transition system and current practical workflow target.
- `qts::TransitionSystem` / `SymTS` in `transition_system_qadd.hpp`: QADD-backed symbolic transition system, useful for experiments and regressions.

Near-term workflow work prioritizes Qiskit parsing, annotation/labelling APIs, convenient `QOperation` construction, and practical debugging over aggressive QADD symbolic optimization.
