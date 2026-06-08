# qreach-tools Current Prompt

## Project Snapshot

This repository currently contains two transition-system implementations for symbolic quantum program analysis and model checking.

- `transition_system.hpp` is the original explicit transition system (`qts_naive::TransitionSystem`). It stores locations and relations explicitly and serves as the semantic baseline.
- `transition_system_qadd.hpp` is the current symbolic transition system (`qts::TransitionSystem`). It uses a QADD whose terminals are `QOperation` objects and whose internal variables encode interleaved source and target location bits.
- `qadd.hpp` implements the QADD node structure, unique table, compute table, terminal canonicalization, and generic `Apply` operators such as `ADD`, `JOIN`, `MEET`, `DIFF`, and `APPLY`.
- `quantum_operation.hpp` implements subspace and operator semantics on top of CFLOBDD-backed quantum terms. This is the terminal algebra layer used by QADD.
- `cl_proposition.hpp` stores classical propositions associated with locations. These propositions are now available in both Python-accessible symbolic and explicit workflows.

The current symbolic post algorithm is no longer the old `Apply -> exists_vars -> rename_vars` pipeline. In the current version, `compute_symbolic_post()` in `transition_system_qadd.hpp` uses a fused relation-side recursion that directly computes the next P-DD while eliminating source bits and projecting target bits into the state-DD shape.

## What Is Already Working

- Minimal SymTS end-to-end post computation is working and matches the explicit reference on the current minimal regression.
- Symbolic control-flow lowering for Qiskit `while_loop` and `if_test` is working in the tested reduced cases.
- Steane-derived bounded regressions are now stabilized in `python_pkg/test_symts_steane_regression.py`.
- The current fused symbolic post path preserves the minimal regression and the bounded Steane regression.
- The Python wrapper exposes the main SymTS functionality needed by the current post-only workflow.

## Current Build And Run Commands

### C++ benchmark build

Run from the repository root:

```bash
cd /home/dac22/stab_dd/qreach-tools
make test
./test_qreach 8
```

Useful profiling form:

```bash
cd /home/dac22/stab_dd/qreach-tools
TS_PROFILE=1 ./test_qreach 8
```

### Python binding build

Run from `python_pkg` inside the `qmc` conda environment:

```bash
source ~/miniconda3/etc/profile.d/conda.sh
conda activate qmc
cd /home/dac22/stab_dd/qreach-tools/python_pkg
export PYTHON_INCLUDE=$(python -c "from sysconfig import get_paths as gp; print(gp()['include'])")
export BOOST_PATH=~/stab_dd/boost_1_81_0
invoke build-qreach && invoke build-pybind11
```

### Current regression commands

Minimal symbolic-vs-explicit check:

```bash
source ~/miniconda3/etc/profile.d/conda.sh
conda activate qmc
cd /home/dac22/stab_dd/qreach-tools/python_pkg
python test_symts_minimal.py
```

Steane bounded regression:

```bash
source ~/miniconda3/etc/profile.d/conda.sh
conda activate qmc
cd /home/dac22/stab_dd/qreach-tools/python_pkg
python test_symts_steane_regression.py
```

Steane profiling commands currently used in diagnosis:

```bash
source ~/miniconda3/etc/profile.d/conda.sh
conda activate qmc
cd /home/dac22/stab_dd/qreach-tools/python_pkg
TS_PROFILE=1 python debug_steane_measure_only.py 1
TS_PROFILE=1 python test_symts_steane_post.py sym 1 full
```

## Important Files And Their Roles

- `qadd.hpp`
  Implements QADD nodes, unique-table and compute-table management, generic `Apply`, terminal canonicalization, and current zero-subtree short-circuits.

- `quantum_operation.hpp`
  Implements the actual quantum terminal semantics. This file defines subspace operations such as disjunction, conjunction, difference, comparison, and `postImage` over CFLOBDD-backed terms.

- `transition_system.hpp`
  Explicit baseline transition system. Use this as the semantic reference when a symbolic result is in doubt.

- `transition_system_qadd.hpp`
  Current symbolic transition system. Encodes locations as interleaved source/target bits, stores symbolic annotation and relation DDs, and contains the current fixed-point post computation.

- `python_pkg/qreach_python_wrapper.cpp`
  Pybind11 wrapper exposing `QOperation`, `ClassicalProposition`, explicit `TransitionSystem`, and symbolic `SymTS` to Python.

- `python_pkg/parse_qiskit.py`
  Lowers Qiskit circuits and control flow into either the explicit TS or SymTS. This is the main Python-side frontend for symbolic construction.

- `python_pkg/qctl.py`
  Labeling and qCTL helper layer. This is the main future integration surface for a full symbolic model-checking workflow, not just post-image fixed points.

- `python_pkg/test_symts_minimal.py`
  Minimal end-to-end regression comparing explicit and symbolic post results.

- `python_pkg/test_symts_control_flow.py`
  Regression coverage for reduced `while_loop` and `if_test` lowering.

- `python_pkg/test_symts_steane_post.py`
  Steane-derived comparison harness that can run either `construct` or `full` stages and compare explicit vs symbolic behavior.

- `python_pkg/debug_steane_measure_only.py`
  Focused profiler for the Steane base case without the full `if_test` comparison layer.

- `python_pkg/test_symts_steane_regression.py`
  Current stable regression entrypoint for bounded Steane checks on the fused post implementation.

## Current Performance Understanding

The current symbolic bottleneck is no longer `exists_vars()` or `rename_vars()`. In the fused post version, those phases have effectively been removed from the hot path.

For the current Steane cases with `measure_count=1`:

- `python debug_steane_measure_only.py 1` takes about 23.4 to 23.7 seconds.
- `python test_symts_steane_post.py sym 1 full` takes about 24.8 to 25.0 seconds.
- Both spend almost all of their time in `apply_ms`.
- Both still require about 38 to 40 fixed-point iterations.

This means the dominant cost is the recursive relation-side traversal inside fused symbolic `Apply(APPLY)`, not terminal algebra and not the old abstraction/renaming passes.

The relation node counts themselves are not huge in absolute size, but the current interleaved relation representation still induces expensive recursive traversal across many fixed-point rounds. So the main issue is better described as relation recursion shape and reuse efficiency, not simply “the relation is too large” in node-count terms.

## Current Constraints And Caveats

- Python-side combined naive and symbolic runs can still be sensitive to repeated global transition-system initialization, so subprocess isolation remains the safe default for comparisons.
- The current symbolic path only covers the post-image fixed-point flow. A full symbolic qCTL model-checking workflow has not been completed yet.
- The current fused post implementation is semantically validated on minimal and bounded Steane regressions, but larger symbolic workloads are still dominated by relation recursion.