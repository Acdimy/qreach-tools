# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project direction

QReach is a quantum model-checking / reachability-analysis tool for quantum Markov chains and Qiskit programs. The original README describes a CFLOBDD-backed symbolic-DD direction, but current development should prioritize practical Qiskit program debugging over aggressive QADD symbolic optimization.

Near-term work should favor:
- Specializing and tooling the Python workflow around Qiskit parsing, labeling, qCTL/SMV generation, and regression harnesses.
- Adding practical Python APIs that make QReach usable as a debugging library.
- Pragmatic construction/runtime optimizations for the explicit `qts_naive::TransitionSystem` path.

Deprioritize the QADD symbolic optimization roadmap from `plan.md` unless the user explicitly asks for it. `transition_system_qadd.hpp`, `qadd.hpp`, and symbolic post-image work remain useful references/regression baselines, but do not assume they are the preferred next implementation target.

## Build and test commands

C++ build from the repository root:

```bash
make all          # builds libqreach.so
make test         # builds test_qreach
./test_qreach 8   # run the C++ benchmark/test executable with an example argument
make clean        # remove libqreach.so and object files
```

Python dependencies and bindings are manual, not packaged through `setup.py`/`pyproject.toml`. Use the `qmc` conda environment if available, then build from `python_pkg`:

```bash
source ~/miniconda3/etc/profile.d/conda.sh
conda activate qmc
cd python_pkg
export PYTHON_INCLUDE=$(python -c "from sysconfig import get_paths as gp; print(gp()['include'])")
export BOOST_PATH=~/stab_dd/boost_1_81_0
invoke build-qreach
invoke build-pybind11
```

`python_pkg/tasks.py` hardcodes `python3.10` for pybind11 include/config commands, so use a Python 3.10 environment or update that task deliberately.

Common Python regression entrypoints, run from `python_pkg` after building `pyqreach`:

```bash
python test_symts_minimal.py              # minimal explicit-vs-symbolic post check
python test_symts_control_flow.py         # reduced if/while Qiskit control-flow checks
python test_symts_steane_regression.py    # stable bounded Steane symbolic regression
python test_RUS.py                        # RUS Qiskit/debugging workflow test
python test_steane_qiskit.py              # older/full Steane Qiskit workflow reference
```

Useful profiling/debug commands:

```bash
TS_PROFILE=1 ./test_qreach 8
cd python_pkg
TS_PROFILE=1 python debug_steane_measure_only.py 1
TS_PROFILE=1 python test_symts_steane_post.py sym 1 full
TS_MAX_POST_ITER=5 python test_symts_minimal.py
```

NuSMV integration in `python_pkg/qctl.py` defaults to `../NuSMV-2.7.0-linux64/bin/NuSMV`; pass `nusmv_path` to `modelChecking(...)` or ensure `NuSMV` is available there/in `PATH`.

## High-level architecture

The C++ layer supplies the semantic engine and exposes two transition-system implementations:

- `transition_system.hpp` defines the explicit baseline namespace `qts_naive`. `qts_naive::TransitionSystem` stores `Location` objects, explicit relation maps, pre/post adjacency, quantum bounds, and classical propositions. Treat this as the semantic reference and the current preferred optimization target for practical workflows.
- `transition_system_qadd.hpp` defines the QADD-backed namespace `qts`. `qts::TransitionSystem` stores location annotations and relations as QADDs with interleaved source/target location bits. It contains the current fused symbolic post-image fixed-point implementation and `TS_PROFILE` counters.
- `qadd.hpp` implements the QADD node/unique-table/compute-table layer and generic `Apply` operators (`ADD`, `JOIN`, `MEET`, `DIFF`, `APPLY`).
- `quantum_operation.hpp` is the terminal algebra layer over the CFLOBDD backend: quantum subspace/operator creation, comparison, conjunction/disjunction/difference, and `postImage` semantics.
- `cl_proposition.hpp` stores classical propositions attached to locations; both explicit and symbolic Python workflows use it.
- `cflobdd/` is the imported/modified CFLOBDD backend used by quantum operations.

The Python layer is the user-facing workflow surface:

- `python_pkg/qreach_python_wrapper.cpp` builds the `pyqreach` extension with pybind11. It exposes `QOperation`, `ClassicalProposition`, explicit `TransitionSystem`/`Location`, and symbolic `SymTS`.
- `python_pkg/parse_qiskit.py` lowers Qiskit `QuantumCircuit` objects into either explicit `pyqreach.TransitionSystem` (`parse_qiskit_cir`) or symbolic `pyqreach.SymTS` (`parse_qiskit_cir_sym`). It handles classical-register bookkeeping, measurements, resets, and Qiskit control flow such as `if_else`, `while_loop`, `for_loop`, and `switch_case`.
- `python_pkg/qctl.py` labels transition-system locations, converts transition systems to dictionaries/NetworkX/SMV, emits CTL formulas for NuSMV, runs NuSMV through `modelChecking`, and parses satisfaction/counterexamples.
- The `python_pkg/test_*.py` files are mostly executable regression scripts rather than pytest-style unit tests. Prefer running the specific script related to the touched workflow.

Typical end-to-end flow:

1. Build `libqreach.so` and `python_pkg/pyqreach...so`.
2. Create or load a Qiskit circuit in Python.
3. Lower it with `parse_qiskit_cir(...)` for explicit TS work, or `parse_qiskit_cir_sym(...)` only when maintaining symbolic regressions.
4. Use `qctl.py` to label states, generate SMV/CTL, and run NuSMV or inspect counterexamples.
5. Compare against explicit `qts_naive::TransitionSystem` behavior when symbolic behavior is uncertain.

## Notes for future changes

- Prefer extending explicit/Qiskit/qCTL workflows before adding new QADD internals.
- When modifying parser behavior, check both explicit and symbolic helper branches in `parse_qiskit.py` only if the symbolic regression surface is affected.
- Python-side explicit-vs-symbolic comparisons can be sensitive to repeated global transition-system initialization; existing docs recommend subprocess isolation for safer comparisons.
- `README.md` is sparse and older than `plan.md`/`prompt.md`; use those files for current context, with the project-direction override above taking precedence.
