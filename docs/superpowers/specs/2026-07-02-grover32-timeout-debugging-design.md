# Grover32 Lazy QASM Timeout Debugging Design

Date: 2026-07-02

## Context

The lazy QASM workflow times out on `single-it-grover32-plus.qasm` and `single-it-grover32-zero.qasm`, while larger Grover64 converted QASM files and Grover32 linear variants finish quickly. The current CCX/padding/split-boundary explanation is only a working hypothesis. The first implementation phase must prioritize localization over fixes.

## Goals

1. Identify the exact instruction or smallest instruction prefix that triggers the Grover32 lazy parse slowdown.
2. Distinguish whether time is spent in Python parser bookkeeping, `QOperation` construction, or lazy post-image propagation.
3. Build a repeatable minimal reproducer suitable for a focused regression after the root cause is known.
4. Preserve default parser, runner, and CSV behavior unless explicit opt-in debug switches are enabled.

## Non-goals

- Do not implement a workflow workaround before localizing the root cause.
- Do not assume CCX is the final root cause until profiling evidence points there.
- Do not change product-state `QOperation(["..."])` syntax.
- Do not add a 300-second regression test.

## Approach

Use a staged diagnostic path.

### 1. Opt-in parser profiling

Add environment-variable controlled profiling inside `python_pkg/parse_qiskit.py` around the main `parse_qiskit_cir(...)` instruction loop.

Initial interface:

- `QREACH_PARSE_PROFILE=1`: enable profiling.
- `QREACH_PARSE_PROFILE_THRESHOLD=<seconds>`: print only instructions whose elapsed time meets the threshold. Default: `0.25`.
- `QREACH_PARSE_PROFILE_VERBOSE=1`: print every instruction regardless of threshold.

Each emitted line should include:

- original instruction index (`pivot + loop_index`),
- operation name,
- qubit indices,
- current location count before and after,
- instruction elapsed seconds,
- total transition-system locations.

Default behavior with profiling disabled must be unchanged.

If instruction-level timing is insufficient, add a second opt-in detail level that separates operation construction from `_add_post_and_propagate(...)`. That deeper probe should still be disabled by default.

### 2. Prefix reproducer and bisect script

Create `python_pkg/workflow_tests/debug_grover32_ccx_pathology.py`.

The script should:

- load a QASM file,
- build a truncated `QuantumCircuit` containing the first N instructions while preserving circuit registers,
- run `parse_qiskit_cir_lazy(...)`,
- optionally run `ts.computingFixedPointPost()`,
- print load, parse, fixed-post, location, and result-location timings,
- support `--bisect` by running prefixes in subprocesses with a timeout,
- report the smallest failing prefix and largest adjacent passing prefix.

The script must not write CSV files or modify benchmark QASM files.

### 3. Regression test after minimization

Create `python_pkg/workflow_tests/test_grover32_ccx_pathology.py` after the minimal reproducer is known.

The test should:

- skip clearly if `pyqreach` or Qiskit is unavailable,
- use a minimized circuit or prefix, not the full 300-second benchmark,
- assert lazy parse completes under a short budget after the fix,
- check basic structural expectations such as `result.lazy is True` and location count matching the prefix shape.

### 4. Root-cause investigation after localization

Only after the profiling and bisect evidence identifies the failing instruction or prefix should C++ fixes be attempted.

Possible investigation targets, depending on evidence:

- `quantum_operation.hpp` CCX concretization if `QOperation("CCX", ...)` construction is slow,
- `postImage(...)` and normalization if lazy propagation is slow,
- `cflobdd/CFLOBDD/matrix1234_node.cpp::MkCCNOTNode(...)` if CFLOBDD Toffoli construction is implicated,
- parser state/location handling if the slowdown is not gate-specific.

The working CCX/padding/split hypothesis should be tested against Grover32 non-linear, Grover32 linear, Grover64, and synthetic prefixes rather than treated as proven.

## Validation Plan

First diagnostic validation:

1. Run Grover64 with parser profiling enabled and confirm it completes quickly with useful timing lines.
2. Run Grover32 with parser profiling enabled under a short external timeout and capture the last completed instruction.
3. Run the prefix/bisect script to identify the smallest slow prefix.
4. Compare the failing prefix against an adjacent passing prefix.

Post-fix validation:

1. Focused Grover32 pathology regression passes.
2. Existing lazy measurement regression passes.
3. `single-it-grover32-plus.qasm` and `single-it-grover32-zero.qasm` no longer timeout in the lazy benchmark runner.
4. Grover64 plus/zero remain fast.
5. Default CSV schema and default parser output remain compatible.
