# Grover32 Converted QASM Timeout Investigation Handoff

> **Status:** Phase 1 complete (normalization fix), Phase 2 diagnosis in progress (CFLOBDD matrix multiply).  
> **Last updated:** 2026-07-03

**Goal:** Explain and fix why `single-it-grover32-plus.qasm` times out in the lazy QASM workflow while `single-it-grover32-zero.qasm`, `single-it-grover64-plus.qasm`, and `single-it-grover32-plus-linear.qasm` all finish quickly.

**Architecture:** The lazy QASM runner parses Qiskit circuits into the explicit `pyqreach.TransitionSystem` using `parse_qiskit_cir_lazy(...)`. During parsing, each added relation propagates a quantum post-image through `QOperation::postImage()` → `SingleVecTerm::applyGate()` → `Matrix1234ComplexFloatBoost::MatrixMultiplyV4WithInfo()` → CFLOBDD DAG traversal.

**Tech Stack:** Python 3 via `uv`/`.venv`, Qiskit 1.4.2, pybind11 extension `pyqreach`, C++20, CFLOBDD backend, explicit transition system `qts_naive::TransitionSystem`.

## Global Constraints

(unchanged from original — see commit history for full list)
- Work in repository: `/Users/ftdac/thu/qreach-tools`.
- Use `.venv` / `uv` workflow.
- Preserve existing workflow APIs.
- Profiling switches must be opt-in (env-var gated).
- Do not commit without user approval.

---

## What We Now Know

### Circuit dimensions

| File | qreg size | real qnum | CFLOBDD qNum | level | gates | initial state |
|------|-----------|-----------|-------------|-------|-------|---------------|
| `single-it-grover32-plus.qasm` | q[63] | 63 | 64 | 7 | 478 | H^⊗63 |
| `single-it-grover32-zero.qasm` | q[63] | 63 | 64 | 7 | 446 | \|0⟩^63 |
| `single-it-grover32-plus-linear.qasm` | q[64] | 64 | 64 | 7 | 482 | H^⊗64 |
| `single-it-grover64-plus.qasm` | q[127] | 127 | 128 | 8 | 958 | H^⊗127 |

- "Grover32" = 32 data qubits + 31 helper qubits = 63 total (non-linear), or 64 total (linear).
- "Grover64" = 64 data qubits + 63 helper qubits = 127 total (non-linear).

### Performance matrix

| Circuit | lazy parse (prefix 102) | full lazy parse | fixed-point post |
|---------|------------------------|-----------------|-----------------|
| Grover32-plus | 22.7s → **1.48s** (fixed) | still >120s | N/A (lazy parse times out first) |
| Grover32-zero | ~0.08s (prefix 179) | — | — |
| Grover32-plus-linear | ~0.02s | — | — |
| Grover64-plus | ~0.03s | ~0.25s | ~0.03s |
| Grover32-plus (no-lazy) | — | parse fast | **>300s** (worse than lazy!) |

### What was done (Phase 1 — normalization fix)

1. **Parser profiling** (`QREACH_PARSE_PROFILE`): Added env-gated per-instruction timing in `python_pkg/parse_qiskit.py`. Confirmed timeout is during lazy parse, with per-instruction time growing from microseconds to multiple seconds.

2. **Prefix bisect script** (`debug_grover32_ccx_pathology.py`): Bisect found `smallest_failing_prefix=102` at 20s timeout. Adjacent passing prefix 101 completes in ~19.3s (already pathological).

3. **C++ gate/post-image profiling** (`QREACH_GATE_PROFILE`): Instrumented `SingleVecTerm::applyGate` and `QOperation::postImage` to separate concretize / matrix multiply / normalization phases.

4. **Normalization skip optimization** (`knownUnitNorm`): For singleton post-image with a basis-string-derived input vector and a unitary gate, skip `normalizeInline()`. This brought Grover32 prefix 102 from 22.7s → 1.48s. **The fix works for early-mid prefixes, but is insufficient for the full circuit.**

5. **No-lazy experiment**: No-lazy mode is even worse — `computingFixedPointPost()` runs for >5 minutes without completing. This is because no-lazy creates many more locations, and fixed-point iteration computes far more post-image calls across all of them.

### Root cause localization (Phase 2 — after normalization fix)

With normalization skipped, `QREACH_GATE_PROFILE` shows the bottleneck is now **`MatrixMultiplyV4WithInfo`** (CFLOBDD matrix-vector multiply):

| instruction | gate | multiply_us | normalize_us | total postImage |
|------------|------|------------|-------------|-----------------|
| idx=89 | CCX [22,23,43] | 39,421 | 0 | 663ms |
| idx=94 | CCX [24,25,44] | 110,251 | 0 | 1,471ms |
| idx=99 | CCX [26,27,45] | 348,382 | 0 | 3,387ms |
| idx=109 | CCX [30,31,47] | ~3,500,000 | 0 | ~4,700ms |
| idx=118 | CCX [44,45,54] | ~4,000,000 | 0 | ~5,200ms |

Multiply time grows exponentially with accumulated entanglement on higher-index qubit pairs.

### Why only Grover32-plus (level=7, H-initialized)?

Static code analysis of `matrix1234_complex_float_boost_top_node.cpp` and `matrix1234_node.cpp` reveals:

- `clearMultMap()` at level≥5 is NOT the cause — Grover64 (level=8) also clears cache but runs fast.
- Grover32-plus and Grover32-zero are both level=7, same qNum=64 — but zero is fast (initial state |0⟩^63 has one nonzero amplitude).
- Grover32-plus-linear is also level=7 — also fast (symmetric linear qubit ordering → CFLOBDD DAG sharing).
- The key differentiator for Grover32-plus is the combination of: **H^⊗63 initial state** + **nonlinear CCX qubit ordering crossing CFLOBDD split boundaries** at level=7.

**Leading hypothesis:** The CFLOBDD `returnMapHandle` (leaf-amplitude array) fragments at level=7 when starting from a uniform superposition and applying crossing-split CCX gates. Each matrix multiply then evaluates `O(N × M)` double-loop products where N = number of distinct output amplitudes and M = average entries per `MatMultMap`. Both N and M blow up because the level=7 CFLOBDD DAG cannot collapse near-identical subtrees when the state is entangled across the 32/32 split with a "dead" padding qubit (q[63] always |0⟩).

At level=8 (Grover64), the 64/64 split and symmetric qubit usage create more opportunities for DAG sharing, keeping returnMapHandle small. At level=7 with zero initial state, the state remains sparse enough to avoid fragmentation.

---

## Files Created/Modified in Phase 1

| File | Status | Purpose |
|------|--------|---------|
| `python_pkg/parse_qiskit.py` | modified | `QREACH_PARSE_PROFILE` env-gated parser profiling |
| `python_pkg/workflow_tests/debug_grover32_ccx_pathology.py` | new | Prefix timing, bisect, timeout worker |
| `python_pkg/workflow_tests/test_grover32_lazy_prefix_performance.py` | new | Focused pytest regression (prefix 102 <5s, controls <1s) |
| `quantum_operation.hpp` | modified | `knownUnitNorm` normalization skip, `QREACH_GATE_PROFILE` profiling, fallback propagation |
| `docs/agent-handoffs/grover32-timeout-localization-notes.md` | new | Evidence from localization experiments |
| `docs/superpowers/specs/2026-07-02-grover32-timeout-debugging-design.md` | new | Design spec |
| `docs/superpowers/plans/2026-07-02-grover32-timeout-diagnostics.md` | new | Diagnostic implementation plan |
| `docs/superpowers/plans/2026-07-02-grover32-lazy-postimage-timeout.md` | new | Fix implementation plan |

### Key commits (on branch `qts-rollback`, not yet pushed)

```
1fc662b chore: propagate knownUnitNorm after normalizeInline fallback
2b0ef00 fix: apply timeout to non-bisect prefix runs too
05fcfef fix: correct import path in lazy prefix performance regression
8d8354e fix: speed up Grover32 lazy post-image path       ← core fix
3382407 chore: guard postimage profiling timing
6272257 chore: add env-gated post-image gate profiling
fe8a02c test: add Grover32 lazy prefix performance regression
ba031f5 plan: target Grover32 lazy post-image timeout
79837a7 Record Grover32 timeout localization evidence
7794890 Add Grover32 QASM prefix debugger
d65a18d Add opt-in parser profiling
```

---

## Phase 2 Investigation Plan: CFLOBDD Matrix Multiply Fragmentation at Level 7

### Task 1: Profile returnMapHandle growth to confirm DAG fragmentation

**Goal:** Measure whether `returnMapHandle.Size()` (the number of distinct leaf amplitudes) grows disproportionately in Grover32-plus vs all fast variants.

**Approach:** Add temporary `c1_returnmap_sz` / `c2_returnmap_sz` / `result_returnmap_sz` fields to the `QREACH_GATE_PROFILE` post-image output. This requires one-line additions in `MatrixMultiplyV4WithInfoTopNode` in `matrix1234_complex_float_boost_top_node.cpp:985-991`.

**Run matrix:**
```
Grover32-plus prefix 102 (slow)  → expected: result_returnmap_sz >> 10
Grover32-zero prefix 102 (fast)  → expected: result_returnmap_sz ≤ 5
Grover32-plus-linear prefix 102  → expected: result_returnmap_sz ≤ 10
Grover64-plus prefix 132 (fast)  → expected: result_returnmap_sz ≤ 20
```

**Acceptance:** If Grover32-plus `returnMapHandle.Size()` is significantly larger than fast variants, DAG fragmentation is confirmed.

### Task 2: Establish a small CCX-only synthetic reproducer

**Goal:** Remove QASM dependency and test CFLOBDD multiply directly from Python.

```python
# Construct H^⊗63 state, then apply CCX[26,27,45] repeatedly
state = pyqreach.QOperation(['+' * 63])  # needs 63-char + string
ccx_op = pyqreach.QOperation("CCX", 63, [26, 27, 45], [])
# Time postImage(state, ccx_op) as function of accumulated state complexity
```

**Acceptance:** A script that reproduces the multi-second multiply time in < 10 lines of Python, suitable for profiling individual `MatrixMultiplyV4WithInfo` calls.

### Task 3: Check whether disabling `clearMultMap()` helps at level 7

**Rationale:** Even though level=8 also clears cache and is fast, the clearing might interact differently with the level=7 fragmented DAG. A quick experiment to rule this out definitively.

**Change:** In `matrix1234_complex_float_boost_top_node.cpp:976-977`, comment out:
```cpp
if (c1->level >= 5)
    clearMultMap();
```
Rebuild, run prefix 102 benchmark. If multiply time drops significantly, caching is a partial contributor.

### Task 4: Inspect `NormFormComplex` precision / `checkFactor` interaction

**Goal:** Understand whether the `round(real*1e30)/1e30` in `NormFormComplex` is failing to merge near-identical amplitudes at level=7 due to `cpp_complex_100` precision.

**Investigate:**
- `NormFormComplex` rounds to 30 decimal digits. With 100-digit `cpp_complex_100`, intermediate results preserve 70 digits of "noise" that prevents merging.
- The `checkFactor = 2^(2^(level-1) - 1)` grows exponentially with level. At level=7, `checkFactor = 2^63 ≈ 9e18`; at level=8, `checkFactor = 2^127 ≈ 1.7e38`.
- The zero-check `abs(real*checkFactor) < 1e-7` has different effective thresholds at different levels.

**Experiment:** Temporarily increase rounding in `NormFormComplex` (e.g., `round(real*1e15)/1e15`) and observe whether returnMapHandle shrinks and multiply time improves.

### Task 5: Based on evidence, design and implement the fix

Options, in order of preference:

1. **If DAG fragmentation confirmed (Task 1):** Optimize the `MatrixMultiplyV4WithInfoTopNode:991-1000` evaluation loop — currently O(N×M) with `std::map` traversal. Options:
   - Pre-compute products of `(c1_returnmap[i] * c2_returnmap[j])` and reuse across similar MatMultMap entries
   - Hash MatMultMap results by (index1, index2) product value to avoid redundant big-integer × complex multiplications

2. **If cache clearing interaction (Task 3):** Gate the `clearMultMap()` on a more precise condition, or add a level-7-specific cache that persists across calls.

3. **If NormFormComplex precision (Task 4):** Adjust rounding based on level to merge physically indistinguishable amplitudes, effectively reducing returnMapHandle size.

4. **If none of the above work:** Consider a Python-level workaround — detect the qNum=63+padding pattern and pre-pad to qNum=64 with an extra explicit qubit in the QASM (opt-in flag).

### Task 6: Validate fix and update full Grover32 benchmark

- Rebuild (C++ + pybind11)
- Run pytest regression: `pytest workflow_tests/test_grover32_lazy_prefix_performance.py -q`
- Run full Grover32-plus benchmark with 120s timeout
- If passes, run Grover32-zero, Grover64-plus, Grover64-zero for regression
- Run `test_lazy_measurement.py`, `test_grover.py`, `test_newapi.py`

---

## Notes for the Implementing Agent

- **Do not commit without user approval.** The user will manually review and commit.
- CFLOBDD rebuild is slow (~minutes). Make minimal changes per iteration.
- The `returnMapHandle` profiling addition is a one-line change — add `.Size()` to existing profiling output. Remove it after confirming the hypothesis.
- The `clearMultMap()` experiment is a one-line comment-out — easy to revert.
- The `NormFormComplex` experiment changes one constant — easy to revert.
- Focus verification on **prefix 102** (fast path for regression tests) and **full Grover32-plus** (the actual timeout case).
- Key files for Phase 2:
  - `cflobdd/CFLOBDD/matrix1234_complex_float_boost_top_node.cpp:974-1040` — `MatrixMultiplyV4WithInfoTopNode`
  - `cflobdd/CFLOBDD/matrix1234_node.cpp:4499-4837` — `MatrixMultiplyV4WithInfoNode`
  - `cflobdd/CFLOBDD/matmult_map.h:96` — `std::map<INT_PAIR, VAL_TYPE>` (MatMultMap storage)
