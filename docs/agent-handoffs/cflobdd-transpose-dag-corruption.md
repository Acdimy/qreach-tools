---
name: cflobdd-transpose-dag-corruption
description: CFLOBDD MatrixTranspose corrupts column→row vector conversion for GramSchmidt-produced DAGs — three failure modes, two mitigated via H*content trick + fallbacks
metadata:
  type: project
  status: mitigated
  branch: qts-rollback
  date: 2026-08-10
---

# CFLOBDD Transpose DAG Corruption Bug

## Summary

`MatrixTranspose` corrupts the column→row vector conversion for CFLOBDD column vectors
with DAG topologies produced by GramSchmidt's PairProduct/MatrixPlus path.  A full
algorithmic fix to `MatrixTransposeNode` has not been found, but the bug is **mitigated**
in both `normalize()` and `dot()` by using an **identity-multiply (H×content) trick**
that transforms the CFLOBDD into a "safe" DAG topology before applying transpose,
plus fallbacks to `[0,0]` entry extraction.

## Three Failure Modes

| Type | Symptom | retMapSz | [0,0] entry | Triggered by | Mitigation |
|------|---------|----------|-------------|-------------|------------|
| **A** | Assertion `retMapSz<=2` | > 2 | ✅ correct | dqc_pe_2 (4q) | H×content in normalize() + [0,0] fallback |
| **B** | SIGSEGV in `MatrixMultiplyV4` | N/A | N/A | pe_7/qft_7 (7-8+q) | H×content in normalize() |
| **C** | `dot()` returns 0 (row 0 destroyed) | == 1 | ❌ 0 | level ≥ 8 (128+q) | H×content in dot() |

**Type C is the most severe**: even direct `EvaluateIteratively` at `[0,0]` returns 0 —
the matrix entry is genuinely destroyed.  No fallback can recover from this.

## Root Cause

`MatrixTransposeNode` at level≥2 uses original B-return-maps to route transposed exits.
After variable swap (voc1↔voc2), the original maps route to semantically wrong parent
exits because the partition of assignments has changed.  The corruption severity
increases with CFLOBDD level (structure size).

See `memory/cflobdd-transpose-dag-corruption-analysis.md` for the detailed DAG trace.

## Fix Applied: H×content Identity-Multiply Trick

### normalize() — H×content revert (REVERT qts-rollback)

The aa52325 commit changed normalize() from:
```cpp
// OLD (safe for pe_7):
c1 = MatrixMultiplyV4WithInfo(H, content);  // I × content
// then conj(transpose(c1)) * c1
```
To:
```cpp
// NEW (SIGSEGV on pe_7):
c_conj = MatrixConjugate(content);
c_conj = MatrixTranspose(c_conj);    // transpose directly on content
// then c_conj * content
```

The change introduced Type B crashes (SIGSEGV on pe_7).  **Reverted** to H×content
approach because the identity-multiply inserts an intermediate CFLOBDD with a
different internal DAG topology that `MatrixTranspose` handles safely.

If even the H×content path produces `retMapSz > 2`, a fallback extracts the `[0,0]`
entry via `EvaluateIteratively`.

Key code location: `quantum_operation.hpp` normalize(), tagged `REVERT(qts-rollback)`.

### dot() — H×content applied (REVERT qts-rollback)

Same trick applied to `dot()` to prevent Type C (level-8 row 0 destruction):
```cpp
// Before:
tmpVec = MatrixTranspose(other.content);           // direct transpose — DESTROYS row 0 at level≥8
// After:
auto other_safe = MatrixMultiplyV4WithInfo(H, other.content);  // I × other.content
tmpVec = MatrixTranspose(other_safe);              // transpose on safe DAG
```

If even the safe path produces `retMapSz > 2`, a fallback extracts the `[0,0]` entry.

Key code location: `quantum_operation.hpp` dot(), tagged `REVERT(qts-rollback)`.

Both `normalize()` and `dot()` now scale the original `content` (not the intermediate
`c1`/`other_safe`) in their return paths, minimizing risk from the I×content DAG wrapper.

## Reproducers

| Script | Circuit | What it tests | Status |
|--------|---------|---------------|--------|
| `python_pkg/repro_minimal.py` | dqc_pe_2 (4q) | Type A — retMapSz assertion | ✅ pass |
| `python_pkg/repro_tslabelling_crash.py` | pe_7 (8q) | Type B — SIGSEGV | ✅ pass |
| `python_pkg/workflow_tests/repro_transpose_bug.py` | synthetic (32-256q) | Type C — level-8 row 0 destruction | ⚠️ Test 1 pass, Test 2 fails |

## Remaining Issue: Level-8 GramSchmidt Pipeline

With the H×content trick in `dot()`, Type C (row 0 destruction) is **prevented** —
dot() now returns correct inner-product values at all levels.  However,
`repro_transpose_bug.py` Test 2 (`|0> ∈ span(|0>, |+>^k|0>^(n-k))` for k=1..64, n≥128)
still fails.  The failure is no longer in `dot()` but in subsequent GramSchmidt
pipeline steps: `projectOnto` (scalar×vector), `MatrixPlus` (`ivec + neg_proj`),
`checkifzero()`, or `normalizeInline()` — one of these operations is affected by
DAG topology issues at level ≥ 8.

## Modified Files

### `quantum_operation.hpp`
- `normalize()`: H×content revert (lines ~1060-1070) + [0,0] fallback + returns scaled `content`
- `dot()`: H×content identity-multiply (lines ~1014-1025) + [0,0] fallback

### `cflobdd/CFLOBDD/matrix1234_node.cpp`
- Line 2951, 2969: `// FIXME(qts-rollback)` — `m1.AddToEnd(v)` → `m1.AddToEnd(return_handle.LookupInv(v))`
  Fixes a real level-1 transpose bug (old exit values used instead of return_handle indices).
  Not the root cause of the current transpose corruption.

## Key Files

- `quantum_operation.hpp` — `normalize()`, `dot()`, `GramSchmidt()` (all fix sites)
- `cflobdd/CFLOBDD/matrix1234_node.cpp` — `MatrixTransposeNode` (root cause)
- `cflobdd/CFLOBDD/matrix1234_complex_float_boost_top_node.cpp` — `MatrixTransposeTop`
- `python_pkg/repro_minimal.py` — Type A reproducer
- `python_pkg/repro_tslabelling_crash.py` — Type B reproducer
- `python_pkg/workflow_tests/repro_transpose_bug.py` — Type C / level-8 reproducer
- `docs/agent-handoffs/cflobdd-level8-transpose-bug.md` — Level-8 specific analysis
- `memory/cflobdd-transpose-dag-corruption-analysis.md` — Deep DAG analysis
- `memory/cflobdd-transpose-dag-corruption-fix.md` — Fix history
