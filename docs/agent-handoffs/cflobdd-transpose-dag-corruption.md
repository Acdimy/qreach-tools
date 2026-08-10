---
name: cflobdd-transpose-dag-corruption
description: CFLOBDD DAG corruption bugs — transpose return-map routing, MatrixMultiplyV4 coefficient overflow, and post_image chain DAG damage. GramSchmidt mitigated via H×content trick; V4 overflow fixed via convert_to<double>; level≥9 post_image remains.
metadata:
  type: project
  status: partially-fixed
  branch: qts-rollback
  date: 2026-08-10
  updated: 2026-08-10
---

# CFLOBDD DAG Corruption Bugs

## Summary

Multiple CFLOBDD bugs prevent correct operation at higher levels (≥7).
Two are **fixed/mitigated**; one remains at level ≥ 9.

| Bug | Operation | Symptom | Level | Status |
|-----|-----------|---------|-------|--------|
| **A/B** (transpose) | `MatrixTransposeNode` | retMapSz>2 / SIGSEGV | ≥3 (4q) | ✅ mitigated (H×content trick) |
| **C** (transpose) | `MatrixTranspose` | dot() row 0 destroyed | ≥8 (128q) | ✅ mitigated (H×content trick in dot()) |
| **D** (V4 coeff overflow) | `MatrixMultiplyV4TopNode` | dot() returns 0 instead of 1.0 | ≥7 (64q) | ✅ **fixed** (`convert_to<double>`) |
| **E** (post_image DAG) | `MatrixMultiplyV4WithInfo` | satisfy(self)=False, wrong dot values | ≥9 (256q) | ❌ remaining |

The most significant remaining issue is **Type E**: at CFLOBDD level ≥ 9,
`MatrixMultiplyV4WithInfo` (used in `post_image` / gate application) produces
states with corrupted amplitudes, causing self-consistency failures and
incorrect inner products.

## Bug D: MatrixMultiplyV4 Coefficient Overflow (FIXED)

### Root Cause

`MatrixMultiplyV4TopNode` evaluates bilinear polynomials by converting
coefficients from `cpp_int` to `unsigned long long int`.  The coefficients
grow doubly-exponentially with CFLOBDD level:

| Level | Coefficient | Fits in uint64? |
|-------|-------------|-----------------|
| 6 | 2^32 | ✅ |
| 7 | 2^64 | ❌ (ULLONG_MAX = 2^64 - 1) |
| 8 | 2^128 | ❌ |
| 9 | 2^256 | ❌ |

At level 7 (64 qubits), the coefficient `2^64` overflows `unsigned long long int`,
wrapping to 0.  This causes `dot(self,self)` to return 0 instead of 1.0.

### Fix

`matrix1234_complex_float_boost_top_node.cpp` line 937:
```cpp
// Before (broken at level ≥7):
auto factor = j.second.convert_to<unsigned long long int>();

// After (fixed up to level 8):
auto factor = j.second.convert_to<double>();
```

`double` exactly represents all power-of-2 coefficients and the final
`factor × amplitude_product` is always O(1).

### Verification

`repro_postimage_corruption.py`: n=32,64,96,128 all pass.
n=256 still fails (Type E).

## Bug E: Level-9 post_image DAG Corruption (REMAINING)

At CFLOBDD level ≥ 9 (256+ qubits), `MatrixMultiplyV4WithInfo`
(gate application via `post_image`) produces states where:

- Amplitudes are wrong (e.g., 5.3e-23 instead of 2^{-128} for |+>^256)
- `dot(self,self)` returns garbage (e.g., 1.8e16 instead of 1.0)
- `satisfy(self)` returns False

This is distinct from the coefficient overflow (Bug D) — the underlying
CFLOBDD DAG structure produced by `MatrixMultiplyV4WithInfo` is incorrect
at level 9.

### Affected Circuits

- grover128-plus-linear (256q): ❌ fail
- grover150-plus-linear (300q → qNum=512, level 10): ❌ fail
- grover128-zero-linear: ✅ pass (final state = |0⟩, trivial DAG)
- repro_postimage_corruption n=256: ❌ fail

### Investigation Trail

1. CFLOBDD DAG node-by-node dump shows correct static structure at level 7
2. `MatrixMultiplyV4Node` (recursive symbolic multiply) produces correct
   `cpp_int` bilinear polynomial coefficients (`2^64` at level 7)
3. The damage is in `MatrixMultiplyV4WithInfoNode`'s DAG construction,
   not in the multiply evaluation
4. "Warmup" (calling Gram-Schmidt on any unrelated state before gate
   application) changes the CFLOBDD global unique-table state, which
   alters `MatrixMultiplyV4WithInfoNode`'s hash-consing behavior and
   produces a different DAG — confirming the root cause is in the
   unique-table / canonicalization interaction at high levels

## Other Fixes Applied

### satisfy() semantics preserved

Investigation confirmed `satisfy()` logic (`lowerBound ⊆ spec ⊆ upperBound`)
is correct.  The test failure was caused by `parse_qiskit_cir_lazy` setting
`upperBound = initial_op` (overwriting the default Identity), which broke
the `spec ⊆ upperBound` check.  Fix: removed the line.

### Reproducer improvements

- `repro_transpose_bug.py`: covers Gram-Schmidt / dot / normalize correctness
- `repro_postimage_corruption.py`: covers post_image chain self-consistency
