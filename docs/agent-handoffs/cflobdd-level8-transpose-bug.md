# CFLOBDD Level 7+ Bugs — Root Cause of grover128/150 Verification Failures

**Date:** 2026-08-10
**Status:** Partially fixed. Two distinct bugs identified: transpose routing (mitigated), MatrixMultiplyV4 coefficient overflow (**fixed** for level≤8). Level≥9 post_image DAG corruption remains.

## Executive Summary

The grover128/150 `-linear` debug-check failures have multiple root causes:

1. **`MatrixTranspose` bug** (mitigated): corrupts row 0 of the transposed vector at level ≥ 8. Mitigated by H×content identity-multiply trick in `normalize()` and `dot()`.

2. **`MatrixMultiplyV4TopNode` coefficient overflow** (**FIXED**): `convert_to<unsigned long long int>()` overflows at level ≥ 7 (coefficient `2^64` > `ULLONG_MAX`). Fixed by `convert_to<double>()`.

3. **Level≥9 post_image DAG corruption** (REMAINING): `MatrixMultiplyV4WithInfo` produces incorrect state amplitudes at CFLOBDD level ≥ 9 (256+ qubits), causing `dot(self,self)` to return garbage values (e.g., 1.8e16 instead of 1.0).

## Updated Impact Matrix

| Circuit | Level | status | Why |
|---------|-------|--------|-----|
| grover32-plus/zero-linear (64q) | 7 | ✅ PASS | V4 overflow fix + warm DAG |
| grover64-plus/zero-linear (128q) | 7 | ✅ PASS | V4 overflow fix |
| grover128-zero-linear (256q) | 8 | ✅ PASS | Final state = \|0⟩, trivial DAG |
| grover128-plus-linear (256q) | 8 | ✅ PASS (after fix) | V4 overflow fix |
| grover128-plus-linear (256q) | 9 | ✗ FAIL | Level≥9 post_image DAG corruption |
| grover150-plus-linear (300q) | 10 | ✗ FAIL | Level≥10 post_image DAG corruption |
| grover150-zero-linear (300q) | 10 | ✅ PASS | Final state = \|0⟩, trivial DAG |
1. **`dot()`** — uses `MatrixTranspose` directly → always hits the bug
2. **`normalize()`** — uses `I × content` (Identity-multiply trick) → **avoids** transpose → works correctly

## Minimal Reproduction

### Test: `dot(|0⟩^n, |0⟩^n)` should be 1.0

```python
import pyqreach
pyqreach.initializeTransitionSystem()
from qiskit import QuantumCircuit
from parse_qiskit import parse_qiskit_cir_lazy

def make_state(n):
    return pyqreach.QOperation(['0' * n])

for n in [64, 128, 256]:
    z1 = make_state(n)
    z2 = make_state(n)
    t = pyqreach.span_qops([z1, z2])  # triggers Gram-Schmidt → dot()
    # Verify: |0> should be in span(|0>, |0>)
    qc = QuantumCircuit(n); qc.x(0)
    pyqreach.initializeTransitionSystem()
    ts = pyqreach.TransitionSystem()
    r = parse_qiskit_cir_lazy(qc, n, ts, initial_state='0'*n, return_metadata=True)
    print(f'n={n}: |0> in span(|0>,|0>) = {ts.Locations[0].satisfy(t)}')
```

**Expected:** All True (trivially, `|0⟩` is one of the spanning vectors).
**Actual:** True for n ≤ 64, depends on transpose corruption for n ≥ 128.

### Test: `dot(|0⟩^n, |+⟩^k|0⟩^(n-k))` is non-zero for k < n

```python
# At n=128, k=64: dot(|0>, |+>^64|0>^64) = 2^{-32} ≈ 2.3e-10 (non-zero!)
z = make_state(n)
p = make_state(n, list(range(k)))  # |+>^k|0>^(n-k)
t = pyqreach.span_qops([z, p])
# |0> should be in span(|0>, |+>^k|0>^(n-k)) since |0> IS one spanning vector
sat = ts.Locations[0].satisfy(t)
# Returns False when transpose bug hits
```

**Observed behavior at n=128, level=8:**

| k | `\|0⟩ ∈ span` | debug output |
|---|---------------|-------------|
| 1 | False | `dot() retMapSz=1, [0,0]=0` |
| 32 | False | `dot() retMapSz=1, [0,0]=0` |
| 64 | False | `dot() retMapSz=1, [0,0]=0` |
| 128 | **True** | `dot() retMapSz > 2` (different code path) |

The symmetric case (k = n/2 = 128) takes a different code path that happens to survive the bug.

## Root Cause Chain

### The Two Code Paths

**`normalize()` — Identity-multiply trick (WORKS):**
```
c1 = I × content          // ApplyGateF(qNum, 0, MkIdRelationInterleaved)
c1_conj = conj(transpose(c1))
mulres = c1_conj × c1      // ⟨c1|c1⟩ = norm²
```
The initial `I × content` creates a matrix-form CFLOBDD where row 0 has a different
internal DAG topology than the original vector.  The subsequent `transpose` operates
on this "safer" topology and row 0 survives.

**`dot()` — direct transpose (BROKEN):**
```
tmpVec = MatrixTranspose(other.content)   // ← BUG HERE
tmpVec = MatrixConjugate(tmpVec)
tmp = tmpVec × this->content              // ⟨other|this⟩
```
`MatrixTranspose` is called directly on the vector CFLOBDD.  At level ≥ 8, the
transpose operation corrupts the **internal group routing** of the DAG, assigning
entries to wrong exits.  Row 0 of the transposed matrix may end up all-zero.

### Why `Checkifzero` and `Normalize` Thresholds Are Not the Issue

The `[0,0]` entry extracted via **direct `EvaluateIteratively`** (which walks the
CFLOBDD tree without going through the corrupted return map) also returns **0**.
This is not a threshold problem — the value is genuinely missing from the DAG.

### Three Failure Modes of `dot()` at Level 8+

| resMap.Size() | What happened | [0,0] status | Fixable by fallback? |
|--------------|---------------|-------------|---------------------|
| `> 2` | Other rows corrupted, row 0 survives | ✅ correct | ✅ (already fixed) |
| `== 2` | No corruption | ✅ correct | N/A (normal path) |
| `== 1` | **Row 0 destroyed** | ❌ also 0 | ❌ **NOT fixable** |

### Propagation into Debug Check Failures

```
dot() returns 0 for |+⟩^k|0⟩^(n-k) · |0⟩^n  (should be 2^{-k/2})
  → Gram-Schmidt projection coefficient = 0  (should be 2^{-k/2})
  → v1 = |+⟩^k|0⟩^(n-k)  (NOT orthogonalized against |0⟩!)
  → span_qops result has non-orthogonal basis vectors
  → satisfy() disjunction → Gram-Schmidt of state against target
  → dot() again returns 0 for projection coefficients
  → residual survives (incorrectly) → dimension mismatch → False
```

### Why `normalize()` in Gram-Schmidt Shows `retMapSz=2` Correctly

In `span_qops([|0⟩, |+⟩^k|0⟩^(n-k)])`:
- v0 = |0⟩ → `normalizeInline()` → Identity-multiply path → `retMapSz=2, amp=1` ✓
- v1 = |+⟩^k|0⟩^(n-k) → after projection (which uses `dot()`, broken), then `normalizeInline()` → Identity-multiply path → correct normalization of the *wrong* vector

The normalize() result is "correct" (amp=1) because it uses the Identity trick. But the vector it's normalizing is already wrong because `dot()` returned 0 for the projection coefficient.

## Impact on Debug Verification

| Circuit | Level | Status | Why |
|---------|-------|--------|-----|
| grover32-plus/zero-linear (64q) | 7 | ✓ PASS | level < 8, transpose works |
| grover64-plus/zero-linear (128q) | 7 | ✓ PASS | level < 8, transpose works |
| grover128-plus-linear (256q) | 8 | ✗ FAIL | dot() broken, span_qops wrong |
| grover128-zero-linear (256q) | 8 | ✓ PASS | final state = \|0⟩ (trivial span check) |
| grover150-plus-linear (300q) | 8 | ✗ FAIL | same as grover128 |
| grover150-zero-linear (300q) | 8 | ✓ PASS | final state = \|0⟩ (trivial span check) |

## Fixes Applied (This Session)

### 1. Precision control: `zeroThreshold()` (Phase 4 in grover32 notes)

Added level-aware threshold to `checkifzero`, `singletonOrthogonalTo`, `isZero`.
Fixes grover64-plus-linear (inner product 2^{-32} now preserved).

### 2. `normalize()` fallback: scale `content`, not `c1` (Phase 4.5)

Changed `normalize()` retMapSz > 2 fallback from `return (1/factor) * c1` to
`return (1/factor) * content`.  `c1` is the Identity-multiplied (potentially
corrupted) vector; `content` is the original uncorrupted vector.

### 3. `dot()` fallback: handle `retMapSz == 1` (Phase 5)

Added `[0,0]` extraction fallback for `dot()` when `retMapSz == 1`.
**Partial fix only** — works when row 0 survives, but row 0 is sometimes
destroyed at level 8+.

## Remaining Work

1. **Fix `MatrixTranspose` at level ≥ 8** — the root cause.  The internal group
   reconstruction in `MatrixTransposeNode` routes assignments to wrong exits.
   See `docs/agent-handoffs/cflobdd-transpose-dag-corruption-analysis.md` and
   `cflobdd-transpose-dag-corruption-fix.md` for prior analysis.

2. **Alternative: make `dot()` use the Identity-multiply trick**, like `normalize()`
   does.  Instead of `transpose(content) → conj → multiply`, use:
   ```
   c1 = I × content     // Identity-multiply to get safe matrix form
   c1_conj = conj(transpose(c1))  // transpose the safe form
   result = c1_conj × other      // inner product
   ```
   This would completely avoid the transpose bug in `dot()`.

3. **Remove debug output** from `normalize()` and `dot()` before committing.

## Related Documents

- `docs/agent-handoffs/grover32-timeout-localization-notes.md` — Grover32 performance investigation
- `docs/agent-handoffs/cflobdd-transpose-dag-corruption-analysis.md` — Prior transpose bug analysis
- `docs/agent-handoffs/cflobdd-transpose-dag-corruption-fix.md` — Prior transpose fix attempt
- `docs/agent-handoffs/backend-replacement-api-contract.md` — Backend refactoring plan
