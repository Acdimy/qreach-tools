---
name: cflobdd-transpose-dag-corruption
description: CFLOBDD DAG corruption bugs — transpose return-map routing, MatrixMultiplyV4 coefficient overflow, and post_image chain DAG damage. GramSchmidt mitigated via H×content trick; V4 overflow fixed via convert_to<double>; level≥9 post_image remains.
metadata:
  type: project
  status: partially-fixed
  branch: qts-rollback
  date: 2026-08-10
  updated: 2026-08-11
---

# CFLOBDD DAG Corruption Bugs

## Summary

Multiple CFLOBDD bugs prevent correct operation at higher levels (≥7).
Two are **fixed/mitigated**; one remains at level ≥ 9.

| Bug | Operation | Symptom | Level | Status |
|-----|-----------|---------|-------|--------|
| **A/B** (transpose) | `MatrixTransposeNode` | retMapSz>2 / SIGSEGV | ≥3 (4q) | ✅ mitigated (H×content trick) |
| **C** (transpose) | `MatrixTranspose` | dot() row 0 destroyed | ≥8 (128q) | ⚠️ partially mitigated — platform-dependent (see §Platform Dependence) |
| **D** (V4 coeff overflow) | `MatrixMultiplyV4TopNode` | dot() returns 0 instead of 1.0 | ≥7 (64q) | ✅ **fixed** (`convert_to<double>`) |
| **E** (post_image DAG) | `MatrixMultiplyV4WithInfo` | B-connection retSz explosion, satisfy(self)=False | ≥9 (256q) | 🔍 analyzed, not fixed |

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

## Bug E: Level>=9 post_image DAG Corruption — B-Connection retSz Explosion (ANALYZED, NOT FIXED)

### Summary

At CFLOBDD level >= 9 (256+ qubits), `MatrixMultiplyV4WithInfo` (used in
`post_image` / gate application) produces states with corrupted amplitudes.
The root cause has been identified: at level 9, the gate's B-connection
decomposition is coarser (2 B-connections instead of 4 at level 8), causing
each B-connection to handle multiple gate exits.  This leads to incorrect
routing of gate amplitude entries, creating sign-split amplitudes that
`NormFormComplex` cannot merge.  The state DAG's returning-value count (retSz)
grows 2 -> 3 -> 8 across consecutive gate applications, eventually producing
total amplitude garbage.

This is distinct from the coefficient overflow (Bug D) -- the underlying
CFLOBDD DAG structure produced by `MatrixMultiplyV4WithInfoNode` is incorrect
at level 9, not just the evaluated values.

### Root Cause: B-Connection Decomposition at Level 9

The Hadamard gate CFLOBDD has 4 non-zero entries:
```
H|0> = +1/sqrt(2) |0> + 1/sqrt(2) |1>    (exit 0)
H|1> = +1/sqrt(2) |0> - 1/sqrt(2) |1>    (exit 1)
```

At CFLOBDD level k, the gate is composed of B-connections connecting
sub-CFLOBDDs.  The number and structure of B-connections depends on the level:

| Level | Qubit range  | B-connections | Qubits per B-connection | Mapping           |
|-------|-------------|---------------|------------------------|-------------------|
| 8     | 128-255     | 4             | each handles 1 exit    | 1-to-1 with exits |
| 9     | 256-511     | 2             | each handles 2+ exits  | coarse grouping   |

At level 9, each B-connection must handle **multiple** gate exits, each with
different amplitude patterns:

- **`c1.B[0]`** has return map `{0,1,2,3}` -- it handles ALL 4 non-zero gate exits
- **`c1.B[1]`** has return map `{4,5,6,7}` -- handles the remaining (mostly zero) exits

The problem: `c1.B[0]` must route 4 distinct sub-exit patterns from the
gate entries through a single B-connection.  The bilinear coefficients
`a[i] * b[j]` become large and complex, and the `NormFormComplex`
canonicalization fails:

1. The `+1/sqrt(2)` and `-1/sqrt(2)` gate entries get incorrectly routed to
   the same state component (|0> instead of the correct |1>)
2. `NormFormComplex` cannot merge these sign-split amplitudes that should
   have been routed to different exits
3. The state DAG accumulates extra returning values: retSz grows from the
   expected 2 (a proper product state) to 3 (precursor corruption) to 8
   (full corruption)

In contrast, at level 8, the 4 non-zero gate exits map 1-to-1 onto 4
B-connections.  Each B-connection handles exactly one amplitude pattern, so
no sign-splitting occurs, and `NormFormComplex` canonicalizes correctly.

### retSz Propagation Chain

The state accumulation proceeds step by step through the post_image chain
(256 qubits, all H gates producing |+>^256):

```
Step 201: H on qubit 201 (0-indexed)
  Input:  State retSz=2, Gate nB=2
  Output: State retSz=3   <-- PRECURSOR CORRUPTION
  (satisfy(self)=True at this point, but DAG structure is already anomalous)

Step 202: H on qubit 202
  Input:  State retSz=3, Gate nB=2
  Output: bb_old.retSz=8   <-- FULL CORRUPTION
  (satisfy(self)=False, amplitudes are garbage)
```

The jump from retSz=2 to retSz=3 at step 201 is the first sign of trouble.
A correct H-on-product-state application should produce retSz=2 exactly
(because |+> = H|0> is a product state, each qubit maintains one
superposition term).  The retSz=3 state already has sign-split amplitudes
that `NormFormComplex` could not merge.  When this damaged state is multiplied
by another H gate at step 202, the bilinear coefficients explode across the
3 * 2 sub-exit combinations, producing retSz=8 and total DAG corruption.

### Investigation Trail

#### Step-by-step post_image chain bisection

- Instrumented `repro_postimage_corruption.py` (n=256) to run the post_image
  chain one gate at a time, testing `satisfy(self)` and `retSz` at each step
- First corruption detected at step ~202 (H on qubit 201 of 256)
- This narrowed the anomaly to a specific `MatrixMultiplyV4WithInfo` call
  at CFLOBDD level 9

#### V4I Diagnostics at the Anomaly Point

- Added retSz tracking to `MatrixMultiplyV4WithInfoTopNode` to trace
  bilinear coefficient evaluation after `NormFormComplex` canonicalization
- At the first corruption point (step 201), the V4I operation processes
  4 B-connections internally (even though the top-level representation
  uses only 2 B-connections for the gate)
- The anomaly originates in the **B-connection processing at level 8**
  within the top-level level 9 call:
  - Gate retVals are identical at both level 8 and level 9:
    `[0.707, 0, -0.707, -0]`
  - But the B-connection **count** differs: 4 at level 8, 2 at level 9
  - At level 8, each of 4 B-connections handles exactly 1 gate exit,
    keeping return-map routing trivially correct
  - At level 9, each of 2 B-connections handles 2 gate exits, and
    `c1.B[0]` with return map `{0,1,2,3}` must route all 4 non-zero
    entries through a single B-connection.  The amplitude sign-splitting
    (+1/sqrt(2) vs -1/sqrt(2)) cannot be resolved by `NormFormComplex`,
    producing the retSz anomaly

#### Warmup Effect Confirmed

- Running Gram-Schmidt (`span_qops`) on an unrelated state **before** the
  post_image chain changes the result at level 9
- This is NOT a fix -- it merely alters the global CFLOBDD unique-table
  state, which changes which internal DAG nodes get hash-consed during
  `MatrixMultiplyV4WithInfoNode` construction
- Different warmup states produce different corruption patterns and shift
  the exact step where corruption first appears
- This confirms the root cause is in the unique-table / hash-consing /
  canonicalization interaction at high B-connection fan-in, not a wrong
  constant, threshold, or conversion

### Why Level 8 Works But Level 9 Does Not

At level 8 (128-255 qubits), the Hadamard gate's 4 non-zero exits map
cleanly onto 4 B-connections (1 exit per B-connection).  Each B-connection
only needs to handle a single amplitude pattern: `+1/sqrt(2)` or
`-1/sqrt(2)` individually.  `NormFormComplex` correctly canonicalizes
these isolated patterns.

At level 9 (256-511 qubits), the gate uses only 2 B-connections for 4
exits.  The B-connections must combine multiple distinct amplitude
patterns through the same return map, and the bilinear coefficient
interaction with `NormFormComplex` produces incorrect state DAG
structures.

At level 10 (512-1023 qubits), the problem may appear at different steps
but the mechanism is the same: coarser B-connection decomposition leads
to incorrect amplitude routing.

### Similarities and Differences from Bug D

| Aspect | Bug D (V4 coeff overflow) | Bug E (B-connection retSz) |
|--------|---------------------------|----------------------------|
| Location | `MatrixMultiplyV4TopNode` (evaluate) | `MatrixMultiplyV4WithInfoTopNode` (construct) |
| Operation | Convert bilinear coeff to float | Build result DAG via `NormFormComplex` |
| Root cause | `cpp_int -> uint64` overflow at 2^64 | B-connection retSz expansion in NormFormComplex |
| Level threshold | >= 7 (64q) | >= 9 (256q) |
| Symptom | dot(self,self) = 0 | satisfy(self) = False, dot = garbage, retSz > 2 |
| Fix complexity | Simple (1-line type change) | Deep (canonicalization / B-connection redesign) |

### Affected Circuits

| Circuit | Qubits | CFLOBDD Level | Status | Root Cause |
|---------|--------|---------------|--------|------------|
| grover32-plus/zero-linear | 64 | 7 | PASS | Levels < 9 unaffected |
| grover64-plus/zero-linear | 128 | 7-8 | PASS | Levels < 9 unaffected |
| grover128-zero-linear | 256 | 8-9 | PASS | Final state = |0>, trivial DAG |
| grover128-plus-linear | 256 | 9 | FAIL | B-connection retSz explosion |
| grover150-plus-linear | 300 -> 512 | 10 | FAIL | Same mechanism at level 10 |
| grover150-zero-linear | 300 -> 512 | 10 | PASS | Final state = |0>, trivial DAG |
| repro_postimage_corruption n=256 | 256 | 9 | FAIL | B-connection retSz explosion |

### Path to Fix

The root cause is architectural: `NormFormComplex` cannot merge sign-split
amplitudes when a single B-connection return map routes multiple gate exits
through the same connection.  Potential approaches:

1. **Refine `NormFormComplex` canonicalization** -- teach it to recognize
   and merge sign-split amplitude pairs (e.g., `+alpha|s>` and `-alpha|s>`)
   that should be routed to different exits but end up on the same component.
   This is the most direct fix but may not cover all edge cases.

2. **Restructure B-connection decomposition at high levels** -- ensure each
   B-connection handles at most 1 gate exit.  At level 9+, this would require
   more B-connections (4 instead of 2 for the Hadamard gate), potentially
   changing the CFLOBDD construction algorithm.

3. **Alternative gate encoding for level 9+** -- encode the Hadamard gate
   differently when the number of distinct non-zero entries exceeds the
   B-connection count, avoiding the problematic routing.

4. **Tactical workaround: force lower level** -- for 256-qubit circuits,
   it may be possible to force CFLOBDD level 8 (which works correctly) by
   adjusting qubit padding.  This is a stopgap, not a fix.

## Platform Dependence: Hash-Consing Sensitivity to Memory Layout

### Discovery (2026-08-11)

The **GHZ benchmark suite** (`benchmark/converted_from_qai_ghz/`) revealed that
Bug C manifests **platform-dependently** at CFLOBDD level 8:

| Platform | ghz002–200 (L≤8) | ghz250–500 (L=8) | ghz550–800 (L≥9) |
|----------|:---:|:---:|:---:|
| macOS (Apple Silicon) | ✅ PASS | ❌ normalize() assert | ❌ float_next<double> overflow |
| Linux (x86-64) | ✅ PASS | ✅ PASS | ❌ float_next<double> overflow |

On macOS, ghz250–500 crash at `quantum_operation.hpp` in `normalize()`:
```
assert(abs(amp.imag()*dimfactor) < 1e-8 && abs(amp.real()*dimfactor) < 1e-8)
```
The Identity-multiply trick produces `retMapSz=1` with a garbled non-zero
amplitude — the transpose corruption destroys row 0 of the intermediate matrix.

On Linux, the same circuits pass.  The state produced by the identical gate
sequence is self-consistent.

### Root Cause

CFLOBDD relies on **hash-consing** (pointer-based hashing of internal DAG
nodes) for canonicalization.  Hash values depend on heap addresses returned
by `malloc`.  macOS and Linux use different `malloc` implementations
(macOS: nanomalloc; Linux: ptmalloc), producing different address layouts
and therefore **different hash-consing decisions**.

This means the *same* CFLOBDD operation can produce **different internal
DAG topologies** on different platforms.  Some topologies happen to route
`MatrixTranspose` assignments to correct exits; others route them to wrong
exits, destroying row 0.

### Connection to the Warmup Effect

This is the same mechanism as the documented **warmup effect** (Bug E
investigation trail, point 4): running Gram-Schmidt on an unrelated state
before the post_image chain changes the unique-table state, which changes
which internal DAG nodes get hash-consed.  Different warmup states produce
different corruption patterns.

Cross-platform memory layouts are simply a different source of unique-table
state variation, leading to the same class of non-deterministic behavior.

### Implication

Bug C is **not fixed** by the Identity-multiply trick.  The trick changes
the DAG topology in a way that *reduces* the probability of hitting the
corruption, but the underlying `MatrixTranspose` routing bug remains.
Whether a given circuit hits the corruption depends on:

1. CFLOBDD level (≥8 required)
2. Physical qubit count and gate sequence (determines DAG structure)
3. Heap memory layout (platform + allocator state)
4. Unique-table fill level (prior computation history)

### GHZ Circuit as a Bug C Detector

The GHZ circuit is a particularly effective detector for Bug C because:

- It has a "star" CNOT structure (all CNOTs target the same qubit),
  creating a distinctive CFLOBDD DAG topology at level 8
- The state after init H gates (`|+⟩^k|0⟩`) is a product state whose
  `normalize()` path exercises the Identity-multiply → transpose code path
- The transition between working (ghz200, k=199) and broken (ghz250, k=249)
  on macOS shows the sensitivity to physical qubit count within the same
  CFLOBDD level

On macOS, the threshold is between 200 and 250 physical qubits (both at
qNum=256, level 8).  On Linux, the threshold is above 500 physical qubits
(moving from level 8 to level 9, where Bug E takes over).

## Other Fixes Applied

### satisfy() semantics preserved

Investigation confirmed `satisfy()` logic (`lowerBound ⊆ spec ⊆ upperBound`)
is correct.  The test failure was caused by `parse_qiskit_cir_lazy` setting
`upperBound = initial_op` (overwriting the default Identity), which broke
the `spec ⊆ upperBound` check.  Fix: removed the line.

### Reproducer improvements

- `repro_transpose_bug.py`: covers Gram-Schmidt / dot / normalize correctness
- `repro_postimage_corruption.py`: covers post_image chain self-consistency
