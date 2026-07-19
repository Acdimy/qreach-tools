---
name: cflobdd-transpose-dag-corruption
description: CFLOBDD transpose produces wrong results for vectors created through GramSchmidt PairProduct path when retMapSz>2
metadata:
  type: project
  status: under-investigation
  branch: qts-rollback
  date: 2026-07-18
---

# CFLOBDD Transpose DAG Corruption Bug

## Symptoms

`dqc_pe` circuits with Pauli error injection crash at:
```
Assertion failed: (resMap.Size() <= 2), function normalize, file quantum_operation.hpp
```

The crash occurs in `loc.satisfy(expected)` → `compare()` → `disjunction()` → `GramSchmidt` → `normalizeInline()` → `normalize()`.

## Reproducer

```bash
cd python_pkg
../../.venv/bin/python -m invoke build-pybind11
../../.venv/bin/python repro_minimal.py
```

Expected: `span_qops` passes (step 1 OK), `satisfy` crashes (step 3 assertion).

The minimal repro uses `dqc_pe_2.qasm` from `benchmark/dqc_pe/`. Attempts to create a smaller circuit-based repro (simple X/H/Z gates) failed because the resulting states have `retMapSz=2` and don't trigger the bug. The problematic vector structure can only be produced through the GramSchmidt `operator+` path.

## Root Cause

`normalize()` computes `<content|content>` via:
```cpp
CFLOBDD_COMPLEX_BIG c_conj = MatrixConjugate(content);
c_conj = MatrixTranspose(c_conj);  // ← corrupts column-vector structure
auto mulres = MatrixMultiplyV4(c_conj, content);
```

`MatrixTranspose` on a column vector (only col 0 non-zero) should produce a row vector (only row 0 non-zero). But when the input has `retMapSz > 2` AND was constructed through GramSchmidt's `ivec->content + neg_proj` (CFLOBDD `MatrixPlus`/`PairProduct`/`ApplyAndReduce`), transpose produces multi-row output.

Example corruption (from `[transpose-test #9]`):
- Input: column vector with entries at rows 2=(1,0.0029) and 10=(-1,-0.0029)
- After conjugate: same structure (column vector) ✓
- After transpose: value at row 2's conjugate leaks to rows 0-7 at col 2 ✗

## Key Findings

1. **Not about bit-span**: vectors with entries at rows differing by bit 2 (rows 2&6) work fine; rows differing by bit 3 (rows 2&10) also work when constructed via simple `make_basis(a) + scalar * make_basis(b)`.

2. **Not about opposite signs**: vectors with opposite-sign entries at rows 2&6 (constructed via basis addition) transpose correctly.

3. **Specific to GramSchmidt-constructed vectors**: The bug only manifests on vectors created through the GramSchmidt orthogonalization path (`ivec->content + neg_proj`), not on identically-valued vectors created by simple basis-vector addition. The internal CFLOBDD DAG structure differs between the two construction paths even when the printed matrix values are identical.

4. **`returnMapHandle.Size()` cannot detect structural corruption**: A vector can have the correct number of distinct leaf values but wrong structural properties (e.g., non-zeros in multiple columns).

5. **Earlier `H * content` fix was related**: `normalize()` originally did `c1 = MatrixMultiplyV4WithInfo(H, content)` where `H = ApplyGateF(..., MkIdRelationInterleaved)` — an identity matrix. This `MatrixMultiplyV4WithInfo` also corrupted column-vector structure in the same way. Fixed by replacing with direct `conj(transpose(content)) * content` pattern (same as `dot()`). This resolved the `span_qops` crash but `satisfy` still crashes.

## Modified Files

### `quantum_operation.hpp` — normalize() fix (committed)

Changed from:
```cpp
auto H = ApplyGateF(..., MkIdRelationInterleaved);
c1 = MatrixMultiplyV4WithInfo(H, content);
// compute norm from c1, return scaled c1
```
To:
```cpp
// Compute <content|content> directly, following dot()'s pattern
c_conj = MatrixConjugate(content);
c_conj = MatrixTranspose(c_conj);
mulres = MatrixMultiplyV4(c_conj, content);
// extract norm, return scaled content
```

### `cflobdd/CFLOBDD/matrix1234_node.cpp` — transpose fix (separate bug, kept, tagged)

Lines 2951, 2969: Changed `m1.AddToEnd(v)` to `m1.AddToEnd(return_handle.LookupInv(v))` in level-1 transpose fork case. Old exit values were used instead of return_handle indices. Masked when retMapSz≤2. Tagged with `// FIXME(qts-rollback):` comments. This is a real bug but does not fix the current crash.

### Debug instrumentation (current state)

`normalize()` in `quantum_operation.hpp` prints `[transpose-test #N]` output when `retMapSz > 2`:
- content (input column vector)
- after conjugate
- after transpose

See `grep 'transpose-test' quantum_operation.hpp` for the debug code locations.

## Next Steps for Investigation

1. **Deep-dive into `MatrixTranspose`**: Trace the interleaved CFLOBDD transpose algorithm (`matrix1234_node.cpp`) with the specific DAG structure produced by GramSchmidt's `MatrixPlus`/`PairProduct`/`ApplyAndReduce`. The corruption source is in the DAG traversal that handles column→row variable swapping.

2. **Compare DAG structures**: Instrument a comparison between two structurally different column vectors that print identically:
   - Vector A: `make_basis(2) + (-1) * make_basis(10)` (transpose works)
   - Vector B: GramSchmidt-produced vector at rows 2,10 (transpose fails)
   - Compare their `entryPointHandle` DAG topology, `returnMapHandle` contents, and internal node structures.

3. **Possible workaround**: `normalize()` could potentially bypass transpose entirely by using a different norm-computation strategy (e.g., computing `dot(*this)` which goes through the existing `dot()` function and may avoid the problematic transpose path).

4. **Investigate `MatrixMultiplyV4WithInfo`**: The original `H * content` corruption suggests `MatrixMultiplyV4WithInfo` may share the same root cause as transpose — both operate on the interleaved CFLOBDD DAG structure and may mishandle certain node configurations created by PairProduct.

## Key Files

- `quantum_operation.hpp:1004-1048` — `normalize()` (with debug)
- `quantum_operation.hpp:1637-1649` — `GramSchmidt()` (the `ivec->content + neg_proj` path)
- `cflobdd/CFLOBDD/matrix1234_node.cpp` — `MatrixTranspose` implementation
- `cflobdd/CFLOBDD/cross_product.cpp` — `PairProduct` implementation
- `cflobdd/CFLOBDD/cflobdd_top_node_t.cpp:329-386` — `ApplyAndReduce` / `MkPlusTopNode`
- `cflobdd/CFLOBDD/vector_complex_float_boost_top_node.cpp:251-291` — `VectorPrintColumnMajor` (debug print utility)
- `python_pkg/repro_minimal.py` — minimal reproducer
- `docs/agent-handoffs/` — directory for handoff documents

## Contact / References

See `docs/agent-handoffs/` for related performance investigation documents (grover32, benchpress).
These are mostly about Reduce fragmentation and amplitude explosion — separate from this DAG corruption bug.
