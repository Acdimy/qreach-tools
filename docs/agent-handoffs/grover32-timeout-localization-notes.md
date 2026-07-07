# Grover32 Timeout Localization Notes

**Last updated:** 2026-07-04 (Phase 2 complete, root cause localized to CFLOBDD Reduce)

## Phase 1: Normalization Fix (Complete)

### Grover64 Reference
- Full lazy parse: `parsed 959 1 0.254s`. No profile lines above 0.05s threshold.

### Grover32 Plus — Original Pathology
- Bisect found `smallest_failing_prefix=102` at 20s timeout.
- Last completed profile line: `idx=100 op=x qubits=[27] elapsed=3.247749s total_locations=102`.
- Next instruction: `idx=101 x [26]`.
- Adjacent passing prefix 101: `time_parse=19.345969249960035`.
- **Fixed:** After `knownUnitNorm` normalization skip: prefix 102 `time_parse=1.48s`. Full circuit still >120s.

### Grover32 Zero
- No failing prefix through 179 instructions. `time_parse=0.08474900003056973`.
- CCX-heavy region (indices 96-104: `ccx [60,61,62]` etc.) parses rapidly.
- Initial state |0⟩^63 = single nonzero amplitude → CFLOBDD stays compact.

### Grover32 Plus Linear
- Prefix 102: `time_parse=0.02s`. Symmetric linear qubit ordering → DAG sharing.

### Grover64 Plus
- Full parse fast. Level=8 (qNum=127→128), different CFLOBDD topology.

## Phase 2: CFLOBDD Matrix Multiply Investigation (Complete)

### Gate profiling evidence (after normalization fix, `QREACH_GATE_PROFILE=1`)

Late-prefix post-image breakdown for Grover32-plus prefix 102:

| idx | gate | concretize | multiply | normalize | total |
|-----|------|-----------|----------|-----------|-------|
| 89 | CCX [22,23,43] | 23μs | 39,421μs | 0μs | 663ms |
| 94 | CCX [24,25,44] | 21μs | 110,251μs | 0μs | 1,471ms |
| 99 | CCX [26,27,45] | 26μs | 348,382μs | 0μs | 3,387ms |
| 100 | X [27] | 52μs | 299,911μs | 0μs | 3,266ms |

Full-circuit late instructions (> idx 109):

| idx | gate | multiply (est.) | total (est.) |
|-----|------|----------------|-------------|
| 109 | CCX [30,31,47] | ~3,500,000μs | ~4,700ms |
| 112 | CCX [32,33,48] | ~4,000,000μs | ~5,100ms |
| 118 | CCX [44,45,54] | ~4,000,000μs | ~5,200ms |

**Key finding:** Normalization is fully skipped (normalize_us=0 for all). The bottleneck is `MatrixMultiplyV4WithInfo` — CFLOBDD matrix-vector multiply.

### No-lazy mode
Even worse than lazy. `computingFixedPointPost()` runs >5 minutes without completing. No-lazy creates many more locations, and fixed-point iteration computes far more post-image calls across all of them.

### Phase 2 Deeper Profiling: `MatrixMultiplyV4WithInfoTopNode` time breakdown

Added temporary instrumentation to split `MatrixMultiplyV4WithInfoTopNode` into three phases:
- **recurse_us**: Time in `MatrixMultiplyV4WithInfoNode` (recursive CFLOBDD DAG traversal)
- **eval_us**: Time in the top-level O(N×M) evaluation loop (lines 991-1027)
- **reduce_us**: Time in `tempHandle.Reduce(reductionMapHandle, v.Size(), true)`

**Grover32-plus prefix 102 (slow) — late CCX gates:**

| Entry | recurse_us | reduce_us | total_us | c1_sz | mm_sz |
|-------|-----------|-----------|---------|-------|-------|
| 79 | 3,382 | 2,344 | 5,734 | 2 | 2 |
| 84 | 6,796 | 6,628 | 13,431 | 2 | 2 |
| 89 | 13,563 | 22,958 | 36,528 | 2 | 2 |
| 94 | 27,032 | 81,121 | 108,160 | 2 | 2 |
| 99 | 59,524 | **298,782** | 358,321 | 2 | 2 |
| 100 | 69,494 | **276,772** | 346,279 | 2 | 2 |

**Control cases (Grover32-zero, Grover32-linear):**

| Metric | Grover32-plus | Control cases | Ratio |
|--------|--------------|---------------|-------|
| Max recurse_us | 69,494μs | 152μs | 457x |
| Max reduce_us | 298,782μs | 58μs | **5,150x** |
| Max total | 358,321μs | 198μs | **1,800x** |
| c1_sz (vector leaf amps) | 2 | 2 | same |
| mm_sz (MatMultMap entries) | 2 | 2 | same |

## Root Cause Analysis

### What IS the bottleneck
The **Reduce** operation (`CFLOBDDNodeHandle::Reduce`, called at the end of `MatrixMultiplyV4WithInfoTopNode`) consumes ~80% of the multiply time for slow calls. The Reduce collapses the intermediate CFLOBDD DAG (built during recursive traversal) into a canonical form by merging nodes that map to the same leaf-amplitude value.

### Mechanism
1. Starting from H^⊗63 (uniform superposition), the vector CFLOBDD has maximal amplitude diversity
2. Applying CCX gates with nonlinear qubit ordering (e.g., [26,27,45]) that cross the CFLOBDD split boundary at level=7 causes the intermediate DAG to become internally fragmented
3. The recursive traversal (`MatrixMultiplyV4WithInfoNode`) builds a **massive intermediate DAG** with many internal nodes
4. Even though the DAG is huge, it eventually maps to only **2 distinct leaf amplitudes** (`c1_sz=2`)
5. The Reduce operation spends ~300ms collapsing this massive-but-redundant DAG into just 2 leaf values

### Why only Grover32-plus (H^⊗63 + nonlinear + level=7)?
- **Grover32-zero** (|0⟩^63): single nonzero amplitude → CFLOBDD stays sparse throughout → Reduce is trivial (58μs max)
- **Grover32-linear**: symmetric qubit ordering (0..63) → CFLOBDD DAG sharing works → Reduce stays fast
- **Grover64-plus**: level=8 with 64/64 split → different CFLOBDD topology with better sharing
- **Grover32-plus**: H^⊗63 + nonlinear ordering + level=7 with 63-in-64 padding = worst-case DAG fragmentation

### Why it's hard to fix
- Reduce is **essential** to CFLOBDD canonicalization — it cannot be skipped
- The fragmentation is a fundamental CFLOBDD DAG property, not a simple algorithmic bug
- The issue only manifests under specific circuit/initial-state combinations
- Any fix would require deep CFLOBDD architecture changes

### Practical workarounds
- Use Grover32-zero when the initial state is |0⟩^n (most common case)
- Use Grover32-plus-linear with explicit q[64] for symmetric ordering
- Use Grover64 (naturally at level=8, fast)
- For circuits that must use H^⊗n with nonlinear qubit ordering at non-power-of-2 qubit counts, expect CFLOBDD performance issues

## Files Modified for Investigation

| File | Status | Purpose |
|------|--------|---------|
| `cflobdd/CFLOBDD/matrix1234_complex_float_boost_top_node.cpp` | **temporary** | Added `<chrono>`, timing instrumentation for recurse/eval/reduce breakdown |
| `quantum_operation.hpp` | **committed** | `knownUnitNorm` flag, `isNormPreservingGate()`, normalization skip, fallback propagation |
| `python_pkg/parse_qiskit.py` | **committed** | `QREACH_PARSE_PROFILE` env-gated parser profiling |
| `python_pkg/workflow_tests/debug_grover32_ccx_pathology.py` | **committed** | Prefix timing, bisect, timeout worker |
| `python_pkg/workflow_tests/test_grover32_lazy_prefix_performance.py` | **committed** | Focused pytest regression |

### Cleanup needed
The temporary profiling instrumentation in `matrix1234_complex_float_boost_top_node.cpp` should be reverted (or guarded by a permanent profiling flag) before merging. It adds `<chrono>` include and timing/profiling blocks that are not needed in production.

## Next Steps

1. **Revert** the temporary CFLOBDD profiling instrumentation in `matrix1234_complex_float_boost_top_node.cpp`
2. **Clean up** the `scripts/patch_multiply_top_node.py` helper script
3. **Document** the pathological pattern: H^⊗n initial state + nonlinear CCX at non-power-of-2 qubit count
4. **Investigate other bad cases** to see if they share the same root cause
