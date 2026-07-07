# Benchpress-Medium Floating-Point Amplitude Explosion

> **Status:** Root cause identified (amplitude explosion from parameterized rotation gates).  
> **Last updated:** 2026-07-04

**Goal:** Explain why circuits with many floating-point parameterized rotation gates (rz/ry) time out even with small qubit counts (~25) and |0⟩ initial state.

**Affected files:**
- `python_pkg/benchmark/benchpress-medium/not-supported/knn_n25.qasm` — 24 ry, 12 swap, 2 h, q[25]
- `python_pkg/benchmark/benchpress-medium/not-supported/knn_n25_transpiled.qasm` — 172 rz, 96 cx, 74 sx, q[25]
- `python_pkg/benchmark/benchpress-medium/not-supported/ising_n26_transpiled.qasm` — 101 rz, 50 cx, 27 s, 26 sx, q[26]
- `python_pkg/benchmark/benchpress-medium/not-supported/swap_test_n25_transpiled.qasm` — 196 rz, 96 cx, 74 sx, q[25]
- `python_pkg/benchmark/benchpress-medium/not-supported/wstate_n27_transpiled.qasm` — 104 sx, 52 rz, 52 cx, 28 s, 1 x, q[27]

## Root Cause

### Mechanism: Exponential Leaf-Amplitude Explosion

Each parameterized rotation gate `ry(θ)` or `rz(λ)` is lowered to a `U3` gate with floating-point parameters. These parameters induce **incommensurate amplitudes** (e.g., `cos(θ/2)` and `sin(θ/2)`) that the CFLOBDD's `NormFormComplex` rounding (30 decimal places) cannot merge.

When U3 gates are applied to **distinct qubits** from a product-state initial vector, each gate **doubles** the number of distinct leaf-amplitude values in the CFLOBDD:

```
State after k U3 gates on distinct qubits:
  distinct_amplitudes ≈ 2^k
```

The `MatrixMultiplyV4WithInfoTopNode` evaluation loop processes `result_mm_sz ≈ 2^k + 1` MatMultMap entries, each with an inner `O(M)` loop. Total work is `O(2^k × M)` per multiply, growing exponentially with each additional U3 gate on a new qubit.

### Evidence: knn_n25.qasm prefix=10 profiling

`QREACH_GATE_PROFILE=1` with `qreach-multmap-breakdown` instrumentation:

| gate | qubit idx | c1_rmap_sz | c2_rmap_sz | result_mm_sz | multiply_us | growth |
|------|-----------|-----------|-----------|-------------|------------|--------|
| U3 | 1 | 5 | 2 | 3 | 197 | — |
| U3 | 2 | 5 | 3 | 5 | 79 | 1.7x |
| U3 | 3 | 5 | 5 | 9 | 100 | 1.8x |
| U3 | 4 | 5 | 9 | 17 | 159 | 1.9x |
| U3 | 5 | 5 | 17 | 33 | 248 | 1.9x |
| U3 | 6 | 5 | 33 | 65 | 450 | 2.0x |
| U3 | 7 | 5 | 65 | 129 | 763 | 1.9x |
| U3 | 8 | 5 | 129 | 257 | 1,376 | 1.9x |
| U3 | 9 | 5 | 257 | 513 | 2,331 | 1.8x |
| U3 | 10 | 5 | 513 | 1025 | 4,869 | 2.0x |

**Pattern:** `result_mm_sz` doubles with each new U3 gate on a distinct qubit (`2^k + 1`). `multiply_us` grows proportionally.

At index 15 (observed before timeout): `multiply_us = 861,603μs` (861ms).

Projected for all 24 U3 gates in knn_n25: `result_mm_sz ≈ 2^24 + 1 ≈ 16.7 million`, with multiply time in the **hours** range.

### Why this differs from Grover32

| Dimension | Grover32-plus | benchpress (knn_n25 etc.) |
|-----------|--------------|---------------------------|
| Initial state | H^⊗63 | \|0⟩^n |
| Gate types | H, X, CCX (no float params) | U3 (rz/ry with float params) |
| Amplitudes created | ±1/√2 (commensurate) | cos(θ/2), sin(θ/2) (irrational) |
| CFLOBDD merging | Can share — amplitudes are ±1/√2 | Cannot merge — all unique |
| c1_sz (leaf values) | Always 2 | Grows with distinct qubit gates |
| result_mm_sz | Always 2 | **2^k exponential growth** |
| Bottleneck | Reduce (80% of time) | Evaluation loop (O(2^k × M)) |
| Fix difficulty | Deep CFLOBDD change | **Need parameterized amplitude handling** |

### Why transpiled versions are also affected

`knn_n25_transpiled.qasm`, `swap_test_n25_transpiled.qasm`, etc. contain `rz(pi/2)` gates. Although `pi/2` is a rational multiple of π, the transpiler decomposes these into hardware-native gates that may reintroduce floating-point parameters. The same exponential accumulation occurs.

### Why `ising_n26` and `wstate_n27` are affected

Both contain `rz` gates with floating-point parameters. The `rz` gates cause the same amplitude diversity accumulation, though the growth rate depends on how many distinct qubits receive rotation gates before CX/SWAP gates entangle them.

## Practical Impact

- Circuits with **parametric rotation gates on distinct qubits** cannot scale beyond ~15 such gates in the current CFLOBDD backend
- Transpiled circuits (with hardware-native basis gates) inherit this limitation from `rz`/`sx` decompositions
- The issue is **independent of qubit count** — even 10-qubit circuits would time out with enough parameterized rotations

## Potential Fix Directions

### Option A: Symbolic/Parameterized CFLOBDD Leaves

Instead of storing floating-point amplitudes as `BIG_COMPLEX_FLOAT` values, store them as symbolic expressions (e.g., `cos(θ₁/2)·sin(θ₂/2)`). The CFLOBDD would share structurally identical subtrees regardless of parameter values.

**Pros:** Handles arbitrary parameterized circuits. General solution.
**Cons:** Major CFLOBDD architecture change. Requires symbolic algebra in the DAG evaluation.

### Option B: Rational Approximation of Parameters

Round floating-point gate parameters to rational approximations before lowering. For example, `ry(0.93346815)` → `ry(π/3.367)` → rational fraction. This would allow `NormFormComplex` to merge more amplitudes.

**Pros:** Simple parser-level change. No CFLOBDD changes.
**Cons:** Changes circuit semantics. May not help with transpiled circuits where parameters are already approximations.

### Option C: Detect Tensor-Product Structure in Parser

Recognize when gates are applied to independent qubits and maintain a tensor-product decomposition rather than a single monolithic CFLOBDD vector.

**Pros:** Handles the common case of rotation gates on distinct qubits.
**Cons:** Complex parser change. Breaks when CX/SWAP entangle qubits.

### Option D: Accept as Known Limitation

Document the limitation and use alternative simulators for parameterized circuits. The CFLOBDD backend is optimized for circuits with discrete gate sets (Clifford + T, Grover, etc.) where amplitudes are commensurate.

**Pros:** No code changes. Honest about tool scope.
**Cons:** Limits applicability of the tool.

## Next Steps

1. Determine which fix direction (if any) is in scope for the current project phase
2. If pursuing Option B or C, start with a Python-level prototype before modifying CFLOBDD internals
3. Add these benchmark files to a known-unsupported list with clear documentation of the limitation
