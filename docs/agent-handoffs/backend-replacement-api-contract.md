# QReach Backend Replacement — Document A: New DD Backend API Contract

> **Audience:** AI agent working on the *new* decision-diagram backend tool.
> **Status:** Draft — API contract specification.
> **Last updated:** 2026-07-08

## 1. Project Context

QReach is a quantum model-checking / reachability-analysis tool. It parses Qiskit circuits into a transition system, propagates quantum states through gate operations, and checks CTL/LTL properties via NuSMV.

The current backend is **CFLOBDD** (Context-Free-Language Ordered Binary Decision Diagrams). The goal is to replace it with an alternative decision-diagram backend while keeping the higher-level quantum semantics layer (`quantum_operation.hpp`) and the Python workflow layer intact.

### Architecture (simplified)

```
Python:  qctl.py  →  parse_qiskit.py  →  pyqreach (pybind11)
                                              │
C++:     transition_system.hpp                 │
              │                                │
         quantum_operation.hpp  ←──────────────┘
              │
         [DD Backend]   ←  THIS IS WHAT YOU REPLACE
```

The DD backend is consumed **only** through `quantum_operation.hpp`. The transition system, parser, and Python layer never touch the DD backend directly.

## Task Summary

Implement a C++ namespace/class that provides the API surface described below. The implementation must be a drop-in replacement for the CFLOBDD backend used by `quantum_operation.hpp`.

---

## Required API Surface

### 1. Core Types

```cpp
// Scalar type for complex amplitudes. Must support:
//   - real(), imag()
//   - operator==, !=
//   - operator* (with scalar)
//   - abs(), norm()
//   - default construction
using DDComplex = /* your complex scalar type */;

// The core DD type. Must support:
//   - operator+ (DD + DD, pointwise addition — used by InitializeWithAmplitudes / span)
//   - operator* (scalar * DD, left scalar multiplication)
//   - operator==
//   - copy construction / assignment (reference-counting semantics preferred)
using DD = /* your DD type */;

// NOTE (2026-08-13, confirmed with LimTDD feasibility review):
//   pointwise DD * DD is NOT required. The QReach semantic layer
//   (quantum_operation.hpp) only ever multiplies a scalar by a DD; it never
//   does element-wise DD*DD. If your backend supports it, fine — but the
//   contract no longer requires it.
```

### 2. Vector Operations (`DDVector` namespace)

```cpp
namespace DDVector {
    // --- Initialization ---
    void Initialize();

    // --- Construction ---
    DD MkBasisVector(unsigned int level, unsigned int index);
    //   Create the computational basis state |index⟩ in a 2^level dimensional space.

    DD MkBasisVector(unsigned int level, std::string bitstring);
    //   Create a basis state from a bitstring (e.g. "0101").

    DD NoDistinctionNode(unsigned int level, DDComplex val);
    //   Create a constant vector where every amplitude equals `val`.

    DD InitializeWithAmplitudes(unsigned int qnum, std::vector<double> amps);
    //   Create a vector from a raw amplitude array. `amps` is interleaved
    //   real/imaginary: [re0, re1, ..., im0, im1, ...].
    //   The vector must be normalized (||v|| = 1).

    // --- Format Conversion ---
    DD VectorToMatrixInterleaved(DD vec);
    //   Convert a vector representation to a matrix representation.
    //   In CFLOBDD this changes the variable ordering from single-vocabulary
    //   to interleaved double-vocabulary.  Your backend may implement this
    //   as a no-op if your matrix multiply works directly on vectors.

    // --- Queries ---
    int GetLevel(DD c);
    //   Return the CFLOBDD level.  Dimension = 2^level.

    bool IsApproximatelyZero(DD c, double threshold = 1e-8);
    //   Return true if the vector is approximately the zero vector.

    DDComplex ExtractSingleAmplitude(DD c);
    //   Assert that the DD has exactly one non-zero leaf value, return it.
    //   Used by dot() and normalize().

    std::vector<std::pair<unsigned int, DDComplex>> GetNonZeroAmplitudes(
        DD c, double threshold = 1e-8);
    //   Return a list of (index, amplitude) for all non-zero entries.

    // --- Normalization ---
    DD Normalize(DD c);
    //   Return a normalized copy of the vector.

    // --- Inner product (added 2026-08-13, after LimTDD feasibility review) ---
    DDComplex InnerProduct(DD a, DD b);
    //   Return <a|b> = conj(a) · b. Backend-native:
    //     - CFLOBDD: Transpose→Conjugate→MatrixMultiply→ExtractSingleAmplitude
    //       (vectors are stored in matrix form).
    //     - LimTDD:  tensor contraction cont(conj(a), b) directly — NO vector
    //       padding to matrix form. This is the point: inner product is a
    //       contraction, not a matrix-multiply on padded vectors.

    // --- Debug ---
    void VectorPrintColumnHead(DD c, std::ostream& out);
    //   Print the vector in a human-readable format.
}
```

### 3. Matrix Operations (`DDMatrix`)

```cpp
namespace DDMatrix {
    // --- Initialization ---
    void Initialize();

    // --- Single-qubit gate matrices (at given level) ---
    DD MkIdRelation(unsigned int level);          // Identity I
    DD MkWalsh(unsigned int level);               // Hadamard H
    DD MkNegation(unsigned int level);            // Pauli X
    DD MkPauliY(unsigned int level);              // Pauli Y
    DD MkPauliZ(unsigned int level);              // Pauli Z
    DD MkSGate(unsigned int level);               // S = diag(1, i) — exact i, NOT cos/sin(π/2)
    DD MkPhaseShift(unsigned int level, double theta);  // P(θ)
    DD MkU3(unsigned int level, std::vector<double> params);  // U3(θ,φ,λ)
    DD MkArbitrary(unsigned int level, std::vector<double> params); // 8-param unitary

    // --- Multi-qubit gate matrices ---
    DD MkCNOT(unsigned int level, unsigned int n, long ctrl, long tgt);
    DD MkCCNOT(unsigned int level, unsigned int n, long c1, long c2, long tgt);
    DD MkSwap(unsigned int level, long i, long j);
    DD MkiSwap(unsigned int level, long i, long j);
    DD MkCP(unsigned int level, long ctrl, long tgt, double theta);

    // --- Single-qubit gate on n-qubit system ---
    DD MkSingleQubitGateOnN(unsigned int n, unsigned int target,
                            DD(*gate1q)(unsigned int));
    DD MkSingleQubitGateOnNWithParam(unsigned int n, unsigned int target,
                            DD(*gate1q)(unsigned int, double), double theta);
    DD MkSingleQubitGateOnNWithParamVec(unsigned int n, unsigned int target,
                            DD(*gate1q)(unsigned int, std::vector<double>),
                            std::vector<double> v);
    //   Given a 1-qubit gate constructor, build the n-qubit matrix
    //   that applies it to `target` and identity to all other qubits.
    //   The WithParam / WithParamVec variants carry the gate's angle/parameter
    //   arguments (needed for PhaseShift / U3 / Arbitrary). These replace
    //   CFLOBDD's ApplyGateF / ApplyGateFWithParam / ApplyGateFWithParamVec.

    // --- Matrix operations ---
    DD KroneckerProduct(DD a, DD b);              // Tensor product ⊗
    DD MatrixMultiply(DD a, DD b);                // Matrix multiplication
    DD MatrixMultiplyWithVector(DD gate, DD vec); // gate|ψ⟩ — HOT PATH
    DD Conjugate(DD c);                           // Complex conjugate
    DD Transpose(DD c);                           // Transpose
}
```

### 4. Key Design Notes

**Level / Dimension convention.** QReach uses `level` where `dimension = 2^level`. The `qNum` in `quantum_operation.hpp` is always a power of 2 (padded from the physical qubit count). For example, 5 physical qubits → `qNum = 8` → `level = 3`.

**VectorToMatrixInterleaved.** This is a CFLOBDD-specific format conversion. If your backend represents vectors and matrices uniformly, this can be a no-op. The caller expects the result to be usable as input to `MatrixMultiplyWithVector`.

**MatrixMultiplyWithVector is the hot path.** This is called for every gate application during lazy parsing and fixed-point iteration. The current CFLOBDD bottleneck is in `MatrixMultiplyV4WithInfoTopNode` (Reduce + evaluation loop). Your implementation should be optimized for this operation.

**Reference counting.** The current `CFLOBDD_T<T>` uses `ref_ptr` for copy semantics. If your DD type uses value semantics, the integration layer will need to adapt.

**Initialization.** Both `DDVector::Initialize()` and `DDMatrix::Initialize()` are called once at startup. Use them to set up any global state (unique tables, compute tables, etc.).

## What NOT to Implement

- **No QOperation / QuantumTerm logic.** These are in `quantum_operation.hpp` and stay unchanged.
- **No transition system logic.** The TS layer is above the DD backend.
- **No Python bindings.** The pybind11 wrapper is maintained separately.
- **No Gram-Schmidt, conjunction, disjunction, pre/post image.** These are all in `quantum_operation.hpp` and operate on `QOperation` objects, not directly on DD objects.

## Verification

After implementing the API surface, the following should work:

1. `make test && ./test_qreach 8` — C++ benchmark/test executable
2. `cd python_pkg && PYTHONPATH=. python workflow_tests/test_newapi.py` — basic API smoke test
3. `PYTHONPATH=. python workflow_tests/test_grover.py` — Grover workflow with model checking
4. `PYTHONPATH=. python workflow_tests/test_grover_wp.py` — weakest-precondition computation

---

## 6. Conventions the backend MUST match (required reading)

> These were reverse-engineered from the CFLOBDD source on 2026-08-13. They are the part of the contract most likely to cause silently-wrong results if you implement them differently. Items marked **[verify]** are best-effort; confirm with the QReach agent if your implementation disagrees with a reference check.

### 6.1 Level / qubit-count / dimension convention (corrected)

The **"dimension = 2^level"** note in §4 is imprecise. The actual convention is:

| Object | `level` L | number of qubits | dimension (entries) |
|---|---|---|---|
| **Vector** | L | `2^L` | `2^(2^L)` |
| **Matrix** | L | `2^(L-1)` | `2^(2^(L-1))` × `2^(2^(L-1))` |

Equivalently, in terms of `qNum` (the number of qubits, **padded to a power of 2** — e.g. 5 physical qubits → `qNum = 8`):

- Vector over `qNum` qubits → `level = log2(qNum)`, dimension `= 2^qNum`.
- Matrix over `qNum` qubits → `level = log2(qNum) + 1`, size `2^qNum × 2^qNum`.

Worked examples:

| `qNum` (qubits) | vector `level` | vector dim | matrix `level` | matrix size |
|---|---|---|---|---|
| 1 | 0 | 2 | 1 | 2×2 |
| 2 | 1 | 4 | 2 | 4×4 |
| 4 | 2 | 16 | 3 | 16×16 |
| 8 | 3 | 256 | 4 | 256×256 |

The matrix at level L consumes `2^L` Boolean variables: `2^(L-1)` row bits + `2^(L-1)` column bits.

### 6.2 Bit endianness and qubit-index semantics (the observable contract)

**The internal variable order is NOT part of the contract.** LimTDD may choose any
internal index ordering it likes, as long as its own gate construction and its own
`MatrixMultiplyWithVector` are mutually consistent. Do **not** try to mimic CFLOBDD's
interleaved "VOC12" order.

What *is* observable (and must match exactly) is the endianness and the qubit-index
semantics:

- `MkBasisVector(level, index)` / `MkBasisVector(level, s)` interpret `index`/`s`
  **big-endian**: bit/character `0` is the most significant qubit (qubit 0, Qiskit
  convention). `MkBasisVector(level, s)` requires `s.length() == 2^level`.
- `MkCNOT(level, n, ctrl, tgt)` / `MkCCNOT` / `MkCP` / `MkSwap`: `ctrl`/`tgt`/`i`/`j`
  are qubit indices under the same convention (qubit 0 = most significant).

Reference (CFLOBDD uses interleaved row/col internally, but you are free to differ —
this is informational only):

```text
row = (x0, x1, ..., x_{N-1})     # x0 = MSB
col = (y0, y1, ..., y_{N-1})     # y0 = MSB
variables: x0, y0, x1, y1, ..., x_{N-1}, y_{N-1}
```

```text
MkBasisVector(level, "10")  ==  MkBasisVector(level, 2)     # s[0]='1' is the MSB
MkBasisVector(level, "010") ==  MkBasisVector(level, 2)     # (for level where width = 3)
```

`MkBasisVector(level, s)` requires `s.length() == 2^level` (the number of qubits). Qubit `0` is the most significant (leftmost) — this matches Qiskit's `|q0 q1 … qn⟩` ordering.

### 6.3 Scalar type — full requirement (extends §1)

`DDComplex` must support, **in addition** to `real()`, `imag()`, `==`, `!=`, `operator*`, `abs()`, `norm()`, default construction:

- `operator+`, `operator-`, `operator/` (the semantic layer computes `dot(a)/dot(b)` and scalar-division factors),
- comparison against integer literals `0` and `1` (`amp != 0`, `amp != 1`),
- construction from `(double re, double im)` and from a single `double`,
- `operator*` / `operator+` with the DD type: `DDComplex * DD` and `DD + DD` are required (`DD` also needs `operator+`, `operator*` with scalar, and `operator==`).

### 6.4 Gate constructor semantics

`level` argument conventions per §6.1. Angle arguments are **in units of π** (implemented with `cos_pi`/`sin_pi`, not radian `cos`/`sin`).

| Function | Matrix / meaning |
|---|---|
| `MkIdRelation(level)` | identity |
| `MkWalsh(level)` | single-qubit Hadamard `[[1,1],[1,-1]]/√2` — the `1/√2` IS baked in (normalized/unitary) |
| `MkNegation(level)` | Pauli X `[[0,1],[1,0]]` |
| `MkPauliY(level)` | Pauli Y `[[0,-i],[i,0]]` |
| `MkPauliZ(level)` | Pauli Z `[[1,0],[0,-1]]` |
| `MkSGate(level)` | S = `diag(1, i)` (exact `i`, not `cos/sin(π/2)`) |
| `MkPhaseShift(level, θ)` | `diag(1, e^{iπθ})` |
| `MkU3(level, [θ,φ,λ])` | Qiskit U3, `[[cos(πθ/2), -e^{iπλ}sin(πθ/2)], [e^{iπφ}sin(πθ/2), e^{iπ(φ+λ)}cos(πθ/2)]]` |
| `MkArbitrary(level, v[0..7])` | 2×2 `[[a,b],[c,d]]`, row-major, complex entries as interleaved real/imag: `a=(v0,v1)`, `b=(v2,v3)`, `c=(v4,v5)`, `d=(v6,v7)` |
| `MkCNOT(level, n, ctrl, tgt)` | CNOT, control `ctrl`, target `tgt`; `n` = number of qubits = `qNum = 2^(level-1)` |
| `MkCCNOT(level, n, c1, c2, tgt)` | Toffoli, controls `c1`, `c2`, target `tgt` |
| `MkSwap(level, i, j)` | SWAP qubits `i`,`j` |
| `MkiSwap(level, i, j)` | iSWAP |
| `MkCP(level, ctrl, tgt, θ)` | controlled-phase: `diag(1,1,1,e^{iπθ})`, control `ctrl`, target `tgt` |

`concretize()` composes these with `MkSingleQubitGateOnN` / `KroneckerProduct` / `MatrixMultiply` — the backend must not re-order control/target roles (the QReach side already handles `controller < controlled` reordering via SWAP conjugation).

### 6.5 Level relationship for the hot path

`MatrixMultiplyWithVector(gate, vec)` is `gate·|vec⟩` where:
- `gate` is a matrix at level `L+1` (`qNum` qubits),
- `vec` is a vector at level `L` (`qNum` qubits; possibly already converted by `VectorToMatrixInterleaved`),
- result is a vector at level `L`.

`VectorToMatrixInterleaved(vec)` promotes a level-L vector to a level-(L+1) matrix form. If your backend multiplies directly against a plain vector, make this a no-op — but the *output level/type* of a subsequent `MatrixMultiplyWithVector` must still be a level-L vector.

### 6.6 Initialization idempotency

`DDVector::Initialize()` and `DDMatrix::Initialize()` are called once at startup **and again in multiple constructors** (the transition system initializes in 3 places). They must be safe to call repeatedly (idempotent).

### 6.7 Exponential operations must be guarded

`InitializeWithAmplitudes(qnum, …)`, `GetNonZeroAmplitudes`, and `ExtractSingleAmplitude` are inherently `O(2^n)`. They are used **only** for initial-state construction and small-`n` verification — never in the hot path. The adapter must:
- add an `n` upper-bound assertion (guard against the historical 32-bit overflow at `n ≥ 32`),
- build `InitializeWithAmplitudes` bottom-up from the backend's own edge/array constructor (do not insert amplitudes one by one).

### 6.8 `Conjugate` / `Transpose` must be defined on vectors (and phase-aware)

`resetall` decomposes as `NoDistinctionNode → VectorToMatrixInterleaved → Conjugate → Transpose`. If `VectorToMatrixInterleaved` is a no-op in your backend, then `Conjugate` and `Transpose` must still produce the correct reset matrix when applied to a *vector-shaped* object. Implement them for both vectors and matrices, or promote the vector to a matrix first inside the adapter. Verify `resetall` with a dedicated dense cross-check.

**`Conjugate` is not just negating edge weights.** If your backend stores phases in a `P^k = e^{ikπ/n}` map, conjugation must also flip `k → -k`. Conjugating only the edge weights is wrong.

---

## 7. Worked examples (self-check your implementation)

Use these to validate a fresh backend before touching QReach's integration:

```text
qNum = 2  →  vector level 1 (4 amplitudes), matrix level 2 (4×4)

MkBasisVector(1, 0)      = [1, 0, 0, 0]ᵀ      # |00⟩
MkBasisVector(1, 1)      = [0, 1, 0, 0]ᵀ      # |01⟩
MkBasisVector(1, 2)      = [0, 0, 1, 0]ᵀ      # |10⟩
MkBasisVector(1, "10")   = [0, 0, 1, 0]ᵀ      # same as index 2 (big-endian)

MkWalsh(1)               = Hadamard (2×2)
MatrixMultiplyWithVector(MkWalsh(1), MkBasisVector(1, 0)) = (|00⟩+|10⟩)/√2   # H on qubit 0

MkCNOT(2, 2, 0, 1)       = 4×4 CNOT, control qubit 0, target qubit 1
MatrixMultiplyWithVector(MkCNOT(2,2,0,1), MkBasisVector(1, 2)) = MkBasisVector(1, 3)   # |10⟩ → |11⟩
```

The QReach side runs a **dense cross-check against Qiskit's `Statevector`** as the authoritative oracle for these small cases — if your backend disagrees on any of the above, the variable-order or endianness convention is wrong.

---

## 8. Open items and decisions

**Resolved (2026-08-13, after LimTDD feasibility review):**

1. **`MkWalsh` normalization** — RESOLVED: the `1/√2` is baked in (normalized Hadamard). Consistent with QReach's `knownUnitNorm` normalization-skip logic.
2. **`MkSGate`** — RESOLVED: added to the contract (§3, §6.4) as `diag(1, i)` with exact values.
3. **`MkSingleQubitGateOnN` param variants** — RESOLVED: `WithParam` / `WithParamVec` added (§3).
4. **Pointwise `DD * DD`** — RESOLVED: dropped from the contract (§1); the semantic layer never uses it.

**Still to decide / confirm:**

- **Scalar precision (P0) — DECIDED (route a).** LimTDD is double-precision; CFLOBDD is 100-digit. Chosen route: **accept double, parameterize `zeroThreshold()` per backend** (behind `dd_backend.hpp`), and rely on the **Qiskit `Statevector` dense cross-check** as the correctness oracle — exposing precision risk explicitly rather than masking it. Fallback if fixed-point saturation regresses: keep CFLOBDD as the exact baseline (route b). QReach's fixed-point terminates on *subspace dimension* saturation (integer), but that dimension is computed by Gram-Schmidt over thresholded `dot()` products — so double-precision noise can prevent exact saturation; the dense cross-check is the guard.
- **Reference counting vs value semantics** — `CFLOBDD_T<T>` uses `ref_ptr`; the semantic layer copies `DD` freely. If LimTDD uses value semantics with deep copies, the hot path may regress; flag it.
- **Feasibility milestone** — before committing further, run a 4–6 qubit Grover/RUS node-count + memory comparison (LimTDD vs CFLOBDD) to confirm LimTDD's map-compression advantage holds on arbitrary entangled reachable states.