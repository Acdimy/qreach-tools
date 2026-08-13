# LimTDD Backend — One-Page Implementation Brief

> **For:** the agent implementing the `DDVector` / `DDMatrix` API inside the LimTDD project.
> **You do not need any knowledge of QReach or CFLOBDD to do this job.** Everything you need is here.
> Full contract: `docs/agent-handoffs/backend-replacement-api-contract.md` (§1–§5 signatures, §6–§8 below).

---

## What you are building

Two C++ namespaces that wrap **your** decision-diagram (or tensor) backend so that a *different* tool (QReach) can use it as a drop-in for its existing CFLOBDD backend. You only implement these functions; QReach's quantum-semantics layer calls them.

- `DDVector` — quantum state vectors.
- `DDMatrix` — quantum gate matrices (unitary operators) and matrix algebra.

Both get an `Initialize()` called once at startup (and may be called again — must be idempotent).

---

## Conventions (the part that must be exact)

### Level ↔ qubits ↔ dimension

`qNum` is the **number of qubits, padded to a power of 2** (e.g. 5 qubits → `qNum = 8`). `level` is a `log2` of qubit count, **not** of dimension:

| Object | `level` L | # qubits | dimension |
|---|---|---|---|
| Vector | `log2(qNum)` | `2^L` | `2^(2^L)` entries |
| Matrix | `log2(qNum)+1` | `2^(L-1)` | `2^(2^(L-1))` × `2^(2^(L-1))` |

Examples:

| `qNum` | vector level | vector dim | matrix level | matrix size |
|---|---|---|---|---|
| 1 | 0 | 2 | 1 | 2×2 |
| 2 | 1 | 4 | 2 | 4×4 |
| 4 | 2 | 16 | 3 | 16×16 |
| 8 | 3 | 256 | 4 | 256×256 |

### Variable order & endianness

- **Internal variable order is YOUR choice.** You do not need to mimic CFLOBDD's interleaved order — just keep your own gate construction and your own `MatrixMultiplyWithVector` mutually consistent.
- **Vectors are big-endian (observable).** In `MkBasisVector(level, "…")`, character 0 is the most significant qubit (qubit 0 = leftmost, Qiskit convention). Same for `MkCNOT`/`MkCCNOT`/`MkCP`/`MkSwap` qubit indices.

```text
MkBasisVector(level, "10") == MkBasisVector(level, 2)     # s[0] is the MSB
MkBasisVector(level, s)  requires  s.length() == 2^level
```

### Scalar type `DDComplex`

Must support: `real()`, `imag()`, `+`, `-`, `*`, `/`, `abs()`, `norm()`, `==`, `!=`, comparison against integer `0`/`1`, construction from `(double re, double im)` and from `double`. The `DD` type must support `DD + DD`, `DDComplex * DD`, and `==`.

### Angles are in units of π

Phase angles are implemented with `cos(π·θ)` / `sin(π·θ)`, **not** radian `cos`/`sin`.

---

## Gate constructor semantics

| Function | Meaning |
|---|---|
| `MkIdRelation(level)` | identity |
| `MkWalsh(level)` | Hadamard `[[1,1],[1,-1]]/√2` (the `1/√2` is baked in) |
| `MkNegation(level)` | X `[[0,1],[1,0]]` |
| `MkPauliY(level)` | Y `[[0,-i],[i,0]]` |
| `MkPauliZ(level)` | Z `[[1,0],[0,-1]]` |
| `MkSGate(level)` | S = `diag(1, i)` (exact `i`) |
| `MkSGate(level)` | S = `diag(1, i)` (exact `i`) |
| `MkPhaseShift(level, θ)` | `diag(1, e^{iπθ})` |
| `MkU3(level, [θ,φ,λ])` | Qiskit U3 (angles in π units) |
| `MkArbitrary(level, v[0..7])` | 2×2 `[[a,b],[c,d]]`, `a=(v0,v1) b=(v2,v3) c=(v4,v5) d=(v6,v7)` (re/imag interleaved, row-major) |
| `MkCNOT(level, n, ctrl, tgt)` | CNOT; `n = qNum = 2^(level-1)` |
| `MkCCNOT(level, n, c1, c2, tgt)` | Toffoli |
| `MkSwap(level, i, j)` | SWAP |
| `MkiSwap(level, i, j)` | iSWAP |
| `MkCP(level, ctrl, tgt, θ)` | controlled-phase `diag(1,1,1,e^{iπθ})` |

Plus: `KroneckerProduct`, `MatrixMultiply`, `MatrixMultiplyWithVector` (gate·vector, the hot path), `Conjugate`, `Transpose`, and `MkSingleQubitGateOnN(n, target, gate1q)` **plus its parameterized variants** `MkSingleQubitGateOnNWithParam` / `MkSingleQubitGateOnNWithParamVec` (build an n-qubit gate = `I⊗…⊗G⊗…⊗I`, carrying the gate's angle/param arguments).

`MatrixMultiplyWithVector(gate_level L+1, vector_level L) → vector_level L`. `VectorToMatrixInterleaved` may be a no-op if you multiply directly on vectors.

**Inner product is a contraction, not matrix-multiply.** Add `DDVector::InnerProduct(a, b) = ⟨a|b⟩ = conj(a)·b`, implemented natively as `cont(conj(a), b)`. Do **not** pad vectors to matrix form and use `MatrixMultiply` — that is a CFLOBDD artifact. (This is the one function the current 44-test LimTDD build is still missing; it is the 13th `DDVector` function.)

**Note:** `Conjugate` must flip the phase *map* (`k → -k`), not just edge weights. `Conjugate`/`Transpose` must work on a vector-shaped object too (`resetall` relies on this).

Pointwise `DD * DD` is **not required** — the semantic layer only does `scalar * DD`.

---

## Self-check (validate before integration)

```text
qNum = 2 → vector level 1 (4 amps), matrix level 2 (4×4)

MkBasisVector(1, 0)    = [1,0,0,0]ᵀ    # |00⟩
MkBasisVector(1, "10") = [0,0,1,0]ᵀ    # |10⟩  (big-endian)
MkCNOT(2, 2, 0, 1) applied to [0,0,1,0]ᵀ = [0,0,0,1]ᵀ   # |10⟩ → |11⟩
MkWalsh(1) applied to [1,0,0,0]ᵀ = (|00⟩+|10⟩)/√2        # H on qubit 0
```

If any of these differ, your variable-order or endianness is wrong (the #1 failure mode).

---

## Open items — confirm with the QReach agent

**Resolved:** `MkWalsh` (√2 baked in) · `MkSGate` (added, exact `diag(1,i)`) · `MkSingleQubitGateOnN` param variants (added) · pointwise `DD*DD` (dropped).

**Still open:**

1. **Scalar precision (P0) — DECIDED.** Route chosen: accept double + parameterize `zeroThreshold()` per backend + Qiskit `Statevector` dense cross-check as oracle. Fallback: CFLOBDD stays as exact baseline.
2. Reference counting vs value semantics for `DD` (deep copies on the hot path could regress).
3. Feasibility: run a 4–6 qubit Grover/RUS node-count + memory comparison vs CFLOBDD before committing further.
