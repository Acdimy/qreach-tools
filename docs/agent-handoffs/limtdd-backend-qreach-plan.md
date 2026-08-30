# QReach × LimTDD Backend Replacement — QReach-Side Plan

> **Audience:** AI agent working in *this* repository (`qreach-tools`), branch `limtdd-backend`.
> **Status:** Phase 0–5 **executed** (see `limtdd-backend-integration-log.md` for results, open items, and remaining LimTDD-side TODOs). Do not merge without user approval.
> **Last updated:** 2026-08-13 (plan); execution results in `limtdd-backend-integration-log.md`.
> **Verified against:** `limtdd-backend` @ `0570b80` (clean tree, in sync with `origin/limtdd-backend`).

---

## 0. Scope and division of labor

The backend replacement is split across two projects / two agents:

| Side | Repository | Responsibilities |
|---|---|---|
| **QReach (this doc)** | `qreach-tools` | Introduce a backend-neutral interface; route `quantum_operation.hpp` + `transition_system.hpp` through it; keep CFLOBDD compiling behind the same interface as the default; make LimTDD a compile-time-switchable drop-in. |
| **LimTDD (other agent)** | LimTDD project | Implement the `DDVector` / `DDMatrix` API surface specified in `docs/agent-handoffs/backend-replacement-api-contract.md`. |

The API contract is **Document A** (`backend-replacement-api-contract.md`). The refactoring sketch is **Document B** (`backend-replacement-qreach-refactoring.md`). This plan supersedes Document B's phase outline with a concrete, line-verified, buildable sequence and flags the places where the contract is currently underspecified.

**Hard boundary:** `quantum_operation.hpp` is the *only* file in the semantic layer that may call the DD backend. `transition_system.hpp`, `qreach_python_wrapper.cpp`, `parse_qiskit.py`, and `qctl.py` must never see a backend type. The Python layer is backend-agnostic and must not change.

---

## 1. Current coupling inventory (verified, with line numbers)

All line numbers refer to `quantum_operation.hpp` unless stated otherwise.

### 1.1 Type-level coupling

```cpp
// quantum_operation.hpp:4-5
#include "cflobdd/CFLOBDD/matrix1234_complex_float_boost.h"
#include "cflobdd/CFLOBDD/vector_complex_float_boost.h"
// quantum_operation.hpp:21
using namespace CFL_OBDD;
```

The two types used everywhere are:

| Symbol | Meaning | ~occurrences |
|---|---|---|
| `CFLOBDD_COMPLEX_BIG` | the DD (vector/matrix) value type | ~95 |
| `BIG_COMPLEX_FLOAT` | the complex scalar (`cpp_complex_100`, 100-digit) | ~12 |

These become `qreach::DD` and `qreach::DDComplex`.

### 1.2 Free functions that use CFLOBDD internals (move into backend namespace)

| Function | Lines | Backend symbols used |
|---|---|---|
| `ApplyGateF` | 23–45 | `MkIdRelationInterleaved`, `KroneckerProduct2Vocs` |
| `ApplyGateFWithParam` | 47–67 | same |
| `ApplyGateFWithParamVec` | 69–89 | same |
| `InitializeWithVector` | 91–125 | `NoDistinctionNode`, `MkBasisVector`, `VectorToMatrixInterleaved`, `BIG_COMPLEX_FLOAT * basisVec`, `res + scaledVec` |

Map to `DDMatrix::MkSingleQubitGateOnN`, `DDMatrix::MkSingleQubitGateOnNWithParam(Vec)`, `DDVector::InitializeWithAmplitudes` (see §3).

### 1.3 Direct CFLOBDD DAG access (must be replaced by query functions)

| Site | Line(s) | Access | Replacement |
|---|---|---|---|
| `checkifzero` | 189 | `c.root->rootConnection.returnMapHandle` | `DDVector::IsApproximatelyZero` / `GetNonZeroAmplitudes` |
| `concretizeInline` | 928 | `content.root->level != 1` | `DDVector::GetLevel(content)` |
| `SingleVecTerm(DD x)` | 980 | `x.root->level` | `DDVector::GetLevel(x)` |
| `dot` | 998–1006 | `tmp.root->rootConnection.returnMapHandle`, `tmp.root->level`, `tmp.root->EvaluateIteratively`, `SH_OBDD::Assignment` | `DDVector::ExtractSingleAmplitude` |
| `normalize` | 1030–1082 | `content.root->level` (×4), `mulres.root->rootConnection.returnMapHandle`, `content.root->EvaluateIteratively`, `SH_OBDD::Assignment` | `DDVector::GetLevel`, `ExtractSingleAmplitude`, `Normalize` |

**Important:** the `dot()` / `normalize()` "FALLBACK" blocks (lines 995–1007, 1041–1068) are *workarounds for a CFLOBDD-specific `MatrixTranspose` corruption bug* (see `cflobdd-level8-transpose-bug.md`, `cflobdd-transpose-dag-corruption*.md`). They must **not** be copied into a LimTDD build path. Once `dot`/`normalize` are expressed as `ExtractSingleAmplitude` + `Normalize`, the fallbacks disappear because the backend implements these primitives correctly.

### 1.4 Backend namespace call sites (complete)

**`Matrix1234ComplexFloatBoost::`** — 19 distinct methods, 95 total call sites:

| Current | New (`DDMatrix`) | Notes |
|---|---|---|
| `MkIdRelationInterleaved` | `MkIdRelation` | identity |
| `MkWalshInterleaved` | `MkWalsh` | H |
| `MkNegationMatrixInterleaved` | `MkNegation` | X |
| `MkPauliYMatrixInterleaved` | `MkPauliY` | Y |
| `MkPauliZMatrixInterleaved` | `MkPauliZ` | Z |
| `MkSGateInterleaved` | **gap** → see §4.1 | S = P(π/2) |
| `MkPhaseShiftGateInterleaved` | `MkPhaseShift` | P(θ); `θ` in units of π |
| `MkU3GateInterleaved` | `MkU3` | U3(θ,φ,λ) |
| `MkArbitraryGateInterleaved` | `MkArbitrary` | 8-param; also meas0/meas1/reset0 |
| `MkCNOT` | `MkCNOT` | |
| `MkCPGate` | `MkCP` | also CZ, CSX |
| `MkSwapGate` | `MkSwap` | |
| `MkiSwapGate` | `MkiSwap` | |
| `MkCCNOT` | `MkCCNOT` | |
| `KroneckerProduct2Vocs` | `KroneckerProduct` | ⊗ (interleaved) |
| `MatrixMultiplyV4` | `MatrixMultiply` | matrix·matrix |
| `MatrixMultiplyV4WithInfo` | `MatrixMultiplyWithVector` | **hot path** |
| `MatrixConjugate` | `Conjugate` | |
| `MatrixTranspose` | `Transpose` | |

**`VectorComplexFloatBoost::`** — 4 distinct methods, 14 total call sites:

| Current | New (`DDVector`) | Notes |
|---|---|---|
| `MkBasisVector(level, index)` | `MkBasisVector(level, index)` | |
| `MkBasisVector(level, string)` | `MkBasisVector(level, string)` | |
| `NoDistinctionNode(level, val)` | `NoDistinctionNode(level, val)` | constant vector |
| `VectorToMatrixInterleaved` | `VectorToMatrixInterleaved` | may be no-op in LimTDD |
| `VectorPrintColumnHead` | `VectorPrintColumnHead` | debug |

### 1.5 Global initialization

`transition_system.hpp:201-203` (and two sibling sites 245–247, 255–257):

```cpp
CFLOBDDNodeHandle::InitNoDistinctionTable();
CFLOBDDNodeHandle::InitAdditionInterleavedTable();
CFLOBDDNodeHandle::InitReduceCache();
```

→ `DDVector::Initialize(); DDMatrix::Initialize();`

### 1.6 Build system

- `Makefile` compiles every `cflobdd/CFLOBDD/*.cpp` into `libqreach.so`. `HEADERS` has a stale `-I.cflobdd/...` (leading dot) that is unused — includes resolve via `-I.` + full relative paths.
- `python_pkg/tasks.py::compile_python_module` compiles `qreach_python_wrapper.cpp` as a single TU against `-I../ -lqreach`. The wrapper itself contains **no** backend symbol usage (only `#include "quantum_operation.hpp"` / `transition_system.hpp`).

---

## 2. Target architecture

Introduce one header that owns the backend selection:

```cpp
// dd_backend.hpp  (new, at repo root)
#pragma once

#ifdef QREACH_USE_LIMTDD
  #include "limtdd/limtdd_adapter.hpp"     // LimTDD-side adapter (other agent)
  namespace qreach {
    using DD        = limtdd::...;         // per adapter
    using DDComplex = limtdd::...;
  }
#else
  #include "cflobdd/CFLOBDD/matrix1234_complex_float_boost.h"
  #include "cflobdd/CFLOBDD/vector_complex_float_boost.h"
  namespace qreach {
    using DD        = CFL_OBDD::CFLOBDD_COMPLEX_BIG;
    using DDComplex = CFL_OBDD::BIG_COMPLEX_FLOAT;
  }
#endif

namespace DDVector { /* forward decls of the 10 functions */ }
namespace DDMatrix { /* forward decls of the 19 functions */ }
```

The `DDVector`/`DDMatrix` namespaces are then **implemented per backend**:

- CFLOBDD impl: `cflobdd/CFLOBDD/dd_backend_cflobdd.h` — thin inline wrappers over `Matrix1234ComplexFloatBoost` / `VectorComplexFloatBoost`, plus the new query functions (`GetLevel`, `ExtractSingleAmplitude`, `Normalize`, `IsApproximatelyZero`, `GetNonZeroAmplitudes`) and the free functions moved out of `quantum_operation.hpp` (`MkSingleQubitGateOnN*`, `InitializeWithAmplitudes`).
- LimTDD impl: provided by the other agent, mapping the same namespaces onto LimTDD.

Both files are header-only (matching the existing "everything is header + inline" style of `quantum_operation.hpp`). Only one is compiled at a time via the `QREACH_USE_LIMTDD` macro.

The key invariant the refactor must hold: **`quantum_operation.hpp` and `transition_system.hpp` contain zero `CFL_OBDD`/`CFLOBDD`/`BIG_COMPLEX`/`.root->` references after the refactor.** Grep this as the completion check (see §6).

---

## 3. Complete API mapping

This is the single source of truth for the QReach side. The LimTDD agent implements the right-hand columns; the QReach agent implements the left-hand columns behind the interface.

### 3.1 Core types (already in contract §1)

```cpp
namespace qreach {
  using DDComplex = /* complex scalar: real(), imag(), ==, !=, *, abs(), norm() */;
  using DD        = /* DD: +, scalar*, pointwise*, ==, copy/assign (refcount or value) */;
}
```

### 3.2 `DDVector` (contract §2 + additions for current code)

| Function | From (current) | Contract status |
|---|---|---|
| `void Initialize()` | `CFLOBDDNodeHandle::InitNoDistinctionTable()` etc. | contract |
| `DD MkBasisVector(level, index)` | `VectorComplexFloatBoost::MkBasisVector` | contract |
| `DD MkBasisVector(level, bitstring)` | `VectorComplexFloatBoost::MkBasisVector` | contract |
| `DD NoDistinctionNode(level, val)` | `VectorComplexFloatBoost::NoDistinctionNode` | contract |
| `DD InitializeWithAmplitudes(qnum, amps)` | free `InitializeWithVector` | contract |
| `DD VectorToMatrixInterleaved(vec)` | `VectorComplexFloatBoost::VectorToMatrixInterleaved` | contract (may be no-op) |
| `int GetLevel(DD)` | `c.root->level` | contract |
| `bool IsApproximatelyZero(DD, threshold)` | `checkifzero` body | contract |
| `DDComplex ExtractSingleAmplitude(DD)` | `dot`/`normalize` `[0,0]` extraction | contract |
| `GetNonZeroAmplitudes(DD, threshold)` | `returnMapHandle` iteration | contract |
| `DD Normalize(DD)` | `normalize()` | contract |
| `void VectorPrintColumnHead(DD, ostream&)` | `VectorComplexFloatBoost::VectorPrintColumnHead` | contract |

### 3.3 `DDMatrix` (contract §3)

| Function | From (current) |
|---|---|
| `void Initialize()` | `CFLOBDDNodeHandle::InitAdditionInterleavedTable()` + `InitReduceCache()` |
| `DD MkIdRelation(level)` | `MkIdRelationInterleaved` |
| `DD MkWalsh(level)` | `MkWalshInterleaved` |
| `DD MkNegation(level)` | `MkNegationMatrixInterleaved` |
| `DD MkPauliY(level)` | `MkPauliYMatrixInterleaved` |
| `DD MkPauliZ(level)` | `MkPauliZMatrixInterleaved` |
| `DD MkPhaseShift(level, theta)` | `MkPhaseShiftGateInterleaved` |
| `DD MkU3(level, params)` | `MkU3GateInterleaved` |
| `DD MkArbitrary(level, params)` | `MkArbitraryGateInterleaved` |
| `DD MkCNOT(level, n, ctrl, tgt)` | `MkCNOT` |
| `DD MkCCNOT(level, n, c1, c2, tgt)` | `MkCCNOT` |
| `DD MkSwap(level, i, j)` | `MkSwapGate` |
| `DD MkiSwap(level, i, j)` | `MkiSwapGate` |
| `DD MkCP(level, ctrl, tgt, theta)` | `MkCPGate` |
| `DD MkSingleQubitGateOnN(n, target, gate1q)` | free `ApplyGateF` |
| `DD KroneckerProduct(a, b)` | `KroneckerProduct2Vocs` |
| `DD MatrixMultiply(a, b)` | `MatrixMultiplyV4` |
| `DD MatrixMultiplyWithVector(gate, vec)` | `MatrixMultiplyV4WithInfo` |
| `DD Conjugate(c)` | `MatrixConjugate` |
| `DD Transpose(c)` | `MatrixTranspose` |

---

## 4. Contract gaps / open questions for the LimTDD agent

Resolve these **before or during** Phase 4, not after. Each is a correctness or semantic-compatibility risk.

### 4.1 `MkSGate` is missing from the contract (P0)

Current code has a dedicated `MkSGateInterleaved` (used by `s`, and by `sx = H·S·H`). The contract only lists `MkPhaseShift`. S = `PhaseShift(π/2)`, i.e. `MkPhaseShift(level, 0.5)` under the existing `θ`-in-units-of-π convention (Sdg = −0.5, T = 0.25).

**Decision needed:** either (a) QReach maps `s` → `MkPhaseShift(level, 0.5)`, or (b) LimTDD adds `MkSGate`. **Risk of (a):** the CFLOBDD `MkSGate` may use exact complex values (0, 1, ±i) rather than `cos/sin(π/2)`, which affects amplitude-merging/equality. **Recommendation:** add `MkSGate` to the contract so the LimTDD agent can choose exact values; QReach keeps the existing gate decomposition unchanged.

### 4.2 Variable-ordering / vocabulary convention (P0 — the #1 correctness risk)

CFLOBDD matrices use an **interleaved** row/col variable order (`VOC12`): `x0,y0,x1,y1,…`. `VectorToMatrixInterleaved`, `KroneckerProduct2Vocs`, and the gate constructors all assume this. `MkSingleQubitGateOnN` must produce a matrix whose variable order matches whatever `MatrixMultiplyWithVector` expects.

The contract must pin down, unambiguously:
1. Matrix variable order (interleaved row-major? `row0,row1` then col?).
2. Vector variable order (single vocabulary, index = basis integer in little-endian bit order).
3. The `level` ↔ `dimension = 2^level` convention, and that a matrix over an n-qubit system is `level = ceil(log2(n)) + 1` with `2^level = 2n` variables.

If LimTDD and CFLOBDD disagree here, `MatrixMultiplyWithVector` will produce silently-wrong results. Add a small dense cross-check (Phase 6, item 2) to catch this.

### 4.3 Scalar precision mismatch (P0 — DECIDED: route a)

CFLOBDD uses `BIG_COMPLEX_FLOAT` = `cpp_complex_100` (100-digit). QReach's `zeroThreshold()` / `checkifzero` / `isZero` / `singletonOrthogonalTo` thresholds were **tuned against 100-digit arithmetic** (see `grover32-timeout-localization-notes.md` Phase 4). A `std::complex<double>` LimTDD backend has a different noise floor.

**Decision (2026-08-13, after LimTDD feasibility review):** accept double precision; **parameterize `zeroThreshold()` per backend** behind `dd_backend.hpp`; rely on the Qiskit `Statevector` dense cross-check (Phase 6) as the oracle. Fallback if fixed-point saturation regresses: route b (CFLOBDD as exact baseline). Fixed-point terminates on *subspace dimension* saturation (integer), but that dimension comes from Gram-Schmidt over thresholded `dot()` — so precision genuinely matters for convergence; the dense cross-check is the guard, not a proof.

### 4.4 Reference counting vs value semantics (P1)

`CFLOBDD_T<T>` uses `ref_ptr` (shared, copy-on-write-agnostic). The contract allows value semantics. The QReach side currently copies `DD` objects freely (`auto c1 = ...; c1 = ...; return c1;`). If LimTDD uses value semantics with deep copies, the hot path (`MatrixMultiplyWithVector`) could regress. This is the LimTDD agent's responsibility, but QReach should avoid introducing any new `DD` copies during refactor (currently none are introduced — this is just a caution).

### 4.5 `QOperation(["+"*n])` / `simpleProductStateAmplitudes` (P2)

The string constructor path (`MkBasisVector(level, string)`) and the `InitializeWithAmplitudes` path must both be provided by LimTDD. Note the 32-bit overflow fix (`commit ecaec8d`) changed `simpleProductStateAmplitudes` to avoid `1 << qNum` overflow — ensure the LimTDD `MkBasisVector(level, string)` and `InitializeWithAmplitudes` do not reintroduce the same overflow for `n ≥ 32`.

### 4.6 `resetall` decomposition uses `Conjugate`+`Transpose` of a vector (P2)

`resetall` (line 760–764) builds `NoDistinctionNode(level,1)` → `VectorToMatrixInterleaved` → `Conjugate` → `Transpose`. Under LimTDD, `VectorToMatrixInterleaved` may be a no-op; verify this composition still yields the correct reset matrix. This is also where the CFLOBDD transpose bug bit before — treat `Transpose` as a correctness-critical primitive in LimTDD.

---

## 5. Phased implementation plan

Each phase leaves the tree building and all CFLOBDD tests green. Do **not** batch phases — the CFLOBDD rebuild is slow, so a broken phase is expensive to bisect.

**Build/test gate after every phase:**
```bash
make test && ./test_qreach 8
cd python_pkg && ../.venv/bin/python -m invoke build-pybind11
../.venv/bin/python workflow_tests/test_newapi.py
../.venv/bin/python workflow_tests/test_grover.py
```

### Phase 0 — Interface header + CFLOBDD passthrough (no behavior change)

1. Create `dd_backend.hpp` with the `qreach::DD`/`qreach::DDComplex` aliases (CFLOBDD branch only for now) and `DDVector`/`DDMatrix` forward declarations.
2. Create `cflobdd/CFLOBDD/dd_backend_cflobdd.h` implementing `DDVector`/`DDMatrix` as inline wrappers that **forward** to the existing `Matrix1234ComplexFloatBoost::` / `VectorComplexFloatBoost::` methods (renamed per §3.3), plus the new query functions:
   - `GetLevel(c) = c.root->level`
   - `ExtractSingleAmplitude(c)` = the `[0,0]` extraction logic currently in `dot()`
   - `Normalize(c)` = the current `normalize()` body (moved verbatim)
   - `IsApproximatelyZero(c, thr)` = the `checkifzero` body
   - `GetNonZeroAmplitudes(c, thr)` = iterate `returnMapHandle`
   - `MkSingleQubitGateOnN*` = the `ApplyGateF*` bodies
   - `InitializeWithAmplitudes` = the `InitializeWithVector` body
3. Add `-DQREACH_USE_CFLOBDD` to the Makefile compile flags for documentation (default branch anyway).
4. Build + run the gate. **No semantic change expected** — this phase only relocates code.

### Phase 1 — Swap type aliases

1. In `quantum_operation.hpp`, replace the two `#include`s and `using namespace CFL_OBDD;` with `#include "dd_backend.hpp"` + `using namespace qreach;`.
2. Mechanical replace: `CFLOBDD_COMPLEX_BIG` → `DD`, `BIG_COMPLEX_FLOAT` → `DDComplex`.
3. `transition_system.hpp`: same include swap; `using namespace CFL_OBDD;` → `using namespace qreach;`.
4. Build + gate.

### Phase 2 — Replace direct DAG access with query functions

Replace each site in §1.3 with the `DDVector::` query:
- `checkifzero` → `DDVector::IsApproximatelyZero`.
- `concretizeInline`, `SingleVecTerm(DD)`, `normalize` → `DDVector::GetLevel`.
- `dot` → `DDVector::ExtractSingleAmplitude` (drops the `returnMapHandle`/`EvaluateIteratively`/`SH_OBDD::Assignment` code and the transpose-corruption fallback).
- `normalize` → `DDVector::Normalize`.

Build + gate. This is the point where the CFLOBDD transpose-corruption fallbacks are removed from the *semantic layer* and become the CFLOBDD backend's responsibility inside `DDVector::Normalize`/`ExtractSingleAmplitude` (they can live in `dd_backend_cflobdd.h`, keeping CFLOBDD correct while LimTDD implements them cleanly).

### Phase 3 — Move free functions into the backend namespace

1. Delete `ApplyGateF` / `ApplyGateFWithParam` / `ApplyGateFWithParamVec` / `InitializeWithVector` from `quantum_operation.hpp`.
2. In `concretize()`, replace their call sites with `DDMatrix::MkSingleQubitGateOnN(..., DDMatrix::MkWalsh)` etc., and `InitializeWithVector` with `DDVector::InitializeWithAmplitudes`.
3. Build + gate.

### Phase 4 — Gate concretization switch

Replace every `Matrix1234ComplexFloatBoost::Mk*` in `concretize()` (lines 655–917) with `DDMatrix::Mk*` per §3.3, and every `VectorComplexFloatBoost::` call with `DDVector::` per §3.2. Resolve §4.1 (`MkSGate`) first.

Build + gate. After this phase, `grep -E "Matrix1234ComplexFloatBoost|VectorComplexFloatBoost|CFLOBDD|BIG_COMPLEX|CFL_OBDD" quantum_operation.hpp transition_system.hpp` must return nothing.

### Phase 5 — Global initialization

Replace the init sequence in `transition_system.hpp` (and any in `qreach_python_wrapper.cpp`) with `DDMatrix::Initialize(); DDVector::Initialize();`.

The full CFLOBDD init sequence (verified in Phase 0) is **seven** calls, not three:

```cpp
CFLOBDDNodeHandle::InitNoDistinctionTable();
CFLOBDDNodeHandle::InitAdditionInterleavedTable();
CFLOBDDNodeHandle::InitReduceCache();
InitPairProductCache();          // needed by KroneckerProduct / pointwise ops
InitTripleProductCache();        // needed by matrix multiply
Matrix1234ComplexFloatBoost::Matrix1234Initializer();
VectorComplexFloatBoost::VectorInitializer();
```

Mapping in `dd_backend_cflobdd.h`: `DDMatrix::Initialize()` does the first six; `DDVector::Initialize()` does `VectorInitializer()`. Call `DDMatrix::Initialize()` first. **The product caches are easy to miss** — omitting them causes a segfault in `Hashtable::Fetch` the first time `KroneckerProduct`/`MatrixMultiply` runs (observed in the Phase 0 smoke test).

### Phase 6 — Full verification + LimTDD integration readiness

1. Full C++ + Python regression suite (§6).
2. Add a small **dense cross-check** test that, for a handful of small circuits (H, CX, CCX, U3, SWAP on 2–5 qubits), compares `pyqreach` post-image results against Qiskit's `Statevector` — this is the oracle that catches variable-ordering mismatches (§4.2) when LimTDD is dropped in.
3. Verify `make clean && make all` still links with no CFLOBDD-only symbols leaking into the semantic layer.
4. Hand the `DDVector`/`DDMatrix` namespaces + `dd_backend.hpp` + this plan to the LimTDD agent as the integration surface.

---

## 6. Verification checklist

**Static (must be clean after Phase 4/5):**
```bash
grep -nE "Matrix1234ComplexFloatBoost|VectorComplexFloatBoost|CFLOBDD|BIG_COMPLEX|CFL_OBDD|CFLOBDDNodeHandle|\.root->|SH_OBDD" \
  quantum_operation.hpp transition_system.hpp cl_proposition.hpp \
  python_pkg/qreach_python_wrapper.cpp
# expected: no matches (except none)
```

**Build:**
```bash
make clean && make test && ./test_qreach 8
cd python_pkg && ../.venv/bin/python -m invoke build-qreach && ../.venv/bin/python -m invoke build-pybind11
```

**Python regression suite (CFLOBDD must stay green through every phase):**
```bash
cd python_pkg
../.venv/bin/python workflow_tests/test_newapi.py
../.venv/bin/python workflow_tests/test_grover.py
../.venv/bin/python workflow_tests/test_RUS.py
../.venv/bin/python workflow_tests/test_lazy_measurement.py
../.venv/bin/python workflow_tests/test_bv_n14.py        # (long, ~60s)
../.venv/bin/python test_symts_minimal.py                 # symbolic regression (shelved, still must pass)
```

**New dense cross-check (Phase 6):** a `pytest` comparing a few small circuits against `qiskit.quantum_info.Statevector`.

---

## 7. Risks

| Risk | Severity | Mitigation |
|---|---|---|
| Variable-order mismatch between CFLOBDD and LimTDD | **Critical** | Pin the convention in the contract (§4.2); add dense cross-check oracle (Phase 6). |
| Scalar precision drift (100-digit → double) | High | Keep QReach thresholds; parameterize `zeroThreshold` per backend if needed (§4.3). |
| Transpose corruption fallbacks leak into LimTDD path | High | Move fallbacks into `dd_backend_cflobdd.h` only; LimTDD implements `Transpose`/`ExtractSingleAmplitude` cleanly (§2, Phase 2). |
| `MkSGate` missing from contract | Medium | Add to contract; keep decomposition unchanged (§4.1). |
| Slow CFLOBDD rebuild obscures phase bisection | Medium | Strict one-phase-at-a-time discipline (§5). |
| Value-semantics DD copies regress hot path | Medium | No new copies introduced on QReach side; LimTDD agent owns (§4.4). |
| Shelved SymTS/QADD regressions | Low | `test_symts_minimal.py` kept green; do not touch `transition_system_qadd.hpp`/`qadd.hpp`. |

---

## 8. Related documents

- `docs/agent-handoffs/backend-replacement-api-contract.md` — Document A (API contract, LimTDD side).
- `docs/agent-handoffs/backend-replacement-qreach-refactoring.md` — Document B (this plan supersedes its phase outline).
- `docs/agent-handoffs/cflobdd-level8-transpose-bug.md`, `cflobdd-transpose-dag-corruption*.md` — why the `dot`/`normalize` fallbacks exist and must not leak.
- `docs/agent-handoffs/grover32-timeout-localization-notes.md` — numerical-threshold history (Phase 3/4).
- `CFLOBDD_THEORY_AND_IMPLEMENTATION_GUIDE.md` — variable-order / vocabulary semantics.
