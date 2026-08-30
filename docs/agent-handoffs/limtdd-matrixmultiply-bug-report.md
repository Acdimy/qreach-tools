# LimTDD `MatrixMultiply` (matrix×matrix) bug — breaks CSX and ctrl>tgt two-qubit gates

**For:** LimTDD backend agent
**From:** QReach side (branch `limtdd-backend`)
**Date:** 2026-08-30
**Severity:** correctness — silently wrong results and spurious `size mismatch` throws
**Found by:** `python_pkg/workflow_tests/test_backend_crosscheck.py` (the Statevector dense cross-check oracle)

## Summary

`DDMatrix::MatrixMultiply` (matrix × matrix) is not yet correct/compactified on
the LimTDD backend. Two-qubit gates that require a **matrix-composition** step
fall apart:

- **CSX** (`√X` controlled) throws `ValueError: MatrixMultiply: size mismatch`.
- **CX / CP / CZ / CCX with `ctrl > tgt`** silently produce the wrong state.

Both symptoms route through the *same* primitive: the QReach semantic layer
(`quantum_operation.hpp` `QuantumGateTerm::concretize`) composes gates with
`MatrixMultiply` whenever it must reorder control/target or build a non-native
gate. The direct path (`ctrl < tgt`) never calls `MatrixMultiply` and works;
the composition path does not.

## Minimal reproductions (QReach + LimTDD, no Qiskit needed to see the failure)

```python
import pyqreach
from qiskit import QuantumCircuit
from qreach.parse_qiskit import parse_qiskit_cir
from qreach.qctl import set_initial_state

# --- 1. CSX (ctrl < tgt, still needs H·C·H composition) ---
qc = QuantumCircuit(2); qc.h(0); qc.csx(0, 1)
ts = pyqreach.TransitionSystem()
parse_qiskit_cir(qc, 2, ts); set_initial_state(ts, "00")
ts.computingFixedPointPost()
# ValueError: MatrixMultiply: size mismatch

# --- 2. CX with ctrl > tgt (needs S·C·S SWAP conjugation) ---
# H(0); CX(0,2); CX(2,1)  ->  Qiskit: (|000>+|111>)/√2
# QReach/LimTDD returns (|000>+|101>)/√2  (the CX(2,1) is silently dropped)
```

## What the semantic layer does (why `MatrixMultiply` is on the hot path)

`quantum_operation.hpp` `concretize()`:

- **CX, `ctrl < tgt`** (line ~579): `res = MkCNOT(ctrl, tgt)` — no `MatrixMultiply`. **Works.**
- **CX, `ctrl > tgt`** (line ~584): `S=MkSwap(ctrl,tgt); C=MkCNOT(ctrl,tgt); res=S·C·S` via two `MatrixMultiply` calls. **Broken** (silent wrong result).
- **CSX** (line ~597): `CSX = H(target)·CP(π/2)(ctrl,tgt)·H(target)` via `MatrixMultiply` even for `ctrl<tgt`. **Throws `size mismatch`.**

So `MatrixMultiply(matrix, matrix)` is the only primitive in play. The CSX case
is the cleanest isolation: it calls `MatrixMultiply(H, C)` where `H =
MkSingleQubitGateOnN(qNum, target, MkWalsh)` and `C = MkCP(level, ctrl, tgt, 0.5)`
are both n-qubit gates of the same dimension, yet it throws `size mismatch`.

## Impact

Under LimTDD, any Qiskit circuit containing:
- `csx`, or
- a two-qubit controlled gate whose control index is **greater** than its target
  index (`cx`, `cp`, `cz`, `ccx` with ctrl>tgt after the parser's qubit mapping),

will either crash or silently produce a wrong reachable set. Standard transpiled
circuits reorder qubits freely, so this is not a corner case.

## Why the existing tests missed it

The passing workflow tests (`test_simulation_grover`, `test_simulation_qft`,
`test_vqss_*`, `test_RUS`) all happen to use controlled gates with `ctrl < tgt`
(or no controlled gates), so they never hit the `MatrixMultiply` composition
path. The new cross-check adds the `ctrl > tgt` and `csx` cases explicitly.

## Where it is tracked

`python_pkg/workflow_tests/test_backend_crosscheck.py` marks these two cases
`XFAIL` under LimTDD (keyed on `not symbolic_available()`). Under CFLOBDD
(`MatrixMultiplyV4` is mature) they are expected to **pass**. When the LimTDD
`MatrixMultiply` is fixed, the XFAIL entries should be removed and the harness
should go fully green on both backends.

## Suggested fix direction (LimTDD side)

`MatrixMultiply` (matrix × matrix) is the remaining "matrix algebra
compactification" item in `limtdd-backend` (see LimTDD `HISTORY.md`). The two
observed modes:

1. `size mismatch` — the dimension/layer check in `MatrixMultiply` rejects two
   same-size operands (likely a `level`/`n` convention leftover after the
   level→n constructor change).
2. silent wrong result — the `cont`-based composition (or `MatrixMultiply`'s
   result renaming) is incorrect for the S·C·S conjugation.

Both are LimTDD-internal; the QReach semantic layer and the CFLOBDD backend are
correct (verified by the cross-check under CFLOBDD conventions).
