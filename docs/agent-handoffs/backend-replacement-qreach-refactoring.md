# QReach Backend Replacement — Document B: Refactoring Plan

**Audience:** AI agent working on the QReach codebase.
**Status:** Planning document — do not implement without user approval.
**Last updated:** 2026-07-08

## Overview

QReach currently uses CFLOBDD as its sole quantum decision-diagram backend. The goal is to enable a **pluggable backend architecture** where alternative DD backends can be substituted with minimal changes to the quantum logic layer.

This document describes the refactoring needed on the QReach side to decouple `quantum_operation.hpp` from CFLOBDD internals and prepare for a second backend.

## Current Architecture

```
Python:  qctl.py  →  parse_qiskit.py  →  pyqreach (pybind11)
                                              │
C++:     transition_system.hpp                 │
              │                                │
         quantum_operation.hpp ←───────────────┘
              │
              ├── QOperation (quantum algebra: Gram-Schmidt, pre/post image, ...)
              ├── QuantumGateTerm (gate concretization → CFLOBDD matrices)
              ├── SingleVecTerm (vector operations → CFLOBDD vectors)
              │
              └── CFLOBDD backend (directly coupled)
                   ├── CFLOBDD_COMPLEX_BIG = CFLOBDD_T<BIG_COMPLEX_FLOAT>
                   ├── Matrix1234ComplexFloatBoost::Mk* (gate constructors)
                   ├── VectorComplexFloatBoost::Mk* (vector constructors)
                   ├── MatrixMultiplyV4WithInfo (hot path)
                   └── Internal access: c.root->level, c.root->rootConnection.returnMapHandle
```

## Coupling Points to Refactor

### Level 1: Type Aliases (easy)

**Current:**
```cpp
// quantum_operation.hpp:21
using namespace CFL_OBDD;
// CFLOBDD_COMPLEX_BIG and BIG_COMPLEX_FLOAT are used everywhere
```

**Plan:** Introduce backend-neutral type aliases in a new header `dd_backend.hpp`:
```cpp
// dd_backend.hpp
namespace qreach {
    using DDComplex = /* backend-specific complex type */;
    using DD = /* backend-specific DD type */;
}
```

`quantum_operation.hpp` would use `DD` and `DDComplex` instead of `CFLOBDD_COMPLEX_BIG` and `BIG_COMPLEX_FLOAT`.

### Level 2: Free Functions Using CFLOBDD Internals (medium)

**Current:** Three free functions in `quantum_operation.hpp` directly use CFLOBDD's Kronecker-product recursion:

| Function | Lines | What it does |
|----------|-------|-------------|
| `ApplyGateF` | 23-45 | Recursively builds n-qubit gate from 1-qubit gate via Kronecker |
| `ApplyGateFWithParam` | 47-67 | Same with a double parameter |
| `ApplyGateFWithParamVec` | 69-89 | Same with a vector parameter |
| `InitializeWithVector` | 91-120 | Builds a vector DD from raw amplitudes |

**Plan:** Move these into the backend namespace:
```cpp
namespace DDMatrix {
    DD MkSingleQubitGateOnN(unsigned int n, unsigned int target,
        DD(*gate1q)(unsigned int));
    DD MkSingleQubitGateOnNWithParam(unsigned int n, unsigned int target,
        DD(*gate1q)(unsigned int, double), double theta);
    DD MkSingleQubitGateOnNWithParamVec(unsigned int n, unsigned int target,
        DD(*gate1q)(unsigned int, std::vector<double>), std::vector<double> v);
}
namespace DDVector {
    DD InitializeWithAmplitudes(unsigned int qnum, std::vector<double> amps);
}
```

### Level 3: Direct CFLOBDD Internal Access (critical)

These are the most problematic coupling points. `quantum_operation.hpp` directly accesses CFLOBDD's internal DAG structure:

| Location | Access Pattern | Purpose |
|----------|---------------|---------|
| `checkifzero()` (L187-203) | `c.root->rootConnection.returnMapHandle` | Check if vector is zero |
| `SingleVecTerm::dot()` (L988-1004) | `tmp.root->rootConnection.returnMapHandle` | Extract overlap amplitude |
| `SingleVecTerm::normalize()` (L1005-1034) | `mulres.root->rootConnection.returnMapHandle` | Extract norm factor |
| `SingleVecTerm::applyGate()` (L1049-1113) | `c.root->level` | Get DD level |
| `SingleVecTerm` constructor (L952-981) | `x.root->level` | Get DD level |
| `QOperation::printFormal()` (L2277-2334) | `VectorPrintColumnHead` | Debug printing |

**Plan:** Replace with backend-neutral query functions:

```cpp
namespace DDVector {
    int GetLevel(DD c);
    bool IsApproximatelyZero(DD c, double threshold = 1e-8);
    DDComplex ExtractSingleAmplitude(DD c);
    std::vector<std::pair<unsigned int, DDComplex>> GetNonZeroAmplitudes(DD c, double threshold);
    DD Normalize(DD c);
    void VectorPrintColumnHead(DD c, std::ostream& out);
}
```

### Level 4: Gate Concretization (large but mechanical)

`QuantumGateTerm::concretize()` (lines 647-924) is a large switch statement that calls `Matrix1234ComplexFloatBoost::Mk*` for each gate type. This is the largest single function to refactor.

**Plan:** The gate constructors are already well-factored — each gate type maps to one `Mk*` call. The refactoring is mechanical: replace `Matrix1234ComplexFloatBoost::MkXxx` with `DDMatrix::MkXxx`.

### Level 5: `QOperation` String Constructors (minor)

`QOperation(std::vector<std::string>)` (lines 1156-1180) and `QOperation(std::vector<double>, unsigned int)` (lines 1181-1206) construct `SingleVecTerm` objects using CFLOBDD-specific constructors.

**Plan:** These already go through `SingleVecTerm` constructors. Once `SingleVecTerm` is decoupled, these follow automatically.

## Proposed Refactoring Strategy

### Phase 1: Introduce Backend Abstraction Header

Create `dd_backend.hpp` with:
- Type aliases: `DD`, `DDComplex`
- Forward declarations of `DDVector` and `DDMatrix` namespaces
- A compile-time backend selector: `#ifdef QREACH_USE_CFLOBDD` / `#elif defined(QREACH_USE_NEW_BACKEND)`

### Phase 2: Extract Internal Access to Query Functions

Add `GetLevel`, `IsApproximatelyZero`, `ExtractSingleAmplitude`, `GetNonZeroAmplitudes`, `Normalize` to the CFLOBDD backend as wrapper functions. Then replace all direct `c.root->...` accesses in `quantum_operation.hpp` with these wrappers.

### Phase 3: Move Free Functions to Backend

Move `ApplyGateF`, `ApplyGateFWithParam`, `ApplyGateFWithParamVec`, `InitializeWithVector` into the CFLOBDD backend namespace. Replace call sites with the new names.

### Phase 4: Gate Concretization Switch

Replace `Matrix1234ComplexFloatBoost::Mk*` calls in `QuantumGateTerm::concretize()` with `DDMatrix::Mk*` calls.

### Phase 5: Integration Test

After each phase, rebuild and run:
```bash
make test && ./test_qreach 8
cd python_pkg && PYTHONPATH=. ../.venv/bin/python workflow_tests/test_newapi.py
```

## Files to Modify

| File | Change |
|------|--------|
| `dd_backend.hpp` | **New** — type aliases, backend selector |
| `quantum_operation.hpp` | Replace CFLOBDD-specific types and internal access |
| `transition_system.hpp` | May need type alias updates |
| `cflobdd/CFLOBDD/matrix1234_complex_float_boost.h` | Add wrapper functions for Phase 2-3 |
| `cflobdd/CFLOBDD/vector_complex_float_boost.h` | Add wrapper functions for Phase 2-3 |
| `python_pkg/qreach_python_wrapper.cpp` | May need type alias updates |
| `Makefile` | Add `USE_CFLOBDD` define |

## Files NOT to Modify

- `qadd.hpp`, `transition_system_qadd.hpp` — shelved
- `qctl.py`, `parse_qiskit.py` — Python layer is backend-agnostic
- `cflobdd/CFLOBDD/matrix1234_node.cpp` — internal CFLOBDD implementation
- `cflobdd/CFLOBDD/cflobdd_node.cpp` — internal CFLOBDD implementation

## Verification Checklist

After full refactoring:

1. `make test && ./test_qreach 8` — C++ benchmark
2. `cd python_pkg && PYTHONPATH=. ../.venv/bin/python workflow_tests/test_newapi.py`
3. `PYTHONPATH=. ../.venv/bin/python workflow_tests/test_grover.py`
4. `PYTHONPATH=. ../.venv/bin/python workflow_tests/test_RUS.py`
5. `PYTHONPATH=. ../.venv/bin/python workflow_tests/test_grover_wp.py`
6. `PYTHONPATH=. ../.venv/bin/python workflow_tests/test_lazy_measurement.py`
7. `PYTHONPATH=. ../.venv/bin/python workflow_tests/test_vqss_correct_lazy.py`
8. `QREACH_GATE_PROFILE=1 PYTHONPATH=. ../.venv/bin/python workflow_tests/test_grover32_lazy_prefix_performance.py` — profiling still works