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
//   - operator+ (pointwise addition)
//   - operator* (scalar * DD, left scalar multiplication)
//   - operator* (DD * DD, pointwise multiplication)
//   - operator==
//   - copy construction / assignment (reference-counting semantics preferred)
using DD = /* your DD type */;
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
    //   Given a 1-qubit gate constructor, build the n-qubit matrix
    //   that applies it to `target` and identity to all other qubits.
    //   This replaces CFLOBDD's ApplyGateF recursive Kronecker-product pattern.

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