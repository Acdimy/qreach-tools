#ifndef DD_BACKEND_HPP
#define DD_BACKEND_HPP

// =============================================================================
// QReach decision-diagram backend interface.
//
// This is the ONLY place a concrete DD backend type is named. The semantic
// layer (`quantum_operation.hpp`) and the transition-system layer
// (`transition_system.hpp`) must use `qreach::DD` / `qreach::DDComplex` /
// `DDVector` / `DDMatrix` and never a CFLOBDD or LimTDD type directly.
//
// Select the backend with the compile-time macro:
//   - undefined or QREACH_USE_CFLOBDD  -> CFLOBDD (default, current)
//   - QREACH_USE_LIMTDD                -> LimTDD (external project at LIMTDD_PATH)
//
// Conventions (level/dimension, endianness, gate semantics) are specified in
// docs/agent-handoffs/backend-replacement-api-contract.md §6–§8.
// =============================================================================

#include <complex>
#include <cstddef>
#include <ostream>
#include <string>
#include <utility>
#include <vector>

#ifdef QREACH_USE_LIMTDD

  // LimTDD provides the full `DDVector` / `DDMatrix` implementation inline, in
  // global namespaces (not `limtdd::`), and the DD/DDComplex types in `limtdd::`.
  // No forward declarations or separate impl header are needed here.
  #include "dd/backend/DDVector.hpp"
  #include "dd/backend/DDMatrix.hpp"
  namespace qreach {
    using DD        = limtdd::DD;
    using DDComplex = limtdd::DDComplex;
  }

#else  // QREACH_USE_CFLOBDD (default)

  #include "cflobdd/CFLOBDD/matrix1234_complex_float_boost.h"
  #include "cflobdd/CFLOBDD/vector_complex_float_boost.h"
  namespace qreach {
    using DD        = CFL_OBDD::CFLOBDD_COMPLEX_BIG;
    using DDComplex = CFL_OBDD::BIG_COMPLEX_FLOAT;
  }

  // ---------------------------------------------------------------------------
  // Vector operations (interface declarations; defined in dd_backend_cflobdd.h)
  // ---------------------------------------------------------------------------
  namespace DDVector {

    void Initialize();

    // Construction
    qreach::DD MkBasisVector(unsigned int level, unsigned int index);
    qreach::DD MkBasisVector(unsigned int level, std::string bitstring);
    qreach::DD NoDistinctionNode(unsigned int level, qreach::DDComplex val);
    qreach::DD InitializeWithAmplitudes(unsigned int qnum, std::vector<double> amps);

    // Format conversion
    qreach::DD VectorToMatrixInterleaved(qreach::DD vec);

    // Queries
    int GetLevel(qreach::DD c);
    bool IsApproximatelyZero(qreach::DD c, double threshold = 1e-8);
    qreach::DDComplex ExtractSingleAmplitude(qreach::DD c);
    std::vector<std::pair<unsigned int, qreach::DDComplex>> GetNonZeroAmplitudes(
        qreach::DD c, double threshold = 1e-8);

    // Normalization
    qreach::DD Normalize(qreach::DD c);

    // Inner product <a|b> = conj(a) · b. Backend-native: CFLOBDD uses
    // Transpose→Conjugate→MatrixMultiply on matrix-form vectors; LimTDD uses
    // tensor contraction cont(conj(a), b) directly (no vector padding).
    qreach::DDComplex InnerProduct(qreach::DD a, qreach::DD b);

    // Debug
    void VectorPrintColumnHead(qreach::DD c, std::ostream& out);

  } // namespace DDVector

  // ---------------------------------------------------------------------------
  // Matrix operations
  // ---------------------------------------------------------------------------
  namespace DDMatrix {

    void Initialize();

    // Single-qubit gate matrices
    qreach::DD MkIdRelation(unsigned int level);
    qreach::DD MkWalsh(unsigned int level);
    qreach::DD MkNegation(unsigned int level);
    qreach::DD MkPauliY(unsigned int level);
    qreach::DD MkPauliZ(unsigned int level);
    qreach::DD MkSGate(unsigned int level);
    qreach::DD MkPhaseShift(unsigned int level, double theta);
    qreach::DD MkU3(unsigned int level, std::vector<double> params);
    qreach::DD MkArbitrary(unsigned int level, std::vector<double> params);

    // Multi-qubit gate matrices
    qreach::DD MkCNOT(unsigned int level, unsigned int n, long ctrl, long tgt);
    qreach::DD MkCCNOT(unsigned int level, unsigned int n, long c1, long c2, long tgt);
    qreach::DD MkSwap(unsigned int level, long i, long j);
    qreach::DD MkiSwap(unsigned int level, long i, long j);
    qreach::DD MkCP(unsigned int level, long ctrl, long tgt, double theta);

    // Single-qubit gate embedded in an n-qubit system (I ⊗ ... ⊗ G ⊗ ... ⊗ I)
    qreach::DD MkSingleQubitGateOnN(unsigned int n, unsigned int target,
                                    qreach::DD(*gate1q)(unsigned int));
    qreach::DD MkSingleQubitGateOnNWithParam(unsigned int n, unsigned int target,
                                             qreach::DD(*gate1q)(unsigned int, double),
                                             double theta);
    qreach::DD MkSingleQubitGateOnNWithParamVec(
        unsigned int n, unsigned int target,
        qreach::DD(*gate1q)(unsigned int, std::vector<double>),
        std::vector<double> v);

    // Matrix operations
    qreach::DD KroneckerProduct(qreach::DD a, qreach::DD b);
    qreach::DD MatrixMultiply(qreach::DD a, qreach::DD b);
    qreach::DD MatrixMultiplyWithVector(qreach::DD gate, qreach::DD vec);  // hot path
    qreach::DD Conjugate(qreach::DD c);
    qreach::DD Transpose(qreach::DD c);

  } // namespace DDMatrix

  // CFLOBDD passthrough implementation of the declarations above.
  #include "cflobdd/CFLOBDD/dd_backend_cflobdd.h"

#endif // QREACH_USE_LIMTDD / QREACH_USE_CFLOBDD

#endif // DD_BACKEND_HPP
