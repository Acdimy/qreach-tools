#ifndef DD_BACKEND_CFLOBDD_H
#define DD_BACKEND_CFLOBDD_H

// =============================================================================
// CFLOBDD implementation of the DDVector / DDMatrix interface.
//
// Included from dd_backend.hpp (which already declares the interface and
// defines qreach::DD / qreach::DDComplex). This header is header-only and
// forwards to the existing CFLOBDD backend, plus the query functions and the
// free functions that were moved out of quantum_operation.hpp.
//
// NOTE: the dot()/normalize() "FALLBACK" blocks below are workarounds for a
// CFLOBDD-specific MatrixTranspose corruption bug at level >= 8. They belong
// HERE (the CFLOBDD backend), NOT in the semantic layer, so that a LimTDD
// build can implement Normalize()/ExtractSingleAmplitude() cleanly.
// =============================================================================

#include <cassert>
#include <cmath>
#include <iostream>

#include "cflobdd/CFLOBDD/assignment.h"   // SH_OBDD::Assignment
#include "cflobdd/CFLOBDD/cross_product.h" // InitPairProductCache / InitTripleProductCache

// Same environment the moved code was written under.
using namespace CFL_OBDD;

namespace DDVector {

    inline void Initialize() {
        VectorComplexFloatBoost::VectorInitializer();
    }

    inline qreach::DD MkBasisVector(unsigned int level, unsigned int index) {
        return VectorComplexFloatBoost::MkBasisVector(level, index);
    }

    inline qreach::DD MkBasisVector(unsigned int level, std::string bitstring) {
        return VectorComplexFloatBoost::MkBasisVector(level, bitstring);
    }

    inline qreach::DD NoDistinctionNode(unsigned int level, qreach::DDComplex val) {
        return VectorComplexFloatBoost::NoDistinctionNode(level, val);
    }

    inline qreach::DD InitializeWithAmplitudes(unsigned int qnum, std::vector<double> amps) {
        // amps is interleaved real/imag: [re0, re1, ..., im0, im1, ...]
        unsigned int n = amps.size();
        assert((n & (n - 1)) == 0 && n != 0);
        assert(2 * (1 << qnum) == n);
        std::vector<std::complex<double>> vec(n / 2);
        for (unsigned int i = 0; i < n / 2; i++) {
            vec[i] = std::complex<double>(amps[i], amps[i + n / 2]);
        }
        // assert vec is normalized
        double norm = 0;
        for (unsigned int i = 0; i < n / 2; i++) {
            norm += std::norm(vec[i]);
        }
        assert(std::abs(norm - 1.0) < 1e-8);
        unsigned int level = ceil(log2(qnum));
        qreach::DD res = VectorComplexFloatBoost::NoDistinctionNode(level, 0);
        for (unsigned int i = 0; i < n / 2; i++) {
            if (std::abs(vec[i]) > 1e-10) {
                auto basisVec = VectorComplexFloatBoost::MkBasisVector(level, i);
                auto scaledVec = qreach::DDComplex(vec[i].real(), vec[i].imag()) * basisVec;
                res = res + scaledVec;
            }
        }
        return VectorComplexFloatBoost::VectorToMatrixInterleaved(res);
    }

    inline qreach::DD VectorToMatrixInterleaved(qreach::DD vec) {
        return VectorComplexFloatBoost::VectorToMatrixInterleaved(vec);
    }

    inline int GetLevel(qreach::DD c) {
        return c.root->level;
    }

    inline bool IsApproximatelyZero(qreach::DD c, double threshold) {
        auto resMap = c.root->rootConnection.returnMapHandle;
        if (resMap.Size() == 0) {
            return true;
        }
        auto sum = abs(resMap[0].real()) + abs(resMap[0].imag());
        for (unsigned int i = 1; i < resMap.Size(); i++) {
            sum += (abs(resMap[i].real()) + abs(resMap[i].imag()));
            if (sum > threshold) {
                return false;
            }
        }
        return true;
    }

    inline qreach::DDComplex ExtractSingleAmplitude(qreach::DD c) {
        auto resMap = c.root->rootConnection.returnMapHandle;

        // FALLBACK (CFLOBDD transpose-corruption workaround, moved from
        // SingleVecTerm::dot): if transpose corrupts the row vector, retMapSz
        // may exceed 2. Extract the [0,0] entry directly.
        if (resMap.Size() > 2) {
            std::cerr << "Warning: ExtractSingleAmplitude() retMapSz=" << resMap.Size()
                      << " > 2, falling back to [0,0] entry." << std::endl;
            unsigned int idxBits = 1 << (c.root->level - 1);
            SH_OBDD::Assignment a(2 * idxBits);
            for (unsigned int k = 0; k < 2 * idxBits; k++) a[k] = false;
            return c.root->EvaluateIteratively(a);
        }
        // end FALLBACK

        if (resMap.Size() == 0) {
            return qreach::DDComplex(0.0, 0.0);
        }
        if (resMap.Size() == 2) {
            return (resMap[0] != 0) ? resMap[0] : resMap[1];
        }
        return resMap[0];
    }

    inline std::vector<std::pair<unsigned int, qreach::DDComplex>> GetNonZeroAmplitudes(
        qreach::DD c, double threshold) {
        std::vector<std::pair<unsigned int, qreach::DDComplex>> out;
        auto resMap = c.root->rootConnection.returnMapHandle;
        for (unsigned int i = 0; i < resMap.Size(); i++) {
            auto mag = abs(resMap[i].real()) + abs(resMap[i].imag());
            if (mag > threshold) {
                out.push_back({i, resMap[i]});
            }
        }
        return out;
    }

    inline qreach::DD Normalize(qreach::DD c) {
        // H*content identity-multiply path (see quantum_operation.hpp history:
        // avoids MatrixTranspose SIGSEGV on Gram-Schmidt-produced DAGs).
        auto H = DDMatrix::MkSingleQubitGateOnN(
            std::pow(2, c.root->level - 1), 0, DDMatrix::MkIdRelation);
        qreach::DD c1 = DDMatrix::MatrixMultiplyWithVector(H, c);
        qreach::DD c1_conj = DDMatrix::Conjugate(c1);
        c1_conj = DDMatrix::Transpose(c1_conj);
        auto mulres = DDMatrix::MatrixMultiply(c1_conj, c1);
        auto resMap = mulres.root->rootConnection.returnMapHandle;

        // Maybe #BUGS here!
        double dimfactor = std::pow(double(2), double(std::pow(2, c.root->level - 1) - 1));

        // FALLBACK: if the H*content path also hits retMapSz>2, use direct
        // row-evaluation (CFLOBDD transpose-corruption workaround).
        if (resMap.Size() > 2) {
            std::cerr << "Warning: DDVector::Normalize() H*content retMapSz="
                      << resMap.Size() << " > 2, using direct row-eval." << std::endl;
            unsigned int level = c.root->level;
            unsigned int indexBits = 1 << (level - 1);
            unsigned int totalBits = 2 * indexBits;
            unsigned long int numRows = 1UL << indexBits;
            qreach::DDComplex normsq = 0;
            SH_OBDD::Assignment a(totalBits);
            for (unsigned long int row = 0; row < numRows; row++) {
                unsigned long int mask = 1UL;
                for (int k = indexBits - 1; k >= 0; k--) {
                    a[2 * k] = (row & mask) ? true : false;
                    mask <<= 1;
                }
                for (unsigned int k = 0; k < indexBits; k++) {
                    a[2 * k + 1] = false;
                }
                qreach::DDComplex val = c.root->EvaluateIteratively(a);
                normsq += val.real() * val.real() + val.imag() * val.imag();
            }
            double factor = double(sqrt(normsq));
            assert(factor > 0);
            if (std::abs(factor - 1.0) < 1e-10) return c;
            return (1.0 / factor) * c;
        }

        assert(resMap.Size() <= 2);
        qreach::DDComplex amp;
        if (resMap.Size() == 2) {
            amp = (resMap[0] != 0) ? resMap[0] : resMap[1];
            assert(abs(amp.imag() * dimfactor) < 1e-8 && amp.real() > 0);
            double factor = double(sqrt(amp.real()));
            c1 = (1 / factor) * c1;
        } else {
            std::cout << "Warning: DDVector::Normalize() has only one factor!" << std::endl;
            amp = resMap[0];
            assert(abs(amp.imag() * dimfactor) < 1e-8 && abs(amp.real() * dimfactor) < 1e-8);
            c1 = NoDistinctionNode(c.root->level, 0);
        }
        return c1;
    }

    inline qreach::DDComplex InnerProduct(qreach::DD a, qreach::DD b) {
        // <a|b> = conj(a)·b. CFLOBDD represents vectors in matrix form, so this
        // is transpose(a) → conjugate → matrix-multiply(b) → extract scalar.
        // (Moved from SingleVecTerm::dot.)
        auto tmpVec = DDMatrix::Transpose(a);
        tmpVec = DDMatrix::Conjugate(tmpVec);
        auto tmp = DDMatrix::MatrixMultiply(tmpVec, b);
        return ExtractSingleAmplitude(tmp);
    }

    inline void VectorPrintColumnHead(qreach::DD c, std::ostream& out) {
        VectorComplexFloatBoost::VectorPrintColumnHead(c, out);
    }

} // namespace DDVector

namespace DDMatrix {

    inline void Initialize() {
        // Full CFLOBDD global init (matches transition_system.hpp /
        // transition_system_qadd.hpp initializeTransitionSystem()). The
        // shared node tables + product caches are needed by both vectors and
        // matrices; the matrix-specific initializer follows. Call
        // DDMatrix::Initialize() before DDVector::Initialize().
        CFLOBDDNodeHandle::InitNoDistinctionTable();
        CFLOBDDNodeHandle::InitAdditionInterleavedTable();
        CFLOBDDNodeHandle::InitReduceCache();
        InitPairProductCache();
        InitTripleProductCache();
        Matrix1234ComplexFloatBoost::Matrix1234Initializer();
    }

    inline qreach::DD MkIdRelation(unsigned int level) {
        return Matrix1234ComplexFloatBoost::MkIdRelationInterleaved(level);
    }

    inline qreach::DD MkWalsh(unsigned int level) {
        return Matrix1234ComplexFloatBoost::MkWalshInterleaved(level);
    }

    inline qreach::DD MkNegation(unsigned int level) {
        return Matrix1234ComplexFloatBoost::MkNegationMatrixInterleaved(level);
    }

    inline qreach::DD MkPauliY(unsigned int level) {
        return Matrix1234ComplexFloatBoost::MkPauliYMatrixInterleaved(level);
    }

    inline qreach::DD MkPauliZ(unsigned int level) {
        return Matrix1234ComplexFloatBoost::MkPauliZMatrixInterleaved(level);
    }

    inline qreach::DD MkSGate(unsigned int level) {
        return Matrix1234ComplexFloatBoost::MkSGateInterleaved(level);
    }

    inline qreach::DD MkPhaseShift(unsigned int level, double theta) {
        return Matrix1234ComplexFloatBoost::MkPhaseShiftGateInterleaved(level, theta);
    }

    inline qreach::DD MkU3(unsigned int level, std::vector<double> params) {
        return Matrix1234ComplexFloatBoost::MkU3GateInterleaved(level, params);
    }

    inline qreach::DD MkArbitrary(unsigned int level, std::vector<double> params) {
        return Matrix1234ComplexFloatBoost::MkArbitraryGateInterleaved(level, params);
    }

    inline qreach::DD MkCNOT(unsigned int level, unsigned int n, long ctrl, long tgt) {
        return Matrix1234ComplexFloatBoost::MkCNOT(level, n, ctrl, tgt);
    }

    inline qreach::DD MkCCNOT(unsigned int level, unsigned int n, long c1, long c2, long tgt) {
        return Matrix1234ComplexFloatBoost::MkCCNOT(level, n, c1, c2, tgt);
    }

    inline qreach::DD MkSwap(unsigned int level, long i, long j) {
        return Matrix1234ComplexFloatBoost::MkSwapGate(level, i, j);
    }

    inline qreach::DD MkiSwap(unsigned int level, long i, long j) {
        return Matrix1234ComplexFloatBoost::MkiSwapGate(level, i, j);
    }

    inline qreach::DD MkCP(unsigned int level, long ctrl, long tgt, double theta) {
        return Matrix1234ComplexFloatBoost::MkCPGate(level, ctrl, tgt, theta);
    }

    inline qreach::DD MkSingleQubitGateOnN(unsigned int n, unsigned int target,
                                           qreach::DD(*gate1q)(unsigned int)) {
        if (n == 1) {
            return gate1q(1);
        }
        int level = ceil(log2(n / 2));
        if (target < n / 2) {
            qreach::DD T = MkIdRelation(level + 1);
            qreach::DD H = MkSingleQubitGateOnN(n / 2, target, gate1q);
            return KroneckerProduct(H, T);
        } else {
            qreach::DD T = MkIdRelation(level + 1);
            return KroneckerProduct(T, MkSingleQubitGateOnN(n / 2, target - n / 2, gate1q));
        }
    }

    inline qreach::DD MkSingleQubitGateOnNWithParam(
        unsigned int n, unsigned int target,
        qreach::DD(*gate1q)(unsigned int, double), double theta) {
        if (n == 1) {
            return gate1q(1, theta);
        }
        int level = ceil(log2(n / 2));
        if (target < n / 2) {
            qreach::DD T = MkIdRelation(level + 1);
            qreach::DD H = MkSingleQubitGateOnNWithParam(n / 2, target, gate1q, theta);
            return KroneckerProduct(H, T);
        } else {
            qreach::DD T = MkIdRelation(level + 1);
            return KroneckerProduct(T, MkSingleQubitGateOnNWithParam(n / 2, target - n / 2, gate1q, theta));
        }
    }

    inline qreach::DD MkSingleQubitGateOnNWithParamVec(
        unsigned int n, unsigned int target,
        qreach::DD(*gate1q)(unsigned int, std::vector<double>), std::vector<double> v) {
        if (n == 1) {
            return gate1q(1, v);
        }
        int level = ceil(log2(n / 2));
        if (target < n / 2) {
            qreach::DD T = MkIdRelation(level + 1);
            qreach::DD H = MkSingleQubitGateOnNWithParamVec(n / 2, target, gate1q, v);
            return KroneckerProduct(H, T);
        } else {
            qreach::DD T = MkIdRelation(level + 1);
            return KroneckerProduct(T, MkSingleQubitGateOnNWithParamVec(n / 2, target - n / 2, gate1q, v));
        }
    }

    inline qreach::DD KroneckerProduct(qreach::DD a, qreach::DD b) {
        return Matrix1234ComplexFloatBoost::KroneckerProduct2Vocs(a, b);
    }

    inline qreach::DD MatrixMultiply(qreach::DD a, qreach::DD b) {
        // matrix × matrix. CFLOBDD's hot path is the "WithInfo" variant.
        return Matrix1234ComplexFloatBoost::MatrixMultiplyV4WithInfo(a, b);
    }

    inline qreach::DD MatrixMultiplyWithVector(qreach::DD gate, qreach::DD vec) {
        return Matrix1234ComplexFloatBoost::MatrixMultiplyV4WithInfo(gate, vec);
    }

    inline qreach::DD Conjugate(qreach::DD c) {
        return Matrix1234ComplexFloatBoost::MatrixConjugate(c);
    }

    inline qreach::DD Transpose(qreach::DD c) {
        return Matrix1234ComplexFloatBoost::MatrixTranspose(c);
    }

} // namespace DDMatrix

#endif // DD_BACKEND_CFLOBDD_H
