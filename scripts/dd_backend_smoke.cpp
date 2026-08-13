// Phase 0 smoke test for the DD backend interface + CFLOBDD passthrough.
//
// Compiles the dd_backend.hpp interface and the dd_backend_cflobdd.h passthrough
// and checks that the wrappers forward to CFLOBDD correctly. This is NOT wired
// into the semantic layer yet (that is Phase 1); it only proves the interface
// header and passthrough are valid.
//
// Build (from repo root, after `make all`):
//   g++ -std=c++2a -I. -I../BOOST/boost_1_81_0 \
//       scripts/dd_backend_smoke.cpp -L. -lqreach -Wl,-rpath,. \
//       -o /tmp/dd_backend_smoke && /tmp/dd_backend_smoke

#include "dd_backend.hpp"
#include <cassert>
#include <iostream>

int main() {
    std::cerr << "step: init matrix" << std::endl;
    DDMatrix::Initialize();
    std::cerr << "step: init vector" << std::endl;
    DDVector::Initialize();

    std::cerr << "step: basis vectors" << std::endl;
    qreach::DD v0 = DDVector::MkBasisVector(1, 0);   // |00>
    qreach::DD v3 = DDVector::MkBasisVector(1, 3);   // |11>
    assert(DDVector::GetLevel(v0) == 1);
    assert(DDVector::GetLevel(v3) == 1);

    std::cerr << "step: zero detection" << std::endl;
    qreach::DD zero = DDVector::NoDistinctionNode(1, 0);
    assert(DDVector::IsApproximatelyZero(zero));
    assert(!DDVector::IsApproximatelyZero(v0));

    std::cerr << "step: extract amplitude" << std::endl;
    qreach::DDComplex a0 = DDVector::ExtractSingleAmplitude(v0);
    qreach::DDComplex a3 = DDVector::ExtractSingleAmplitude(v3);
    assert(a0 == qreach::DDComplex(1.0, 0.0));
    assert(a3 == qreach::DDComplex(1.0, 0.0));

    std::cerr << "step: hash-consing" << std::endl;
    assert(DDMatrix::MkIdRelation(1) == DDMatrix::MkIdRelation(1));
    assert(DDMatrix::MkWalsh(1) == DDMatrix::MkWalsh(1));

    std::cerr << "step: gate embedding" << std::endl;
    qreach::DD h2 = DDMatrix::MkSingleQubitGateOnN(2, 1, DDMatrix::MkWalsh);
    (void)h2;

    std::cerr << "dd_backend smoke: OK" << std::endl;
    return 0;
}
