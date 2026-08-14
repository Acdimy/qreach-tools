// Minimal explicit-path semantic-layer smoke test against the LimTDD backend.
//
// Verifies quantum_operation.hpp + transition_system.hpp (qts_naive) compile
// and run correctly when the backend is LimTDD. Same shape as test_newapi.py.
//
// Build (from repo root):
//   g++ -std=c++2a -w -DQREACH_USE_LIMTDD \
//       -I. -I../LimTDDexpr/LimTDD/DDPackage -I../LimTDDexpr/include \
//       scripts/limtdd_explicit_smoke.cpp \
//       ../LimTDDexpr/LimTDD/DDPackage/dd/Edge.cpp \
//       ../LimTDDexpr/LimTDD/DDPackage/dd/Maps.cpp \
//       ../LimTDDexpr/LimTDD/DDPackage/dd/Node.cpp \
//       -o /tmp/limtdd_explicit_smoke

#include "transition_system.hpp"
#include <cassert>
#include <iostream>

int main() {
    qts_naive::initializeTransitionSystem();

    int qNum = 2;
    // Bell pair: H(0); CX(0,1) on |00>  ->  (|00> + |11>)/sqrt(2)
    QOperation init(std::vector<std::string>{std::string(qNum, '0')});
    QOperation h("H", qNum, std::vector<unsigned int>{0}, std::vector<double>{});
    QOperation cx("cx", qNum, std::vector<unsigned int>{0, 1}, std::vector<double>{});

    qts_naive::TransitionSystem ts(qNum);
    ts.addLocation(qts_naive::Location(qNum));  // loc 0
    ts.addLocation(qts_naive::Location(qNum));  // loc 1
    ts.addLocation(qts_naive::Location(qNum));  // loc 2
    ts.setAnnotation(std::vector<std::tuple<unsigned int, QOperation>>{{0, init}});
    ts.setInitLocation(0);
    ts.addRelation(0, 1, h);
    ts.addRelation(1, 2, cx);

    ts.computingFixedPointPost();

    // Bell state (|00>+|11>)/sqrt(2) is a SINGLE state -> lowerDim == 1,
    // upperDim == 2^qNum == 4. (Dimension check only — see NOTE below.)
    auto dims = ts.printDims(2);
    std::cout << "loc2 dims = [" << dims.first << ", " << dims.second << "]" << std::endl;
    assert(dims.first == 4);
    assert(dims.second == 1);

    // NOTE: the Gram-Schmidt / satisfy / disjunction path (which computes
    // inner products of ORTHOGONAL states, i.e. a ZERO result) currently
    // throws "ExtractSingleAmplitude: not exactly one non-zero amplitude" on
    // the LimTDD backend. LimTDD's InnerProduct/ExtractSingleAmplitude must
    // return 0 for the zero-inner-product case instead of asserting. Reported
    // to the LimTDD agent. The post-image path above (H then CX on |00>) is
    // unaffected and produces the correct Bell-state dimension.

    std::cout << "explicit semantic layer with LimTDD: OK (Bell-state post-image)" << std::endl;
    return 0;
}
