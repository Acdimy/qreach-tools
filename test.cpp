#include <iostream>
#include <memory>
#include "transition_system_qadd.hpp"
using namespace qts;

// Test the basic functionality of transition_system
int main() {
    initializeTransitionSystem();
    int qNum = 16;
    // Create a all-zero string with length qNum
    std::vector<std::string> terms = {std::string(qNum, '0')};
    QOperation op(terms);
    TransitionSystem ts(qNum);
    int loc0 = ts.addLocation();
    ts.setAnnotation(loc0, op);
    int loc1 = ts.addLocation();
    QOperation oph(std::string("H"), qNum, std::vector<unsigned int>{0}, std::vector<double>{});
    ts.addRelation(loc0, loc1, oph);
    std::vector<int> currLoc = {loc1};
    int totalLocs = 2;
    for(int i = 0; i <= qNum - 1; i++) {
        // Apply a measurement on the i-th qubit, and add a relation from loc0 to loc1 with the projective operations meas0 and meas1.
        QOperation meas0(std::string("meas0"), qNum, std::vector<unsigned int>{static_cast<unsigned int>(i)}, std::vector<double>{});
        QOperation meas1(std::string("meas1"), qNum, std::vector<unsigned int>{static_cast<unsigned int>(i)}, std::vector<double>{});
        // For each location in currLoc, create two new locations for the two post-locations of the projective measurement.
        std::vector<int> newLocs;
        for (int loc : currLoc) {
            int newLoc0 = ts.addLocation();
            int newLoc1 = ts.addLocation();
            ts.addRelation(loc, newLoc0, meas0);
            ts.addRelation(loc, newLoc1, meas1);
            newLocs.push_back(newLoc0);
            newLocs.push_back(newLoc1);
            totalLocs += 2;
        }
        currLoc = newLocs;
    }
    // Record time comsumption
    auto start = std::chrono::high_resolution_clock::now();
    ts.postConditions();
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    std::cout << "Time consumed: " << duration.count() << " ms" << std::endl;
    std::cout << "Total locations: " << totalLocs << std::endl;
    // ts.printAnnotation();
    // ts.printRelation();
    std::cout << "Complete" << std::endl;
    return 0;
}
