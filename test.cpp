#include <iostream>
#include <memory>
#include "transition_system.hpp"

// Test the basic functionality of transition_system
int main() {
    TransitionSystem ts;
    ts.initLocation = 0; // Initialize the first location
    std::vector<std::string> strings =  {std::string("00000000")};
    QOperation op(strings);
    QOperation oph(std::string("H"), 8, std::vector<unsigned int>{0}, std::vector<double>{});
    QOperation opi(std::string("I"), 8, std::vector<unsigned int>{0}, std::vector<double>{});
    QOperation opx(std::string("X"), 8, std::vector<unsigned int>{0}, std::vector<double>{});
    QOperation opz(std::string("Z"), 8, std::vector<unsigned int>{0}, std::vector<double>{});
    QOperation opy(std::string("Y"), 8, std::vector<unsigned int>{0}, std::vector<double>{});
    QOperation opm0(std::string("meas0"), 8, std::vector<unsigned int>{0}, std::vector<double>{});
    QOperation opm1(std::string("meas1"), 8, std::vector<unsigned int>{0}, std::vector<double>{});
    Location loc(8,0); // Create a location with 8 qubits
    Location loc1(8,1);
    Location loc2(8,2);
    Location loc3(8,3);
    Location loc4(8,4);
    ts.addLocation(loc);
    ts.addLocation(loc1);
    ts.addLocation(loc2);
    ts.addLocation(loc3);
    ts.addLocation(loc4);
    loc.appendPostLocation(&loc1); // Add a post-location relation
    loc1.appendPostLocation(&loc2);
    loc1.appendPostLocation(&loc3);
    loc2.appendPostLocation(&loc4);
    loc3.appendPostLocation(&loc4);
    ts.setInitLocation(0); // Set the initial location to 0
    ts.addRelation(0, 1, oph); // Add a self-loop relation with the operation
    ts.addRelation(1, 2, opm0);
    ts.addRelation(1, 3, opm1);
    ts.addRelation(2, 4, opi);
    ts.addRelation(3, 4, opi); // Add a self-loop relation for location 3
    ts.addRelation(4, 4, opi);
    // ts.addRelation(4, 0, opi);
    ts.setAnnotation({{0, op}}); // Set the annotation for the location
    ComputingFixedPointPost(ts); // Compute the fixed point of the post-conditions
    // The result should be zero space!
    ts.printDims(0);
    ts.printDims(1);
    ts.printDims(2);
    ts.printDims(3);
    ts.printDims(4);
    ts.printSupp(4);

    // QOperation test = CreateIdentityQO(8);
    // QOperation test2 = CreateZeroQO(8);
    // std::cout << "Test QOperation: " << test.isIdentity << " " << test.qNum << std::endl;
    // std::cout << "Test2 QOperation: " << test2.isIdentity << " " << test2.qNum << std::endl;

    return 0;
}
