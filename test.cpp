#include <iostream>
#include <memory>
#include "transition_system.hpp"

// Test the basic functionality of transition_system
int main() {
    TransitionSystem ts;
    ts.initLocation = 0; // Initialize the first location
    std::vector<std::string> strings =  {std::string("00000000")};
    QOperation op(strings);
    QOperation op2(std::string("H"), 8, std::vector<unsigned int>{0}, std::vector<double>{});
    Location loc(8,0); // Create a location with 8 qubits
    Location loc2(8,1);
    ts.addLocation(loc);
    ts.addLocation(loc2);
    loc.appendPostLocation(&loc2); // Add a post-location relation
    loc2.appendPreLocation(&loc); // Add a pre-location relation
    ts.setInitLocation(0); // Set the initial location to 0
    ts.addRelation(0, 1, op2); // Add a self-loop relation with the operation
    ts.setAnnotation({{0, op}}); // Set the annotation for the location
    ts.postConditionInit(); // Initialize the post-condition
    ts.postConditionOneStep(0); // Compute the post-condition for location 0

    // QOperation test = CreateIdentityQO(8);
    // QOperation test2 = CreateZeroQO(8);
    // std::cout << "Test QOperation: " << test.isIdentity << " " << test.qNum << std::endl;
    // std::cout << "Test2 QOperation: " << test2.isIdentity << " " << test2.qNum << std::endl;

    return 0;
}
