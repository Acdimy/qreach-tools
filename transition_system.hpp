#include "quantum_operation.hpp"
#include <functional>
#include <deque>
#ifndef TRANSITION
#define TRANSITION

//TODO: check the unique table and the computing table
// It seems that CFLOBDD does not provide the two tables.

using namespace CFL_OBDD;

// Provide a hash specialization for std::tuple<unsigned int, unsigned int, unsigned int, unsigned int>
namespace std {
    template <>
    struct hash<std::tuple<unsigned int, unsigned int, unsigned int, unsigned int>> {
        std::size_t operator()(const std::tuple<unsigned int, unsigned int, unsigned int, unsigned int>& t) const {
            std::size_t h1 = std::hash<unsigned int>()(std::get<0>(t));
            std::size_t h2 = std::hash<unsigned int>()(std::get<1>(t));
            std::size_t h3 = std::hash<unsigned int>()(std::get<2>(t));
            std::size_t h4 = std::hash<unsigned int>()(std::get<3>(t));
            // Combine the hashes
            return (((h1 ^ (h2 << 1)) >> 1) ^ (h3 << 1)) ^ (h4 << 1);
        }
    };
}

class Location
{
public:
    unsigned int idx;
    QOperation upperBound;
    QOperation lowerBound;
    // QOperation annotation;
    std::vector<Location*> preLocations;
    std::vector<Location*> postLocations;
public:
    Location(/* args */) {};
    Location(const int qNum);
    Location(const int qNum, const unsigned int idx);
    ~Location();
    // void setAnnotation(QOperation annotation);
    void appendPreLocation(Location* loc);
    void appendPostLocation(Location* loc);
};

Location::Location(const int qNum)
{
    this->upperBound = CreateIdentityQO(qNum);
    this->lowerBound = CreateZeroQO(qNum);
}

Location::Location(const int qNum, const unsigned int idx)
{
    this->idx = idx;
    this->upperBound = CreateIdentityQO(qNum);
    this->lowerBound = CreateZeroQO(qNum);
}

Location::~Location()
{
}

// void Location::setAnnotation(QOperation annotation)
// {
//     // this->annotation = annotation;
// }
void Location::appendPreLocation(Location* loc)
{
    // Make sure loc is not in this->preLocations
    if (std::find(this->preLocations.begin(), this->preLocations.end(), loc) != this->preLocations.end()) {
        return; // Already exists
    }
    this->preLocations.push_back(loc);
}
void Location::appendPostLocation(Location* loc)
{
    // Make sure loc is not in this->postLocations
    if (std::find(this->postLocations.begin(), this->postLocations.end(), loc) != this->postLocations.end()) {
        return; // Already exists
    }
    this->postLocations.push_back(loc);
    loc->appendPreLocation(this); // Ensure the reverse relation is also established
}

/***
 * Transition system for model checking algorithms.
 ***/

class TransitionSystem
{
public:
    std::deque<int> currPreLocs;
    std::deque<int> currPostLocs;
    std::vector<bool> visitedPre;
    std::vector<bool> visitedPost;
public:
    std::vector<Location> Locations;
    unsigned int initLocation;
    // I need to build a map from relation names to the quantum operations.
    std::map<std::tuple<unsigned int, unsigned int>, QOperation> relations;

    
    // A computed table: tuple(from_loc, from_dim, to_loc, to_dim), indicating the relation(from_loc, to_loc) with the dimensions is computed.
    std::unordered_set<std::tuple<unsigned int, unsigned int, unsigned int, unsigned int>> computedTablePre;
    std::unordered_set<std::tuple<unsigned int, unsigned int, unsigned int, unsigned int>> computedTablePost;
    
public:
    TransitionSystem() {
        CFLOBDDNodeHandle::InitNoDistinctionTable();
        CFLOBDDNodeHandle::InitAdditionInterleavedTable();
        CFLOBDDNodeHandle::InitReduceCache();
        InitPairProductCache();
        InitTripleProductCache();
        Matrix1234ComplexFloatBoost::Matrix1234Initializer();
        VectorComplexFloatBoost::VectorInitializer();
    };
    ~TransitionSystem();
    // void initialization();
    void addLocation(Location loc);
    void addRelation(unsigned int from, unsigned int to, QOperation op);
    void setAnnotation(std::vector<std::tuple<unsigned int, QOperation>> annotations);
    void preConditionInit();
    void preConditionOneStep(unsigned int loc);
    void preConditions();
    void postConditionInit();
    void postConditionOneStep(unsigned int loc);
    void postConditions();
    void createAdd();
    bool satisfy();
    void setInitLocation(unsigned int loc);
    void printDims(unsigned int loc);
    void printSupp(unsigned int loc);
};

void TransitionSystem::addLocation(Location loc)
{
    this->Locations.push_back(loc);
    loc.idx = static_cast<int>(this->Locations.size()) - 1;
}

void TransitionSystem::addRelation(unsigned int from, unsigned int to, QOperation op)
{
    // The problem is how to represent a projective operation' support vectors.
    this->relations[std::make_tuple(from, to)] = op;
    this->Locations[from].appendPostLocation(&this->Locations[to]);
    this->Locations[to].appendPreLocation(&this->Locations[from]);
}

void TransitionSystem::setAnnotation(std::vector<std::tuple<unsigned int, QOperation>> annotations)
{
    // We don't allow the same location to have different annotations yet.
    for (const auto& annotation : annotations) {
        unsigned int loc = std::get<0>(annotation);
        QOperation op = std::get<1>(annotation);
        // this->Locations[loc].setAnnotation(op);
        this->currPreLocs.push_back(loc);
        this->currPostLocs.push_back(loc);
        this->Locations[loc].upperBound = op; // Must be copy assignment!
        this->Locations[loc].lowerBound = op;
    }
    // For any location that is not in annotations, set the upperBound to zero and lowerBound to identity. (Already implemented in the constructor of Location)
    // for (unsigned int i = 0; i < this->Locations.size(); i++) {
    //     if (std::find(this->currPreLocs.begin(), this->currPreLocs.end(), i) == this->currPreLocs.end()) {
    //         this->Locations[i].upperBound = CreateZeroQO(this->Locations[i].upperBound.qNum);
    //         this->Locations[i].lowerBound = CreateIdentityQO(this->Locations[i].lowerBound.qNum);
    //     }
    // }
}

void TransitionSystem::setInitLocation(unsigned int loc)
{
    this->initLocation = loc;
}

TransitionSystem::~TransitionSystem()
{
}

void TransitionSystem::preConditionInit() {
    /*
    Initialize the pre-condition of the transition system.
    Select all locations whose postLocations are only themselves. Append them to currPreLocs.
    */
   // Initialize the visitedPre vectors
    this->visitedPre.resize(this->Locations.size(), false);
    this->computedTablePre.clear(); // Clear the computed table for pre-conditions
    for (unsigned int i = 0; i < this->Locations.size(); i++) {
        if ((this->Locations[i].postLocations.size() == 1 && this->Locations[i].postLocations[0] == &this->Locations[i]) || this->Locations[i].postLocations.size() == 0) {
            // If i not in currPreLocs, add it. (TODO: Use a set to avoid duplicates)
            if (std::find(this->currPreLocs.begin(), this->currPreLocs.end(), i) == this->currPreLocs.end()) {
                this->currPreLocs.push_back(i);
                this->visitedPre[i] = true; // Mark this location as visited
            }
        }
    }
    
}

void TransitionSystem::preConditionOneStep(unsigned int loc) {
    /*
    For the location loc, compute the pre-condition of the upperBound in loc.
    Tranverse all the preLocations of loc, and compute the pre-image of the upperBound in loc.
    Use the method QOperation::preImage
    Do the conjunction with the existed upperBound of preLocations (Use QOperation.conjunction).
    */
    std::cout << this->Locations[loc].preLocations.size() << " pre locations for location " << loc << std::endl;
    for (unsigned int i = 0; i < this->Locations[loc].preLocations.size(); i++) {
        Location* preLoc = this->Locations[loc].preLocations[i];
        // If it is a self-loop and the relation is identity, skip it without appending currPreLoc.
        if (preLoc->idx == loc && this->relations[std::make_tuple(preLoc->idx, loc)].isIdentity) {
            std::cout << "Skip self-loop for location " << loc << std::endl;
            continue;
        }
        if (this->computedTablePre.find(std::make_tuple(loc, this->Locations[loc].upperBound.oplist.size(), preLoc->idx, preLoc->upperBound.oplist.size())) == this->computedTablePre.end()) {
            std::cout << this->relations[std::make_tuple(preLoc->idx, loc)].oplist.size() << " operations in relation from " << loc << " to " << preLoc->idx << std::endl;
            // std::cout << this->relations[std::make_tuple(preLoc->idx, loc)].type << std::endl;
            QOperation preImage = this->Locations[loc].upperBound.preImage(this->relations[std::make_tuple(preLoc->idx, loc)]);
            computedTablePre.insert(std::make_tuple(loc, this->Locations[loc].upperBound.oplist.size(), preLoc->idx, preLoc->upperBound.oplist.size()));
            int dimBefore = preLoc->upperBound.oplist.size(); // Upper bound: dimension 2^n as default
            preLoc->upperBound = preLoc->upperBound.conjunction_simp(preImage); // TODO: Conjunction inline
            if (visitedPre[preLoc->idx] == false) {
                std::cout << "Visit a new pre location " << preLoc->idx << std::endl;
                this->currPreLocs.push_back(preLoc->idx);
                visitedPre[preLoc->idx] = true;
            } else if (preLoc->upperBound.oplist.size() < dimBefore && preLoc->upperBound.oplist.size() > 0) {
                // If the dimension of the upperBound is reduced, we need to recheck the pre-condition.
                if (std::find(this->currPreLocs.begin(), this->currPreLocs.end(), preLoc->idx) == this->currPreLocs.end()) {
                    std::cout << "Pre condition for location " << preLoc->idx << " is updated from " << dimBefore << " to " << preLoc->lowerBound.oplist.size() << std::endl;
                    this->currPreLocs.push_back(preLoc->idx);
                }
            } else {
                std::cout << "Pre condition for location " << preLoc->idx << " is not updated." << std::endl;
            }
        } else {}
    }
    std::cout << "Pre condition for location " << loc << " computed." << std::endl;
}

void TransitionSystem::preConditions() {
    /*
    Compute the weakest pre-condition of the annotation in all currPreLocs.
    Use the method QOperation::preImage
    For different locations having the same predecessor, the result of pre-image should be conjuncted.
    */
    while (!this->currPreLocs.empty()) {
        unsigned int loc = this->currPreLocs.front();
        this->currPreLocs.pop_front();
        this->preConditionOneStep(loc);
    }
}

void TransitionSystem::postConditionInit() {
    /*
    Initialize the post-condition of the transition system.
    Select all locations whose preLocations are only themselves. Append them to currPostLocs.
    */
    // Initialize the visitedPost vectors
    this->visitedPost.resize(this->Locations.size(), false);
    this->computedTablePost.clear(); // Clear the computed table for post-conditions
    for (unsigned int i = 0; i < this->Locations.size(); i++) {
        if ((this->Locations[i].preLocations.size() == 1 && this->Locations[i].preLocations[0] == &this->Locations[i]) || this->Locations[i].preLocations.size() == 0) {
            // If i not in currPostLocs, add it.
            if (std::find(this->currPostLocs.begin(), this->currPostLocs.end(), i) == this->currPostLocs.end()) {
                this->currPostLocs.push_back(i);
                this->visitedPost[i] = true; // Mark this location as visited
            }
        }
    }
}

void TransitionSystem::postConditionOneStep(unsigned int loc) {
    /*
    For the location loc, compute the post-condition of the lowerBound in loc.
    Tranverse all the postLocations of loc, and compute the post-image of the lowerBound in loc.
    Use the method QOperation::postImage
    Do the disjunction with the existed lowerBound of postLocations (Use QOperation.disjunction).
    */
    std::cout << this->Locations[loc].postLocations.size() << " post locations for location " << loc << std::endl;
    for (unsigned int i = 0; i < this->Locations[loc].postLocations.size(); i++) {
        Location* postLoc = this->Locations[loc].postLocations[i];
        // If it is a self-loop and the relation is identity, skip it without appending currPostLocs.
        if (postLoc->idx == loc && this->relations[std::make_tuple(loc, postLoc->idx)].isIdentity) {
            std::cout << "Skip self-loop for location " << loc << std::endl;
            continue;
        }
        // If (loc, locDim, postLoc, postLocDim) is not computed, compute it.
        if (this->computedTablePost.find(std::make_tuple(loc, this->Locations[loc].lowerBound.oplist.size(), postLoc->idx, postLoc->lowerBound.oplist.size())) == this->computedTablePost.end()) {
            std::cout << this->relations[std::make_tuple(loc, postLoc->idx)].oplist.size() << " operations in relation from " << loc << " to " << postLoc->idx << std::endl;
            QOperation postImage = this->Locations[loc].lowerBound.postImage(this->relations[std::make_tuple(loc, postLoc->idx)]);
            computedTablePost.insert(std::make_tuple(loc, this->Locations[loc].lowerBound.oplist.size(), postLoc->idx, postLoc->lowerBound.oplist.size()));
            int dimBefore = postLoc->lowerBound.oplist.size(); // Lower bound: dimension 0 as default
            postLoc->lowerBound = postLoc->lowerBound.disjunction(postImage); // TODO: Disjunction inline
            if (visitedPost[postLoc->idx] == false) {
                std::cout << "Visit a new post location " << postLoc->idx << std::endl;
                this->currPostLocs.push_back(postLoc->idx);
                visitedPost[postLoc->idx] = true;
            } else if (postLoc->lowerBound.oplist.size() > dimBefore && postLoc->lowerBound.oplist.size() < std::pow(2, postLoc->lowerBound.qNum)) {
                // If the dimension of the lowerBound is reduced, we need to recheck the post-condition.
                // If postLoc->idx is not in currPostLocs, append it.
                if (std::find(this->currPostLocs.begin(), this->currPostLocs.end(), postLoc->idx) == this->currPostLocs.end()) {
                    std::cout << "Post condition for location " << postLoc->idx << " is updated from " << dimBefore << " to " << postLoc->lowerBound.oplist.size() << std::endl;
                    this->currPostLocs.push_back(postLoc->idx);
                }
            } else {
                std::cout << "Post condition for location " << postLoc->idx << " is not updated." << std::endl;
            }
        } else {}
    }
    std::cout << "Post condition for location " << loc << " computed." << std::endl;
}

void TransitionSystem::postConditions() {
    /*
    Compute the weakest post-condition of the annotation in all currPostLocs.
    Use the method QOperation::postImage
    For different locations having the same successor, the result of post-image should be conjuncted.
    */
    while (!this->currPostLocs.empty()) {
        unsigned int loc = this->currPostLocs.front();
        this->currPostLocs.pop_front();
        this->postConditionOneStep(loc);
    }
}

void TransitionSystem::printDims(unsigned int loc) {
    /*
    Print the dimensions of the upperBound and lowerBound of the location.
    */
    std::cout << "Location " << loc << ": upperBound dimension = " << this->Locations[loc].upperBound.oplist.size() 
              << ", lowerBound dimension = " << this->Locations[loc].lowerBound.oplist.size() << std::endl;
}

void TransitionSystem::printSupp(unsigned int loc) {
    std::cout << "Location " << loc << ": upperBound support = \n";
    this->Locations[loc].upperBound.printFormal();
    std::cout << "Location " << loc << ": lowerBound support = \n";
    this->Locations[loc].lowerBound.printFormal();
}

void ComputingFixedPointPre(TransitionSystem& ts) {
    /*
    Compute the fixed point of the pre-conditions.
    The fixed point is the status when preConditions converge.
    */
    ts.preConditionInit();
    ts.preConditions();
}

void ComputingFixedPointPost(TransitionSystem& ts) {
    /*
    Compute the fixed point of the post-conditions.
    The fixed point is the status when postConditions converge.
    */
    ts.postConditionInit();
    ts.postConditions();
}

void ComputingFixedPoint(TransitionSystem& ts) {
    /*
    Compute the fixed point of the transition system.
    The fixed point is the status when preConditions and postConditions converge.
    */
   ComputingFixedPointPre(ts);
   ComputingFixedPointPost(ts);
}

void fromProgramToTransitionSystem(std::string filename, TransitionSystem& ts) {
    /*
    Read the program from the file and convert it to a transition system.
    The program is in the form of a list of locations and transitions.
    Each location has an annotation and a set of pre- and post-conditions.
    */
   // Need to write a parser!
}

/***
 * Location: The location of a transition system.
 * It contains the state and the annotation of the location.
 ***/

#endif
