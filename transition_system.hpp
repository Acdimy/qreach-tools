#include "quantum_operation.hpp"
#ifndef TRANSITION
#define TRANSITION

//TODO: check the unique table and the computing table

using namespace CFL_OBDD;

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

Location::~Location()
{
}

// void Location::setAnnotation(QOperation annotation)
// {
//     // this->annotation = annotation;
// }
void Location::appendPreLocation(Location* loc)
{
    this->preLocations.push_back(loc);
}
void Location::appendPostLocation(Location* loc)
{
    this->postLocations.push_back(loc);
}

/***
 * Transition system for model checking algorithms.
 ***/

class TransitionSystem
{
private:
    std::vector<int> currPreLocs;
    std::vector<int> currPostLocs;
    std::vector<bool> visitedPre;
    std::vector<bool> visitedPost;
public:
    std::vector<Location> Locations;
    unsigned int initLocation;
    std::map<std::tuple<unsigned int, unsigned int>, QOperation> relations;
    // I need to build a map from relation names to the quantum operations.
public:
    TransitionSystem() {};
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
};

void TransitionSystem::addLocation(Location loc)
{
    this->Locations.push_back(loc);
    loc.idx = this->Locations.size() - 1;
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
}

TransitionSystem::~TransitionSystem()
{
}

void TransitionSystem::preConditionInit() {
    /*
    Initialize the pre-condition of the transition system.
    Select all locations whose postLocations are only themselves. Append them to currPreLocs.
    */
    for (unsigned int i = 0; i < this->Locations.size(); i++) {
        if (this->Locations[i].postLocations.size() == 1 && this->Locations[i].postLocations[0] == &this->Locations[i]) {
            this->currPreLocs.push_back(i);
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
    for (unsigned int i = 0; i < this->Locations[loc].preLocations.size(); i++) {
        Location* preLoc = this->Locations[loc].preLocations[i];
        QOperation preImage = this->Locations[loc].upperBound.preImage(this->relations[std::make_tuple(loc, preLoc->idx)]);
        preLoc->upperBound = preLoc->upperBound.conjunction_simp(preImage); // TODO: Conjunction inline
    }
}

void TransitionSystem::preConditions() {
    /*
    Compute the weakest pre-condition of the annotation in all currPreLocs.
    Use the method QOperation::preImage
    For different locations having the same predecessor, the result of pre-image should be conjuncted.
    */
   for (unsigned int i = 0; i < this->currPreLocs.size(); i++) {
        unsigned int loc = this->currPreLocs[i];
        this->preConditionOneStep(loc);
    }
}

void TransitionSystem::postConditionInit() {
    /*
    Initialize the post-condition of the transition system.
    Select all locations whose preLocations are only themselves. Append them to currPostLocs.
    */
    for (unsigned int i = 0; i < this->Locations.size(); i++) {
        if (this->Locations[i].preLocations.size() == 1 && this->Locations[i].preLocations[0] == &this->Locations[i]) {
            this->currPostLocs.push_back(i);
        }
    }
}

void TransitionSystem::postConditionOneStep(unsigned int loc) {
    /*
    For the location loc, compute the post-condition of the lowerBound in loc.
    Tranverse all the postLocations of loc, and compute the post-image of the lowerBound in loc.
    Use the method QOperation::postImage
    Do the conjunction with the existed lowerBound of postLocations (Use QOperation.conjunction).
    */
    for (unsigned int i = 0; i < this->Locations[loc].postLocations.size(); i++) {
        Location* postLoc = this->Locations[loc].postLocations[i];
        QOperation postImage = this->Locations[loc].lowerBound.postImage(this->relations[std::make_tuple(loc, postLoc->idx)]);
        postLoc->lowerBound = postLoc->lowerBound.disjunction(postImage); // TODO: Conjunction inline
    }
}

void TransitionSystem::postConditions() {
    /*
    Compute the weakest post-condition of the annotation in all currPostLocs.
    Use the method QOperation::postImage
    For different locations having the same successor, the result of post-image should be conjuncted.
    */
   for (unsigned int i = 0; i < this->currPostLocs.size(); i++) {
        unsigned int loc = this->currPostLocs[i];
        this->postConditionOneStep(loc);
    }
}

void ComputingFixedPoint(TransitionSystem& ts) {
    /*
    Compute the fixed point of the transition system.
    The fixed point is the status when preConditions and postConditions converge.
    */
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
