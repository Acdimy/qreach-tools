#include "quantum_operation.hpp"


/***
 * Location: The location of a transition system.
 * It contains the state and the annotation of the location.
 ***/
class Location
{
public:
    QOperation upperBound;
    QOperation lowerBound;
    QOperation annotation;
    std::vector<Location*> preLocations;
    std::vector<Location*> postLocations;
public:
    Location(/* args */) {};
    ~Location();
    void setAnnotation(QOperation annotation);
    void appendPreLocation(Location* loc);
    void appendPostLocation(Location* loc);
};

Location::~Location()
{
}

void Location::setAnnotation(QOperation annotation)
{
    // this->annotation = annotation;
}
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
    std::vector<std::tuple<unsigned int, unsigned int, std::string>> relations;
public:
    TransitionSystem() {};
    ~TransitionSystem();
    // void initialization();
    void addLocation(Location loc);
    void addRelation(unsigned int from, unsigned int to, std::string name);
    void setAnnotation(std::vector<std::tuple<unsigned int, QOperation>> annotations);
    void preConditionInit();
    void preConditionOneStep(int loc);
    void preConditions();
    void postConditionInit();
    void postConditionOneStep(int loc);
    void postConditions();
    void createAdd();
    bool satisfy();
    void setInitLocation(unsigned int loc);
};

void TransitionSystem::addLocation(Location loc)
{
    this->Locations.push_back(loc);
}

void TransitionSystem::addRelation(unsigned int from, unsigned int to, std::string name)
{
    this->relations.push_back(std::make_tuple(from, to, name));
}

void TransitionSystem::setAnnotation(std::vector<std::tuple<unsigned int, QOperation>> annotations)
{
    for (const auto& annotation : annotations) {
        unsigned int loc = std::get<0>(annotation);
        QOperation op = std::get<1>(annotation);
        this->Locations[loc].setAnnotation(op);
        this->currPreLocs.push_back(loc);
        this->currPostLocs.push_back(loc);
    }
}

TransitionSystem::~TransitionSystem()
{
}

void TransitionSystem::preConditions() {
    /*
    Compute the weakest pre-condition of the annotation in all currPreLocs.
    Use the method QOperation::preImage
    For different locations having the same predecessor, the result of pre-image should be conjuncted.
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

