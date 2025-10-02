#include "quantum_operation.hpp"
#include "cl_proposition.hpp"
#include "QASM_parser.hpp"
#include <functional>
#include <deque>
#include <unordered_map>
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
    unsigned int qNum;
    /* Flag for some delicate settings
     * Only use in pre-condition computation
     * flag == -1: a normal location
     * flag in {0,1,2}: a branching location whose semantic is a binary projective measurement
     * flag == 0: a measurement location when non of the post-locations reached here
     * flag == 1: a measurement location when one of the post-locations reached here
     * flag == 2: a measurement location when all of the post-locations reached here
     */
    int flag = -1;
    std::string identifier;
    QOperation upperBound;
    QOperation lowerBound;
    // QOperation annotation;
    std::vector<unsigned int> preLocations;
    std::vector<unsigned int> postLocations;
    std::vector<QOperation> tempOperations;
    ClassicalProposition cp;
    std::vector<std::string> APs;
public:
    Location(/* args */) {};
    Location(const int qNum);
    Location(const int qNum, const unsigned int idx);
    Location(const int qNum, const unsigned int idx, const unsigned int cNum);
    ~Location();
    // void setAnnotation(QOperation annotation);
    void resetBounds();
    void appendPreLocation(unsigned int loc);
    void appendPostLocation(unsigned int loc);
    void appendClassicalAP(std::string ap);
    void copyClassicalAP(const Location& other) {
        // Copy the classical proposition from another location
        this->cp = other.cp;
    }
    bool satisfy(QOperation spec);
    bool satisfyDefault();
    bool find(std::string ap);
    bool equalAP(const Location& other) const {
        // Check if the classical propositions are equal
        return this->cp == other.cp;
    }
    void setClassicalValue(unsigned int index, bool value) {
        // Set the value of the classical proposition at a specific index
        this->cp.setValue(index, value);
    }
    std::vector<std::string> satisfyBit(std::vector<unsigned int> indexs, std::vector<bool> values);
    std::vector<std::string> unsatisfyBit(std::vector<unsigned int> indexs, std::vector<bool> values);
    int termNum() const {
        return static_cast<int>(this->cp.terms.size());
    }
    void setLabel(std::string label) {
        this->APs.push_back(label);
    }
    std::vector<std::string> getLabels() const {
        return this->APs;
    }
    void setIdentifier(std::string id) {
        this->identifier = id;
    }
    std::string getIdentifier() const {
        return this->identifier;
    }
};

Location::Location(const int qNum)
{
    this->qNum = qNum;
    this->upperBound = CreateIdentityQO(qNum);
    this->lowerBound = CreateZeroQO(qNum);
}

Location::Location(const int qNum, const unsigned int idx)
{
    this->qNum = qNum;
    this->idx = idx;
    this->upperBound = CreateIdentityQO(qNum);
    this->lowerBound = CreateZeroQO(qNum);
}

Location::Location(const int qNum, const unsigned int idx, const unsigned int cNum)
{
    this->qNum = qNum;
    this->idx = idx;
    this->upperBound = CreateIdentityQO(qNum);
    this->lowerBound = CreateZeroQO(qNum);
    this->cp = ClassicalProposition(cNum);
}

Location::~Location()
{
}

void Location::resetBounds()
{
    this->upperBound = CreateIdentityQO(this->qNum);
    this->lowerBound = CreateZeroQO(this->qNum);
}

// void Location::setAnnotation(QOperation annotation)
// {
//     // this->annotation = annotation;
// }
void Location::appendPreLocation(unsigned int loc)
{
    // Make sure loc is not in this->preLocations
    if (std::find(this->preLocations.begin(), this->preLocations.end(), loc) != this->preLocations.end()) {
        return; // Already exists
    }
    this->preLocations.push_back(loc);
}
void Location::appendPostLocation(unsigned int loc)
{
    // Make sure loc is not in this->postLocations
    if (std::find(this->postLocations.begin(), this->postLocations.end(), loc) != this->postLocations.end()) {
        return; // Already exists
    }
    this->postLocations.push_back(loc);
    // loc->appendPreLocation(this->idx); // Ensure the reverse relation is also established
}

void Location::appendClassicalAP(std::string ap)
{
    this->cp.addTerm(ap);
}

bool Location::satisfy(QOperation spec)
{
    assert(spec.isProj < 0); // Only projective operations are supported for now
    // std::cout << this->lowerBound.compare(spec) << " " << spec.compare(this->upperBound) << std::endl;
    return (this->lowerBound.compare(spec) % 4 == 0 && spec.compare(this->upperBound) % 4 == 0);
}

bool Location::satisfyDefault()
{
    return this->lowerBound.compare(this->upperBound) % 4 == 0;
}

bool Location::find(std::string ap)
{
    // Check if the classical proposition satisfies the given AP
    return this->cp.find(ap);
}

std::vector<std::string> Location::satisfyBit(std::vector<unsigned int> indexs, std::vector<bool> values)
{
    // Returns a list of terms that satisfy the proposition with the given bit assignment.
    return this->cp.satisfyBit(indexs, values);
}
std::vector<std::string> Location::unsatisfyBit(std::vector<unsigned int> indexs, std::vector<bool> values)
{
    // Returns a list of terms that do not satisfy the proposition with the given bit assignment.
    return this->cp.unsatisfyBit(indexs, values);
}

void initializeTransitionSystem()
{
    CFLOBDDNodeHandle::InitNoDistinctionTable();
    CFLOBDDNodeHandle::InitAdditionInterleavedTable();
    CFLOBDDNodeHandle::InitReduceCache();
    InitPairProductCache();
    InitTripleProductCache();
    Matrix1234ComplexFloatBoost::Matrix1234Initializer();
    VectorComplexFloatBoost::VectorInitializer();
}

/***
 * Transition system for model checking algorithms.
 ***/

class TransitionSystem
{
private:
    std::unordered_map<int, int> locationBuf;
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
    TransitionSystem(bool init) {
        if (init) {
            CFLOBDDNodeHandle::InitNoDistinctionTable();
            CFLOBDDNodeHandle::InitAdditionInterleavedTable();
            CFLOBDDNodeHandle::InitReduceCache();
            InitPairProductCache();
            InitTripleProductCache();
            Matrix1234ComplexFloatBoost::Matrix1234Initializer();
            VectorComplexFloatBoost::VectorInitializer();
        }
    };
    ~TransitionSystem();
    // void initialization();
    void addLocation(Location loc);
    void addRelation(unsigned int from, unsigned int to, QOperation op);
    void setAnnotation(std::vector<std::tuple<unsigned int, QOperation>> annotations);
    void resetLocationBounds();
    void preConditionInit();
    void preConditionOneStep(unsigned int loc);
    void preConditions();
    void postConditionInit();
    void postConditionOneStep(unsigned int loc);
    void postConditions();
    void createAdd();
    bool satisfy(unsigned int loc, QOperation spec);
    void setInitLocation(unsigned int loc);
    unsigned int getInitLocation() const {
        return this->initLocation;
    }
    void computingFixedPointPre();
    void computingFixedPointPost();
    std::pair<int, int> printDims(unsigned int loc);
    void printSupp(unsigned int loc);
    unsigned int getLocationNum() const {
        return static_cast<unsigned int>(this->Locations.size());
    }
    std::string getRelationName(unsigned int from, unsigned int to) const {
        auto it = this->relations.find(std::make_tuple(from, to));
        if (it != this->relations.end()) {
            return it->second.getName();
        }
        return "";
    }
    void setLabel(unsigned int loc, std::string label) {
        // Set the label for the location
        this->Locations[loc].setLabel(label);
    }
    std::vector<std::string> getLabels(unsigned int loc) const {
        // Get the labels for the location
        return this->Locations[loc].getLabels();
    }
    bool isLeafLoc(unsigned int loc) const {
        // Check if the location is a leaf location (no post-locations)
        return this->Locations[loc].postLocations.empty() || 
               (this->Locations[loc].postLocations.size() == 1 && this->Locations[loc].postLocations[0] == loc);
    }
};

void TransitionSystem::addLocation(Location loc)
{
    this->Locations.push_back(loc);
    this->Locations.back().idx = static_cast<unsigned int>(this->Locations.size()) - 1;
}

void TransitionSystem::addRelation(unsigned int from, unsigned int to, QOperation op)
{
    // The problem is how to represent a projective operation's support vectors.
    this->relations[std::make_tuple(from, to)] = op;
    // The bi-directional relation is established here.
    this->Locations[from].appendPostLocation(to);
    this->Locations[to].appendPreLocation(from);
    // std::cout << "Add relation from " << &this->Locations[from] << " " << this->Locations[from].idx << " to " << &this->Locations[to] << " " << this->Locations[to].idx << " with operation\n";
    // this->Locations[to].appendPreLocation(&this->Locations[from]);
    // The operation is a projective operation, and from has two post-locations.
    // Haven't finished yet.
    if (op.isProj >= 0) {
        // Assert that the operation is a binary projective operation or a resetting.
        assert(op.oplist.size() <= 2);
        assert(op.oplist[0]->getType() == true);
        auto name = dynamic_cast<QuantumGateTerm*>(op.oplist[0].get())->name;
        assert(name == "meas0" || name == "meas1");
        this->Locations[from].flag = 0; // Set the flag to 0 for projective operations
    }
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

void TransitionSystem::resetLocationBounds()
{
    for (auto& loc : this->Locations) {
        loc.resetBounds();
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
   // Initialize the visitedPre vectors
    this->visitedPre.resize(this->Locations.size(), false);
    this->computedTablePre.clear(); // Clear the computed table for pre-conditions
    for (unsigned int i = 0; i < this->Locations.size(); i++) {
        if ((this->Locations[i].postLocations.size() == 1 && this->Locations[i].postLocations[0] == i) || this->Locations[i].postLocations.size() == 0) {
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
        unsigned int preLocIdx = this->Locations[loc].preLocations[i];
        Location* preLoc = this->Locations.data() + preLocIdx; // Get the pointer to the preLocation
        assert(preLoc->idx == preLocIdx); // Ensure the index is correct
        // If it is a self-loop and the relation is identity, skip it without appending currPreLoc.
        if (preLoc->idx == loc && this->relations[std::make_tuple(preLoc->idx, loc)].isIdentity) {
            std::cout << "Skip self-loop for location " << loc << std::endl;
            continue;
        }
        if (this->computedTablePre.find(std::make_tuple(loc, this->Locations[loc].upperBound.oplist.size(), preLoc->idx, preLoc->upperBound.oplist.size())) == this->computedTablePre.end()) {
            std::cout << this->relations[std::make_tuple(preLoc->idx, loc)].oplist.size() << " operations in relation from " << loc << " to " << preLoc->idx << std::endl;
            // std::cout << this->relations[std::make_tuple(preLoc->idx, loc)].type << std::endl;
            int dimBefore;
            if (preLoc->flag < 0) {
                QOperation preImage = this->Locations[loc].upperBound.preImage(this->relations[std::make_tuple(preLoc->idx, loc)]);
                computedTablePre.insert(std::make_tuple(loc, this->Locations[loc].upperBound.oplist.size(), preLoc->idx, preLoc->upperBound.oplist.size()));
                dimBefore = preLoc->upperBound.oplist.size(); // Upper bound: dimension 2^n as default
                preLoc->upperBound = preLoc->upperBound.conjunction_simp(preImage); // TODO: Conjunction inline
            } else if (preLoc->flag == 0) {
                std::cout << "Flag is 0 for location " << preLoc->idx << std::endl;
                QOperation tempConjunction = this->Locations[loc].upperBound.conjunction_simp(this->relations[std::make_tuple(preLoc->idx, loc)]);
                preLoc->tempOperations.push_back(tempConjunction);
                preLoc->flag = 1;
                this->locationBuf[preLoc->idx] = loc;
                continue; // Skip the conjunction for flag == 0, don't put a new location into list, until the other branch reaches here.
            } else if (preLoc->flag == 1) {
                // The other branch has reached here, we can do the conjunction.
                std::cout << "Flag is 1 for location " << preLoc->idx << std::endl;
                QOperation tempConjunction2 = this->Locations[loc].upperBound.conjunction_simp(this->relations[std::make_tuple(preLoc->idx, loc)]);
                QOperation tempDisjunction = preLoc->tempOperations.back().disjunction(tempConjunction2);
                preLoc->tempOperations.pop_back(); // Remove the last operation
                preLoc->upperBound = preLoc->upperBound.conjunction_simp(tempDisjunction);
                preLoc->flag = 0; // Set the flag to 0, indicating all branches have reached here.
                if (this->locationBuf.count(preLoc->idx)) {
                    this->locationBuf.erase(preLoc->idx);
                }
            } else if (preLoc->flag == -2) {
                
            }
            if (visitedPre[preLoc->idx] == false) {
                std::cout << "Visit a new pre location " << preLoc->idx << std::endl;
                this->currPreLocs.push_back(preLoc->idx);
                visitedPre[preLoc->idx] = true;
            } else if (preLoc->upperBound.oplist.size() < dimBefore && dimBefore > 0) {
                // If the dimension of the upperBound is reduced, we need to recheck the pre-condition.
                // TODO: Check it carefully!!
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
    while (true) {
        if (!this->currPreLocs.empty()) {
            unsigned int loc = this->currPreLocs.front();
            this->currPreLocs.pop_front();
            this->preConditionOneStep(loc);
        } else if (!this->locationBuf.empty()) {
            // If there are locations in locationBuf, we need to process them.
            // Arbitrarily select a location from locationBuf.
            std::cout << "Processing location from locationBuf." << std::endl;
            auto it = this->locationBuf.begin();
            int loc = it->first;
            int postLoc = it->second;
            this->locationBuf.erase(it);
            this->Locations[loc].flag = -1; // TODO: Here is a tricky part, we view this location as a normal one and never try a special treatment.
            assert(this->relations[std::make_tuple(loc, postLoc)].isProj >= 0);
            QOperation negOther = this->relations[std::make_tuple(loc, postLoc)].negation();
            negOther.genProjMeasSpace();
            QOperation preImage = this->Locations[loc].tempOperations.back().disjunction(negOther);
            this->Locations[loc].tempOperations.pop_back();
            computedTablePre.insert(std::make_tuple(postLoc, this->Locations[postLoc].upperBound.oplist.size(), loc, this->Locations[loc].upperBound.oplist.size()));
            this->Locations[loc].upperBound = this->Locations[loc].upperBound.conjunction_simp(preImage); // TODO: Conjunction inline
            // We are sure the dimension has changed here, so don't need to check it.
            this->visitedPre[loc] = true;
            this->preConditionOneStep(loc);
        } else {
            // No more locations to process
            break;
        }
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
        if ((this->Locations[i].preLocations.size() == 1 && this->Locations[i].preLocations[0] == i) || this->Locations[i].preLocations.size() == 0) {
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
    // std::cout << this->Locations[loc].postLocations.size() << " post locations for location " << loc << std::endl;
    for (unsigned int i = 0; i < this->Locations[loc].postLocations.size(); i++) {
        unsigned int postLocIdx = this->Locations[loc].postLocations[i];
        Location* postLoc = this->Locations.data() + postLocIdx; // Get the pointer to the postLocation
        if (postLoc->idx != postLocIdx) {
            std::cout << "Warning: postLoc->idx (" << postLoc->idx << ") != postLocIdx (" << postLocIdx << ")" << std::endl;
        }
        assert(postLoc->idx == postLocIdx); // Ensure the index is correct
        // If it is a self-loop and the relation is identity, skip it without appending currPostLocs.
        if (postLoc->idx == loc && this->relations[std::make_tuple(loc, postLoc->idx)].isIdentity) {
            std::cout << "Skip self-loop for location " << loc << std::endl;
            continue;
        }
        // If (loc, locDim, postLoc, postLocDim) is not computed, compute it.
        if (this->computedTablePost.find(std::make_tuple(loc, this->Locations[loc].lowerBound.oplist.size(), postLoc->idx, postLoc->lowerBound.oplist.size())) == this->computedTablePost.end()) {
            // std::cout << this->relations[std::make_tuple(loc, postLoc->idx)].oplist.size() << " operations in relation from " << loc << " to " << postLoc->idx << std::endl;
            QOperation postImage = this->Locations[loc].lowerBound.postImage(this->relations[std::make_tuple(loc, postLoc->idx)]);
            computedTablePost.insert(std::make_tuple(loc, this->Locations[loc].lowerBound.oplist.size(), postLoc->idx, postLoc->lowerBound.oplist.size()));
            int dimBefore = postLoc->lowerBound.oplist.size(); // Lower bound: dimension 0 as default
            postLoc->lowerBound = postLoc->lowerBound.disjunction(postImage); // TODO: Disjunction inline
            if (visitedPost[postLoc->idx] == false) {
                // std::cout << "Visit a new post location " << postLoc->idx << std::endl;
                this->currPostLocs.push_back(postLoc->idx);
                visitedPost[postLoc->idx] = true;
            } else if (postLoc->lowerBound.oplist.size() > dimBefore && dimBefore < std::pow(2, postLoc->lowerBound.qNum)) {
                // If the dimension of the lowerBound is reduced, we need to recheck the post-condition.
                // If postLoc->idx is not in currPostLocs, append it.
                if (std::find(this->currPostLocs.begin(), this->currPostLocs.end(), postLoc->idx) == this->currPostLocs.end()) {
                    // std::cout << "Post condition for location " << postLoc->idx << " is updated from " << dimBefore << " to " << postLoc->lowerBound.oplist.size() << std::endl;
                    this->currPostLocs.push_back(postLoc->idx);
                }
            } else {
                // std::cout << "Post condition for location " << postLoc->idx << " is not updated." << std::endl;
            }
        } else {}
    }
    // std::cout << "Post condition for location " << loc << " computed." << std::endl;
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

void TransitionSystem::computingFixedPointPre() {
    this->preConditionInit();
    this->preConditions();
}

void TransitionSystem::computingFixedPointPost() {
    this->postConditionInit();
    this->postConditions();
}

bool TransitionSystem::satisfy(unsigned int loc, QOperation spec) {
    /*
    Check if the location loc satisfies the specification spec.
    The specification is a projective operation.
    */
    assert(spec.isProj < 0); // Only projective operations are supported for now
    return this->Locations[loc].satisfy(spec);
}

std::pair<int, int> TransitionSystem::printDims(unsigned int loc) {
    /*
    Print the dimensions of the upperBound and lowerBound of the location.
    */
    // std::cout << "Location " << loc << ": upperBound dimension = " << (this->Locations[loc].upperBound.isIdentity ? (std::pow(2, this->Locations[loc].upperBound.qNum)) : this->Locations[loc].upperBound.oplist.size()) 
    //           << ", lowerBound dimension = " << this->Locations[loc].lowerBound.oplist.size() << std::endl;
    int upperDim = this->Locations[loc].upperBound.isIdentity ? (std::pow(2, this->Locations[loc].upperBound.qNum)) : this->Locations[loc].upperBound.oplist.size();
    // exclude the zero operators!!
    int lowerDim = 0;
    for (const auto& op : this->Locations[loc].lowerBound.oplist) {
        if (!checkifzero(op->content)) {
            lowerDim++;
        }
    }
    return std::make_pair(upperDim, lowerDim);
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

TransitionSystem fromProgramToTransitionSystem(const std::string& filename, TransitionSystem& ts) {
    /*
    Read the program from the file and convert it to a transition system.
    The program is in the form of a list of locations and transitions.
    Each location has an annotation and a set of pre- and post-conditions.
    */
   // Need to write a parser!
    QASMProgram program = parse_qasm_file(filename);
    // Confirm the number of locations
    int numLocs = program.operations.size();
   // Confirm the number and type of relations
}

/***
 * Location: The location of a transition system.
 * It contains the state and the annotation of the location.
 ***/

#endif
