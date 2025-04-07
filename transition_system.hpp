#include "quantum_operation.hpp"

class Location
{
public:
    QOperation state;
    QOperation annotation;
    std::vector<Location*> preLocations;
    std::vector<Location*> postLocations;
public:
    Location(/* args */);
    ~Location();
    void setAnnotation();
};

Location::Location(/* args */)
{
}

Location::~Location()
{
}

/***
 * Doesn't consider the representation of a CTL formula now for simplicity
 ***/

class TransitionSystem
{
public:
    std::vector<Location> Locations;
    unsigned int initLocation;
    std::vector<std::tuple<unsigned int, unsigned int, std::string>> relations;
public:
    TransitionSystem(/* args */);
    ~TransitionSystem();
    void initialization();
    void setAnnotation();
    void preConditions();
    void postConditions();
    bool satisfy();
};

TransitionSystem::TransitionSystem(/* args */)
{
}

TransitionSystem::~TransitionSystem()
{
}

