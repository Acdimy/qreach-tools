#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <vector>
#include "../quantum_operation.hpp"
#include "transition_system.hpp"
#include "transition_system_qadd.hpp"

namespace py = pybind11;

PYBIND11_MODULE(pyqreach, m) {
    m.doc() = "python wrapper for Quantum Simulation"; // Optional module docstring

    py::class_<QOperation>(m, "QOperation")
        .def(py::init<>())
        .def(py::init<std::vector<std::string>>())
        .def(py::init<std::vector<double>, unsigned int>())
        .def(py::init<std::string, unsigned int, std::vector<unsigned int>, std::vector<double>>())
        .def_readonly("type", &QOperation::type)
        .def_readonly("normalized", &QOperation::normalized)
        .def_readonly("qNum", &QOperation::qNum)
        .def_readonly("isIdentity", &QOperation::isIdentity)
        .def_readonly("isProj", &QOperation::isProj);
    
    m.def("CreateIdentityQO", &CreateIdentityQO, "Create an identity quantum operation");
    m.def("CreateZeroQO", &CreateZeroQO, "Create a zero quantum operation");

    py::class_<ClassicalProposition>(m, "ClassicalProposition")
        .def(py::init<>())
        .def("toString", &ClassicalProposition::toString, "toString")
        .def("print", &ClassicalProposition::print, "print")
        .def_readonly("terms", &ClassicalProposition::terms);

    py::class_<qts_naive::Location>(m, "Location")
        .def(py::init<>())
        .def(py::init<int>())
        .def(py::init<int, unsigned int>())
        .def_readonly("qNum", &qts_naive::Location::qNum)
        .def("appendPreLocation", &qts_naive::Location::appendPreLocation, "appendPreLocation")
        .def("appendPostLocation", &qts_naive::Location::appendPostLocation, "appendPostLocation")
        .def("appendClassicalAP", &qts_naive::Location::appendClassicalAP, "appendClassicalAP")
        .def("copyClassicalAP", &qts_naive::Location::copyClassicalAP, "copyClassicalAP")
        .def("satisfyBit", &qts_naive::Location::satisfyBit, "satisfyBit")
        .def("unsatisfyBit", &qts_naive::Location::unsatisfyBit, "unsatisfyBit")
        .def("termNum", &qts_naive::Location::termNum, "termNum")
        .def("equalAP", &qts_naive::Location::equalAP, "equalAP")
        .def("setClassicalValue", &qts_naive::Location::setClassicalValue, "setClassicalValue")
        .def("find", &qts_naive::Location::find, "find")
        .def("satisfy", &qts_naive::Location::satisfy, "satisfy with QOperation")
        .def("satisfyDefault", &qts_naive::Location::satisfyDefault, "satisfy without argument, check if lowerBound <= upperBound")
        .def("setIdentifier", &qts_naive::Location::setIdentifier, "setIdentifier")
        .def("getIdentifier", &qts_naive::Location::getIdentifier, "getIdentifier")
        .def_readwrite("idx", &qts_naive::Location::idx)
        .def_readwrite("flag", &qts_naive::Location::flag)
        .def_readwrite("upperBound", &qts_naive::Location::upperBound)
        .def_readwrite("lowerBound", &qts_naive::Location::lowerBound)
        .def_readwrite("cp", &qts_naive::Location::cp)
        .def_readonly("postLocations", &qts_naive::Location::postLocations);

    m.def("initializeTransitionSystem", &qts_naive::initializeTransitionSystem, "initializeTransitionSystem");

    py::class_<qts_naive::TransitionSystem>(m, "TransitionSystem")
        .def(py::init<>())
        .def(py::init<bool>())
        .def("addLocation", &qts_naive::TransitionSystem::addLocation, "addLocation")
        .def("addRelation", &qts_naive::TransitionSystem::addRelation, "addRelation")
        .def("setInitLocation", &qts_naive::TransitionSystem::setInitLocation, "setInitLocation")
        .def("getInitLocation", &qts_naive::TransitionSystem::getInitLocation, "getInitLocation")
        .def("setAnnotation", &qts_naive::TransitionSystem::setAnnotation, "setAnnotation")
        .def("resetLocationBounds", &qts_naive::TransitionSystem::resetLocationBounds, "resetLocationBounds")
        .def("preConditionInit", &qts_naive::TransitionSystem::preConditionInit, "preConditionInit")
        .def("preConditionOneStep", &qts_naive::TransitionSystem::preConditionOneStep, "preConditionOneStep")
        .def("preConditions", &qts_naive::TransitionSystem::preConditions, "preConditions")
        .def("postConditionInit", &qts_naive::TransitionSystem::postConditionInit, "postConditionInit")
        .def("postConditionOneStep", &qts_naive::TransitionSystem::postConditionOneStep, "postConditionOneStep")
        .def("postConditions", &qts_naive::TransitionSystem::postConditions, "postConditions")
        .def("computingFixedPointPre", &qts_naive::TransitionSystem::computingFixedPointPre, "computingFixedPointPre")
        .def("computingFixedPointPost", &qts_naive::TransitionSystem::computingFixedPointPost, "computingFixedPointPost")
        .def("satisfy", &qts_naive::TransitionSystem::satisfy, "satisfy")
        .def("getLocationNum", &qts_naive::TransitionSystem::getLocationNum, "getLocationNum")
        .def("printDims", &qts_naive::TransitionSystem::printDims, "printDims")
        .def("printSupp", &qts_naive::TransitionSystem::printSupp, "printSupp")
        .def("getRelationName", &qts_naive::TransitionSystem::getRelationName, "getRelationName")
        .def("setLabel", &qts_naive::TransitionSystem::setLabel, "setLabel")
        .def("getLabels", &qts_naive::TransitionSystem::getLabels, "getLabels")
        .def("isLeafLoc", &qts_naive::TransitionSystem::isLeafLoc, "isLeafLoc")
        .def_readonly("relations", &qts_naive::TransitionSystem::relations)
        .def_readonly("Locations", &qts_naive::TransitionSystem::Locations);

    py::class_<qts::TransitionSystem>(m, "SymTS")
        .def(py::init<>())
        .def(py::init<int, int>())
        .def("addLocation", &qts::TransitionSystem::addLocation, "addLocation")
        .def("addRelation", &qts::TransitionSystem::addRelation, "addRelation")
        .def("setAnnotation", &qts::TransitionSystem::setAnnotation, "setAnnotation")
        .def("postConditions", &qts::TransitionSystem::postConditions, "postConditions");
    
}
