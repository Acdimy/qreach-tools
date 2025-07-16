#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <vector>
#include "../quantum_operation.hpp"
#include "transition_system.hpp"

namespace py = pybind11;

PYBIND11_MODULE(pyqreach, m) {
    m.doc() = "python wrapper for Quantum Simulation"; // Optional module docstring

    py::class_<QOperation>(m, "QOperation")
        .def(py::init<>())
        .def(py::init<std::vector<std::string>>())
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

    py::class_<Location>(m, "Location")
        .def(py::init<>())
        .def(py::init<int>())
        .def(py::init<int, unsigned int>())
        .def_readonly("qNum", &Location::qNum)
        .def("appendPreLocation", &Location::appendPreLocation, "appendPreLocation")
        .def("appendPostLocation", &Location::appendPostLocation, "appendPostLocation")
        .def("appendClassicalAP", &Location::appendClassicalAP, "appendClassicalAP")
        .def("copyClassicalAP", &Location::copyClassicalAP, "copyClassicalAP")
        .def("satisfyBit", &Location::satisfyBit, "satisfyBit")
        .def("unsatisfyBit", &Location::unsatisfyBit, "unsatisfyBit")
        .def("termNum", &Location::termNum, "termNum")
        .def("equalAP", &Location::equalAP, "equalAP")
        .def("setClassicalValue", &Location::setClassicalValue, "setClassicalValue")
        .def("find", &Location::find, "find")
        // .def("satisfy", py::overload_cast<QOperation>(&Location::satisfy), "satisfy with QOperation")
        .def_readwrite("idx", &Location::idx)
        .def_readwrite("flag", &Location::flag)
        .def_readwrite("upperBound", &Location::upperBound)
        .def_readwrite("lowerBound", &Location::lowerBound)
        .def_readwrite("cp", &Location::cp)
        .def_readonly("postLocations", &Location::postLocations);

    py::class_<TransitionSystem>(m, "TransitionSystem")
        .def(py::init<>())
        .def("addLocation", &TransitionSystem::addLocation, "addLocation")
        .def("addRelation", &TransitionSystem::addRelation, "addRelation")
        .def("setInitLocation", &TransitionSystem::setInitLocation, "setInitLocation")
        .def("setAnnotation", &TransitionSystem::setAnnotation, "setAnnotation")
        .def("preConditionInit", &TransitionSystem::preConditionInit, "preConditionInit")
        .def("preConditionOneStep", &TransitionSystem::preConditionOneStep, "preConditionOneStep")
        .def("preConditions", &TransitionSystem::preConditions, "preConditions")
        .def("postConditionInit", &TransitionSystem::postConditionInit, "postConditionInit")
        .def("postConditionOneStep", &TransitionSystem::postConditionOneStep, "postConditionOneStep")
        .def("postConditions", &TransitionSystem::postConditions, "postConditions")
        .def("computingFixedPointPre", &TransitionSystem::computingFixedPointPre, "computingFixedPointPre")
        .def("computingFixedPointPost", &TransitionSystem::computingFixedPointPost, "computingFixedPointPost")
        .def("satisfy", &TransitionSystem::satisfy, "satisfy")
        .def("getLocationNum", &TransitionSystem::getLocationNum, "getLocationNum")
        .def("printDims", &TransitionSystem::printDims, "printDims")
        .def("printSupp", &TransitionSystem::printSupp, "printSupp")
        .def("getRelationName", &TransitionSystem::getRelationName, "getRelationName")
        .def_readonly("Locations", &TransitionSystem::Locations);
    
}
