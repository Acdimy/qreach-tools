#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <vector>
#include <unordered_map>
#include "../quantum_operation.hpp"
#include "transition_system.hpp"
#include "transition_system_qadd.hpp"

namespace py = pybind11;

namespace {

py::dict export_qadd_graph(qadd::QADDNode* root) {
    py::list nodes;
    py::list edges;
    std::unordered_map<qadd::QADDNode*, int> node_ids;
    std::vector<qadd::QADDNode*> stack;

    if (!root) {
        py::dict result;
        result["nodes"] = nodes;
        result["edges"] = edges;
        return result;
    }

    stack.push_back(root);
    int next_id = 0;
    while (!stack.empty()) {
        qadd::QADDNode* node = stack.back();
        stack.pop_back();

        if (node_ids.find(node) != node_ids.end()) {
            continue;
        }

        int node_id = next_id++;
        node_ids[node] = node_id;

        py::dict node_info;
        node_info["id"] = node_id;
        if (qadd::is_terminal(node)) {
            node_info["kind"] = "terminal";
            node_info["label"] = node->val.getName();
            node_info["is_operation"] = node->val.type;
            node_info["is_identity"] = node->val.isIdentity;
            node_info["dimension"] = node->val.type ? -1 : qadd::get_dimension(node->val);
        } else {
            node_info["kind"] = "internal";
            node_info["label"] = "X" + std::to_string(node->var);
            node_info["var"] = node->var;
            stack.push_back(node->low);
            stack.push_back(node->high);
        }
        nodes.append(node_info);
    }

    for (const auto& entry : node_ids) {
        qadd::QADDNode* node = entry.first;
        int src_id = entry.second;
        if (qadd::is_terminal(node)) {
            continue;
        }

        py::dict low_edge;
        low_edge["source"] = src_id;
        low_edge["target"] = node_ids.at(node->low);
        low_edge["branch"] = 0;
        edges.append(low_edge);

        py::dict high_edge;
        high_edge["source"] = src_id;
        high_edge["target"] = node_ids.at(node->high);
        high_edge["branch"] = 1;
        edges.append(high_edge);
    }

    py::dict result;
    result["nodes"] = nodes;
    result["edges"] = edges;
    result["node_count"] = static_cast<int>(py::len(nodes));
    result["edge_count"] = static_cast<int>(py::len(edges));
    return result;
}

} // namespace

PYBIND11_MODULE(pyqreach, m) {
    m.doc() = "python wrapper for Quantum Simulation"; // Optional module docstring

    py::class_<QOperation>(m, "QOperation")
        .def(py::init<>())
        .def(py::init<std::vector<std::string>>())
        .def(py::init<std::vector<double>, unsigned int>())
        .def(py::init<std::string, unsigned int, std::vector<unsigned int>, std::vector<double>>())
        .def("getName", &QOperation::getName, "getName")
        .def("printFormal", &QOperation::printFormal, py::arg("print") = true, "printFormal")
        .def("disjunction", &QOperation::disjunction, py::arg("other"), "Return the span/disjunction of two subspace QOperations")
        .def_readonly("type", &QOperation::type)
        .def_readonly("normalized", &QOperation::normalized)
        .def_readonly("qNum", &QOperation::qNum)
        .def_readonly("isIdentity", &QOperation::isIdentity)
        .def_readonly("isProj", &QOperation::isProj);
    
    m.def("CreateIdentityQO", &CreateIdentityQO, "Create an identity quantum operation");
    m.def("CreateZeroQO", &CreateZeroQO, "Create a zero quantum operation");
    m.def("span_qops", &SpanQOperations, py::arg("ops"), "Construct the normalized span of subspace QOperations");

    py::class_<ClassicalProposition>(m, "ClassicalProposition")
        .def(py::init<>())
        .def(py::init<unsigned int>())
        .def(py::init<unsigned int, std::vector<std::string>>())
        .def("addTerm", &ClassicalProposition::addTerm, "addTerm")
        .def("removeTerm", &ClassicalProposition::removeTerm, "removeTerm")
        .def("setValue", &ClassicalProposition::setValue, "setValue")
        .def("find", &ClassicalProposition::find, "find")
        .def("satisfyBit", &ClassicalProposition::satisfyBit, "satisfyBit")
        .def("unsatisfyBit", &ClassicalProposition::unsatisfyBit, "unsatisfyBit")
        .def("toString", &ClassicalProposition::toString, "toString")
        .def("print", &ClassicalProposition::print, "print")
        .def("termNum", [](const ClassicalProposition& cp) { return static_cast<int>(cp.terms.size()); }, "termNum")
        .def_readonly("n", &ClassicalProposition::n)
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
    m.def("initializeSymTransitionSystem", &qts::initializeTransitionSystem, "initializeSymTransitionSystem");

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
           .def(py::init<int, int>(), py::arg("num_qubits"), py::arg("max_locations") = 0)
              .def("addLocation",
                    [](qts::TransitionSystem& ts,
                      const ClassicalProposition& cp,
                      const std::string& identifier) {
                          return ts.addLocation(cp, identifier);
                    },
                    py::arg("cp") = ClassicalProposition(),
                    py::arg("identifier") = "")
        .def("addRelation", &qts::TransitionSystem::addRelation, "addRelation")
           .def("setAnnotation",
               py::overload_cast<int, const QOperation&>(&qts::TransitionSystem::setAnnotation),
               "setAnnotation")
              .def("appendClassicalAP", &qts::TransitionSystem::appendClassicalAP, "appendClassicalAP")
              .def("setClassicalValue", &qts::TransitionSystem::setClassicalValue, "setClassicalValue")
              .def("find", &qts::TransitionSystem::find, "find")
              .def("satisfyBit", &qts::TransitionSystem::satisfyBit, "satisfyBit")
              .def("unsatisfyBit", &qts::TransitionSystem::unsatisfyBit, "unsatisfyBit")
              .def("termNum", &qts::TransitionSystem::termNum, "termNum")
              .def("getClassicalProposition", &qts::TransitionSystem::getClassicalProposition, "getClassicalProposition")
              .def("setClassicalProposition", &qts::TransitionSystem::setClassicalProposition, "setClassicalProposition")
              .def("setIdentifier", &qts::TransitionSystem::setIdentifier, "setIdentifier")
              .def("getIdentifier", &qts::TransitionSystem::getIdentifier, "getIdentifier")
              .def("setLabel", &qts::TransitionSystem::setLabel, "setLabel")
              .def("getLabels", &qts::TransitionSystem::getLabels, "getLabels")
              .def("setInitLocation", &qts::TransitionSystem::setInitLocation, "setInitLocation")
              .def("getInitLocation", &qts::TransitionSystem::getInitLocation, "getInitLocation")
              .def("getLocationNum", &qts::TransitionSystem::getNumLocations, "getLocationNum")
              .def("getLocationIDs", &qts::TransitionSystem::getLocationIDs, "getLocationIDs")
              .def("getPostLocations", &qts::TransitionSystem::getPostLocations, "getPostLocations")
              .def("getRelationName", &qts::TransitionSystem::getRelationName, "getRelationName")
              .def("getLocationAnnotation", &qts::TransitionSystem::getLocationAnnotation, "getLocationAnnotation")
              .def("getLocationDimension", &qts::TransitionSystem::getLocationDimension, "getLocationDimension")
              .def("locationHasNonZeroAnnotation", &qts::TransitionSystem::locationHasNonZeroAnnotation, "locationHasNonZeroAnnotation")
              .def("filterReachableLocations", &qts::TransitionSystem::filterReachableLocations, "filterReachableLocations")
              .def("printDims", &qts::TransitionSystem::printDims, "printDims")
              .def("satisfy", &qts::TransitionSystem::satisfy, "satisfy")
              .def("isLeafLoc", &qts::TransitionSystem::isLeafLoc, "isLeafLoc")
              .def("getAnnotationNodeCount", &qts::TransitionSystem::getAnnotationNodeCount, "getAnnotationNodeCount")
              .def("getRelationNodeCount", &qts::TransitionSystem::getRelationNodeCount, "getRelationNodeCount")
              .def("getTotalUniqueNodeCount", &qts::TransitionSystem::getTotalUniqueNodeCount, "getTotalUniqueNodeCount")
              .def("exportAnnotationQADD",
                   [](const qts::TransitionSystem& ts) {
                       return export_qadd_graph(ts.getAnnotation());
                   },
                   "exportAnnotationQADD")
              .def("exportRelationQADD",
                   [](const qts::TransitionSystem& ts) {
                       return export_qadd_graph(ts.getRelation());
                   },
                   "exportRelationQADD")
           .def("postConditions", &qts::TransitionSystem::postConditions, "postConditions")
           .def("computingFixedPointPost", &qts::TransitionSystem::postConditions, "computingFixedPointPost");
    
}
