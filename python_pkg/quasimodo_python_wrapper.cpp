#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <vector>
#include "../quantum_circuit.h"

namespace py = pybind11;

PYBIND11_MODULE(pyquasimodo, m) {
    m.doc() = "python wrapper for Quantum Simulation"; // Optional module docstring

    py::class_<QuantumCircuit>(m, "QuantumCircuit");
        // .def(py::init<>())
        // .def(py::init<unsigned int>());
    
    py::class_<CFLOBDDQuantumCircuit, QuantumCircuit>(m, "CFLOBDDQuantumCircuit")
        .def(py::init<>())
        .def(py::init<unsigned int, int>())
        .def("setNumQubits", &CFLOBDDQuantumCircuit::setNumQubits, "setNumQubits")
        .def("i", &CFLOBDDQuantumCircuit::ApplyIdentityGate, "ApplyIdentityGate")
        .def("h", &CFLOBDDQuantumCircuit::ApplyHadamardGate, "ApplyHadamardGate")
        .def("x", &CFLOBDDQuantumCircuit::ApplyNOTGate, "ApplyNOTGate")
        .def("swap", &CFLOBDDQuantumCircuit::ApplySwapGate, "ApplySwapGate")
        .def("iswap", &CFLOBDDQuantumCircuit::ApplyiSwapGate, "ApplyiSwapGate")
        .def("prob", &CFLOBDDQuantumCircuit::GetProbability, "GetProbability")
        .def("measure", &CFLOBDDQuantumCircuit::Measure, "Measure")
        .def("measurement_counts", &CFLOBDDQuantumCircuit::GetPathCount, "MeasurementCount")
        .def("y", &CFLOBDDQuantumCircuit::ApplyPauliYGate, "ApplyPauliYGate")
        .def("z", &CFLOBDDQuantumCircuit::ApplyPauliZGate, "ApplyPauliZGate")
        .def("p", &CFLOBDDQuantumCircuit::ApplyPhaseShiftGate, "ApplyPhaseShiftGate")
        .def("s", &CFLOBDDQuantumCircuit::ApplySGate, "ApplySGate")
        .def("t", &CFLOBDDQuantumCircuit::ApplyTGate, "ApplyTGate") 
        .def("cz", &CFLOBDDQuantumCircuit::ApplyCZGate, "ApplyCZGate")
        .def("cx", &CFLOBDDQuantumCircuit::ApplyCNOTGate, "ApplyCNOTGate")
        .def("ccx", &CFLOBDDQuantumCircuit::ApplyCCNOTGate, "ApplyCCNOTGate")
        .def("gp", &CFLOBDDQuantumCircuit::ApplyGlobalPhase, "ApplyGlobalPhase")
        .def("cp", &CFLOBDDQuantumCircuit::ApplyCPGate, "ApplyCPGate")
        .def("cs", &CFLOBDDQuantumCircuit::ApplyCSGate, "ApplyCSGate")
        .def("cswap", &CFLOBDDQuantumCircuit::ApplyCSwapGate, "ApplyCSwapGate")
        .def("u3", &CFLOBDDQuantumCircuit::ApplyU3Gate, "ApplyU3Gate")
        .def("u", &CFLOBDDQuantumCircuit::ApplyArbitraryGate, "ApplyArbitraryGate")
        .def("setState", &CFLOBDDQuantumCircuit::setBasicStateVector, "setBasicStateVector")
        .def("setInitGate", &CFLOBDDQuantumCircuit::setInitGate, "setInitGate")
        .def("setProjector", &CFLOBDDQuantumCircuit::setProjector, "setProjector")
        .def("setProjectorFS", &CFLOBDDQuantumCircuit::setProjectorFromS, "setProjectorFS")
        .def("ApplyProjToEnt", &CFLOBDDQuantumCircuit::ApplyProjectorToEntangle, "ApplyProjectorToEntangle")
        .def("pushState2Cache", &CFLOBDDQuantumCircuit::pushStateToCache, "pushStateToCache")
        .def("popCache2State", &CFLOBDDQuantumCircuit::addCacheToState, "popCache2State")
        .def("appendGateSeries", &CFLOBDDQuantumCircuit::appendGateSeries, "appendGateSeries")
        .def("ApplyGateSeries", &CFLOBDDQuantumCircuit::ApplyGateSeries, "ApplyGateSeries")
        .def("reachability", &CFLOBDDQuantumCircuit::reachability, "reachability")
        .def("test", &CFLOBDDQuantumCircuit::test, "test")
        .def("print", &CFLOBDDQuantumCircuit::print, "print")
        .def("printRV", &CFLOBDDQuantumCircuit::printRV, "printRV")
        .def("printSize", &CFLOBDDQuantumCircuit::printSize, "printSize")
        .def("printYield", &CFLOBDDQuantumCircuit::printYield, "printYield")
        .def("printCol", &CFLOBDDQuantumCircuit::printStateColMajor, "printStateColMajor")
        .def("printColHead", &CFLOBDDQuantumCircuit::printStateColHead, "printStateColHead")
        .def("printProjector", &CFLOBDDQuantumCircuit::printProjector, "printProjector");
        // .def("u2", &CFLOBDDQuantumCircuit::ApplyU2Gate, "ApplyU2Gate");

    
    py::class_<WeightedBDDQuantumCircuit, QuantumCircuit>(m, "WeightedBDDQuantumCircuit")
        .def(py::init<>())
        .def(py::init<unsigned int, int>())
        .def("setNumQubits", &WeightedBDDQuantumCircuit::setNumQubits, "setNumQubits")
        .def("i", &WeightedBDDQuantumCircuit::ApplyIdentityGate, "ApplyIdentityGate")
        .def("h", &WeightedBDDQuantumCircuit::ApplyHadamardGate, "ApplyHadamardGate")
        .def("x", &WeightedBDDQuantumCircuit::ApplyNOTGate, "ApplyNOTGate")
        .def("swap", &WeightedBDDQuantumCircuit::ApplySwapGate, "ApplySwapGate")
        .def("iswap", &WeightedBDDQuantumCircuit::ApplyiSwapGate, "ApplyiSwapGate")
        .def("prob", &WeightedBDDQuantumCircuit::GetProbability, "GetProbability")
        .def("measure", &WeightedBDDQuantumCircuit::Measure, "Measure")
        .def("measurement_counts", &WeightedBDDQuantumCircuit::GetPathCount, "MeasurementCount")
        .def("y", &WeightedBDDQuantumCircuit::ApplyPauliYGate, "ApplyPauliYGate")
        .def("z", &WeightedBDDQuantumCircuit::ApplyPauliZGate, "ApplyPauliZGate")
        .def("p", &WeightedBDDQuantumCircuit::ApplyPhaseShiftGate, "ApplyPhaseShiftGate")
        .def("s", &WeightedBDDQuantumCircuit::ApplySGate, "ApplySGate")
        .def("t", &WeightedBDDQuantumCircuit::ApplyTGate, "ApplyTGate") 
        .def("cz", &WeightedBDDQuantumCircuit::ApplyCZGate, "ApplyCZGate")
        .def("cx", &WeightedBDDQuantumCircuit::ApplyCNOTGate, "ApplyCNOTGate")
        .def("ccx", &WeightedBDDQuantumCircuit::ApplyCCNOTGate, "ApplyCCNOTGate")
        .def("gp", &WeightedBDDQuantumCircuit::ApplyGlobalPhase, "ApplyGlobalPhase")
        .def("cp", &WeightedBDDQuantumCircuit::ApplyCPGate, "ApplyCPGate")
        .def("cs", &WeightedBDDQuantumCircuit::ApplyCSGate, "ApplyCSGate")
        .def("cswap", &WeightedBDDQuantumCircuit::ApplyCSwapGate, "ApplyCSwapGate");

    // py::class_<QuantumState>(m, "QuantumState");
    
    // py::class_<CFLOBDDQuantumState, QuantumState>(m, "CFLOBDDQuantumState")
    //     .def(py::init<>())
    //     .def("print", &CFLOBDDQuantumState::Print, "Print");
}
