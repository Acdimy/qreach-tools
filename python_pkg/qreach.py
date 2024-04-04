import pyqreach


def QuantumCircuit(model_str = "CFLOBDD", numQubits = 0, seed = 0):
    if model_str == "CFLOBDD":
        return pyqreach.CFLOBDDQuantumCircuit(numQubits, seed)
    elif model_str == "BDD":
        return pyqreach.BDDQuantumCircuit(numQubits, seed)
    elif model_str == "WBDD":
        return pyqreach.WeightedBDDQuantumCircuit(numQubits, seed)
    else:
        print ('Unsupported backend: ', model_str)
        exit(1)