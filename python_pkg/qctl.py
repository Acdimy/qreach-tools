import pyqreach
from qiskit import QuantumCircuit, QuantumRegister, ClassicalRegister

def tsLabelling(ts, op: pyqreach.QOperation, label: str):
    for loc in range(ts.getLocationNum()):
        if ts.Locations[loc].satisfy(op):
            ts.setLabel(loc, label)
            # print(f"Location {loc} labelled with {label}")

