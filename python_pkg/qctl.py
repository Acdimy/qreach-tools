import pyqreach
from qiskit import QuantumCircuit, QuantumRegister, ClassicalRegister

def tsLabelling(ts, op: pyqreach.QOperation, label: str):
    for loc in range(ts.getLocationNum()):
        if ts.Locations[loc].satisfy(op):
            ts.setLabel(loc, label)
            # print(f"Location {loc} labelled with {label}")

def tsLabellingClRegList(ts, clRegList: list, label: str):
    """
    Label locations in the transition system based on classical register values.
    :param ts: Transition system
    :param clRegList: List of strings as classical registers
    :label: Label to assign to the locations that satisfy the classical register values
    """
    for loc in range(ts.getLocationNum()):
        for clReg in clRegList:
            # convert clReg to a binary list
            clRegBin = [int(bit) for bit in clReg]
            if ts.Locations[loc].satisfyBit(list(range(len(clRegBin))), clRegBin):
                ts.setLabel(loc, label)
                break
