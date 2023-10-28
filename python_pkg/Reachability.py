import sys
import quasimodo
import time
from utils import *
from qiskit import QuantumCircuit
from math import pi, log2, ceil

# class KrausOperator:
#     def __init__(self) -> None:
#         self.ops = []
#     def append(filename:str):
#         cir=QuantumCircuit.from_qasm_file(filename)
    # def concat():
    #     # Ei\tensor Ei*
    #     pass

def clean_pi(x):
    return round(x/pi, 6)

def str_padding(s):
    l = 2**ceil(log2(len(s)))
    return s+("0"*(l-len(s)))

def applyQiskitGates(cir, qc, isConj=False, l=0, r=100000):
    gates = cir.data
    qnum = cir.num_qubits
    # print(len(gates))
    for gate in gates[l:min(r, len(gates))]:
        idx = [cir.find_bit(bit).index for bit in gate.qubits]
        if isConj:
            idx = [i + qnum for i in idx]
        name = gate[0].name
        # print(name)
        if name == 'x':
            qc.x(idx[0])
        elif name == 'y':
            if not isConj:
                qc.y(idx[0])
            else:
                qc.y(idx[0])
                qc.gp(1)
        elif name == 'z':
            qc.z(idx[0])
        elif name == 'h':
            qc.h(idx[0])
        elif name == 's':
            if not isConj:
                qc.s(idx[0])
            else:
                qc.u3(idx[0], 0, 0, -0.5)
        elif name == 'cx':
            qc.cx(idx[0], idx[1])
        elif name == 'u1':
            if not isConj:
                qc.u3(idx[0], 0, 0, clean_pi(gate[0].params[0]))
            else:
                qc.u3(idx[0], 0, 0, -clean_pi(gate[0].params[0]))
        elif name == 'u2':
            if not isConj:
                qc.u3(idx[0], 1/2, clean_pi(gate[0].params[0]), clean_pi(gate[0].params[1]))
            else:
                qc.u3(idx[0], 1/2, -clean_pi(gate[0].params[0]), -clean_pi(gate[0].params[1]))
        elif name == 'u3':
            if not isConj:
                qc.u3(idx[0], clean_pi(gate[0].params[0]), clean_pi(gate[0].params[1]), clean_pi(gate[0].params[2]))
            else:
                qc.u3(idx[0], clean_pi(gate[0].params[0]), -clean_pi(gate[0].params[1]), -clean_pi(gate[0].params[2]))
        else:
            print("Not Supported gate")
    return qc


def generateCir(qnum):
    qc = quasimodo.QuantumCircuit("CFLOBDD", 2**ceil(log2(qnum)))
    return qc

def readFile(path:str, filename:str, init_state:str):
    cir=QuantumCircuit.from_qasm_file(path+filename)
    # print(get_real_qubit_num(cir))
    # qc = quasimodo.QuantumCircuit("CFLOBDD", 2**ceil(log2(cir.num_qubits)))
    qc = generateCir(cir.num_qubits)
    qc.setState(init_state)
    qc = applyQiskitGates(cir, qc)
    return cir, qc

def fromMatrix(mat):
    pass

# Number 0, make sure everything is right!!
# First, store the circuit in python, and make cpp just a calculation tool
# Then try to partite circuit in CFLOBDD in cpp, try to store different operators in cpp quantum circuit.
def imageComputation(path:str, filename:str, init_state:str):
    cir, qc = readFile(path, filename, init_state)
    return cir, qc

def applyMatRep(cir_list, qc):
    for cir in cir_list:
        qc.pushState2Cache()
        applyQiskitGates(cir, qc)
        applyQiskitGates(cir, qc, True)
        # qc.popCache2State()

def spaceUnion() -> None:
    pass

def getReachableSp() -> None:
    pass

def checkReachablility() -> None:
    pass

