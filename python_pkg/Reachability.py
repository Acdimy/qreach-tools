import sys
import quasimodo
import time
from utils import *
from qiskit import QuantumCircuit
from math import pi, log2, ceil

class QuantumError:
    def __init__(self, err_type='', err_pos=[-1,-1], err_channel=-1, err_params=[]) -> None:
        self.err_type = err_type
        self.err_pos = err_pos
        self.err_channel = err_channel
        self.err_params = err_params

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
        elif name == 'ccx':
            qc.ccx(idx[0], idx[1], idx[2])
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

def loadQiskitGates(cir, qc, e_list:[]):
    """
    err_pos: [position in gate series, qubit index]
    """
    gates = list(cir.data)
    e_list = sorted(e_list, key=lambda x:x.err_pos[0], reverse=True)
    for e in e_list:
        gates.insert(e.err_pos[0], e)
    for i,gate in enumerate(gates):
        if isinstance(gate, QuantumError):
            err_idx = [gate.err_pos[1], gate.err_channel]
            qc.appendGateSeries(gate.err_type, err_idx, gate.err_params, i==0)
        else:
            idx = [cir.find_bit(bit).index for bit in gate.qubits]
            name = gate[0].name
            params = gate[0].params
            params = [clean_pi(x) for x in params]
            qc.appendGateSeries(name, idx, params, i==0)
    return qc

# def loadQiskitGates(cir, qc, e:QuantumError):
#     """
#     err_pos: [position in gate series, qubit index]
#     """
#     gates = cir.data
#     for i,gate in enumerate(gates):
#         if i == e.err_pos[0] and i != 0:
#             err_idx = [e.err_pos[1], e.err_channel]
#             qc.appendGateSeries(e.err_type, err_idx, e.err_params, False)
#         idx = [cir.find_bit(bit).index for bit in gate.qubits]
#         name = gate[0].name
#         params = gate[0].params
#         params = [clean_pi(x) for x in params]
#         # print(name)
#         if i == 0:
#             if i == e.err_pos[0]:
#                 err_idx = [e.err_pos[1], e.err_channel]
#                 qc.appendGateSeries(e.err_type, err_idx, e.err_params, True)
#             qc.appendGateSeries(name, idx, params, True)
#         else:
#             qc.appendGateSeries(name, idx, params, False)
#     if e.err_pos[0] == len(gates):
#         err_idx = [e.err_pos[1], e.err_channel]
#         qc.appendGateSeries(e.err_type, err_idx, e.err_params, False)
#     return qc

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

