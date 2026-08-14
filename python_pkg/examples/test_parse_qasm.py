# This is a test file to run DQC_PE and Grover in the paper.

import pyqreach
### Verifiable quantum secret sharing
from qiskit import QuantumCircuit, QuantumRegister, ClassicalRegister, transpile
from qiskit.circuit import CircuitInstruction
from qiskit.circuit.library import XGate, YGate, ZGate
from qiskit.quantum_info import Statevector
from qiskit_aer import Aer
import numpy as np
from math import pi
import random
from time import time
from qreach.parse_qiskit import *
from qreach.qctl import *
from qreach.circ_utils import *
import pandas as pd

def insert_random_pauli(qc):
    """Returns a new QuantumCircuit with a random Pauli gate (X, Y, or Z) inserted at a random position."""
    import copy
    new_qc = copy.deepcopy(qc)
    num_qubits = new_qc.num_qubits
    num_ops = len(new_qc.data)
    pos = random.randint(0, num_ops)
    q = random.randint(0, num_qubits - 1)
    gate_name = random.choice(["x", "y", "z"])
    gate_map = {"x": XGate(), "y": YGate(), "z": ZGate()}
    instr = CircuitInstruction(gate_map[gate_name], [new_qc.qubits[q]], [])
    new_qc.data.insert(pos, instr)
    return new_qc, (pos, q, gate_name)

def simulate_circuit(qc:QuantumCircuit) -> pyqreach.QOperation:
    ts = pyqreach.TransitionSystem(False)
    resultList = parse_qiskit_cir(qc, qc.num_qubits, ts)
    op00 = pyqreach.QOperation(["0"*qc.num_qubits])
    ts.setAnnotation([[0, op00]])
    ts.computingFixedPointPost()
    final_op = ts.Locations[resultList[-1]].lowerBound
    return final_op

def generate_debug_info(num_qubits:int, filename:str, ts:pyqreach.TransitionSystem, resultList:list, qc_init:QuantumCircuit=None) -> dict:
    # if "grover" in filename:
    if "grover" in filename:
        assert num_qubits % 2 == 1, "Grover benchmark should have odd number of qubits"
        work_qubits = int((num_qubits+1)/2)
        grover_init = ts.Locations[work_qubits].lowerBound # superposition |++...+00...01>
        grover_good = pyqreach.QOperation(["1"*work_qubits + "0"*(num_qubits - work_qubits - 1) + "1"]) # good state |11...100...01>
        ts_temp = pyqreach.TransitionSystem(False)
        loc0, loc1, loc2 = pyqreach.Location(num_qubits,0), pyqreach.Location(num_qubits,1), pyqreach.Location(num_qubits,2)
        ts_temp.addLocation(loc0)
        ts_temp.addLocation(loc1)
        ts_temp.addLocation(loc2)
        ts_temp.addRelation(0, 2, pyqreach.QOperation("I", num_qubits, [0], []))
        ts_temp.addRelation(1, 2, pyqreach.QOperation("I", num_qubits, [0], []))
        ts_temp.setAnnotation([[0, grover_init], [1, grover_good]])
        ts_temp.computingFixedPointPost()
        grover_final = ts_temp.Locations[2].lowerBound
        debug_dict = {}
        debug_res = ts.Locations[resultList[-1]].satisfy(grover_final)
        debug_dict["satisfied"] = debug_res
        return debug_dict
    elif "dqc_pe" in filename:
        debug_dict = {}
        if num_qubits < 11:
            # Cannot get an exact phase estimation
            # Property: Globally has some states in different measurement results
            tsLabellingDefault(ts, "reached")
            # debug_dict = modelChecking(ts, "AG reached")
            pe_pattern = "1"+"0"*(num_qubits-2)
            tsLabellingClRegList(ts, [pe_pattern], "pe_success", locList=resultList)
            debug_dict = modelChecking(ts, "AG ((pe_success -> reached))")
            print(debug_dict["satisfied"])
        else:
            tsLabellingDefault(ts, "reached", resultList)
            pe_pattern = (num_qubits-11)*"0"+"1000000000"
            tsLabellingClRegList(ts, [pe_pattern], "pe_success", locList=resultList)
            debug_dict = modelChecking(ts, "AG ((pe_success -> reached))")
            print(debug_dict["satisfied"])
        return debug_dict
    else:
        debug_op = simulate_circuit(qc_init)
        ts.setLabel(resultList[-1], "final")
        tsLabelling(ts, debug_op, "debug_op", locList=resultList)
        debug_dict = modelChecking(ts, "AG (debug_op <-> final)")
        print(debug_dict["satisfied"])
        return debug_dict

def run_single_test(filename:str, savefile:str=None, init_state:str=None, error_injection:bool=False, debug:bool=False) -> dict:
    qc = QuantumCircuit.from_qasm_file(filename)
    qc_init = None
    if error_injection:
        # Store the original circuit
        import copy
        qc_init = copy.deepcopy(qc)
        qc, injection_info = insert_random_pauli(qc)
        # print(f"Inserted random Pauli gate {injection_info[2].upper()} on qubit {injection_info[1]} at position {injection_info[0]}")
    ts = pyqreach.TransitionSystem(False)
    start_time_prepare = time()
    resultList = parse_qiskit_cir(qc, qc.num_qubits, ts)
    end_time_prepare = time()
    print(f"Time taken for constructing transition system: {end_time_prepare - start_time_prepare:.2f} seconds")
    print("Transition System Locations:", ts.getLocationNum())
    print("Result List size:", len(resultList))
    op00 = pyqreach.QOperation(["0"*qc.num_qubits if init_state is None else init_state])
    ts.setAnnotation([[0, op00]])
    start_time = time()
    ts.computingFixedPointPost()
    end_time = time()
    print(f"Time taken for computing fixed point post: {end_time - start_time:.2f} seconds")
    verification_begin = time()
    if debug:
        debug_info = generate_debug_info(qc.num_qubits, filename, ts, resultList, qc_init=qc_init)
        # print("Debug Info:", debug_info)
    verification_end = time()
    result_dict = {
        "filename": filename,
        "num_qubits": qc.num_qubits,
        "num_gates": len(qc.data),
        "num_locations": ts.getLocationNum(),
        "num_result_list": len(resultList),
        "time_fixed_point_post": end_time - start_time,
        "time_prepare_ts": end_time_prepare - start_time_prepare,
        # "init_state": init_state if init_state is not None else "0",
        "verification_time": verification_end - verification_begin if debug else "None",
        "error_injection": injection_info if error_injection else "None",
        "debug_info": debug_info["satisfied"] if debug else "None"
    }
    if savefile is not None:
        df = pd.DataFrame([result_dict])
        try:
            df_existing = pd.read_csv(savefile)
            df = pd.concat([df_existing, df], ignore_index=True)
        except FileNotFoundError:
            pass
        df.to_csv(savefile, index=False)
    return result_dict

GROVER_LIST = list(range(3,100,2))
DQCQFT_LIST = list(range(2,13))
RUS_LIST = list(range(1,4))
QFT_LIST = list(range(2,13))
PE_LIST = list(range(2,13))
TYPE_DICT = {"grover": GROVER_LIST, "dqc_qft": DQCQFT_LIST, "rus": RUS_LIST, "qft": QFT_LIST, "dqc_pe": PE_LIST}

def gen_random_basis_state(num_qubits:int) -> str:
    state = ''.join(random.choice(['0', '1']) for _ in range(num_qubits))
    # state = "1"*num_qubits
    return state

def run_type(type_name:str="grover", debug:bool=False, error_injection:bool=False, rep=1):
    if type_name not in TYPE_DICT:
        print(f"Type {type_name} not recognized. Available types: {list(TYPE_DICT.keys())}")
        return
    init_state = None
    for val in TYPE_DICT[type_name]:
        if type_name == "grover":
            filename = f"benchmark/grover/grover_{val}.qasm"
            init_state = "0"*(val-1)+"1"
        elif type_name == "dqc_qft":
            filename = f"benchmark/dqc_qft/dqc_qft_{val}.qasm"
            init_state = gen_random_basis_state(val)
        elif type_name == "rus":
            filename = f"benchmark/rus_{val}.qasm"
        elif type_name == "qft":
            filename = f"benchmark/qft/qft_{val}.qasm"
            init_state = gen_random_basis_state(val)
        elif type_name == "dqc_pe":
            filename = f"benchmark/dqc_pe/dqc_pe_{val}.qasm"
            init_state = "0"*val+"1" if not error_injection else "0"*(val+1)
        else:
            continue
        print(f"Running test for {filename}")
        injection_info = "inj" if error_injection else "ori"
        error_inj = False if type_name == "dqc_pe" else error_injection
        result_dict_list = []
        for _i in range(rep):
            result_dict = run_single_test(filename, savefile=None, init_state=init_state, debug=debug, error_injection=error_inj)
            result_dict_list.append(result_dict)
        
        final_result_dict = {}
        final_result_dict["filename"] = filename
        final_result_dict["num_qubits"] = result_dict_list[0]["num_qubits"]
        # gates, locations, result_list, time: average over rep
        for key in ["num_gates", "num_locations", "num_result_list", "time_fixed_point_post", "time_prepare_ts", "verification_time"]:
            final_result_dict[key] = sum(d[key] for d in result_dict_list) / rep
        final_result_dict["error_injection"] = [d["error_injection"] for d in result_dict_list]
        final_result_dict["debug_info"] = [d["debug_info"] for d in result_dict_list]
        debug_sat_percent = sum(1 for d in result_dict_list if d["debug_info"]==True) / rep * 100 if debug else "N/A"
        final_result_dict["debug_satisfaction_percent"] = debug_sat_percent
        
        savefile = f"eval/scale_debug/{type_name}_results_{injection_info}_0218.csv"
        df = pd.DataFrame([final_result_dict])
        try:
            df_existing = pd.read_csv(savefile)
            df = pd.concat([df_existing, df], ignore_index=True)
        except FileNotFoundError:
            pass
        df.to_csv(savefile, index=False)
        print("-"*40)

if __name__ == "__main__":
    pyqreach.initializeTransitionSystem()
    random.seed(42)
    # run_single_test("benchmark/grover/grover_55.qasm", init_state="0"*54+"1", debug=True, error_injection=True)
    # test_pe_num = 8
    # run_single_test(f"benchmark/dqc_pe/dqc_pe_{test_pe_num}.qasm", init_state="0"*test_pe_num+"0", debug=True, error_injection=False)
    
    # run_type("grover", debug=True, error_injection=True, rep=50)
    run_type("dqc_pe", debug=True, error_injection=True)
