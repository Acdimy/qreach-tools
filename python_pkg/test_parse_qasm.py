import pyqreach
### Verifiable quantum secret sharing
from qiskit import QuantumCircuit, QuantumRegister, ClassicalRegister, transpile
from qiskit.quantum_info import Statevector
from qiskit_aer import Aer
import numpy as np
from math import pi
import random
from time import time
from parse_qiskit import *
from qctl import *
from circ_utils import *
import pandas as pd

def run_single_test(filename:str, savefile:str=None, init_state:str=None):
    # filename = "benchmark/grover_32.qasm"
    qc = QuantumCircuit.from_qasm_file(filename)
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
    result_dict = {
        "filename": filename,
        "num_qubits": qc.num_qubits,
        "num_gates": len(qc.data),
        "num_locations": ts.getLocationNum(),
        "num_result_list": len(resultList),
        "time_fixed_point_post": end_time - start_time,
        "time_prepare_ts": end_time_prepare - start_time_prepare,
        "init_state": init_state if init_state is not None else "0"
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
DQCQFT_LIST = list(range(12,13))
RUS_LIST = list(range(1,4))
QFT_LIST = list(range(2,16))
PE_LIST = list(range(2,13))
TYPE_DICT = {"grover": GROVER_LIST, "dqc_qft": DQCQFT_LIST, "rus": RUS_LIST, "qft": QFT_LIST, "pe": PE_LIST}

def gen_random_basis_state(num_qubits:int) -> str:
    state = ''.join(random.choice(['0', '1']) for _ in range(num_qubits))
    # state = "1"*num_qubits
    return state

def run_type(type_name:str="grover"):
    if type_name not in TYPE_DICT:
        print(f"Type {type_name} not recognized. Available types: {list(TYPE_DICT.keys())}")
        return
    init_state = None
    for val in TYPE_DICT[type_name]:
        if type_name == "grover":
            filename = f"benchmark/grover/grover_{val}.qasm"
        elif type_name == "dqc_qft":
            filename = f"benchmark/dqc_qft/dqc_qft_{val}.qasm"
            init_state = gen_random_basis_state(val)
        elif type_name == "rus":
            filename = f"benchmark/rus_{val}.qasm"
        elif type_name == "qft":
            filename = f"benchmark/qft/qft_{val}.qasm"
            init_state = gen_random_basis_state(val)
        elif type_name == "pe":
            filename = f"benchmark/pe/dqc_pe_{val}.qasm"
            init_state = gen_random_basis_state(val+1)
        else:
            continue
        print(f"Running test for {filename}")
        run_single_test(filename, f"eval/scale/{type_name}_results.csv", init_state)
        print("-"*40)

if __name__ == "__main__":
    pyqreach.initializeTransitionSystem()
    # run_single_test("benchmark/dqc_qft/dqc_qft_10.qasm", init_state=gen_random_basis_state(10))
    run_type("grover")
    # run_type("dqc_qft")
    # run_type("pe")
    run_type("qft")
    # run_type("rus")
