import importlib.util
import csv
import sys
import os
import pandas as pd
from unittest.mock import patch
import gc
from qiskit.quantum_info import random_statevector, Statevector
import pyqreach
from qiskit import QuantumCircuit, QuantumRegister, ClassicalRegister, transpile
from qiskit.quantum_info import Statevector
from qiskit_aer import Aer
from qiskit.quantum_info import Pauli
import numpy as np
from math import pi
import random
from time import time
from parse_qiskit import *
from qctl import *
from circ_utils import *
from contextlib import contextmanager
import random
import importlib

PATH = "benchmark"
PROPERTY_CHECKS = {}

def register_property(algorithm_name):
    def wrapper(func):
        PROPERTY_CHECKS[algorithm_name] = func
        return func
    return wrapper

def import_function(module_name, path, function_name):
    """强制从指定路径加载一个函数"""
    spec = importlib.util.spec_from_file_location(module_name, path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[module_name] = module
    spec.loader.exec_module(module)
    return getattr(module, function_name)

@contextmanager
def use_function(target_module, function_name, new_function):
    """
    临时替换某个模块中的函数。
    例如：use_function("case_studies.quantum_teleportation.quantum_teleportation", "quantum_teleportation", mutant_func)
    """
    module = importlib.import_module(target_module)
    old_func = getattr(module, function_name)
    setattr(module, function_name, new_function)
    try:
        yield
    finally:
        setattr(module, function_name, old_func)

def generate_random_state(num_qubits, seed):
    """Generate a random quantum state for a given number of qubits."""
    dim = 2 ** num_qubits
    state = random_statevector(dim, seed).data
    return state

@register_property("quantum_teleportation")
def test_quantum_teleportation(circ: QuantumCircuit, num_tests=10):
    ts = pyqreach.TransitionSystem(False)
    satisfied_num = 0
    qc = circ.copy()
    qc_check = QuantumCircuit(3)
    qc_check.reset(0)
    qc_check.reset(1)
    qc_check.swap(0, 2)
    qc = qc.compose(qc_check)
    prepare_start_time = time()
    resultList = parse_qiskit_cir(qc, 3, ts)
    prepare_end_time = time()
    model_start_time = time()
    for t in range(num_tests):
        ts.resetLocationBounds()
        seed = random.randint(0, 2147483647)
        initial_state = generate_random_state(1, seed)
        initial_state = expand_amplitude(initial_state, [0], 3)
        initial_state = [x.real for x in initial_state] + [x.imag for x in initial_state]
        prop_init = pyqreach.QOperation(initial_state, 3)
        ts.setAnnotation([[0, prop_init]])
        ts.computingFixedPointPost()
        # qc = QuantumCircuit(3, 3)
        # qc.initialize(initial_state, 0)
        # qc = qc.compose(circ)
        # resultList = parse_qiskit_cir(qc, 3, ts)
        # prop_init = ts.Locations[2].lowerBound  # after init
        # visualize_transition_system(ts, f"quantum_teleportation_mutation_test")
        satisfied = True
        for loc in resultList:
            if ts.Locations[loc].satisfy(prop_init) == False:
                satisfied = False
                break
        if satisfied:
            satisfied_num += 1
    model_end_time = time()
    resultDict = {}
    resultDict["num_tests"] = num_tests
    resultDict["satisfied"] = (satisfied_num == num_tests)
    resultDict["satisfiedPercent"] = satisfied_num / num_tests
    resultDict["prepare_time"] = prepare_end_time - prepare_start_time
    resultDict["model_time"] = model_end_time - model_start_time
    return resultDict

@register_property("superdense_coding")
def test_superdense_coding(circ: QuantumCircuit, num_tests=10):
    # num_tests is not used in checking superdense coding
    ts = pyqreach.TransitionSystem(False)
    prepare_start_time = time()
    resultList = parse_qiskit_cir(circ, 3, ts)
    op0 = pyqreach.QOperation(["000"])
    ts.setAnnotation([[0, op0]])
    ts.computingFixedPointPost()
    for loc in resultList:
        ts.setLabel(loc, "leaf")
    # op00 = pyqreach.QOperation(["000", "001"])
    # op01 = pyqreach.QOperation(["010", "011"])
    # op10 = pyqreach.QOperation(["100", "101"])
    # op11 = pyqreach.QOperation(["110", "111"])
    # tsLabelling(ts, op00, "p00")
    # tsLabelling(ts, op01, "p01")
    # tsLabelling(ts, op10, "p10")
    # tsLabelling(ts, op11, "p11")
    for i in range(4):
        # Switch i to bit strings
        bit_str = format(i, '02b')
        op = pyqreach.QOperation([f"{bit_str}0", f"{bit_str}1"])
        tsLabelling(ts, op, f"q{bit_str}")
        tsLabellingClRegList(ts, [bit_str], f"c{bit_str}")
    prepare_end_time = time()
    model_start_time = time()
    result = modelChecking(ts, 'AG (leaf -> ((q00 <-> c00) & (q01 <-> c01) & (q10 <-> c10) & (q11 <-> c11)))')
    model_end_time = time()
    result["prepare_time"] = prepare_end_time - prepare_start_time
    result["model_time"] = model_end_time - model_start_time
    return result

def run_single_experiment(algorithm_name, mutant_file=None, num_tests=10):
    """
    运行一次实验：可以是原始算法，也可以是 mutant
    """
    if mutant_file:
        mutant_name = mutant_file[:-3]  # 去掉 .py
        func = import_function(mutant_name,
                               f"{PATH}/{algorithm_name}/mutants/{mutant_file}",
                               algorithm_name)
    else:
        func = import_function(algorithm_name,
                               f"{PATH}/{algorithm_name}/{algorithm_name}.py",
                               algorithm_name)

    print(f"Running {algorithm_name} mutant={mutant_file or 'original'}")
    

    # 这里临时替换原始函数
    with use_function(f"benchmark.{algorithm_name}.{algorithm_name}", algorithm_name, func):
        qc = func()   # 调用替换后的函数，得到 QuantumCircuit
        # TODO: 在这里调用 property checks
        result = PROPERTY_CHECKS[algorithm_name](qc, num_tests=num_tests)
        return result
    return None
        

def run_all_experiments(algorithm_name, num_tests=10):
    """
    遍历某个算法下的所有 mutants 并运行实验
    """    
    result_list = []
    mutant_folder = f"{PATH}/{algorithm_name}/mutants"
    mutants = [f for f in os.listdir(mutant_folder) if f.endswith(".py")]
    pyqreach.initializeTransitionSystem()
    # 先运行原始算法
    print("Running original algorithm")
    original_result = run_single_experiment(algorithm_name, num_tests=num_tests)
    # append the experiment name to the result
    original_result['mutant'] = 'original'
    result_list.append(original_result)    

    # 再运行 mutants
    for mf in mutants:
        print(f"Running mutant {mf}")
        mf_result = run_single_experiment(algorithm_name, mf, num_tests=num_tests)
        mf_result['mutant'] = mf[:-3]  # 去掉 .py
        result_list.append(mf_result)
    # dump the results into a csv file, path: eval/mutation_results/{algorithm_name}_results.csv
    # The keys of the result dict: satisfied, counterexample, output, prepare_time, model_time
    keys = result_list[0].keys()
    # exclude output
    keys = [k for k in keys if k != 'output' and k != 'counterexample']
    result_df = pd.DataFrame(result_list, columns=keys)
    os.makedirs("eval/mutation_results", exist_ok=True)
    result_df.to_csv(f"eval/mutation_results/{algorithm_name}_results.csv", index=False)
    return result_df


if __name__ == "__main__":
    random.seed(random.randint(0, 2147483647))
    # pyqreach.initializeTransitionSystem()
    # print(run_single_experiment("quantum_teleportation", num_tests=1))
    # result = run_single_experiment("superdense_coding", mutant_file="superdense_coding_m3.py")
    # print(result['output'])
    # run_all_experiments("superdense_coding")
    run_all_experiments("quantum_teleportation", num_tests=64)

