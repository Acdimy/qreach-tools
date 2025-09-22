import importlib.util
import csv
import time
import sys
import os
import pandas as pd
from unittest.mock import patch
import gc

PATH = "./../qucheck-main"

def import_function(module_str, path, function_name):
    spec = importlib.util.spec_from_file_location(module_str, path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[module_str] = module
    spec.loader.exec_module(module)
    return getattr(module, function_name)

def reload_classes(folder_path):
    sys.path.insert(0, folder_path)
    for file in os.listdir(folder_path):
        if file.endswith('.py'):
            module = importlib.import_module(file[:-3])
            importlib.reload(module)
    sys.path.pop(0)

def run_single_test(algorithm_name, num_inputs, measurements, mutant_type, index, number_of_properties, run_optimization=True,
                    csvwriter=None):
    mutant_name = f"{algorithm_name}_{mutant_type}{index}"

    circuit_function = import_function(mutant_name,
                                       f"{PATH}/case_studies/{algorithm_name}/mutants/{mutant_name}.py",
                                       algorithm_name)
    print(f"Testing {mutant_name}")

    with patch(f"case_studies.{algorithm_name}.{algorithm_name}.{algorithm_name}", circuit_function):
        reload_classes(f"{PATH}/case_studies/{algorithm_name}")
