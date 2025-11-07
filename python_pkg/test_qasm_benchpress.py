from test_parse_qasm import *
import os
import shutil

BENCHPRESS_PATH = "../../benchpress/benchpress/qasm/qasmbench-medium/"
TARGET_PATH = "benchmark/benchpress-medium/supported/"

if __name__ == "__main__":
    # traverse all folders in BENCHPRESS_PATH and find all .qasm files in these folders
    qasm_files = []
    for root, dirs, files in os.walk(TARGET_PATH):
        for file in files:
            if file.endswith(".qasm"):
                src_path = os.path.join(root, file)
                qasm_files.append(src_path)
                # Copy the file to path benchmark/benchpress-medium/
                # os.makedirs(TARGET_PATH, exist_ok=True)

                # # Construct destination file path
                # dest_path = os.path.join(TARGET_PATH, file)

                # # Copy the file (preserving metadata)
                # shutil.copy2(src_path, dest_path)

                # print(f"Copied {src_path} -> {dest_path}")
                
    pyqreach.initializeTransitionSystem()
    random.seed(42)
    qasm_files = ["benchmark/benchpress-medium/supported/bv_n14_transpiled.qasm"]
    for qasm_file in qasm_files:
        print(f"Running test for {qasm_file}")
        run_single_test(qasm_file, "eval/scale_debug/benchpress_medium_inj_1106.csv", debug=True, error_injection=True)
        print("-"*40)
