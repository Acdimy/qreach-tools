from checker import *
import sys

info = "Expected arguments: experiment type ['qrw', 'rus', 'grover']; qubits; quantum errors which are set as 'amplitude dumping'"
assert len(sys.argv) >= 3, info
expr_name = sys.argv[1]
qubits = sys.argv[2]
if len(sys.argv) >= 4:
    err_list = sys.argv[3:] # Not available
cir = "benchmark/"+expr_name+"_"+qubits+".qasm"


t_start=time.time()
if expr_name == "qrw":
    if len(sys.argv) >= 4:
        err_model = [Noise("ad", [-1,0], [0.5])]
    else:
        err_model = []
    qmc = QuantumMarkovChain(cir=cir, err_model=err_model)
elif expr_name == "rus":
    pass
elif expr_name == "grover":
    pass
else:
    raise RuntimeError(info)
qchecker = fromMarkovModel(qmc)
qchecker = initWithStr(qchecker, str_list=["0","1"])

reachable_dim = qchecker.reachability()
# qchecker.printProjector()

print(f"use {time.time()-t_start} seconds")

print("Reachable dimensions: ", reachable_dim)
qchecker.printSize("state")
qchecker.printSize("projector")
