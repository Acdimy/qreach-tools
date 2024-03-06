from checker import *

t_start=time.time()
qmc = QuantumMarkovChain(cir="benchmark/qrw_2.qasm", err_model=[Noise("ad", [-1,0], [0.5])])
qchecker = fromMarkovModel(qmc)
qchecker = initWithStr(qchecker, str_list=["001"])

reachable_dim = qchecker.reachability()
# qchecker.printProjector()

print(f"use {time.time()-t_start} seconds")

print(reachable_dim)
qchecker.printSize("state")
qchecker.printSize("projector")
