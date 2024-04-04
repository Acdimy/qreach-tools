from checker import *

err_model = [Noise("ad", [-1,0], [0.5])]
qmc = QuantumMarkovChain(cir="benchmark/qrw_2.qasm", err_model=err_model)
qchecker = fromMarkovModel(qmc)
qchecker = initWithStr(qchecker, str_list=["0"])
reachable_dim = qchecker.reachability()
print(reachable_dim)
