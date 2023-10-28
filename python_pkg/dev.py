import quasimodo
import time
from math import pi
from Reachability import *
from math import sqrt

t_start=time.time()
# cir, qc = imageComputation("benchmark/", "random_walk_14.qasm", str_padding("0000000000000000"))

cir=QuantumCircuit.from_qasm_file("benchmark/random_walk_10.qasm")
# Reachability
# cir_err = QuantumCircuit.from_qasm_file("benchmark/test.qasm")
# qc = generateCir(cir.num_qubits)
# qc.setState("00000000")
# qc.setProjectorFS()
# qc.setState("10000000")
# qc.setProjectorFS()
# qc.ApplyProjToEnt()
# for i in range(16):
#     print("turn: ", i)
#     applyMatRep([cir], qc)
#     qc.printSize("state")

# representation
qc = generateCir(cir.num_qubits)
qc.setInitGate()
for i in range(0, len(cir.data), 1000):
    qc = applyQiskitGates(cir, qc, l=i, r=i+1000)
    cnt = qc.printSize("state")
    if cnt > 1000:
        qc.setProjectorFS()
        qc.setInitGate()

print(f"use {time.time()-t_start} seconds")



# qc = quasimodo.QuantumCircuit("CFLOBDD", 2, seed=round(time.time()))
# qc.setState("00")

# qc.u3(0,0,0,-1/4)
# qc.x(0)
# print(qc.measure())
# qc.print()
# qc.printColInterleaved()

# qc.test()

# qc.printRV("kraus")

# qc.printSize("kraus")

# print(str_padding("01000000"))

# qc = generateCir(8)
# qc.setEntangle()
qc.printRV("state")
cnt = qc.printSize("state")
print(cnt)
# qc.print()
# qc.test()