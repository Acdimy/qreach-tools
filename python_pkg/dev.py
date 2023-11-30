import quasimodo
import time
from math import pi
from Reachability import *
from math import sqrt

t_start=time.time()
# cir, qc = imageComputation("benchmark/", "random_walk_14.qasm", str_padding("0000000000000000"))

cir=QuantumCircuit.from_qasm_file("benchmark/test.qasm")
# cir_err = QuantumCircuit.from_qasm_file("benchmark/test.qasm")
# Reachability, fomula method
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

# representation, apply state in a bunch
# qc = generateCir(cir.num_qubits)
# qc.setInitGate()
# for i in range(0, len(cir.data), 1000):
#     qc = applyQiskitGates(cir, qc, l=i, r=i+1000)
#     cnt = qc.printSize("state")
#     if cnt > 1000:
#         qc.setProjectorFS()
#         qc.setInitGate()

# # arbitrary u test
# qc = generateCir(2)
# qc.setState("10")
# qc.u(0, [0,0,0,0,0,1,1,0])
# qc.printColHead()

# real reachability testcase 1
# qc = generateCir(cir.num_qubits)
# e = QuantumError('ad', [len(cir), 0], 1, [0.5])
# qc = loadQiskitGates(cir, qc, e)
# e.err_channel = 2
# qc = loadQiskitGates(cir, qc, e)
# qc.setState("00100000")
# qc.setProjectorFS()
# qc.setState("01000000")
# qc.setProjectorFS()
# reachable_dim = qc.reachability()

qc = generateCir(cir.num_qubits)
e1 = QuantumError('measure', [len(cir), 0], 1, [])
e2 = QuantumError('measure', [len(cir), 0], 2, [])
e3 = QuantumError('measure', [len(cir), 1], 1, [])
e4 = QuantumError('measure', [len(cir), 1], 2, [])
qc = loadQiskitGates(cir, qc, [e1,e3])
qc = loadQiskitGates(cir, qc, [e1,e4])
qc.appendGateSeries("x", [1], [], False)
qc = loadQiskitGates(cir, qc, [e2,e3])
qc.appendGateSeries("x", [0], [], False)
qc = loadQiskitGates(cir, qc, [e2,e4])
qc.appendGateSeries("x", [0], [], False)
qc.appendGateSeries("x", [1], [], False)
qc.setState(str_padding("000"))
qc.setProjectorFS()
reachable_dim = qc.reachability()

print(f"use {time.time()-t_start} seconds")

print(reachable_dim)

# qc = quasimodo.QuantumCircuit("CFLOBDD", 2, seed=round(time.time()))
# qc.setState("00")

## print
# qc.print()
# qc.printRV("state")
# cnt = qc.printSize("state")
# print(cnt)
qc.printProjector()
# qc.test()