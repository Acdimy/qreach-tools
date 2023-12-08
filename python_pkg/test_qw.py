import quasimodo
import time
from math import pi
from Reachability import *
from math import sqrt

t_start=time.time()
cir=QuantumCircuit.from_qasm_file("benchmark/qrw_2.qasm")
qc = generateCir(cir.num_qubits)
e = QuantumError('ad', [len(cir), 0], 1, [0.5])
qc = loadQiskitGates(cir, qc, [e])
e.err_channel = 2
qc = loadQiskitGates(cir, qc, [e])
qc.setState(str_padding("0010"))
qc.setProjectorFS()
qc.setState(str_padding("0001"))
qc.setProjectorFS()
reachable_dim = qc.reachability()

print(f"use {time.time()-t_start} seconds")

print(reachable_dim)
qc.printProjector()
