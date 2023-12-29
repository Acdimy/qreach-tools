import quasimodo
import time
from math import pi
from Reachability import *
from math import sqrt

t_start=time.time()
cir=QuantumCircuit.from_qasm_file("benchmark/qrw_6.qasm")
qc = generateCir(cir.num_qubits)
qc.setRealQubits(cir.num_qubits)
e = QuantumError('ad', [len(cir), 0], 1, [0.5])
# modified to "loadGates" with QMarkov module qc.fromQiskitGates() fromMarkovModel
qc = loadQiskitGates(cir, qc, [e])
e.err_channel = 2
qc = loadQiskitGates(cir, qc, [e])

# setInitialStates([])
qc.setState(str_padding("00100000"))
qc.setProjectorFS()
qc.setState(str_padding("00010000"))
qc.setProjectorFS()
reachable_dim = qc.reachability()

print(f"use {time.time()-t_start} seconds")

print(reachable_dim)
if cir.num_qubits <= 4:
    qc.printProjector()
