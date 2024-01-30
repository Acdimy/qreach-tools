# Amplitude amplification

import qreach
import time
from math import pi
from Reachability import *
from math import sqrt
import sys

n = int(sys.argv[1])
t_start=time.time()
cir=QuantumCircuit.from_qasm_file("benchmark/grover_b"+str(n)+".qasm")
# cir=QuantumCircuit.from_qasm_file("benchmark/test.qasm")
qc = generateCir(cir.num_qubits)
qc.setRealQubits(cir.num_qubits)
qc = loadQiskitGates(cir, qc, [])

qc.setState(str_padding("0"*(2*n-2)+"1"))
for i in range(n):
    qc.h(i)
qc.h(2*n-2)
qc.setProjectorFS()
reachable_dim = qc.reachability()
print(f"use {time.time()-t_start} seconds")
print(reachable_dim)
qc.printSize("state")
qc.printSize("projector")
# qc.printProjector()
