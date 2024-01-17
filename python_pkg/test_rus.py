import quasimodo
import time
from math import pi
from Reachability import *
from math import sqrt

n = int(sys.argv[1])
t_start=time.time()
cir=QuantumCircuit.from_qasm_file("benchmark/rus_"+str(n)+".qasm")

qc = generateCir(cir.num_qubits)
qc.setRealQubits(cir.num_qubits)
e1 = QuantumError('measure', [len(cir), 0], 1, [])
e2 = QuantumError('measure', [len(cir), 0], 2, [])
if n == 1:
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
else:
    qc = loadQiskitGates(cir, qc, [e1])
    qc = loadQiskitGates(cir, qc, [e2])
    qc.appendGateSeries("x", [0], [], False)
    qc.setState("01")
    qc.setProjectorFS()
    reachable_dim = qc.reachability()

print(f"use {time.time()-t_start} seconds")

print(reachable_dim)
qc.printSize("state")
qc.printSize("projector")
qc.printProjector()
