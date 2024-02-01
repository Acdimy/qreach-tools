import qreach
import time
from math import pi
from Reachability import *
from math import sqrt
import sys

n = int(sys.argv[1])
assert len(sys.argv) <= 3
if len(sys.argv) == 3:
    err_type = sys.argv[2]

t_start=time.time()
cir=QuantumCircuit.from_qasm_file("benchmark/qrw_"+str(n)+".qasm")
#Integration
qc = generateCir(cir.num_qubits)
qc.setRealQubits(cir.num_qubits)

if len(sys.argv) == 3:
    e = ChannelMode(err_type, [len(cir), 0], 1, [0.5])
    # modified to "loadGates" with QMarkov module qc.fromQiskitGates(cir, [e]) qc.fromMarkovModel(cir, qmc)
    qc = loadQiskitGates(cir, qc, [e])
    e.err_channel = 2
    qc = loadQiskitGates(cir, qc, [e])
else:
    qc = loadQiskitGates(cir, qc, [])



'''
from checker import *

QMC = QMarkov(cir_body="path", channels=[Qchannel])
qchecker = fromMarkovModel(QMC)
qchecker = initWithStr(qchecker, str_list=[])
reachable_dim = qchecker.reachability()
qchecker.printProjector()
'''

qc.setState(str_padding("0"+"0"*(cir.num_qubits-1)))
qc.setProjectorFS()
qc.setState(str_padding("1"+"0"*(cir.num_qubits-1)))
qc.setProjectorFS()

reachable_dim = qc.reachability()

print(f"use {time.time()-t_start} seconds")

print(reachable_dim)
qc.printSize("state")
qc.printSize("projector")
if cir.num_qubits <= 4:
    qc.printProjector()

# For qrw10
# use 1910.7084205150604 seconds
# 4097