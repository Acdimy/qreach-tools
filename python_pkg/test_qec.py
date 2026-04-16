import pyqreach
### Verifiable quantum secret sharing
from qiskit import QuantumCircuit, QuantumRegister, ClassicalRegister, transpile
from qiskit.quantum_info import Statevector, state_fidelity, partial_trace
# from qiskit_aer import AerSimulator
import numpy as np
from time import time
from parse_qiskit import *
from qctl import *
from circ_utils import *
from qiskit_aer import AerSimulator

# ---------- 构造电路 ----------
data = QuantumRegister(5, 'data')
anc  = QuantumRegister(4, 'anc')
creg = ClassicalRegister(4, 'c')

circ = QuantumCircuit(data, anc, creg)

# ---------- Step 1: 准备输入态 ----------
# |ψ> = (|0> + i|1>)/sqrt(2)
circ.h(data[0])
circ.s(data[0])

# 保存理想态
ideal_single = Statevector.from_instruction(circ)

prepare_5perfect_code(circ, data)

# 插入错误
circ.x(data[2])


perfect_code_syndrome_and_correct(
    circ,
    data_idx=data,
    anc_idx=anc,
    creg=creg
)

decode_5perfect_code(circ, data)

# final_state = Statevector.from_instruction(circ.remove_final_measurements(inplace=False))

# # 只取第一个逻辑比特
# # 把其余比特 trace 掉
# final_logical = final_state.partial_trace([1,2,3,4,5,6,7,8])  # 只保留 data[0]

# # ---------- Step 7: Fidelity ----------
# fidelity = state_fidelity(final_logical, ideal_single)

# print("Fidelity =", fidelity)

