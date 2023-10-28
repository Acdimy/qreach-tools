from qiskit import QuantumCircuit
import time

cir=QuantumCircuit.from_qasm_file("benchmark/random_walk_4.qasm")
# # Import Aer
# from qiskit import Aer

# # Run the quantum circuit on a statevector simulator backend
# backend = Aer.get_backend('statevector_simulator')
# print(type(backend))
# # Create a Quantum Program for execution
# job = backend.run(cir)

# result = job.result()
# outputstate = result.get_statevector(cir, decimals=3)
# print(outputstate)

# from qiskit import Aer
# backend = Aer.get_backend('unitary_simulator')
# job = backend.run(cir)
# result = job.result()

# # Show the results
# print(result.get_unitary(cir, decimals=3))

from qiskit.quantum_info import Statevector
# Set the initial state of the simulator to the ground state using from_int
state = Statevector.from_int(0, 2**6)
# Evolve the state by the quantum circuit
t_start=time.time()
state = state.evolve(cir)
print(f"use {time.time()-t_start} seconds")

state = state.data

import numpy as np
state = [np.round(s,3) for s in state]
val_dic = list(set(state))


print(val_dic)

