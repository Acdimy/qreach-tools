import sys
import qreach
from qchannel import *

from qiskit import QuantumCircuit

class QuantumMarkovChain:
    def __init__(self, cir, err_model, bound=1e7) -> None:
        if isinstance(cir, str):
            try:
                cir=QuantumCircuit.from_qasm_file(cir)
            except IOerror as e:
                print(e)
        assert(isinstance(cir, QuantumCircuit))
        self.cir = cir
        # check err_model
        self.err_model = err_model
        self.bound = bound