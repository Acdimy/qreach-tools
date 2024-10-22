import sys
import qreach
from qchannel import *
import json

from qiskit import QuantumCircuit

class QuantumMarkovChain:
    def __init__(self, cir, err_model, bound=1e7) -> None:
        if isinstance(cir, str):
            try:
                cir=QuantumCircuit.from_qasm_file(cir)
            except Exception as e:
                print(e)
        assert(isinstance(cir, QuantumCircuit))
        self.cir = cir
        # check err_model
        self.err_model = err_model
        for err in self.err_model:
            if err.pos[0] == -1:
                # Maybe several position -1!
                err.pos[0] = len(cir)
        self.bound = bound

class AtomicProposition:
    """A tuple <M,I>, M is Projector, and I is the interval of the probability under measurement M"""
    def __init__(self, subspace, iterval) -> None:
        self.subspace = subspace
        self.iterval = iterval

class QuantumTransitionSystem:
    def __init__(self) -> None:
        pass
    def fromCircuit(self, cir, err_model=None):
        pass
    def fromSeqCircuit(self, seqCir):
        pass
    def fromKripkeStr(self, kripkeStr):
        pass
    def fromKripkeFile(self, filename):
        pass
    def toQuantumKripke(self):
        pass

"""
    char *content = "{\n"
              "  \"states\": [\n"
              "    { \"name\": \"s0\", \"labels\": [ \"a\", \"bar\" ] },\n"
              "    { \"name\": \"s1\", \"labels\": [ \"a\" ] },\n"
              "    { \"name\": \"s2\", \"labels\": [ \"bar\" ] },\n"
              "    { \"name\": \"s3\", \"labels\": [ \"a\", \"bar\" ] }\n"
              "  ],\n"
              "  \"initialStates\": [ \"s0\" ],\n"
              "  \"relations\": [ [ \"s0\", \"s1\" ], [ \"s1\", \"s2\" ], [ \"s1\", \"s0\" ], [ \"s2\", \"s3\" ], [ \"s2\", \"s0\" ], [ \"s3\", \"s3\" ] ]\n"
              "}";
"""
