import pyqreach
### Verifiable quantum secret sharing
from qiskit import QuantumCircuit, QuantumRegister, ClassicalRegister, transpile
from qiskit.quantum_info import Statevector
from qiskit_aer import Aer
import numpy as np
import math
import random
from time import time

from parse_qiskit import *

def q_k(k):
    comb = math.comb(2 * k, k)  # C(2k, k)
    denom = (2 * k - 1) * (2 ** (2 * k))
    return comb / denom

def build_cdf(threshold=1e-10, max_k=100):
    probs = []
    cdf = []
    total = 0.0
    for k in range(1,max_k):
        prob = q_k(k)
        probs.append(prob)
        total += prob
        cdf.append(total)
        if 1 - total < threshold:
            break
    cdf = [x / total for x in cdf]
    return probs, cdf

def sample_from_distribution(cdf):
    r = random.random()
    for k, p in enumerate(cdf):
        if r < p:
            return k+1
    return len(cdf) - 1  # fallback in case of numerical issue


def p_coin(circ: QuantumCircuit, p: float, k: int, clbits: ClassicalRegister) -> bool:
    theta = 2 * np.arccos(np.sqrt(p))
    circ.reset(0)
    circ.reset(1)
    for i in range(k):
        with circ.if_test((clbits, 0b00)):
            circ.ry(theta, 0)
            circ.ry(theta, 1)
            circ.cx(0, 1)
            circ.h(0)
            circ.measure(0,0)
            circ.measure(1,1)
        with circ.switch(clbits) as case:
            with case(0b00, 0b11):
                circ.reset(0)
                circ.reset(1)
            with case(0b10):
                pass
            with case(0b11):
                circ.reset(0)
                circ.reset(1)
    return True

def p_coin_qiskit(circ: QuantumCircuit, p: float, k: int, clbits: ClassicalRegister):
    theta = 2 * np.arccos(np.sqrt(p))
    circ.reset(0)
    circ.reset(1)
    circ.reset(2)
    for i in range(k):
        with circ.if_test((clbits[:2], 0b00)):
            circ.ry(theta, 0)
            circ.ry(theta, 1)
            circ.cx(0, 1)
            circ.h(0)
            circ.measure(0,0)
            circ.measure(1,1)
        with circ.if_test((clbits[:2], 0b00)):
            circ.x(2)
        with circ.if_test((clbits[:2], 0b11)):
            circ.x(2)
        circ.measure(2,2)
        with circ.while_loop((clbits[2], 0b1)):
            circ.reset(0)
            circ.reset(1)
            circ.reset(2)
            circ.ry(theta, 0)
            circ.ry(theta, 1)
            circ.cx(0, 1)
            circ.h(0)
            circ.measure(0,0)
            circ.measure(1,1)
            with circ.if_test((clbits[:2], 0b00)):
                circ.x(2)
            with circ.if_test((clbits[:2], 0b11)):
                circ.x(2)
            circ.measure(2,2)
        with circ.switch(clbits[:2]) as case:
            with case(0b10):
                pass
            with case(0b11):
                circ.reset(0)
                circ.reset(1)
    return True