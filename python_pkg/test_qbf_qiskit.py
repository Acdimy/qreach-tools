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

def p_coin(p: float, k: int) -> bool:
    backend = Aer.get_backend('aer_simulator')
    theta = 2 * np.arccos(np.sqrt(p))

    i = 0
    while i < k:
        # 每次创建新的电路
        circ = QuantumCircuit(2, 2)
        circ.initialize(0,0)
        circ.ry(theta, 0)
        circ.ry(theta, 1)
        circ.cx(0, 1)
        circ.h(0)
        circ.measure([0,1], [0,1])

        result = backend.run(circ, shots=1).result()
        counts = result.get_counts()
        measured = list(counts.keys())[0]  # e.g., '01', '10', '00', '11'

        # Qiskit qubit顺序：从右到左，bit[0]是右边
        bit0 = int(measured[-1])  # qubit 0
        bit1 = int(measured[-2])  # qubit 1

        if bit0 == bit1:
            continue
        elif bit0 == 1 and bit1 == 0:
            return False
        elif bit0 == 0 and bit1 == 1:
            i += 1
    return True

def p_coin_qiskit(p: float, k: int) -> bool:
    """
    Seems to have a bug, the figure of estimated probability is not correct.
    """
    qubits = QuantumRegister(3)
    midBits = ClassicalRegister(2, name='mid')
    ancBits = ClassicalRegister(1, name='anc')
    circ = QuantumCircuit(qubits, midBits, ancBits)
    theta = 2 * np.arccos(np.sqrt(p))
    circ.reset(0)
    circ.reset(1)
    circ.reset(2)
    for i in range(k):
        with circ.if_test((midBits, 0b00)):
            circ.ry(theta, 0)
            circ.ry(theta, 1)
            circ.cx(0, 1)
            circ.h(0)
            circ.measure(0,0)
            circ.measure(1,1)
        with circ.if_test((midBits, 0b00)):
            circ.x(2)
            circ.measure(2,2)
        with circ.if_test((midBits, 0b11)):
            circ.x(2)
            circ.measure(2,2)
        with circ.while_loop((ancBits, 0b1)):
            circ.reset(0)
            circ.reset(1)
            circ.reset(2)
            
            # Necessary to re-measure qubit 2 to reset the ancBits
            circ.measure(2,2)

            circ.ry(theta, 0)
            circ.ry(theta, 1)
            circ.cx(0, 1)
            circ.h(0)
            circ.measure(0,0)
            circ.measure(1,1)
            with circ.if_test((midBits, 0b00)):
                circ.x(2)
            with circ.if_test((midBits, 0b11)):
                circ.x(2)
            circ.measure(2,2)
        with circ.if_test((midBits, 0b01)):
            circ.reset(0)
            circ.reset(1)
            circ.reset(2)
            circ.measure(0,0)
            circ.measure(1,1)
    backend = Aer.get_backend('aer_simulator')
    result = backend.run(circ, shots=1).result()
    counts = result.get_counts()
    measured = list(counts.keys())[0]  # e.g., '01', '10', '00', '11'
    bit0 = int(measured[-1])  # qubit 0
    bit1 = int(measured[-2])  # qubit 1
    if bit0 == 1 and bit1 == 0:
        return False
    return True

# for p in [0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9], repeat the tests to get sufficient samples, and estimate the probability
total_results = []
for p in [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]:
    print(f"Testing for p = {p}")
    repeats = 1000
    prob_results = []
    for _ in range(repeats):
        probs, cdf = build_cdf()
        k = sample_from_distribution(cdf)
        result = p_coin_qiskit(p, k)
        prob_results.append(result)
    estimated_prob = sum(prob_results) / repeats
    total_results.append((p, estimated_prob))
    print(f"Estimated probability for p={p}: {estimated_prob:.4f}")
# draw a figure to show the results
import matplotlib.pyplot as plt
p_values, estimated_probs = zip(*total_results)
plt.plot(p_values, estimated_probs, marker='o')
plt.xlabel('p (Probability of Heads)')
plt.ylabel('Estimated Probability of Heads')
plt.title('Estimated Probability of Heads vs p')
plt.grid()
# save the figure
plt.savefig('estimated_probability_vs_p_qiskit1.png')

# probs, cdf = build_cdf()
# k = sample_from_distribution(cdf)
# result = p_coin(0.5, k)


print(f"Sampled k = {k}, Bernoulli result: {result}")
