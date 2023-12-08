from qiskit import QuantumCircuit, BasicAer, transpile
import numpy as np
import random

gen_cir = {}
gen_sp = {}
def register_function(func_dict, name):
    def decorator(func):
        func_dict[name] = func
        return func
    return decorator

def get_qasm(qc:QuantumCircuit)->QuantumCircuit:
    backend = BasicAer.get_backend('qasm_simulator')
    qc = transpile(qc,backend)

    return qc

def mct_ancilla_m(qc:QuantumCircuit,work_qubits:list,ancilla_qubits:list,target_qubit:int)->QuantumCircuit:
    n = len(work_qubits)
    if n == 1:
        qc.cx(work_qubits[0],target_qubit)
        return qc
    if  n == 2:
        qc.ccx(work_qubits[0],work_qubits[1],target_qubit)
        return qc

    ancilla_qubits.append(target_qubit)
    j = len(ancilla_qubits)-n+1
    sequence = [[work_qubits[i],ancilla_qubits[i-2+j],ancilla_qubits[i-1+j]] for i in range(2,n)]
    ancilla_qubits.pop()

    for c1,c2,t in sequence[::-1]: qc.ccx(c1,c2,t)
    qc.ccx(work_qubits[0],work_qubits[1],ancilla_qubits[j])
    for c1,c2,t in sequence: qc.ccx(c1,c2,t)

    for c1,c2,t in sequence[-2::-1]: qc.ccx(c1,c2,t)
    qc.ccx(work_qubits[0],work_qubits[1],ancilla_qubits[j])
    for c1,c2,t in sequence[:-1:]: qc.ccx(c1,c2,t)
    return qc

def mct_ancilla_1(qc:QuantumCircuit,controlled_qubits:list,ancilla_qubit:int,target_qubit:int)->QuantumCircuit:
    n = len(controlled_qubits)
    if n == 0:
        qc.x(target_qubit)
        return qc
    if n == 1:
        qc.cx(controlled_qubits[0],target_qubit)
        return qc
    elif n == 2:
        qc.ccx(controlled_qubits[0],controlled_qubits[1],target_qubit)
        return qc
    n = n//2
    part_a = controlled_qubits[:n+1]
    part_b = controlled_qubits[n+1:]

    qc = mct_ancilla_m(qc,part_a,part_b+[target_qubit],ancilla_qubit)
    qc = mct_ancilla_m(qc,part_b+[ancilla_qubit],part_a,target_qubit)
    qc = mct_ancilla_m(qc,part_a,part_b+[target_qubit],ancilla_qubit)
    qc = mct_ancilla_m(qc,part_b+[ancilla_qubit],part_a,target_qubit)

    return qc

@register_function(gen_cir, "rand_clif")
def rand_clif(qubit_num:int, size=5000)->QuantumCircuit:
    clifqc = QuantumCircuit(qubit_num)
    gate_type = random.randint(1,6)
    for _ in range(size):
        idx0 = random.randint(0, qubit_num-1)
        idx1 = random.randint(0, qubit_num-1)
        if idx0 == idx1:
            idx1 = (idx1+1) % qubit_num
        if gate_type == 1: # cx
            clifqc.cx(idx0, idx1)
        elif gate_type == 2: # h
            clifqc.h(idx0)
        elif gate_type == 3: # s
            clifqc.s(idx0)
        elif gate_type == 4: # x
            clifqc.x(idx0)
        elif gate_type == 5: # y
            clifqc.y(idx0)
        elif gate_type == 6: # z
            clifqc.z(idx0)
    return clifqc

@register_function(gen_cir, "grover")
def grover(qubit_num:int)->QuantumCircuit:
    """Generate a simple circuit of Grover's algorithm

    Args:
        qubit_num (int): The size of the search space is 2^(qubit_num-1). the second-to-last qubit(qubit_num - 1) is oracle workspace and the last qubit (qubit_num) is an ancilla qubit.

    Returns:
        QuantumCircuit: A quantum circuit implementing Grover's algorithm for the specified number of qubits.

    Here the oracle operator we used is f(|111..>) = 1, that is an mct gate with an ancilla.

    This function uses a custom multi-controlled Toffoli gate implementation called `mct_ancilla_1`, which takes in a list of control qubits, an ancilla qubit, and a target qubit.
    """
    search_space = list(range(qubit_num - 1))
    oracle_work = qubit_num - 1

    if qubit_num < 3: raise(ValueError("search space too small"))
    elif qubit_num == 3 :
        ancilla_qubit = None
        grover_qc = QuantumCircuit(qubit_num)
    else:
        ancilla_qubit = qubit_num
        grover_qc = QuantumCircuit(qubit_num + 1)

    # oracle operator
    grover_qc = mct_ancilla_1(grover_qc,search_space,ancilla_qubit,oracle_work)

    # diffusion operator
    grover_qc.h(search_space)
    grover_qc.x(search_space)

    #  mcz
    grover_qc.h(search_space[-1])
    grover_qc = mct_ancilla_1(grover_qc,search_space[:-1],ancilla_qubit,search_space[-1])
    grover_qc.h(search_space[-1])

    grover_qc.x(search_space)
    grover_qc.h(search_space)

    return grover_qc

@register_function(gen_sp, "grover_main")
def grover_0_main(qubit_num:int)->QuantumCircuit:
    """
    Generate a simple circuit of Grover's algorithm
    But qubit 0 is |0>
    """
    search_space = list(range(qubit_num - 1))
    oracle_work = qubit_num - 1

    if qubit_num < 3: raise(ValueError("search space too small"))
    elif qubit_num == 3:
        ancilla_qubit = None
        grover_qc = QuantumCircuit(qubit_num)
    else:
        ancilla_qubit = qubit_num
        grover_qc = QuantumCircuit(qubit_num + 1)

    # delete oracle operator
    # diffusion operator
    grover_qc.h(search_space)
    grover_qc.x(search_space)

    #  mcz
    grover_qc.h(search_space[-1])
    grover_qc = mct_ancilla_1(grover_qc,search_space[:-1],ancilla_qubit,search_space[-1])
    grover_qc.h(search_space[-1])

    grover_qc.x(search_space)
    grover_qc.h(search_space)

    return grover_qc

@register_function(gen_sp, "grover_add")
def grover_0_add(qubit_num:int)->QuantumCircuit:
    """
    Generate a simple circuit of Grover's algorithm
    But qubit 0 is |1>
    """

    search_space = list(range(qubit_num - 1))
    oracle_work = qubit_num - 1

    if qubit_num < 3: raise(ValueError("search space too small"))
    elif qubit_num == 3:
        ancilla_qubit = None
        grover_qc = QuantumCircuit(qubit_num)
    else:
        ancilla_qubit = qubit_num
        grover_qc = QuantumCircuit(qubit_num + 1)

    # oracle operator
    grover_qc = mct_ancilla_1(grover_qc,search_space[1:],ancilla_qubit,oracle_work)

    return grover_qc

@register_function(gen_cir, "qft")
def qft(qubit_num:int)->QuantumCircuit:
    """
    Generate a circuit for the quantum Fourier transform (QFT) on a given number of qubits.

    Args:
        qubit_num (int): The number of qubits to apply the QFT to.

    Returns:
        QuantumCircuit: A quantum circuit implementing the QFT on the specified number of qubits.

    The QFT applies a series of Hadamard and controlled rotation gates to the input qubits in order to transform them into their frequency-domain representation. The circuit generated by this function applies the following operations:
    1. Apply a Hadamard gate to each qubit.
    2. For each qubit, apply a controlled rotation gate to each target qubit with a higher index, where the rotation angle is pi/2^(tar-qubit).

    Note that this implementation uses the cp gate for the controlled rotations, which can lead to issues with numerical precision for certain angles. The code includes commented-out lines for using the Decimal library to address this issue if desired.
    """
    qft_qc=QuantumCircuit(qubit_num)
    for qubit in range(qubit_num):
        qft_qc.h(qubit)
        for target_qubit in range(qubit+1,qubit_num):
            theta=np.pi/2**(target_qubit-qubit)
            # theta=Decimal(theta).quantize(Decimal(1)/Decimal(10**16))
            if theta==0:
                print(theta)
            qft_qc.cp(theta,qubit,target_qubit)

    return qft_qc

@register_function(gen_sp, "qft_main")
def qft_0_main(qubit_num:int)->QuantumCircuit:
    """
    Generate a circuit for the quantum Fourier transform (QFT) on a given number of qubits.
    But qubit 0 is 0.
    """
    qft_qc=QuantumCircuit(qubit_num)

    for qubit in range(1,qubit_num):
        qft_qc.h(qubit)
        for target_qubit in range(qubit+1,qubit_num):
            theta= np.pi/2**(target_qubit-qubit)
            # theta=Decimal(theta).quantize(Decimal(1)/Decimal(10**16))
            if theta == 0:
                print(theta)
            qft_qc.cp(theta,qubit,target_qubit)

    return qft_qc

@register_function(gen_sp, "qft_add")
def qft_0_add(qubit_num:int)->QuantumCircuit:
    """
    Generate a circuit for the quantum Fourier transform (QFT) on a given number of qubits.
    But qubit 0 is 1.
    """
    qft_qc=QuantumCircuit(qubit_num)

    for target_qubit in range(1,qubit_num):
        theta=np.pi/2**(target_qubit+1)
            # theta=Decimal(theta).quantize(Decimal(1)/Decimal(10**16))
        if theta==0:
            print(theta)
        qft_qc.p(theta,target_qubit)
    return qft_qc

def generate_a_string(qubit_num:int, prob=1.0)->str :
    """
        generate a random binary string as a hidden bit string
    """
    answer = ["1" if random.random() <= prob else "0" for _ in range(qubit_num)]

    return "".join(answer)

@register_function(gen_cir, "bv")
def bv(qubit_num:int, hiddenString="")->QuantumCircuit:
    """
    Generate a quantum circuit for the Bernstein-Vazirani algorithm.

    Args:
        qubit_num (int): The number of qubits in the circuit.
        hiddenString (str): A binary string from function 'generate_a_string' consisting of "0" and "1", representing the hidden bit string that the algorithm is trying to determine.

    Returns:
        QuantumCircuit: A quantum circuit implementing the Bernstein-Vazirani algorithm on the specified number of qubits.

    The Bernstein-Vazirani algorithm is a quantum algorithm that finds the hidden bit string in a black box function that computes f(x) = s*x for some bit string s. The algorithm has the following steps:
    1. Initialize the input qubits to the |0> state.
    2. Apply a Hadamard gate to each qubit except the last one.
    3. Apply an X gate to the last qubit and then a Hadamard gate to it.
    4. Query the black box function by applying a series of controlled-not (CX) gates with the last qubit as the control and each other qubit as the target if the corresponding bit in the hidden string is 1.
    5. Apply a Hadamard gate to each qubit except the last one.

    This function generates a quantum circuit that implements the Bernstein-Vazirani algorithm for a given number of qubits and a given hidden bit string. The circuit is returned as a QuantumCircuit object.
    """
    bv_qc = QuantumCircuit(qubit_num)
    cont_space= list(range(qubit_num-1))
    target_qubit = qubit_num-1

    if not hiddenString:
        hiddenString = generate_a_string(qubit_num-1)
    bv_qc.h(cont_space)

    bv_qc.x(target_qubit)
    bv_qc.h(target_qubit)

    hiddenString = hiddenString[::-1]
    for qubit in range(len(hiddenString)):
        if hiddenString[qubit] == "1":
            bv_qc.cx(qubit, target_qubit)
    hiddenString = hiddenString[::-1]

    # Apply Hadamard gates after querying the oracle
    bv_qc.h(cont_space)

    return bv_qc

@register_function(gen_sp, "bv_main")
def bv_0_main(qubit_num:int, hiddenString="")->QuantumCircuit:
    """
    Generate a quantum circuit for the Bernstein-Vazirani algorithm.
    but qubit 0 is 0
    """
    if not hiddenString :
        hiddenString = generate_a_string(qubit_num-1)
    bv_qc = QuantumCircuit(qubit_num)
    cont_space= list(range(1,qubit_num-1))
    target_qubit = qubit_num-1

    bv_qc.h(cont_space)

    bv_qc.x(target_qubit)
    bv_qc.h(target_qubit)

    hiddenString = hiddenString[::-1]
    for qubit in range(len(hiddenString[1:])):
        if hiddenString[qubit+1] == "1":
            bv_qc.cx(qubit, target_qubit)
    hiddenString = hiddenString[::-1]

    # Apply Hadamard gates after querying the oracle
    bv_qc.h(cont_space)

    return bv_qc

@register_function(gen_sp, "bv_add")
def bv_0_add(qubit_num:int, hiddenString="")->QuantumCircuit:
    """
    Generate a quantum circuit for the Bernstein-Vazirani algorithm.
    But qubit 0 is 1
    """
    if not hiddenString :
        hiddenString = generate_a_string(qubit_num-1)

    bv_qc = QuantumCircuit(qubit_num)
    cont_space= list(range(1,qubit_num-1))
    target_qubit = qubit_num-1


    bv_qc.x(target_qubit)
    bv_qc.h(target_qubit)

    hiddenString = hiddenString[::-1]
    if hiddenString[0] == "1":
        bv_qc.x(target_qubit)
    hiddenString = hiddenString[::-1]

    # Apply Hadamard gates after querying the oracle

    return bv_qc

@register_function(gen_cir, "ghz")
def ghz(qubit_num:int)->QuantumCircuit:
    """
    Returns a QuantumCircuit that prepares a qubit state in the GHZ (Greenberger-Horne-Zeilinger) state, which is an entangled state of `qubit_num` qubits where each qubit is in a superposition of |0> and |1>, and all qubits
    are correlated with each other.

    Args:
        qubit_num (int): The number of qubits to include in the GHZ state

    Returns:
        QuantumCircuit: A QuantumCircuit object that prepares a GHZ state with the specified number of qubits.
    """
    ghz_qc = QuantumCircuit(qubit_num)
    ghz_qc.h(0)
    for qubit in range(1,qubit_num):
        ghz_qc.cx(0,qubit)
    return ghz_qc

@register_function(gen_sp, "ghz_main")
def ghz_0_main(qubit_num:int)->QuantumCircuit:
    """
    Returns a QuantumCircuit that prepares a qubit state in the GHZ (Greenberger-Horne-Zeilinger) state.
    But qubit 0 is 0
    """
    ghz_qc = QuantumCircuit(qubit_num)
    return ghz_qc

@register_function(gen_sp, "ghz_add")
def ghz_0_add(qubit_num:int)->QuantumCircuit:
    """
    Returns a QuantumCircuit that prepares a qubit state in the GHZ (Greenberger-Horne-Zeilinger) state.
    But qubit 0 is 1
    """
    ghz_qc = QuantumCircuit(qubit_num)
    ghz_qc.x(range(1,qubit_num))
    return ghz_qc

@register_function(gen_cir, "qrw")
def qrw(qubit_num:int)->QuantumCircuit:
    """
    Returns a QuantumCircuit that implements a quantum random walk with `qubit_num` qubits.In this function, the random walk is in a (2^qubit_num) cycle.

    Args:
        qubit_num (int): The number of qubits to use for the random walk. Must be greater than or equal to 2.

    Raises:
        ValueError: If `qubit_num` is less than 2.

    Returns:
        QuantumCircuit: A QuantumCircuit object that implements the quantum random walk with the specified number of qubits.
    """

    if qubit_num < 2:
        raise(ValueError("qubit not enough"))
    elif qubit_num == 2:
        qrw_qc = QuantumCircuit(qubit_num+1)
    else:
        qrw_qc = QuantumCircuit(qubit_num+2)

    walk_space = list(range(qubit_num))
    coin_qubit = qubit_num
    ancilla_qubit = qubit_num + 1

    qrw_qc.h(coin_qubit)

    # inc
    for target_qubit in walk_space[:-1]:
        qrw_qc=mct_ancilla_1(qrw_qc,walk_space[target_qubit+1:]+[coin_qubit],ancilla_qubit,target_qubit)
    qrw_qc.x(walk_space[-1])
    # dec
    qrw_qc.x(coin_qubit)
    for target_qubit in walk_space[-2::-1]:
        qrw_qc=mct_ancilla_1(qrw_qc,walk_space[target_qubit+1:]+[coin_qubit],ancilla_qubit,target_qubit)
    qrw_qc.x(coin_qubit)
    return qrw_qc

@register_function(gen_cir, "qrw_noise")
def qrw_noise(qubit_num:int)->QuantumCircuit:
    """
    Returns a QuantumCircuit that implements a quantum random walk with `qubit_num` qubits.In this function, the random walk is in a (2^qubit_num) cycle.

    Args:
        qubit_num (int): The number of qubits to use for the random walk. Must be greater than or equal to 2.

    Raises:
        ValueError: If `qubit_num` is less than 2.

    Returns:
        QuantumCircuit: A QuantumCircuit object that implements the quantum random walk with the specified number of qubits.
    """

    if qubit_num < 2:
        raise(ValueError("qubit not enough"))
    elif qubit_num == 2:
        qrw_qc = QuantumCircuit(qubit_num+1)
    else:
        qrw_qc = QuantumCircuit(qubit_num+2)

    walk_space = list(range(qubit_num))
    coin_qubit = qubit_num
    ancilla_qubit = qubit_num + 1

    # noise
    # qrw_qc.h(coin_qubit)
    # qrw_qc.x(coin_qubit)

    # inc
    for target_qubit in walk_space[:-1]:
        qrw_qc=mct_ancilla_1(qrw_qc,walk_space[target_qubit+1:]+[coin_qubit],ancilla_qubit,target_qubit)
    qrw_qc.x(walk_space[-1])
    # dec
    qrw_qc.x(coin_qubit)
    for target_qubit in walk_space[-2::-1]:
        qrw_qc=mct_ancilla_1(qrw_qc,walk_space[target_qubit+1:]+[coin_qubit],ancilla_qubit,target_qubit)
    qrw_qc.x(coin_qubit)
    return qrw_qc

@register_function(gen_sp, "qrw_0")
def qrw_coin_0(qubit_num:int)->QuantumCircuit:
    """
    Returns a QuantumCircuit that implements a quantum random walk with `qubit_num` qubits.In this function, the random walk is in a (2^qubit_num) cycle.
    But coin qubit is 0
    """

    if qubit_num < 2:
        raise(ValueError("qubit not enough"))
    elif qubit_num == 2:
        qc = QuantumCircuit(qubit_num+1)
    else:
        qc = QuantumCircuit(qubit_num+2)

    walk_space = list(range(qubit_num))
    coin_qubit = qubit_num
    ancilla_qubit = qubit_num + 1

    # dec
    qc.x(walk_space[-1])
    for target_qubit in walk_space[-2::-1]:
        qc=mct_ancilla_1(qc,walk_space[target_qubit+1:],ancilla_qubit,target_qubit)
    return qc

@register_function(gen_sp, "qrw_1")
def qrw_coin_1(qubit_num:int)->QuantumCircuit:
    """
    Returns a QuantumCircuit that implements a quantum random walk with `qubit_num` qubits.In this function, the random walk is in a (2^qubit_num) cycle.
    But coin qubit 1 is 1
    """

    if qubit_num < 2:
        raise(ValueError("qubit not enough"))
    elif qubit_num == 2:
        qc = QuantumCircuit(qubit_num+1)
    else:
        qc = QuantumCircuit(qubit_num+2)

    walk_space = list(range(qubit_num))
    coin_qubit = qubit_num
    ancilla_qubit = qubit_num + 1

    # inc
    for target_qubit in walk_space[:-1]:
        qc=mct_ancilla_1(qc,walk_space[target_qubit+1:],ancilla_qubit,target_qubit)
    qc.x(walk_space[-1])

    return qc