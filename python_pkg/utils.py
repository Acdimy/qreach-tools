
def get_real_qubit_num(cir):
    """Calculate the real number of qubits of a circuit"""
    return cir.num_qubits
    # gates=cir.data
    # q=0
    # for k in range(len(gates)):
    #     q=max(q,max([qbit.index for qbit in gates[k][1]]))
    # return q+1
