from checker import *
import sys

if __name__ == '__main__':
    info = "Expected arguments: experiment type ['qrw', 'rus', 'grover']; qubits; #input dimensions (1 or 2);\
    quantum errors which are set initially as 'amplitude dumping'"
    assert len(sys.argv) >= 3, info
    expr_name = sys.argv[1]
    num = sys.argv[2]

    input_dim = int(sys.argv[3])
    if input_dim == 1:
        input_str = ["0"]
    elif input_dim == 2:
        input_str = ["0","1"]
    else:
        raise RuntimeError(info)

    if len(sys.argv) >= 5:
        err_list = sys.argv[4] # Not available
    cir = "benchmark/"+expr_name+"_"+num+".qasm"


    t_start=time.time()
    if expr_name == "qrw":
        if len(sys.argv) >= 5:
            err_model = [Noise("ad", [-1,0], [0.5])]
        else:
            err_model = []
        qmc = QuantumMarkovChain(cir=cir, err_model=err_model)
        s_list = input_str
        qchecker = fromMarkovModel(qmc)
        qchecker = initWithStr(qchecker, str_list=s_list)
    elif expr_name == "rus":
        if num == "1":
            err_model = [Measurement("measure", [-1,0]), Measurement("measure", [-1,1])]
        else:
            err_model = [Measurement("measure", [-1,0])]
        qmc = QuantumMarkovChain(cir=cir, err_model=err_model)
        s_list = input_str
        qchecker = fromMarkovModel(qmc)
        qchecker = initWithStr(qchecker, str_list=s_list)
    elif expr_name == "grover":
        qmc = QuantumMarkovChain(cir=cir, err_model=[])
        # manually appoint the initial state because the grover walk's initialization is complicate
        s = "0"*(2*int(num)-2)+"1"
        qchecker = fromMarkovModel(qmc)
        qchecker.setState(str_padding(s))
        for i in range(int(num)):
            qchecker.h(i)
        qchecker.h(2*int(num)-2)
        qchecker.setProjectorFS()
    else:
        raise RuntimeError(info)

    reachable_dim = qchecker.reachability()

    print(f"use {time.time()-t_start} seconds")

    # if expr_name == "rus":
    #     channel = "Measure"
    # elif expr_name == "grover":
    #     channel = "Unitary"
    # elif len(sys.argv) < 5:
    #     channel = "Unitary"
    # else:
    #     channel = "Noise"
    # print("---Eval info: ", expr_name, "  #Qubit: ", qmc.cir.num_qubits, "  Input dim: ", input_dim, "  Channel: ", channel, " ---")

    # if qmc.cir.num_qubits <= 4:
    #     print("Output basis of reachable space:\n")
    #     qchecker.printProjector()

    print("Reachable dimensions: ", reachable_dim)
    # qchecker.printSize("state")
    qchecker.printSize("projector")
