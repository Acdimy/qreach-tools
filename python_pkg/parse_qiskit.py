import pyqreach
from qiskit import QuantumCircuit
from graphviz import Digraph
from math import floor
from math import ceil, log2

# Get the binary representation of a number
def getBinary(num, length):
    return format(num, '0' + str(length) + 'b')

def prepare_5perfect_code(circ, idx):
    assert(len(idx) == 5)
    circ.z(idx[0])

    circ.h(idx[1])
    circ.cz(idx[1],idx[2])
    circ.cz(idx[1],idx[4])
    circ.cx(idx[1],idx[0])

    circ.h(idx[4])
    circ.cz(idx[4],idx[3])
    circ.cz(idx[4],idx[1])
    circ.cx(idx[4],idx[0])

    circ.h(idx[3])
    circ.cz(idx[3],idx[2])
    circ.cz(idx[3],idx[0])
    circ.cx(idx[3],idx[4])

    circ.h(idx[2])
    circ.cz(idx[2],idx[1])
    circ.cz(idx[2],idx[4])
    circ.cx(idx[2],idx[3])

def decode_5perfect_code(circ, idx):
    pass

def prepare_steane_code(circ, idx):
    """
    Prepare a Steane code in the given quantum circuit.
    
    Args:
        circ (QuantumCircuit): The quantum circuit to modify.
        idx (list): List of qubit indices to apply the Steane code.
    """
    # Apply the Steane code preparation steps
    assert(len(idx) == 7)
    circ.h(idx[0])
    circ.h(idx[1])
    circ.h(idx[2])
    circ.cx(idx[0], idx[3])
    circ.cx(idx[6], idx[4])
    circ.cx(idx[6], idx[5])
    circ.cx(idx[1], idx[3])
    circ.cx(idx[0], idx[5])
    circ.cx(idx[2], idx[3])
    circ.cx(idx[1], idx[4])
    circ.cx(idx[0], idx[6])
    circ.cx(idx[2], idx[4])
    circ.cx(idx[1], idx[6])
    circ.cx(idx[2], idx[5])

def decode_steane_code(circ, idx):
    pass

def parse_qiskit(qc: QuantumCircuit, loc0_idx=0) -> pyqreach.TransitionSystem:
    """
    Parse a Qiskit QuantumCircuit into a pyqreach TransitionSystem.
    
    Args:
        qc (QuantumCircuit): The quantum circuit to parse.
        
    Returns:
        pyqreach.TransitionSystem: The corresponding transition system.
    """
    qnum = qc.num_qubits
    ts = pyqreach.TransitionSystem()
    ts.setInitLocation(0)

    # Add locations based on the number of qubits
    loc_list = []
    loc0 = pyqreach.Location(qnum, loc0_idx)
    ts.addLocation(loc0)
    loc_list.append(loc0)

    # Define operations based on the gates in the circuit
    loc_idx = loc0_idx
    for _,gate in enumerate(qc.data):
        op_name = gate[0].name
        qubits = [q._index for q in gate[1]]
        cbits = [c._index for c in gate[2]] if gate[2] else []
        loc = pyqreach.Location(qnum, loc_idx+1)
        ts.addLocation(loc)
        loc_idx += 1
        if op_name == 'h':
            op = pyqreach.QOperation("H", qnum, qubits, [])
        elif op_name == 'id':
            op = pyqreach.QOperation("I", qnum, qubits, [])
        elif op_name == 'x':
            op = pyqreach.QOperation("X", qnum, qubits, [])
        elif op_name == 'y':
            op = pyqreach.QOperation("Y", qnum, qubits, [])
        elif op_name == 'z':
            op = pyqreach.QOperation("Z", qnum, qubits, [])
        elif op_name == 's':
            op = pyqreach.QOperation("S", qnum, qubits, [])
        elif op_name == 't':
            op = pyqreach.QOperation("T", qnum, qubits, [])
        elif op_name == 'cx':
            op = pyqreach.QOperation("CX", qnum, qubits, [])
        elif op_name == 'cz':
            op = pyqreach.QOperation("CZ", qnum, qubits, [])
        elif op_name == 'measure':
            print("Measurement operation detected")
            opm0 = pyqreach.QOperation("meas0", qnum, qubits, [])
            ts.addRelation(loc_idx-1, loc_idx, opm0)
            loc_m1 = pyqreach.Location(qnum, loc_idx+1)
            ts.addLocation(loc_m1)
            loc_idx += 1
            opm1 = pyqreach.QOperation("meas1", qnum, qubits, [])
            ts.addRelation(loc_idx-1, loc_idx, opm1)
            loc_merge = pyqreach.Location(qnum, loc_idx+1)
            op = pyqreach.QOperation("I", qnum, qubits, [])
            ts.addLocation(loc_merge)
            loc_idx += 1
            ts.addRelation(loc_idx-2, loc_idx, op)
        else:
            print(f"Unsupported gate: {op_name}")
        # Add the operation to the transition system
        ts.addRelation(loc_idx-1, loc_idx, op)  # Simplified relation for demonstration
    return ts

def visualize_transition_system(ts: pyqreach.TransitionSystem, filename='transition_system'):
    """
    Visualize the transition system using Graphviz.
    
    Args:
        ts (pyqreach.TransitionSystem): The transition system to visualize.
        filename (str): The name of the output file.
    """
    dot = Digraph(comment='Transition System')
    
    for loc in ts.Locations:
        dot.node(str(loc.idx), label=str(loc.idx))
    
        for post in loc.postLocations:
            rel = ts.getRelationName(loc.idx, post)
            dot.edge(str(loc.idx), str(post), label=str(rel))
    
    dot.render(filename, format='png', cleanup=True)

def applyFinalMeasurement(ts: pyqreach.TransitionSystem, PauliString: str, qlist: list, qnum: int) -> list:
    """
    Apply a final measurement to the transition system.
    
    Args:
        ts (pyqreach.TransitionSystem): The transition system to modify.
        PauliString (str): The Pauli string representing the measurement.
        Each Pauli is a local measurement on a qubit.
        qlist (list): List of locations to apply the measurement to.
    """
    currLocNum = ts.getLocationNum()
    nontrivialPauli = 0
    numMeasuredQubits = sum(1 for p in PauliString if p != 'I')
    indexInPauli = [i for i, p in enumerate(PauliString) if p != 'I']
    totalNewLocation = 2**(numMeasuredQubits + 1) - 2
    # for pauli in PauliString:
    #     if pauli != 'I':
    #         nontrivialPauli += 1
    #         if pauli == 'X':
    #             for i in range(2**(nontrivialPauli - 1)):
    #                 loc = pyqreach.Location(qnum, 0)
    #                 ts.addLocation(loc)
    for i in range(totalNewLocation):
        loc = pyqreach.Location(qnum, 0)
        ts.addLocation(loc)
    for layer in range(numMeasuredQubits): # traverse each parent layer
        currLayerSize = 2**(layer)
        parentStart = currLocNum - 1 + 2**layer - 1
        childLayerStart = currLocNum - 1 + 2**(layer+1) - 1
        pauliIndex = indexInPauli[layer]
        pauli = PauliString[pauliIndex]
        qubit = qlist[pauliIndex]
        for p in range(currLayerSize):  # traverse each parent node
            if pauli == 'X':
                oph = pyqreach.QOperation("H", qnum, [qubit], [])
                ts.addLocation(pyqreach.Location(qnum, 0))  # add a new location for Hadamard
                tmpHloc = ts.getLocationNum() - 1
                ts.addRelation(parentStart + p, tmpHloc, oph)
            for outcomeBit in [0, 1]:
                parentLoc = parentStart + p
                childLoc = childLayerStart + 2*p
                if pauli == 'Z':
                    measOp = pyqreach.QOperation("meas0" if outcomeBit == 0 else "meas1", qnum, [qubit], [])
                    ts.addRelation(parentLoc, childLoc + outcomeBit, measOp)
                elif pauli == 'X':
                    measOp = pyqreach.QOperation("meas0" if outcomeBit == 0 else "meas1", qnum, [qubit], [])
                    ts.addRelation(tmpHloc, childLoc + outcomeBit, measOp)

    resultList = []
    for j in range(2**(numMeasuredQubits)):
        resultList.append(currLocNum + 2**(numMeasuredQubits) - 2 + j)
    return resultList


def applyMeasureAndReset(ts: pyqreach.TransitionSystem, PauliString: str, qlist: list, qnum: int) -> list:
    """
    Apply a measurement and reset operation to the transition system.
    
    Args:
        ts (pyqreach.TransitionSystem): The transition system to modify.
        PauliString (str): The Pauli string representing the measurement.
        qlist (list): List of qubit indices to apply the measurement to.
        qnum (int): Total number of qubits.
    
    Returns:
        list: List of terminal locations corresponding to each measurement outcome.
    """
    resultList = applyFinalMeasurement(ts, PauliString, qlist, qnum)
    # From resultList to original locations. Note that PauliString may contains 'I', which means no measurement on that qubit. Just change qlist
    qlist_meas = [qlist[i] for i, p in enumerate(PauliString) if p != 'I']

    for i, loc in enumerate(resultList):
        resultStr = getBinary(i, len(qlist_meas))
        currLoc = loc
        for j, bit in enumerate(resultStr):
            if bit == '1':
                resetOp = pyqreach.QOperation("X", qnum, [qlist_meas[j]], [])
                ts.addLocation(pyqreach.Location(qnum, 0))
                ts.addRelation(currLoc, ts.getLocationNum() - 1, resetOp)
                currLoc = ts.getLocationNum() - 1
            resultList[i] = currLoc
    return resultList


def analyze_pruned_binary_tree(keepBitstrings: list[str]) -> tuple[int, list[list[str]]]:
    """
    Given a set of binary strings representing success paths,
    build a minimal pruned binary decision tree and return:
        - total node count
        - each layer's node list (by bit prefix)
    """
    keepSet = set(keepBitstrings)
    maxDepth = len(keepBitstrings[0])
    assert all(len(k) == maxDepth for k in keepBitstrings)

    from collections import deque, defaultdict

    visited = set()
    layers = defaultdict(list)

    queue = deque()
    queue.append("")  # root prefix

    while queue:
        prefix = queue.popleft()
        if prefix in visited:
            continue
        visited.add(prefix)
        layers[len(prefix)].append(prefix)

        if len(prefix) < maxDepth:
            left = prefix + "0"
            right = prefix + "1"
            if any(k.startswith(left) for k in keepSet):
                queue.append(left)
            if any(k.startswith(right) for k in keepSet):
                queue.append(right)

    totalNodes = len(visited)
    layerList = [sum([len(layers[i]) for i in range(d+1)]) for d in range(maxDepth + 1)]

    return totalNodes, layerList


def applySelectiveFinalMeasurement(
    ts: pyqreach.TransitionSystem,
    PauliString: str,
    qlist: list,
    qnum: int,
    entryLocList: list,
    keepBitstrings: list[str]
) -> tuple[list[list[int]], int]:
    """
    Apply a shared measurement tree from multiple entry locations, pruning paths
    based on allowed outcome bitstrings. Invalid paths go to a common error location.

    Args:
        ts (pyqreach.TransitionSystem): Transition system to modify.
        PauliString (str): Pauli measurement string (e.g., 'XZI').
        qlist (list[int]): List of qubit indices, same length as PauliString.
        qnum (int): Total number of qubits.
        entryLocList (list[int]): List of entry location indices to attach the measurement tree.
        keepBitstrings (list[str]): List of allowed outcome bitstrings (e.g., ['00001']).

    Returns:
        tuple: (List of list of valid final location indices per entry, error location index)
    """
    assert len(PauliString) == len(qlist)
    depth = sum(1 for p in PauliString if p != 'I')
    assert all(len(s) == depth for s in keepBitstrings)

    # Step 1: 构造测量子树结构（只构造一次）
    baseLocIndex = ts.getLocationNum() - 1

    # totalTreeNodes = 0
    # layerStartIndices = [baseLocIndex-1]
    # for d in range(1, depth + 1):
    #     layerStart = totalTreeNodes
    #     layerStartIndices.append(layerStart)
    #     totalTreeNodes += 2 ** d
    # for _ in range(totalTreeNodes*len(entryLocList)):
    #     ts.addLocation(pyqreach.Location(qnum, 0))

    totalTreeNodes, layerStartIndices = analyze_pruned_binary_tree(keepBitstrings)
    assert(len(layerStartIndices) == depth + 1)
    totalTreeNodes -= 1
    for _ in range(totalTreeNodes*len(entryLocList)):
        ts.addLocation(pyqreach.Location(qnum, 0))

    # Step 2: 构造 error location
    errorLoc = ts.getLocationNum()
    ts.addLocation(pyqreach.Location(qnum, 0))  # dummy error location

    # Step 3: 构造精简版测量树：只连符合 keepBitstrings 的路径
    indexInPauli = [i for i, p in enumerate(PauliString) if p != 'I']
    keepSet = set(keepBitstrings)
    

    for e, entry in enumerate(entryLocList):
        bitPrefixList = [['']]
        for d in range(depth):  # 每一层测量
            bitPrefixLayer = []
            # print(f"Processing layer {d} for entry {entry}...")
            pauliIndex = indexInPauli[d]
            qubit = qlist[pauliIndex]
            pauli = PauliString[pauliIndex]

            parentStart = (layerStartIndices[d-1] + baseLocIndex + e * totalTreeNodes) if d != 0 else entry
            childStart = layerStartIndices[d] + baseLocIndex + e * totalTreeNodes
            parentTrace = 0
            childTrace = 0

            for i in range(layerStartIndices[d] - layerStartIndices[d-1] if d != 0 else layerStartIndices[d]):
                bitPrefix = bitPrefixList[d][i]
                parentTrace += 1

                for outcomeBit in [0, 1]:
                    newBitstring = bitPrefix + str(outcomeBit)
                    parentLoc = parentStart + i
                    childLoc = childStart + childTrace

                    if any(s.startswith(newBitstring) for s in keepSet):
                        # 构造合法测量路径
                        if pauli == 'Z':
                            meas = pyqreach.QOperation("meas0" if outcomeBit == 0 else "meas1", qnum, [qubit], [])
                            ts.addRelation(parentLoc, childLoc, meas)
                        elif pauli == 'X':
                            had = pyqreach.QOperation("H", qnum, [qubit], [])
                            meas = pyqreach.QOperation("meas0" if outcomeBit == 0 else "meas1", qnum, [qubit], [])
                            tempLoc = ts.getLocationNum()
                            ts.addLocation(pyqreach.Location(qnum, 0))
                            ts.addRelation(parentLoc, tempLoc-1, had)
                            ts.addRelation(tempLoc-1, childLoc, meas)
                        childTrace += 1
                        bitPrefixLayer.append(newBitstring)
                    else:
                        # 剪枝路径连接到 error location，并带有与合法路径相反的测量操作
                        wrongMeas = pyqreach.QOperation("meas0" if outcomeBit == 0 else "meas1", qnum, [qubit], [])
                        ts.addRelation(parentLoc, errorLoc, wrongMeas)
                bitPrefixList.append(bitPrefixLayer)

    # Step 4: 将所有 entry location 连接到测量树的根节点
    validResultsPerEntry = []
    finalLayerStart = layerStartIndices[-1]
    for e, entry in enumerate(entryLocList):
        entryResults = []
        for j in range(2 ** depth):
            resultLoc = finalLayerStart + j + e * totalTreeNodes
            entryResults.append(resultLoc)
        validResultsPerEntry.append(entryResults)

    return validResultsPerEntry, errorLoc

# TODO: Randomly inject errors into the transition system
def injectRandomErrors(ts: pyqreach.TransitionSystem, errorRate: float):
    pass
