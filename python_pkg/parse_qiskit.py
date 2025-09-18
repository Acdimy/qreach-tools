import pyqreach
from qiskit import QuantumCircuit, QuantumRegister, ClassicalRegister
from qiskit.circuit import Clbit
from graphviz import Digraph
from math import floor
from math import ceil, log2

# Get the binary representation of a number
def getBinary(num, length):
    return format(num, '0' + str(length) + 'b')



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
    loc0 = pyqreach.Location(qnum, loc0_idx)
    ts.addLocation(loc0)

    # Define operations based on the gates in the circuit
    loc_idx = loc0_idx
    curr_loc_list = [loc0_idx]
    for _,gate in enumerate(qc.data):
        op_name = gate[0].name
        qubits = [q._index for q in gate[1]]
        cbits = [c._index for c in gate[2]] if gate[2] else []
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
        elif op_name == 'reset':
            # The case if reset on all qubits
            pass
        elif op_name == 'if_else':
            # Two instruction blocks, one for if and one for else, double the curr_loc_list
            pass
        elif op_name == 'switch_case':
            pass
        elif op_name == 'for_loop':
            pass
        elif op_name == 'while_loop':
            pass
        elif op_name == 'measure':
            pass
        else:
            print(f"Unsupported gate: {op_name}")
        
        # Make new locations
        loc = pyqreach.Location(qnum, loc_idx+1)
        ts.addLocation(loc)
        loc_idx += 1
        # Add the operation to the transition system
        ts.addRelation(loc_idx-1, loc_idx, op)  # Simplified relation for demonstration
    return ts

# Global variables for mapping classical registers to indices
CLREG_ORDER = []

def init_parse(qc: QuantumCircuit):
    """
    Initialize the classical register mapping for parsing.
    
    Args:
        qc (QuantumCircuit): The quantum circuit to parse.
    """
    global CLREG_ORDER
    CLREG_ORDER = []
    for creg in qc.cregs:
        CLREG_ORDER.append(creg)
    # print(CLREG_ORDER, "Classical registers in the circuit")

def get_condition_info(cregs: list, condition: tuple) -> str:
    """
    cregs: list of classical registers (including bit number and name)
    condition: (ClassicalRegister, value)
    Returns: tuple(list of clbits indexes, list of clbits values in condition)
    """
    # Use global CLREG_ORDER instead of passing cregs
    global CLREG_ORDER
    clreg, value = condition
    clbits_idx, clbits_vals = [], []
    if isinstance(clreg, ClassicalRegister):
        clreg_name = clreg.name
        clreg_size = clreg.size
        clreg_reg_idx = next((i for i, reg in enumerate(CLREG_ORDER) if reg.name == clreg_name), None)
        reg_start_idx = sum([reg.size for reg in CLREG_ORDER[:clreg_reg_idx]])
        clbits_idx = [reg_start_idx + i for i in range(clreg_size)]
        # need reverse? Yes!
        clbits_vals = [int(bit) for bit in format(value, f'0{clreg_size}b')][::-1]
    elif isinstance(clreg, Clbit):
        clreg_name = clreg._register.name
        clreg_bit_idx = clreg._index
        clreg_reg_idx = next((i for i, reg in enumerate(CLREG_ORDER) if reg.name == clreg_name), None)
        reg_start_idx = sum([reg.size for reg in CLREG_ORDER[:clreg_reg_idx]]) + clreg_bit_idx
        clbits_idx = [reg_start_idx]
        clbits_vals = [int(value)]  # value is a single bit, so just convert it to int
    else:
        raise ValueError("Unsupported classical register type in condition")
    return clbits_idx, clbits_vals

def groupby_classical_aps(ts: pyqreach.TransitionSystem, currLocs: list) -> dict:
    """
    Group current locations by their classical APs.

    ts: Transition system (including all locations)
    currLocs: List of indexes of current locations. Each of them has a classical APs.
    Returns: A dictionary where keys are tuples of classical APs and values are lists of location indexes.
    """
    grouped = {}
    for loc_idx in currLocs:
        loc = ts.Locations[loc_idx]
        # classical proposition (cp) is originally ordered in the construction phase, so we directly use its string representation as the key.
        key = loc.cp.toString() if loc.cp else "empty"
        if key not in grouped:
            grouped[key] = []
        grouped[key].append(loc_idx)
    return grouped

def merge_locations(ts: pyqreach.TransitionSystem, currLocs: list, toMergeLocs: list, identifier: str="") -> list:
    """
    ts: Transition system to modify
    currLocs: List of current locations (including locations to be merged and not to be merged)
    toMergeLocs: List of locations to be merged. Merge them into a new location.
    Returns: The list of index of all new current locations. (replacing toMergeLocs with the new location)
    """
    assert len(toMergeLocs) > 0, "No locations to merge."
    assert len(currLocs) > 0, "Current locations list is empty."
    assert max(currLocs) < ts.getLocationNum(), "Current locations exceed the number of locations in the transition system."
    qNum = ts.Locations[currLocs[0]].qNum
    mg_location = pyqreach.Location(qNum)
    mg_location.setIdentifier(identifier)
    ts.addLocation(mg_location)
    newLocIdx = ts.getLocationNum() - 1
    newCurrLocs = [loc for loc in currLocs if loc not in toMergeLocs]
    newCurrLocs.append(newLocIdx)
    for l in toMergeLocs:
        # Merge classical APs
        loc = ts.Locations[l]
        for term in loc.cp.terms:
            if not ts.Locations[newLocIdx].find(term):
                ts.Locations[newLocIdx].appendClassicalAP(term)
        # Add relations from the merged locations to the new location
        ts.addRelation(l, newLocIdx, pyqreach.QOperation("I", qNum, [0], []))
    return newCurrLocs

def build_while_loop(qc: QuantumCircuit, qnum: int, clbits_idx, clbits_vals, whileStarter: list, startIdx: int, ts: pyqreach.TransitionSystem, identifier="" , abstractLevel=1) -> list:
    """
    Build a while loop in the transition system based on the Qiskit QuantumCircuit.
    
    Args:
        qc (QuantumCircuit): The loop body of the while loop.
        qnum (int): Total number of qubits.
        whileStarter (list): List of starting locations for the while loop.
        startIdx (int): The index to start parsing the while loop.
        ts (pyqreach.TransitionSystem): The transition system to modify.
        abstractLevel (int): Level of abstraction for merging locations.
        
    Returns:
        list: List of resulting locations (exiting locations) after processing the while loop.
    """
    # print("Calling build_while_loop")
    exitList = []
    afterLoopBodyList = []
    satisfyTerms = ts.Locations[startIdx].satisfyBit(clbits_idx, clbits_vals)
    unsatisfyTerms = ts.Locations[startIdx].unsatisfyBit(clbits_idx, clbits_vals)
    assert(len(satisfyTerms) + len(unsatisfyTerms) == ts.Locations[startIdx].termNum()), "Condition terms do not match the location's terms."
    if identifier != "" and not identifier.endswith('.'):
        identifier += "."
    if len(satisfyTerms) != 0:
        # Create a new location for the while loop
        whileLocation = pyqreach.Location(qnum)
        for term in satisfyTerms:
            whileLocation.appendClassicalAP(term)
        whileLocation.setIdentifier(identifier + "W")
        ts.addLocation(whileLocation)
        ts.addRelation(startIdx, ts.getLocationNum()-1, pyqreach.QOperation("I", qnum, [0], []))
        # Parse the while loop body
        afterLoopBodyList = parse_qiskit_cir(qc, qnum, ts, [ts.getLocationNum()-1], identifier+"W")
        # For each afterLoopBody location, recursively call the while loop
        for loc in afterLoopBodyList:
            # Check whether the AP of the loc satisfies one of the whileStarter locations
            findEqLoc = False
            for sloc in whileStarter:
                if ts.Locations[loc].equalAP(ts.Locations[sloc]):
                    findEqLoc = True
                    # A back edge to the while loop starter
                    ts.addRelation(loc, sloc, pyqreach.QOperation("I", qnum, [0], []))
                    break
            # print(len(whileStarter), "whileStarter locations", "findEqLoc:", findEqLoc)
            if not findEqLoc:
                newStarter = pyqreach.Location(qnum)
                newStarter.copyClassicalAP(ts.Locations[loc])
                # Seems don't need this
                newStarter.setIdentifier(identifier + "W")
                ts.addLocation(newStarter)  # Create a new location for the after loop body
                ts.addRelation(loc, ts.getLocationNum()-1, pyqreach.QOperation("I", qnum, [0], []))
                whileStarter.append(ts.getLocationNum()-1)  # Add the new location to the whileStarter list
                # Add the while loop location to the exit list DEBUG 0904!!!
                exitList.extend(build_while_loop(qc, qnum, clbits_idx, clbits_vals, whileStarter, ts.getLocationNum()-1, ts, identifier, abstractLevel))
    if len(unsatisfyTerms) != 0:
        # Create a new location for the exit of the while loop
        exitLocation = pyqreach.Location(qnum)
        for term in unsatisfyTerms:
            exitLocation.appendClassicalAP(term)
        exitLocation.setIdentifier(identifier + "EW")
        ts.addLocation(exitLocation)
        ts.addRelation(startIdx, ts.getLocationNum()-1, pyqreach.QOperation("I", qnum, [0], []))
        exitList.append(ts.getLocationNum()-1)
    # update exitList by grouping by classical APs (Default abstractLevel is 1)
    grouped = groupby_classical_aps(ts, exitList)
    for key, locs in grouped.items():
        if len(locs) > 1:
            exitList = merge_locations(ts, whileStarter, locs, identifier+"W.M")
        else:
            exitList.append(locs[0])
    # print("Exit locations after while loop:", exitList)
    return exitList

def simplify_gates(instruction: list, qnum: int) -> list:
    """
    This function merges single reset gates into a single resetAll gate if all qubits are reset.
    instruction: List of Qiskit instructions (gates) to simplify.
    qnum: Total number of qubits in the circuit.
    Returns: A list of simplified Qiskit instructions.
    """
    pass

def parse_qiskit_cir(qc: QuantumCircuit, qnum: int, ts: pyqreach.TransitionSystem, startNodes: list=[], identifier: str="", pivot: int=0, pivotend: int=1000000, abstractLevel: int=1) -> list:
    """
    Parse a Qiskit QuantumCircuit into a pyqreach TransitionSystem, starting from specified nodes.
    qc: QuantumCircuit to parse
    startNodes: List of starting locations (indexes) to begin parsing from
    ts: Transition system to modify
    pivot: Index to start parsing from (default is 0)
    pivotend: Index to end parsing (default is 1000000, meaning until the end of the circuit)
    abstractLevel: Level of abstraction for merging locations (default is 1)
    Returns: List of resulting locations after parsing the circuit
    """
    # Assert start nodes don't exceed numLocations of ts
    for node in startNodes:
        assert node < ts.getLocationNum(), f"Start node {node} exceeds the number of locations in the transition system."
    if startNodes == []:
        # Append an initial location to the transition system
        loc0 = pyqreach.Location(qnum, 0)
        # Initialize the Clasical APs with all zero by the number of clbits of qc
        loc0.appendClassicalAP('0' * qc.num_clbits)
        loc0.setIdentifier("S0")
        ts.addLocation(loc0)
        ts.setInitLocation(0)
        startNodes = [0]  # Start from the initial location
        init_parse(qc)  # Initialize the classical register mapping
    instructions = qc.data[pivot:pivotend] if pivotend != 1000000 else qc.data[pivot:]
    currLoc = startNodes
    resultLocs = []
    pruning_resets = False
    if identifier != "" and not identifier.endswith('.'):
        identifier += "."
    for _,gate in enumerate(instructions):
        # Assume each Locs in currLoc has different classical APs (In the current abstractlevel==1)
        op_name = gate.operation.name
        if op_name != 'reset':
            pruning_resets = False
        if pruning_resets:
            # If we are pruning resets, skip the reset gates
            if op_name == 'reset':
                # pass
                continue
        # Note: Assume there is a single quantum register in the circuit!!!
        qubits = [q._index for q in gate.qubits]
        cbits = [c._index for c in gate.clbits] if gate.clbits else []
        if op_name == 'if_else':
            # Two instruction blocks, one for if and one for else, double the curr_loc_list
            # 1. For each current location, create a branch for ITE (in case part of the clVars satisfy if and part satisfy else). Otherwise, create a single postLoc.
            # 2. Call parse_qiskit_cir recursively for each branch.
            # 3. For each branch, do a heuristic merge.
            if_block_cir = gate.operation.params[0]
            else_block_cir = gate.operation.params[1] if len(gate.operation.params) > 1 else [] # Maybe None
            condition = gate.operation.condition
            clbits_idx, clbits_vals = get_condition_info(qc.cregs, condition)
            # print(clbits_idx, clbits_vals, "Condition info for if_else operation.")
            tempNewCurrLoc = []
            for cLoc in currLoc:
                # First judge if cLoc satisfies the condition
                satisfyTerms = ts.Locations[cLoc].satisfyBit(clbits_idx, clbits_vals)
                # if cLoc == 2606:
                #     print([reg.name for reg in qc.cregs], "Classical registers in the circuit")
                #     print("Debugging location 2606")
                #     print(clbits_idx, clbits_vals, "Condition info for if_else operation.")
                #     print(ts.Locations[cLoc].getIdentifier(), "Identifier of location")
                #     print(satisfyTerms, "Satisfy terms for location 2606")
                #     print(ts.Locations[cLoc].satisfyBit([0,1], [1,0]), "Satisfy terms for condition 10")
                unsatisfyTerms = ts.Locations[cLoc].unsatisfyBit(clbits_idx, clbits_vals)
                assert(len(satisfyTerms) + len(unsatisfyTerms) == ts.Locations[cLoc].termNum()), "Condition terms do not match the location's terms."
                if_result_locs, else_result_locs = [], []
                if len(satisfyTerms) != 0:
                    # Create a new location for the if block
                    ifLocation = pyqreach.Location(qnum)
                    for term in satisfyTerms:
                        ifLocation.appendClassicalAP(term)
                    ifLocation.setIdentifier(identifier + "S" + str(pivot + _ + 1) + ".I")
                    ts.addLocation(ifLocation)
                    ts.addRelation(cLoc, ts.getLocationNum()-1, pyqreach.QOperation("I", qnum, [0], []))
                    # Parse the if block
                    if_result_locs = parse_qiskit_cir(if_block_cir, qnum, ts, [ts.getLocationNum()-1], identifier+"S"+str(pivot + _ + 1)+".I")
                if len(unsatisfyTerms) != 0:
                    # Create a new location for the else block
                    elseLocation = pyqreach.Location(qnum)
                    for term in unsatisfyTerms:
                        elseLocation.appendClassicalAP(term)
                    elseLocation.setIdentifier(identifier + "S" + str(pivot + _ + 1) + ".E")
                    ts.addLocation(elseLocation)
                    ts.addRelation(cLoc, ts.getLocationNum()-1, pyqreach.QOperation("I", qnum, [0], []))
                    # if cLoc == 2606:
                    #     print("Debugging location 2606 - else branch")
                    #     print(unsatisfyTerms, "Unsatisfy terms for location 2606")
                    #     print(ts.Locations[cLoc].getIdentifier(), "Identifier of location")
                    #     print("New else location idx:", ts.getLocationNum()-1)
                    #     print(ts.Locations[ts.getLocationNum()-1].getIdentifier(), "Identifier of else location")
                    # Parse the else block
                    if else_block_cir is not None:
                        else_result_locs = parse_qiskit_cir(else_block_cir, qnum, ts, [ts.getLocationNum()-1], identifier+"S"+str(pivot + _ + 1)+".E")
                    else:
                        else_result_locs = [ts.getLocationNum()-1]
                tempNewCurrLoc.extend(if_result_locs)
                tempNewCurrLoc.extend(else_result_locs)
            # Merge locations if needed
            currLoc = tempNewCurrLoc
            if abstractLevel == 1:
                # merge locations in currLoc when they share the same classical APs. (Don't merge measured locations)
                grouped = groupby_classical_aps(ts, tempNewCurrLoc)
                if 2636 in currLoc:
                    print("Grouped locations by classical APs:", grouped)
                tempNewCurrLoc = []
                for key, locs in grouped.items():
                    if len(locs) > 1:
                        currLoc = merge_locations(ts, currLoc, locs, identifier+"I.M")
                    else:
                        currLoc.append(locs[0])
            else:
                raise ValueError("Unsupported merge level for if_else operation.")
        elif op_name == 'while_loop':
            # print("While loop operation detected")
            while_block_cir = gate.operation.params[0]
            condition = gate.operation.condition
            clbits_idx, clbits_vals = get_condition_info(qc.cregs, condition)
            outLoopLocs = []
            for cLoc in currLoc:
                exitLocs = build_while_loop(while_block_cir, qnum, clbits_idx, clbits_vals, [cLoc], cLoc, ts, identifier+"S"+str(pivot+_+1), abstractLevel)
                outLoopLocs.extend(exitLocs)
            currLoc = outLoopLocs
            # Merge locations if needed
            if abstractLevel == 1:
                # merge locations in currLoc when they share the same classical APs. (Don't merge measured locations)
                grouped = groupby_classical_aps(ts, outLoopLocs)
                # print("Grouped locations by classical APs:", grouped)
                for key, locs in grouped.items():
                    if len(locs) > 1:
                        currLoc = merge_locations(ts, currLoc, locs, identifier+"W.M")
                    else:
                        currLoc.append(locs[0])
            else:
                raise ValueError("Unsupported merge level for while_loop operation.")
                
        elif op_name == 'measure':
            # Another operation that split the locations. The only operation that can modify classical bits.
            for i,cl in enumerate(currLoc):
                for j in range(i+1, len(currLoc)):
                    assert not ts.Locations[cl].equalAP(ts.Locations[currLoc[j]]), "Measure operation cannot be applied to locations with equal classical APs."
            assert len(qubits) == 1, "Measure operation can only be applied to one qubit at a time."
            measuredLocDict = {}
            tempNewCurrLoc = []
            for cLoc in currLoc:
                # create two locations with meas0 and meas1, update their classical APs, let their string representation as keys, if in measuredLocDict,
                # then use the existing location, otherwise, append the new location to the transition system and update the measuredLocDict.
                loc_meas0 = pyqreach.Location(qnum, 0)
                loc_meas1 = pyqreach.Location(qnum, 0)
                loc_meas0.copyClassicalAP(ts.Locations[cLoc])  # Copy classical APs from the current location
                loc_meas1.copyClassicalAP(ts.Locations[cLoc])
                loc_meas0.setClassicalValue(qubits[0], 0)
                loc_meas1.setClassicalValue(qubits[0], 1)
                meas0_key = loc_meas0.cp.toString()
                meas1_key = loc_meas1.cp.toString()
                if meas0_key not in measuredLocDict:
                    ts.addLocation(loc_meas0)
                    # if _ + 1 < len(instructions):
                    #     if instructions[_+1].operation.name == 'measure':
                    #         if loc_meas0.satisfyBit([0,1],[1,0]):
                    #             print("Debugging after two consecutive measures - meas0 branch")
                    #             print("Identifier: ", ts.Locations[cLoc].getIdentifier(), "Location idx:", ts.getLocationNum()-1)
                    #             print("cLoc: ", cLoc)
                    # if _ - 1 >= 0:
                    #     if instructions[_-1].operation.name == 'measure':
                    #         if loc_meas0.satisfyBit([0,1],[1,1]):
                    #             print("Debugging after two consecutive measures - meas0 branch")
                    #             print("Identifier: ", ts.Locations[cLoc].getIdentifier(), "Location idx:", ts.getLocationNum()-1)
                    #             print("cLoc: ", cLoc)
                    measuredLocDict[meas0_key] = ts.getLocationNum() - 1
                    tempNewCurrLoc.append(ts.getLocationNum() - 1)
                    ts.addRelation(cLoc, ts.getLocationNum() - 1, pyqreach.QOperation("meas0", qnum, qubits, []))
                else:
                    ts.addRelation(cLoc, measuredLocDict[meas0_key], pyqreach.QOperation("meas0", qnum, qubits, []))
                if meas1_key not in measuredLocDict:
                    ts.addLocation(loc_meas1)
                    # if _ + 1 < len(instructions):
                    #     if instructions[_+1].operation.name == 'measure':
                    #         if loc_meas1.satisfyBit([0,1],[1,0]):
                    #             print("Debugging after the first consecutive measures - meas1 branch")
                    #             print("Identifier: ", ts.Locations[cLoc].getIdentifier(), "Location idx:", ts.getLocationNum()-1)
                    #             print("cLoc: ", cLoc)
                    # if _ - 1 >= 0:
                    #     if instructions[_-1].operation.name == 'measure' and ts.Locations[cLoc].satisfyBit([0,1],[1,0]):
                    #         if loc_meas1.satisfyBit([0,1],[1,1]):
                    #             print("Debugging after the second consecutive measures - meas1 branch")
                    #             print("Identifier: ", ts.Locations[cLoc].getIdentifier(), "Location idx:", ts.getLocationNum()-1)
                    #             print("cLoc: ", cLoc)
                    measuredLocDict[meas1_key] = ts.getLocationNum() - 1
                    tempNewCurrLoc.append(ts.getLocationNum() - 1)
                    ts.addRelation(cLoc, ts.getLocationNum() - 1, pyqreach.QOperation("meas1", qnum, qubits, []))
                else:
                    ts.addRelation(cLoc, measuredLocDict[meas1_key], pyqreach.QOperation("meas1", qnum, qubits, []))
            # Update the current locations
            currLoc = tempNewCurrLoc
            for l in currLoc:
                ts.Locations[l].setIdentifier(identifier + "S" + str(pivot + _ + 1))
        else:
            op = None
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
            elif op_name == 'u':
                theta, phi, lam = gate.operation.params
                op = pyqreach.QOperation("U3", qnum, qubits, [theta, phi, lam])
            elif op_name == 'cx':
                op = pyqreach.QOperation("CX", qnum, qubits, [])
            elif op_name == 'cz':
                op = pyqreach.QOperation("CZ", qnum, qubits, [])
            elif op_name == 'cp':
                op = pyqreach.QOperation("CP", qnum, qubits, [gate.operation.params[0]])
            elif op_name == 'swap':
                op = pyqreach.QOperation("SWAP", qnum, qubits, [])
            elif op_name == 'reset':
                # Reset operation, we assume it resets all qubits to |0>, using resetAll QOperation, or reset a single qubit to |0>,
                # resulting in a mixed state, we use reset QOperation.
                # Check if the continuous qnum gates are all resets, if so, we set pruning_resets to True.
                doResetAll = True if _ + qnum <= len(instructions) and all(instructions[i][0].name == 'reset' for i in range(_, _ + qnum-1)) else False
                # If there exists a qubit that is not reset, we set doResetAll to False.
                recordResetSet = set()
                if doResetAll:
                    for gidx in range(_, _+qnum-1):
                        resetBit = instructions[gidx].qubits[0]._index
                        if resetBit not in recordResetSet:
                            recordResetSet.add(resetBit)
                        else:
                            doResetAll = False
                if not doResetAll:
                    op = pyqreach.QOperation("reset", qnum, qubits, [])
                else:
                    op = pyqreach.QOperation("resetAll", qnum, qubits, [])
                    pruning_resets = True  # Set pruning_resets to True to skip the reset gates in the next iterations
            else:
                raise ValueError(f"Unsupported gate: {op_name}")
            tempNewCurrLoc = []
            for cLoc in currLoc:
                # Make new locations
                newloc = pyqreach.Location(qnum, 0)
                newloc.copyClassicalAP(ts.Locations[cLoc])  # Copy classical APs from the current location!!
                newloc.setIdentifier(identifier + "S" + str(pivot + _ + 1))
                ts.addLocation(newloc)
                # Add the operation to the transition system
                ts.addRelation(cLoc, ts.getLocationNum()-1, op)
                tempNewCurrLoc.append(ts.getLocationNum()-1)
            # Update the current locations
            currLoc = tempNewCurrLoc
    resultLocs = currLoc
    return resultLocs

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

def applyFinalMeasurement(ts: pyqreach.TransitionSystem, PauliString: str, qlist: list, qnum: int, entryNode=-1) -> list:
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


def applyMeasureAndReset(ts: pyqreach.TransitionSystem, PauliString: str, qlist: list, qnum: int, entryNode=-1) -> list:
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
    if len(keepBitstrings) > 0:
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

