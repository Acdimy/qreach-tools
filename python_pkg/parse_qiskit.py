import pyqreach
from qiskit import QuantumCircuit, QuantumRegister, ClassicalRegister
from qiskit.circuit import Clbit
from graphviz import Digraph
from math import floor
from math import ceil, log2, pi
import numpy as np


def _is_symbolic_ts(ts) -> bool:
    return hasattr(ts, "getLocationIDs") and hasattr(ts, "getClassicalProposition")


def _make_cp(terms: list[str] | None = None):
    cp = pyqreach.ClassicalProposition()
    if terms is not None:
        for term in terms:
            cp.addTerm(term)
    return cp


def _add_location(ts, qnum: int, identifier: str = "", terms: list[str] | None = None, copy_from: int | None = None) -> int:
    if _is_symbolic_ts(ts):
        if copy_from is not None:
            cp = ts.getClassicalProposition(copy_from)
        else:
            cp = _make_cp(terms)
        return ts.addLocation(cp, identifier)

    loc = pyqreach.Location(qnum)
    if copy_from is not None:
        loc.copyClassicalAP(ts.Locations[copy_from])
    else:
        for term in terms or []:
            loc.appendClassicalAP(term)
    loc.setIdentifier(identifier)
    ts.addLocation(loc)
    return ts.getLocationNum() - 1


def _append_classical_ap(ts, loc: int, term: str):
    if _is_symbolic_ts(ts):
        ts.appendClassicalAP(loc, term)
    else:
        ts.Locations[loc].appendClassicalAP(term)


def _set_classical_value(ts, loc: int, index: int, value: bool):
    if _is_symbolic_ts(ts):
        ts.setClassicalValue(loc, index, value)
    else:
        ts.Locations[loc].setClassicalValue(index, value)


def _set_identifier(ts, loc: int, identifier: str):
    if _is_symbolic_ts(ts):
        ts.setIdentifier(loc, identifier)
    else:
        ts.Locations[loc].setIdentifier(identifier)


def _get_cp(ts, loc: int):
    if _is_symbolic_ts(ts):
        return ts.getClassicalProposition(loc)
    return ts.Locations[loc].cp


def _get_cp_terms(ts, loc: int) -> list:
    return list(_get_cp(ts, loc).terms)


def _get_cp_string(ts, loc: int) -> str:
    return _get_cp(ts, loc).toString()


def _term_num(ts, loc: int) -> int:
    if _is_symbolic_ts(ts):
        return ts.termNum(loc)
    return ts.Locations[loc].termNum()


def _satisfy_bit(ts, loc: int, clbits_idx: list, clbits_vals: list) -> list:
    if _is_symbolic_ts(ts):
        return ts.satisfyBit(loc, clbits_idx, clbits_vals)
    return ts.Locations[loc].satisfyBit(clbits_idx, clbits_vals)


def _unsatisfy_bit(ts, loc: int, clbits_idx: list, clbits_vals: list) -> list:
    if _is_symbolic_ts(ts):
        return ts.unsatisfyBit(loc, clbits_idx, clbits_vals)
    return ts.Locations[loc].unsatisfyBit(clbits_idx, clbits_vals)


def _equal_ap(ts, loc_a: int, loc_b: int) -> bool:
    return _get_cp_string(ts, loc_a) == _get_cp_string(ts, loc_b)


def _get_identifier(ts, loc: int) -> str:
    if _is_symbolic_ts(ts):
        return ts.getIdentifier(loc)
    return ts.Locations[loc].getIdentifier()


def _get_post_locations(ts, loc: int) -> list:
    if _is_symbolic_ts(ts):
        return ts.getPostLocations(loc)
    return list(ts.Locations[loc].postLocations)


def _get_location_ids(ts) -> list:
    if _is_symbolic_ts(ts):
        return ts.getLocationIDs()
    return [loc.idx for loc in ts.Locations]


def _get_qnum(ts, loc: int | None = None) -> int:
    if _is_symbolic_ts(ts):
        if loc is not None:
            return ts.getLocationAnnotation(loc).qNum
        if ts.getLocationNum() > 0:
            return ts.getLocationAnnotation(0).qNum
        return 0
    if loc is not None:
        return ts.Locations[loc].qNum
    if ts.getLocationNum() > 0:
        return ts.Locations[0].qNum
    return 0


def _add_location_with_cp(ts, qnum: int, cp, identifier: str = "") -> int:
    if _is_symbolic_ts(ts):
        return ts.addLocation(cp, identifier)

    loc = pyqreach.Location(qnum)
    for term in cp.terms:
        loc.appendClassicalAP(term)
    loc.setIdentifier(identifier)
    ts.addLocation(loc)
    return ts.getLocationNum() - 1


def _copy_cp_with_value(ts, loc: int, index: int, value: bool):
    cp = pyqreach.ClassicalProposition()
    for term in _get_cp_terms(ts, loc):
        cp.addTerm(term)
    cp.setValue(index, value)
    return cp

# Get the binary representation of a number
def getBinary(num, length):
    return format(num, '0' + str(length) + 'b')

import numpy as np

def expand_amplitude(p: np.ndarray, idx: list[int], n: int) -> np.ndarray:
    """
    Expand a partial amplitude vector `p` (on qubits `idx`) into a full
    n-qubit amplitude vector, assuming all other qubits are |0>.

    Args:
        p   : np.ndarray, shape (2**k,), complex amplitudes
        idx : list of int, length k, distinct qubit indices in [0, n-1]
        n   : total number of qubits

    Returns:
        np.ndarray, shape (2**n,), complex amplitudes
    """
    k = len(idx)
    assert p.size == 2**k, "Length of p must be 2**len(idx)"
    
    # reshape p into a tensor with k qubits
    tensor = p.reshape([2]*k)

    # prepare full tensor with all qubits
    full_tensor = np.zeros([2]*n, dtype=complex)

    # insert p into the positions specified by idx, others fixed to |0>
    # we do this by using advanced indexing
    index = [0]*n
    for basis in np.ndindex(*([2]*k)):
        for i, qubit in enumerate(idx):
            index[qubit] = basis[i]
        full_tensor[tuple(index)] = tensor[basis]
        # reset unused positions back to 0 (since only one slice is filled)
        for qubit in set(range(n)) - set(idx):
            index[qubit] = 0

    # flatten to vector
    return full_tensor.reshape(-1)


# Global variables for mapping classical registers to indices
CLREG_ORDER = []
# Global variables for mapping quantum registers to indices
QREG_ORDER = []

def init_parse(qc: QuantumCircuit):
    """
    Initialize the classical register mapping for parsing.
    
    Args:
        qc (QuantumCircuit): The quantum circuit to parse.
    """
    global CLREG_ORDER
    global QREG_ORDER
    CLREG_ORDER = []
    for creg in qc.cregs:
        CLREG_ORDER.append(creg)
    # print(CLREG_ORDER, "Classical registers in the circuit")
    QREG_ORDER = []
    for qreg in qc.qregs:
        QREG_ORDER.append(qreg)

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

def get_global_cl_index(clbits) -> list:
    """
    Get the global indices of the classical bit of a gate.
    
    clbits: list of Clbit (from a Qiskit instruction)
    Returns: list of global indices of the classical bits
    """
    global CLREG_ORDER
    clbits_idx = []
    for c in clbits:
        clreg_name = c._register.name
        clreg_bit_idx = c._index
        clreg_reg_idx = next((i for i, reg in enumerate(CLREG_ORDER) if reg.name == clreg_name), None)
        reg_start_idx = sum([reg.size for reg in CLREG_ORDER[:clreg_reg_idx]]) + clreg_bit_idx
        clbits_idx.append(reg_start_idx)
    return clbits_idx

def get_global_qb_index(qubits) -> list:
    """
    Get the global indices of the quantum bit of a gate.
    
    qubits: list of Qubit (from a Qiskit instruction)
    Returns: list of global indices of the quantum bits
    """
    global QREG_ORDER
    qbits_idx = []
    for q in qubits:
        qreg_name = q._register.name
        qreg_bit_idx = q._index
        qreg_reg_idx = next((i for i, reg in enumerate(QREG_ORDER) if reg.name == qreg_name), None)
        reg_start_idx = sum([reg.size for reg in QREG_ORDER[:qreg_reg_idx]]) + qreg_bit_idx
        qbits_idx.append(reg_start_idx)
    return qbits_idx

def groupby_classical_aps(ts: pyqreach.TransitionSystem, currLocs: list) -> dict:
    """
    Group current locations by their classical APs.

    ts: Transition system (including all locations)
    currLocs: List of indexes of current locations. Each of them has a classical APs.
    Returns: A dictionary where keys are tuples of classical APs and values are lists of location indexes.
    """
    grouped = {}
    for loc_idx in currLocs:
        key = _get_cp_string(ts, loc_idx)
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
    qnum = _get_qnum(ts, currLocs[0])
    newLocIdx = _add_location(ts, qnum, identifier)
    newCurrLocs = [loc for loc in currLocs if loc not in toMergeLocs]
    newCurrLocs.append(newLocIdx)
    # if newLocIdx == 1332:
    #     print(toMergeLocs, "Merging locations into new location 1332")
    for l in toMergeLocs:
        # Merge classical APs
        for term in _get_cp_terms(ts, l):
            if term not in _get_cp_terms(ts, newLocIdx):
                _append_classical_ap(ts, newLocIdx, term)
        # Add relations from the merged locations to the new location
        ts.addRelation(l, newLocIdx, pyqreach.QOperation("I", qnum, [0], []))
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
    satisfyTerms = _satisfy_bit(ts, startIdx, clbits_idx, clbits_vals)
    unsatisfyTerms = _unsatisfy_bit(ts, startIdx, clbits_idx, clbits_vals)
    # print(ts.Locations[startIdx].cp.toString())
    # print(len(satisfyTerms), "satisfy terms, ", len(unsatisfyTerms), "unsatisfy terms at location", startIdx)
    assert(len(satisfyTerms) + len(unsatisfyTerms) == _term_num(ts, startIdx)), "Condition terms do not match the location's terms."
    if identifier != "" and not identifier.endswith('.'):
        identifier += "."
    if len(satisfyTerms) != 0:
        # Create a new location for the while loop
        while_loc = _add_location(ts, qnum, identifier + "W", terms=satisfyTerms)
        ts.addRelation(startIdx, while_loc, pyqreach.QOperation("I", qnum, [0], []))
        # Parse the while loop body
        afterLoopBodyList = parse_qiskit_cir(qc, qnum, ts, [while_loc], identifier+"W")
        # For each afterLoopBody location, recursively call the while loop
        for loc in afterLoopBodyList:
            # Check whether the AP of the loc satisfies one of the whileStarter locations
            findEqLoc = False
            for sloc in whileStarter:
                if _equal_ap(ts, loc, sloc):
                    findEqLoc = True
                    # A back edge to the while loop starter
                    ts.addRelation(loc, sloc, pyqreach.QOperation("I", qnum, [0], []))
                    break
            # print(len(whileStarter), "whileStarter locations", "findEqLoc:", findEqLoc)
            if not findEqLoc:
                new_starter = _add_location(ts, qnum, identifier + "WN", copy_from=loc)
                ts.addRelation(loc, new_starter, pyqreach.QOperation("I", qnum, [0], []))
                whileStarter.append(new_starter)  # Add the new location to the whileStarter list
                # Add the while loop location to the exit list DEBUG 0904!!!
                exitList.extend(build_while_loop(qc, qnum, clbits_idx, clbits_vals, whileStarter, new_starter, ts, identifier, abstractLevel))
    if len(unsatisfyTerms) != 0:
        # Create a new location for the exit of the while loop
        exit_loc = _add_location(ts, qnum, identifier + "EW", terms=unsatisfyTerms)
        ts.addRelation(startIdx, exit_loc, pyqreach.QOperation("I", qnum, [0], []))
        exitList.append(exit_loc)
    # update exitList by grouping by classical APs (Default abstractLevel is 1)
    grouped = groupby_classical_aps(ts, exitList)
    merged_exit_list = list(exitList)
    for key, locs in grouped.items():
        if len(locs) > 1:
            merged_exit_list = merge_locations(ts, merged_exit_list, locs, identifier+"W.M")
    exitList = merged_exit_list
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
        loc0 = _add_location(ts, qnum, "S0", terms=['0' * qc.num_clbits])
        ts.setInitLocation(loc0)
        startNodes = [loc0]  # Start from the initial location
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
        # qubits = [q._index for q in gate.qubits]
        qubits = get_global_qb_index(gate.qubits) if gate.qubits else []
        cbits = get_global_cl_index(gate.clbits) if gate.clbits else []
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
                satisfyTerms = _satisfy_bit(ts, cLoc, clbits_idx, clbits_vals)
                unsatisfyTerms = _unsatisfy_bit(ts, cLoc, clbits_idx, clbits_vals)
                assert(len(satisfyTerms) + len(unsatisfyTerms) == _term_num(ts, cLoc)), "Condition terms do not match the location's terms."
                if_result_locs, else_result_locs = [], []
                if len(satisfyTerms) != 0:
                    # Create a new location for the if block
                    if_loc = _add_location(ts, qnum, identifier + "S" + str(pivot + _ + 1) + ".I", terms=satisfyTerms)
                    ts.addRelation(cLoc, if_loc, pyqreach.QOperation("I", qnum, [0], []))
                    # Parse the if block
                    if_result_locs = parse_qiskit_cir(if_block_cir, qnum, ts, [if_loc], identifier+"S"+str(pivot + _ + 1)+".I")
                if len(unsatisfyTerms) != 0:
                    # Create a new location for the else block
                    else_loc = _add_location(ts, qnum, identifier + "S" + str(pivot + _ + 1) + ".E", terms=unsatisfyTerms)
                    ts.addRelation(cLoc, else_loc, pyqreach.QOperation("I", qnum, [0], []))
                    # Parse the else block
                    if else_block_cir is not None:
                        else_result_locs = parse_qiskit_cir(else_block_cir, qnum, ts, [else_loc], identifier+"S"+str(pivot + _ + 1)+".E")
                    else:
                        else_result_locs = [else_loc]
                tempNewCurrLoc.extend(if_result_locs)
                tempNewCurrLoc.extend(else_result_locs)
            # Merge locations if needed
            currLoc = tempNewCurrLoc
            if abstractLevel == 1:
                # merge locations in currLoc when they share the same classical APs. (Don't merge measured locations)
                grouped = groupby_classical_aps(ts, tempNewCurrLoc)
                merged_curr_locs = list(currLoc)
                for key, locs in grouped.items():
                    if len(locs) > 1:
                        merged_curr_locs = merge_locations(ts, merged_curr_locs, locs, identifier+"I.M")
                currLoc = merged_curr_locs
            else:
                raise ValueError("Unsupported merge level for if_else operation.")
        elif op_name == 'while_loop':
            # print("While loop operation detected")
            while_block_cir = gate.operation.params[0]
            condition = gate.operation.condition
            clbits_idx, clbits_vals = get_condition_info(qc.cregs, condition)
            # print(clbits_idx, clbits_vals, "Condition info for while_loop operation.")
            outLoopLocs = []
            for cLoc in currLoc:
                exitLocs = build_while_loop(while_block_cir, qnum, clbits_idx, clbits_vals, [cLoc], cLoc, ts, identifier+"S"+str(pivot+_+1), abstractLevel)
                outLoopLocs.extend(exitLocs)
            currLoc = outLoopLocs
            # Merge locations if needed
            if abstractLevel == 1:
                # merge locations in currLoc when they share the same classical APs. (Don't merge measured locations)
                grouped = groupby_classical_aps(ts, outLoopLocs)
                merged_curr_locs = list(currLoc)
                for key, locs in grouped.items():
                    if len(locs) > 1:
                        merged_curr_locs = merge_locations(ts, merged_curr_locs, locs, identifier+"W.M")
                currLoc = merged_curr_locs
            else:
                raise ValueError("Unsupported merge level for while_loop operation.")
        elif op_name == 'for_loop':
            # print("For loop operation detected")
            for_block_cir = gate.operation.params[0]
            loop_range = gate.operation.params[1]
            outLoopLocs = []
            for cLoc in currLoc:
                loopStarter = [cLoc]
                for _ in range(loop_range):
                    loopResultLocs = parse_qiskit_cir(for_block_cir, qnum, ts, loopStarter, identifier+"S"+str(pivot+_+1)+".F")
                    loopStarter = loopResultLocs
                outLoopLocs.extend(loopStarter)
            currLoc = outLoopLocs
            # Merge locations if needed
            if abstractLevel == 1:
                # merge locations in currLoc when they share the same classical APs. (Don't merge measured locations)
                grouped = groupby_classical_aps(ts, outLoopLocs)
                merged_curr_locs = list(currLoc)
                for key, locs in grouped.items():
                    if len(locs) > 1:
                        merged_curr_locs = merge_locations(ts, merged_curr_locs, locs, identifier+"F.M")
                currLoc = merged_curr_locs
            else:
                raise ValueError("Unsupported merge level for for_loop operation.")
        elif op_name == 'switch_case':
            # Multiple instruction blocks, one for each case, double the curr_loc_list
            # 1. For each current location, create a branch for each case (in case part of the clVars satisfy one case and part satisfy another). Otherwise, create a single postLoc.
            # 2. Call parse_qiskit_cir recursively for each branch.
            # 3. For each branch, do a heuristic merge.
            case_block_cirs = gate.operation.params[0]  # A list of QuantumCircuits for each case
            condition = gate.operation.condition
            clbits_idx, clbits_vals = get_condition_info(qc.cregs, condition)
            # print(clbits_idx, clbits_vals, "Condition info for switch_case operation.")
            tempNewCurrLoc = []
            for cLoc in currLoc:
                # First judge if cLoc satisfies the condition
                satisfyTerms = _satisfy_bit(ts, cLoc, clbits_idx, clbits_vals)
                unsatisfyTerms = _unsatisfy_bit(ts, cLoc, clbits_idx, clbits_vals)
                assert(len(satisfyTerms) + len(unsatisfyTerms) == _term_num(ts, cLoc)), "Condition terms do not match the location's terms."
                case_result_locs = []
                if len(satisfyTerms) != 0:
                    # Create a new location for the switch_case block
                    switch_loc = _add_location(ts, qnum, identifier + "S" + str(pivot + _ + 1) + ".C", terms=satisfyTerms)
                    ts.addRelation(cLoc, switch_loc, pyqreach.QOperation("I", qnum, [0], []))
                    # Parse each case block
                    for idx,case_cir in enumerate(case_block_cirs):
                        case_result_locs = parse_qiskit_cir(case_cir, qnum, ts, [switch_loc], identifier+"S"+str(pivot + _ + 1)+".C"+str(idx))
                        tempNewCurrLoc.extend(case_result_locs)
                if len(unsatisfyTerms) != 0:
                    # Create a new location for the default block (not satisfying any case)
                    default_loc = _add_location(ts, qnum, identifier + "S" + str(pivot + _ + 1) + ".D", terms=unsatisfyTerms)
                    ts.addRelation(cLoc, default_loc, pyqreach.QOperation("I", qnum, [0], []))
                    # Parse the default block
                    default_result_locs = parse_qiskit_cir(gate.operation.default, qnum, ts, [default_loc], identifier+"S"+str(pivot + _ + 1)+".D") if gate.operation.default is not None else [default_loc]
                    tempNewCurrLoc.extend(default_result_locs)
            # Merge locations if needed
            currLoc = tempNewCurrLoc
            if abstractLevel == 1:
                # merge locations in currLoc when they share the same classical APs. (Don't merge measured locations)
                grouped = groupby_classical_aps(ts, tempNewCurrLoc)
                merged_curr_locs = list(currLoc)
                for key, locs in grouped.items():
                    if len(locs) > 1:
                        merged_curr_locs = merge_locations(ts, merged_curr_locs, locs, identifier+"C.M")
                currLoc = merged_curr_locs
            else:
                raise ValueError("Unsupported merge level for switch_case operation.")
                
        elif op_name == 'measure':
            # Another operation that split the locations. The only operation that can modify classical bits.
            # Assert single qubit measurement!
            # for i_meas,cl in enumerate(currLoc):
            #     for j in range(i_meas+1, len(currLoc)):
            #         assert not ts.Locations[cl].equalAP(ts.Locations[currLoc[j]]), "Measure operation cannot be applied to locations with equal classical APs."
            assert len(qubits) == 1, "Measure operation can only be applied to one qubit at a time."
            measuredLocDict = {}
            tempNewCurrLoc = []
            for cLoc in currLoc:
                cp_meas0 = _copy_cp_with_value(ts, cLoc, cbits[0], False)
                cp_meas1 = _copy_cp_with_value(ts, cLoc, cbits[0], True)
                meas0_key = cp_meas0.toString()
                meas1_key = cp_meas1.toString()
                if meas0_key not in measuredLocDict:
                    meas0_loc = _add_location_with_cp(ts, qnum, cp_meas0)
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
                    measuredLocDict[meas0_key] = meas0_loc
                    tempNewCurrLoc.append(meas0_loc)
                    ts.addRelation(cLoc, meas0_loc, pyqreach.QOperation("meas0", qnum, qubits, []))
                else:
                    ts.addRelation(cLoc, measuredLocDict[meas0_key], pyqreach.QOperation("meas0", qnum, qubits, []))
                if meas1_key not in measuredLocDict:
                    meas1_loc = _add_location_with_cp(ts, qnum, cp_meas1)
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
                    measuredLocDict[meas1_key] = meas1_loc
                    tempNewCurrLoc.append(meas1_loc)
                    # print(ts.Locations[ts.getLocationNum()-1].cp.toString(), "Classical APs of the new measured location at ", cbits[0])
                    ts.addRelation(cLoc, meas1_loc, pyqreach.QOperation("meas1", qnum, qubits, []))
                else:
                    ts.addRelation(cLoc, measuredLocDict[meas1_key], pyqreach.QOperation("meas1", qnum, qubits, []))
            # Update the current locations
            currLoc = tempNewCurrLoc
            for l in currLoc:
                _set_identifier(ts, l, identifier + "S" + str(pivot + _ + 1))
        elif op_name == 'initialize':
            # 1. Apply reset to the qubit indexes;
            # 2. Apply an init gate to the qubit indexes;
            # Assert the indexes are sequential ordered
            assert all(qubits[i] + 1 == qubits[i + 1] for i in range(len(qubits) - 1)), "For now, initialize operation can only be applied to sequential qubits."
            tempNewCurrLoc = []
            for cLoc in currLoc:
                indexNum = len(qubits)
                # For each index, apply reset gates seperately, and append new locations in a sequence
                prevLoc = cLoc
                for i in range(indexNum):
                    reset_loc = _add_location(ts, qnum, identifier + "S" + str(pivot + _ + 1) + ".R" + str(i+1), copy_from=prevLoc)
                    ts.addRelation(prevLoc, reset_loc, pyqreach.QOperation("reset", qnum, [qubits[i]], []))
                    prevLoc = reset_loc
                # Then apply the init gate
                init_loc = _add_location(ts, qnum, identifier + "S" + str(pivot + _ + 1) + ".N", copy_from=prevLoc)
                params_real = [param.real for param in gate.operation.params]
                params_imag = [param.imag for param in gate.operation.params]
                ts.addRelation(prevLoc, init_loc, pyqreach.QOperation("init", qnum, qubits, params_real + params_imag))
                tempNewCurrLoc.append(init_loc)
            # Update the current locations
            currLoc = tempNewCurrLoc
        elif op_name == 'dcx':
            # Slightly difficult to implement in C++
            # Just apply cx[q0,q1], cx[q1,q0]
            assert len(qubits) == 2, "DCX operation must be applied to two qubits."
            op1 = pyqreach.QOperation("CX", qnum, [qubits[0], qubits[1]], [])
            op2 = pyqreach.QOperation("CX", qnum, [qubits[1], qubits[0]], [])
            tempNewCurrLoc = []
            for cLoc in currLoc:
                # Make new locations
                loc1 = _add_location(ts, qnum, identifier + "S" + str(pivot + _ + 1) + ".1", copy_from=cLoc)
                ts.addRelation(cLoc, loc1, op1)
                loc2 = _add_location(ts, qnum, identifier + "S" + str(pivot + _ + 1), copy_from=loc1)
                ts.addRelation(loc1, loc2, op2)
                tempNewCurrLoc.append(loc2)
            # Update the current locations
            currLoc = tempNewCurrLoc
        elif op_name == 'barrier':
            # Barrier operation, do nothing
            # for l in currLoc:
            #     ts.Locations[l].setIdentifier(identifier + "S" + str(pivot + _ + 1))
            continue
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
            elif op_name == 'tdg':
                op = pyqreach.QOperation("U3", qnum, qubits, [0, 0, -1/4])
            elif op_name == 'u':
                theta, phi, lam = gate.operation.params
                op = pyqreach.QOperation("U3", qnum, qubits, [theta/pi, phi/pi, lam/pi])
            elif op_name == 'u1':
                lam = gate.operation.params[0]
                op = pyqreach.QOperation("U3", qnum, qubits, [0, 0, lam/pi])
            elif op_name == 'rx':
                pass
            elif op_name == 'ry':
                theta = gate.operation.params[0]
                op = pyqreach.QOperation("U3", qnum, qubits, [theta/pi, 0, 0])
            elif op_name == 'rz':
                # Lack a global phase
                lam = gate.operation.params[0]
                op = pyqreach.QOperation("U3", qnum, qubits, [0, 0, lam/pi])
            elif op_name == 'p':
                lam = gate.operation.params[0]
                op = pyqreach.QOperation("U3", qnum, qubits, [0, 0, lam/pi])
            elif op_name == 'sx':
                op = pyqreach.QOperation("SX", qnum, qubits, [])
            elif op_name == 'iX':
                op = pyqreach.QOperation("U3", qnum, qubits, [1, 1/2, -1/2])
            elif op_name == 'iY':
                op = pyqreach.QOperation("U3", qnum, qubits, [1, 1, 1])
            elif op_name == 'iZ':
                op = pyqreach.QOperation("arb", qnum, qubits, [0,1,0,0,0,0,0,-1])
            elif op_name == '-iX':
                op = pyqreach.QOperation("U3", qnum, qubits, [1, -1/2, 1/2])
            elif op_name == '-iY':
                op = pyqreach.QOperation("U3", qnum, qubits, [1, 0, 0])
            elif op_name == '-iZ':
                op = pyqreach.QOperation("arb", qnum, qubits, [0,-1,0,0,0,0,0,1])
            elif op_name == 'cx':
                op = pyqreach.QOperation("CX", qnum, qubits, [])
            elif op_name == 'cz':
                op = pyqreach.QOperation("CZ", qnum, qubits, [])
            elif op_name == 'cp':
                op = pyqreach.QOperation("CP", qnum, qubits, [gate.operation.params[0]/pi])
            elif op_name == 'cu1':
                lam = gate.operation.params[0]/pi
                op = pyqreach.QOperation("CP", qnum, qubits, [lam/pi])
            elif op_name == 'csx':
                op = pyqreach.QOperation("CSX", qnum, qubits, [])
            elif op_name == 'swap':
                op = pyqreach.QOperation("SWAP", qnum, qubits, [])
            elif op_name == 'iswap':
                op = pyqreach.QOperation("iSWAP", qnum, qubits, [])
            elif op_name == 'ccx':
                op = pyqreach.QOperation("CCX", qnum, qubits, [])
            elif op_name == 'reset':
                # Reset operation, we assume it resets all qubits to |0>, using resetAll QOperation, or reset a single qubit to |0>,
                # resulting in a mixed state, we use reset QOperation.
                # Check if the continuous qnum gates are all resets, if so, we set pruning_resets to True.
                # TODO: Here is a trick! we assume all reset gates on different qubits are grouped together!
                doResetAll = True if _ + qnum <= len(instructions) and all(instructions[i][0].name == 'reset' for i in range(_, _ + qnum)) else False
                # If there exists a qubit that is not reset, we set doResetAll to False.
                recordResetSet = set()
                if doResetAll:
                    # print("Into resetAll pruning at instruction index:", _)
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
            gate_condition = gate.condition
            tempNewCurrLoc = []
            if gate_condition is None:
                for cLoc in currLoc:
                    # Make new locations
                    new_loc = _add_location(ts, qnum, identifier + "S" + str(pivot + _ + 1), copy_from=cLoc)
                    # Add the operation to the transition system
                    ts.addRelation(cLoc, new_loc, op)
                    tempNewCurrLoc.append(new_loc)
            else:
                clbits_idx, clbits_vals = get_condition_info(qc.cregs, gate_condition)
                for cLoc in currLoc:
                    if _satisfy_bit(ts, cLoc, clbits_idx, clbits_vals):
                        # Make new locations
                        new_loc = _add_location(ts, qnum, identifier + "S" + str(pivot + _ + 1), copy_from=cLoc)
                        # Add the operation to the transition system
                        ts.addRelation(cLoc, new_loc, op)
                        tempNewCurrLoc.append(new_loc)
                    else:
                        # Apply identity operation
                        new_loc = _add_location(ts, qnum, identifier + "S" + str(pivot + _ + 1), copy_from=cLoc)
                        # Add the identity operation to the transition system
                        ts.addRelation(cLoc, new_loc, pyqreach.QOperation("I", qnum, qubits, []))
                        tempNewCurrLoc.append(new_loc)
            # Update the current locations
            currLoc = tempNewCurrLoc
    resultLocs = currLoc
    return resultLocs


def parse_qiskit_cir_sym(qc: QuantumCircuit, qnum: int, ts: pyqreach.SymTS, startNodes: list=None, identifier: str="", pivot: int=0, pivotend: int=1000000, abstractLevel: int=1) -> list:
    """
    Symbolic alias of parse_qiskit_cir for explicit SymTS call sites.
    """
    if startNodes is None:
        startNodes = []
    return parse_qiskit_cir(qc, qnum, ts, startNodes, identifier, pivot, pivotend, abstractLevel)


def build_sym_ts_from_qiskit(qc: QuantumCircuit, qnum: int | None = None, max_locations: int = 0, identifier: str = "", abstractLevel: int = 1):
    """
    Construct a symbolic transition system from a Qiskit circuit and return `(ts, end_locations)`.
    """
    if qnum is None:
        qnum = qc.num_qubits
    pyqreach.initializeSymTransitionSystem()
    ts = pyqreach.SymTS(qnum, max_locations)
    end_locs = parse_qiskit_cir_sym(qc, qnum, ts, identifier=identifier, abstractLevel=abstractLevel)
    return ts, end_locs

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

def visualize_transition_system(ts: pyqreach.TransitionSystem, filename='transition_system'):
    """
    Visualize the transition system using Graphviz.
    
    Args:
        ts (pyqreach.TransitionSystem): The transition system to visualize.
        filename (str): The name of the output file.
    """
    dot = Digraph(comment='Transition System')
    
    for loc in _get_location_ids(ts):
        dot.node(str(loc), label=str(loc))
    
        for post in _get_post_locations(ts, loc):
            rel = ts.getRelationName(loc, post)
            dot.edge(str(loc), str(post), label=str(rel))
    
    dot.render(filename, format='png', cleanup=True)


