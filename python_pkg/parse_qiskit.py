import pyqreach
from qiskit import QuantumCircuit



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
