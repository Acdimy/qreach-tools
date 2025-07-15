from qiskit import QuantumCircuit, QuantumRegister, ClassicalRegister
import pyqreach

### Duplicated: seems useless

class ModChecker:
    """
    A singleton class for the total model checking process.
    """
    def __init__(self):
        self.classicalRegisterList = []
        self.ts = None
        self.qnum = 0
    # Target: when build a quantum circuit, we should record the rough structure of the transition system manully.
    # Esspecially the branch, merge and loop structure.
    def parse_qiskit_cir(self):
        pass
    
"""
from qiskit import QuantumCircuit, ClassicalRegister, QuantumRegister
from qiskit.circuit.controlflow import IfElseOp, WhileLoopOp, ForLoopOp, SwitchCaseOp
from typing import List, Dict, Any

class CFGNode:
    def __init__(self, id: int, instructions: List[Any] = None):
        self.id = id
        self.instructions = instructions if instructions else []
        self.successors = []  # List of CFGNode

    def add_successor(self, node):
        self.successors.append(node)

    def __repr__(self):
        return f"Node({self.id}, instr={len(self.instructions)}, succ={[n.id for n in self.successors]})"

class CFGBuilder:
    def __init__(self, circuit: QuantumCircuit):
        self.circuit = circuit
        self.nodes = []
        self.node_id_counter = 0

    def new_node(self, instructions=None):
        node = CFGNode(self.node_id_counter, instructions)
        self.node_id_counter += 1
        self.nodes.append(node)
        return node

    def build(self):
        entry = self.new_node()
        self._build_block(self.circuit.data, entry)
        return entry, self.nodes

    def _build_block(self, instructions, current_node):
        i = 0
        while i < len(instructions):
            instr, qargs, cargs = instructions[i]

            if isinstance(instr, IfElseOp):
                then_node = self.new_node()
                else_node = self.new_node()
                after_node = self.new_node()

                self._build_block(instr.blocks[0], then_node)
                then_node.add_successor(after_node)

                if len(instr.blocks) > 1:
                    self._build_block(instr.blocks[1], else_node)
                else:
                    else_node.instructions.append(('pass',))
                else_node.add_successor(after_node)

                current_node.add_successor(then_node)
                current_node.add_successor(else_node)
                current_node = after_node

            elif isinstance(instr, WhileLoopOp):
                loop_cond_node = self.new_node()
                loop_body_node = self.new_node()
                after_loop_node = self.new_node()

                self._build_block(instr.blocks[0], loop_body_node)
                loop_body_node.add_successor(loop_cond_node)

                current_node.add_successor(loop_cond_node)
                loop_cond_node.add_successor(loop_body_node)
                loop_cond_node.add_successor(after_loop_node)

                current_node = after_loop_node

            elif isinstance(instr, ForLoopOp):
                loop_header_node = self.new_node()
                loop_body_node = self.new_node()
                after_loop_node = self.new_node()

                self._build_block(instr.blocks[0], loop_body_node)
                loop_body_node.add_successor(loop_header_node)

                current_node.add_successor(loop_header_node)
                loop_header_node.add_successor(loop_body_node)
                loop_header_node.add_successor(after_loop_node)

                current_node = after_loop_node

            elif isinstance(instr, SwitchCaseOp):
                case_nodes = []
                after_switch_node = self.new_node()

                for block in instr.blocks:
                    case_node = self.new_node()
                    self._build_block(block, case_node)
                    case_node.add_successor(after_switch_node)
                    case_nodes.append(case_node)

                for cn in case_nodes:
                    current_node.add_successor(cn)

                current_node = after_switch_node

            else:
                current_node.instructions.append((instr, qargs, cargs))

            i += 1

        return current_node

# Example usage:
qr = QuantumRegister(1)
cr = ClassicalRegister(1)
circ = QuantumCircuit(qr, cr)

# Dummy example circuit (no real control flow added here)
circ.h(qr[0])
circ.measure(qr[0], cr[0])

cfg_builder = CFGBuilder(circ)
entry_node, all_nodes = cfg_builder.build()

for node in all_nodes:
    print(node)  # For debugging
"""