import pyqreach
from qiskit import QuantumCircuit, QuantumRegister, ClassicalRegister

class Proposition:
    def __init__(self, name: str, content=None, condition=None):
        self.name = name
        self.content = content
        self.condition = condition

def tsLabelling(ts, op: pyqreach.QOperation, label: str):
    for loc in range(ts.getLocationNum()):
        if ts.Locations[loc].satisfy(op):
            ts.setLabel(loc, label)
            # print(f"Location {loc} labelled with {label}")

def tsLabellingDefault(ts, label: str):
    for loc in range(ts.getLocationNum()):
        # if ts.Locations[loc].satisfyDefault():
        if ts.printDims(loc)[1] > 0:
            ts.setLabel(loc, label)

def tsLabellingClRegList(ts, clRegList: list, label: str):
    """
    Label locations in the transition system based on classical register values: |= BigVee clRegList.
    :param ts: Transition system
    :param clRegList: List of strings as classical registers
    :label: Label to assign to the locations that satisfy the classical register values
    """
    for loc in range(ts.getLocationNum()):
        for clReg in clRegList:
            # convert clReg to a binary list
            clRegBin = [int(bit) for bit in clReg]
            if ts.Locations[loc].satisfyBit(list(range(len(clRegBin))), clRegBin):
                ts.setLabel(loc, label)
                break

def labelling(ts: pyqreach.TransitionSystem, propositions: list):
    pass

def ts2Dict(ts: pyqreach.TransitionSystem) -> dict:
    """
    Convert the transition system to a dictionary representation.
    
    Args:
        ts (pyqreach.TransitionSystem): The transition system to convert.
    
    Returns:
        dict: Dictionary representation of the transition system.
    """
    labelDict = {}
    for loc in ts.Locations:
        loclabels = ts.getLabels(loc.idx)
        # print(f"Location {loc.idx} labels: {loclabels}")
        labelDict[str(loc.idx)] = loclabels if loclabels else []
    ts_dict = {
        'locationsTuple': [(loc.idx, ts.printDims(loc.idx)[1]>0) for loc in ts.Locations],
        'locations': [str(loc.idx) for loc in ts.Locations],
        'relations': {f"{rel[0]}->{rel[1]}": ts.getRelationName(rel[0], rel[1]) for rel, op in ts.relations.items()},
        'init_location': str(ts.getInitLocation()),
        'num_locations': str(ts.getLocationNum()),
        'labels': labelDict
    }
    return ts_dict

# def ts2Dict_simp(ts: pyqreach.TransitionSystem) -> dict:
#     """
#     Simplified version of ts2Dict, only includes locations and relations.
    
#     Args:
#         ts (pyqreach.TransitionSystem): The transition system to convert.
    
#     Returns:
#         dict: Simplified dictionary representation of the transition system.
#     """
#     ts_dict = {
#         'locations': [str(loc.idx) for loc in ts.Locations],
#         'relations': {f"{rel[0]}->{rel[1]}": ts.getRelationName(rel[0], rel[1]) for rel, op in ts.relations.items()},
#         'init_location': str(ts.getInitLocation()),
#         'num_locations': str(ts.getLocationNum())
#     }
#     return ts_dict

import networkx as nx
import matplotlib.pyplot as plt

def dict2NX(dts: dict) -> nx.DiGraph:
    """
    Convert a dictionary representation of a transition system to a NetworkX directed graph.
    
    Args:
        dts (dict): Dictionary representation of the transition system.
    
    Returns:
        nx.DiGraph: NetworkX directed graph representation of the transition system.
    """
    G = nx.DiGraph()
    
    # Add nodes
    for loc,highlight in dts['locationsTuple']:
        G.add_node(loc, label=str(loc), highlight=highlight)
    
    # Add edges
    for rel, op in dts['relations'].items():
        src, dst = map(int, rel.split('->'))
        G.add_edge(src, dst, label=op)
    
    return G

def nx2Graph(G: nx.DiGraph, filename='transition_system', layout='spring'):
    """
    Visualize a NetworkX directed graph using matplotlib.

    Args:
        G (nx.DiGraph): The directed graph to visualize.
        filename (str): The name of the output file.
        layout (str): Layout type: 'spring', 'circular', 'shell', 'kamada_kawai', 'spectral'.
    """
    # 选择布局
    if layout == 'spring':
        pos = nx.spring_layout(G)
    elif layout == 'shell':
        pos = nx.shell_layout(G)
    elif layout == 'kamada_kawai':
        pos = nx.kamada_kawai_layout(G)
    elif layout == 'spectral':
        pos = nx.spectral_layout(G)
    else:
        raise ValueError(f"Unsupported layout: {layout}")

    # 边标签
    labels = nx.get_edge_attributes(G, 'label')

    plt.figure(figsize=(10, 6))
    nx.draw(
        G, pos,
        with_labels=False,         # 不显示节点标签
        node_size=10,             # 节点尺寸
        node_color='grey',    # 节点颜色
        font_size=10,              # 字体大小
        font_color='black',        # 节点字体颜色
        arrows=False
        # edgecolors='black',        # 节点轮廓颜色
        # linewidths=1               # 节点轮廓宽度
    )
    # nx.draw_networkx_edge_labels(G, pos, edge_labels=labels, font_size=9)

    plt.axis('off')
    # plt.tight_layout()
    plt.savefig(filename + '.png', dpi=300)
    plt.close()

from networkx.drawing.nx_agraph import graphviz_layout
def nx2Graph_hierarchical(G, filename="tree_layout"):
    pos = graphviz_layout(G, prog="dot")  # 分层布局
    labels = nx.get_edge_attributes(G, 'label')
    node_colors = ['orange' if G.nodes[n].get('highlight', False) else 'grey' for n in G.nodes()]
    plt.figure(figsize=(8, 6))
    nx.draw(G, pos,
            with_labels=True, # 显示节点标签
            node_size=20,
            # node label size small
            font_size=5,
            node_color=node_colors,
            edgecolors="black",
            arrows=False,
            linewidths=0.3)
    # nx.draw_networkx_edge_labels(G, pos, edge_labels=labels, font_size=9)
    plt.axis("off")
    # plt.tight_layout()
    plt.savefig(filename + ".png", dpi=300)
    plt.close()

def dict2SMV(dts, ctl_formula):
    smv = "MODULE main\nVAR\n  state: {" + ", ".join(dts['locations']) + "};\n"
    smv += "ASSIGN\n  init(state) := " + dts['init_location'] + ";\n"
    smv += "  next(state) := case\n"
    for s in dts['locations']:
        next_states = [t.split('->')[1] for t in dts['relations'] if t.split('->')[0] == s]
        if next_states:
            smv += f"    state = {s} : {{ {', '.join(next_states)} }};\n"
    smv += "    TRUE : state;\n  esac;\n"
    prop_states = {}
    for s, props in dts['labels'].items():
        for p in props:
            prop_states.setdefault(p, []).append(s)
    smv += "DEFINE\n"
    for p, states in prop_states.items():
        cond = " | ".join(f"state = {st}" for st in states)
        smv += f"  {p} := {cond};\n"
    smv += f"SPEC\n  {ctl_formula};\n"
    return smv

import subprocess
import re

def ts2SMV(ts: pyqreach.TransitionSystem, ctl_formula: str):
    """
    Convert a transition system to SMV format and save it to a file.
    
    Args:
        ts (pyqreach.TransitionSystem): The transition system to convert.
        ctl_formula (str): The CTL formula to include in the SMV file.
        filename (str): The name of the output file (without extension).
    
    Returns:
        str: The SMV content as a string.
    """
    dts = ts2Dict(ts)
    smv_content = dict2SMV(dts, ctl_formula)
    return smv_content

import tempfile
import os

def modelChecking(ts: pyqreach.TransitionSystem, ctl_formula: str, nusmv_path='../NuSMV-2.7.0-linux64/bin/NuSMV'):
    smv_code = ts2SMV(ts, ctl_formula)
    nusmv_cmd = nusmv_path if nusmv_path else 'NuSMV'
    with tempfile.NamedTemporaryFile(mode='w', suffix='.smv', delete=False) as temp_file:
        temp_file.write(smv_code)
        temp_file_name = temp_file.name
    try:
        result = subprocess.run([nusmv_cmd, temp_file_name], capture_output=True, text=True, timeout=60)
        output = result.stdout + result.stderr
    except subprocess.TimeoutExpired:
        output = 'Timeout: NuSMV took too long to respond.'
        return {'satisfied': None, 'counterexample': None, 'output': output}
    except FileNotFoundError:
        output = f'Error: {nusmv_cmd} not found. Please verify the path or ensure NuSMV is in PATH.'
        return {'satisfied': None, 'counterexample': None, 'output': output}
    except Exception as e:
        output = f'Error running NuSMV: {str(e)}'
        return {'satisfied': None, 'counterexample': None, 'output': output}
    finally:
        os.unlink(temp_file_name)
    
    # Parse the output to check if the specification is satisfied
    match = re.search(r"-- specification (.+) is (true|false)", output, re.MULTILINE | re.DOTALL)
    if match:
        satisfied = match.group(2) == 'true'
        counterexample = None
        if not satisfied:
            # Extract counterexample if the specification is false
            cex_start = output.find("-- as demonstrated by the following execution sequence")
            if cex_start != -1:
                cex_end = output.find("********", cex_start)  # counterexample ends with a line of asterisks
                counterexample = output[cex_start:cex_end].strip() if cex_end != -1 else output[cex_start:].strip()
        return {'satisfied': satisfied, 'counterexample': counterexample, 'output': output}
    else:
        return {'satisfied': None, 'counterexample': None, 'output': f'Unexpected output: {output}'}

