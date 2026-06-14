import pyqreach
from qiskit import QuantumCircuit, QuantumRegister, ClassicalRegister
import re
from annotations import (
    AnnotationRegistry,
    annotate,
    annotate_classical,
    annotate_identifier,
    annotate_where,
    default_registry,
)

class Proposition:
    def __init__(self, name: str, content=None, condition=None):
        self.name = name
        self.content = content
        self.condition = condition


def _is_symbolic_ts(ts) -> bool:
    return hasattr(ts, "getLocationIDs") and hasattr(ts, "getLocationAnnotation")


def _get_location_ids(ts):
    if _is_symbolic_ts(ts):
        return ts.getLocationIDs()
    return [loc.idx for loc in ts.Locations]


def _get_post_locations(ts, loc):
    if _is_symbolic_ts(ts):
        return ts.getPostLocations(loc)
    return list(ts.Locations[loc].postLocations)


def _satisfy_quantum(ts, loc, op):
    if _is_symbolic_ts(ts):
        return ts.satisfy(loc, op)
    return ts.Locations[loc].satisfy(op)


def _satisfy_bit(ts, loc, idxs, vals):
    if _is_symbolic_ts(ts):
        return ts.satisfyBit(loc, idxs, vals)
    return ts.Locations[loc].satisfyBit(idxs, vals)

def tsLabelling(ts, op: pyqreach.QOperation, label: str, locList: list=None):
    iterList = range(ts.getLocationNum()) if locList is None else locList
    for loc in iterList:
        if _satisfy_quantum(ts, loc, op):
            ts.setLabel(loc, label)
            # print(f"Location {loc} labelled with {label}")

def tsLabellingDefault(ts, label: str, locList: list=None):
    iterList = range(ts.getLocationNum()) if locList is None else locList
    for loc in iterList:
        # if ts.Locations[loc].satisfyDefault():
        if ts.printDims(loc)[1] > 0:
            ts.setLabel(loc, label)

def tsLabellingClRegList(ts, clRegList: list, label: str, locList:list=None):
    """
    Label locations in the transition system based on classical register values: |= BigVee clRegList.
    :param ts: Transition system
    :param clRegList: List of strings as classical registers
    :label: Label to assign to the locations that satisfy the classical register values
    """
    clRegBins = [[int(bit) for bit in clReg] for clReg in clRegList]
    iterList = range(ts.getLocationNum()) if locList is None else locList
    for loc in iterList:
        # Once one of the clReg in clRegList is satisfied, label the location and break
        for clRegBin in clRegBins:
            # convert clReg to a binary list
            if len(_satisfy_bit(ts, loc, list(range(len(clRegBin))), clRegBin)) != 0:
                ts.setLabel(loc, label)
                break

def BellProposition(num_qubits: int, qubit_indices: list, Bell_idx: int) -> pyqreach.QOperation:
    """
    Construct a proposition for Bell state on specified qubits.
    Bell_idx: 0 for |Φ+>, 1 for |Φ->, 2 for |Ψ+>, 3 for |Ψ->
    """
    assert num_qubits == 3 and qubit_indices == [0, 1], "BellProposition currently only supports 2 qubits at indices [0, 1] in a 3-qubit system."
    if len(qubit_indices) != 2:
        raise ValueError("BellProposition requires exactly 2 qubits.")
    if Bell_idx not in {0, 1, 2, 3}:
        raise ValueError("Bell_idx must be in {0, 1, 2, 3}.")
    ts_temp = pyqreach.TransitionSystem(False)
    loc0, loc1, loc2 = pyqreach.Location(num_qubits,0), pyqreach.Location(num_qubits,1), pyqreach.Location(num_qubits,2)
    ts_temp.addLocation(loc0)
    ts_temp.addLocation(loc1)
    ts_temp.addLocation(loc2)
    ts_temp.addRelation(0, 1, pyqreach.QOperation("H", num_qubits, [qubit_indices[0]], []))
    ts_temp.addRelation(1, 2, pyqreach.QOperation("CX", num_qubits, [qubit_indices[0], qubit_indices[1]], []))
    if Bell_idx == 0:
        ts_temp.setAnnotation([[0, pyqreach.QOperation(["000"])]])
    elif Bell_idx == 1:
        ts_temp.setAnnotation([[0, pyqreach.QOperation(["100"])]])
    elif Bell_idx == 2:
        ts_temp.setAnnotation([[0, pyqreach.QOperation(["010"])]])
    elif Bell_idx == 3:
        ts_temp.setAnnotation([[0, pyqreach.QOperation(["110"])]])
    ts_temp.computingFixedPointPost()
    bell_op = ts_temp.Locations[2].lowerBound
    return bell_op

def labelling(ts: pyqreach.TransitionSystem, propositions: list):
    pass

# 可扩展的保留关键字集合（大小写不敏感）
CTL_NUSMV_RESERVED = {
    # CTL / temporal operators (常见写法)
    "A","E","AX","EX","AF","EF","AG","EG","AU","EU",
    "X","F","G","U","R",
    # 布尔 / 逻辑 / 常量
    "NOT","AND","OR","XOR","IMPLIES","TRUE","FALSE",
    # NuSMV / SMV 中常见关键字（预防被误识别）
    "MODULE","VAR","DEFINE","SPEC","ASSIGN","INIT","NEXT","CASE","ESAC",
    # 你在模型中可能会用到但不应当被当成原子命题的名字
    "state"
}

# 小写化的集合，便于不区分大小写匹配
_reserved_lc = {w.lower() for w in CTL_NUSMV_RESERVED}

def extract_atoms_from_ctl(ctl_formula: str, extra_reserved: set = None) -> set:
    """
    从 ctl_formula 中提取可能的原子命题标识符（不包含保留关键字）。
    返回一个标识符集合（原样大小写保留）。
    """
    if extra_reserved is None:
        extra_reserved = set()
    extra_reserved_lc = {w.lower() for w in extra_reserved}

    # 匹配标识符的正则（以字母或下划线开头，后面字母数字或下划线）
    tokens = set(re.findall(r'\b[A-Za-z_][A-Za-z0-9_]*\b', ctl_formula))

    atoms = set()
    for t in tokens:
        tl = t.lower()
        # 过滤掉保留关键字和额外排除项
        if tl in _reserved_lc or tl in extra_reserved_lc:
            continue
        # 过滤掉布尔常量（重复保险）
        if tl in ("true", "false"):
            continue
        # 现在剩下的基本都是用户定义的原子命题或变量名
        atoms.add(t)
    return atoms

def ts2Dict(ts: pyqreach.TransitionSystem) -> dict:
    """
    Convert the transition system to a dictionary representation.
    
    Args:
        ts (pyqreach.TransitionSystem): The transition system to convert.
    
    Returns:
        dict: Dictionary representation of the transition system.
    """
    labelDict = {}
    loc_ids = _get_location_ids(ts)
    for loc in loc_ids:
        loclabels = ts.getLabels(loc)
        # print(f"Location {loc.idx} labels: {loclabels}")
        labelDict[str(loc)] = loclabels if loclabels else []
    ts_dict = {
        'locationsTuple': [(loc, ts.printDims(loc)[1]>0) for loc in loc_ids],
        'locations': [str(loc) for loc in loc_ids],
        'relations': {
            f"{src}->{dst}": ts.getRelationName(src, dst)
            for src in loc_ids
            for dst in _get_post_locations(ts, src)
        },
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
            with_labels=False, # 显示节点标签
            node_size=20,
            # node label size small
            font_size=5,
            node_color=node_colors,
            edgecolors="lightgrey",
            edge_color="grey",
            width=0.2,
            alpha=0.9,
            arrows=False,
            linewidths=0.1)
    # nx.draw_networkx_edge_labels(G, pos, edge_labels=labels, font_size=9)
    plt.axis("off")
    # plt.tight_layout()
    # save as pdf
    plt.savefig(filename + '.pdf', dpi=300)    
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
    # Collect labels that occured in the transition system
    prop_states = {}
    for s, props in dts['labels'].items():
        for p in props:
            prop_states.setdefault(p, []).append(s)
    # Collect labels that never occured in TS but are in the ctl_formula
    atoms_in_formula = extract_atoms_from_ctl(ctl_formula, extra_reserved={'state', 'init', 'next'})
    smv += "DEFINE\n"
    for p, states in prop_states.items():
        cond = " | ".join(f"state = {st}" for st in states)
        smv += f"  {p} := {cond};\n"
    # Add false definitions for atoms not in the TS
    missing_atoms = atoms_in_formula - set(prop_states.keys())
    for p in sorted(missing_atoms):
        smv += f"  {p} := FALSE;\n"
    smv += f"SPEC\n  {ctl_formula};\n"
    return smv

import subprocess

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

def modelChecking(ts: pyqreach.TransitionSystem, ctl_formula: str, nusmv_path='../../NuSMV-2.7.0-macos-universal/bin/NuSMV'):
    # ../NuSMV-2.7.0-linux64/bin/NuSMV
    smv_code = ts2SMV(ts, ctl_formula)
    nusmv_cmd = nusmv_path if nusmv_path else 'NuSMV'
    with tempfile.NamedTemporaryFile(mode='w', suffix='.smv', delete=False) as temp_file:
        temp_file.write(smv_code)
        temp_file_name = temp_file.name
    try:
        result = subprocess.run([nusmv_cmd, temp_file_name], capture_output=True, text=True, timeout=240)
        output = result.stdout + result.stderr
    except subprocess.TimeoutExpired:
        output = 'Timeout: NuSMV took too long to respond.'
        print(output)
        return {'satisfied': None, 'counterexample': None, 'output': output}
    except FileNotFoundError:
        output = f'Error: {nusmv_cmd} not found. Please verify the path or ensure NuSMV is in PATH.'
        print(output)
        return {'satisfied': None, 'counterexample': None, 'output': output}
    except Exception as e:
        output = f'Error running NuSMV: {str(e)}'
        print(output)
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
        print(f"Unexpected NuSMV output:\n{output}")
        return {'satisfied': None, 'counterexample': None, 'output': f'Unexpected output: {output}'}

