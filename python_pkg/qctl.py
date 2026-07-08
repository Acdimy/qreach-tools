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
from inline_annotations import QReachCircuit, mark

class Proposition:
    def __init__(self, name: str, content=None, condition=None):
        self.name = name
        self.content = content
        self.condition = condition


def quantum_state(bitstring: str) -> pyqreach.QOperation:
    """Create a QOperation representing one simple product state.

    `bitstring` may contain 0/1 computational-basis symbols and +/- Hadamard-basis
    symbols.  '+' means H|0> and '-' means H|1>.
    """
    return pyqreach.QOperation([bitstring])


def zero_subspace(qnum: int) -> pyqreach.QOperation:
    """Create a QOperation representing the zero-dimensional quantum subspace."""
    return pyqreach.CreateZeroQO(qnum, False)

def whole_subspace(qnum: int) -> pyqreach.QOperation:
    """Create a QOperation representing the whole quantum subspace."""
    return pyqreach.CreateIdentityQO(qnum, False)

def span_qops(ops) -> pyqreach.QOperation:
    """Return the normalized span of several subspace QOperations.

    This is the direct API for constructing multi-dimensional quantum
    propositions from one-dimensional states, replacing the older temporary-TS
    workaround used to force Gram-Schmidt through computingFixedPointPost.
    """
    ops = list(ops)
    if not ops:
        raise ValueError("Cannot construct the span of an empty QOperation list")
    return pyqreach.span_qops(ops)


def span_states(states) -> pyqreach.QOperation:
    """Return the normalized span of simple product-state strings."""
    unique_states = []
    for state in states:
        if state not in unique_states:
            unique_states.append(state)
    return span_qops(quantum_state(state) for state in unique_states)


def _infer_qnum(ts, loc: int) -> int:
    if hasattr(ts, "getLocationAnnotation"):
        return ts.getLocationAnnotation(loc).qNum
    return ts.Locations[loc].qNum


def _set_quantum_annotations(ts, annotations: list[tuple[int, pyqreach.QOperation]]) -> list[int]:
    """Set quantum annotations on transition-system locations.

    Explicit TransitionSystem expects a batch ``[[loc, op], ...]`` while SymTS
    exposes ``setAnnotation(loc, op)``.  This helper keeps qctl-level APIs
    independent of that binding difference.
    """
    if not annotations:
        return []
    if _is_symbolic_ts(ts):
        for loc, op in annotations:
            ts.setAnnotation(loc, op)
    else:
        ts.setAnnotation([[loc, op] for loc, op in annotations])
    return [loc for loc, _ in annotations]


def set_initial_state(ts, bitstring: str, loc: int | None = None) -> pyqreach.QOperation:
    """Set a transition system location's initial quantum-state annotation.

    Defaults to the transition system's init location.  Returns the created
    QOperation so callers can reuse it for debugging or additional labelling.
    """
    op = quantum_state(bitstring)
    set_initial_operation(ts, op, loc=loc)
    return op


def set_initial_operation(ts, op: pyqreach.QOperation, loc: int | None = None) -> pyqreach.QOperation:
    """Set a transition system location's initial quantum-operation annotation.

    This is the operation-level dual of ``set_initial_state``.  It is useful when
    the desired initial annotation has already been built with ``span_qops``,
    ``span_states``, ``snapshot_operation``, or another QOperation builder.
    """
    if loc is None:
        loc = ts.getInitLocation()
    _set_quantum_annotations(ts, [(loc, op)])
    return op


def set_zero_initial_state(ts, qnum: int | None = None, loc: int | None = None) -> pyqreach.QOperation:
    """Set the init location to the all-zero computational-basis state."""
    if loc is None:
        loc = ts.getInitLocation()
    if qnum is None:
        qnum = _infer_qnum(ts, loc)
    return set_initial_state(ts, "0" * qnum, loc=loc)


def leaf_locations(ts, loc_list=None) -> list[int]:
    """Return the leaf locations in a transition system.

    If ``loc_list`` is given, only those locations are examined; otherwise all
    locations are scanned.
    """
    locs = list(_get_location_ids(ts) if loc_list is None else loc_list)
    return [loc for loc in locs if ts.isLeafLoc(loc)]


def annotate_leaf_operation(ts, op: pyqreach.QOperation, *, loc_list=None) -> list[int]:
    """Set the same quantum-operation annotation on all leaf locations.

    Returns the list of leaf locations that were annotated.  ``loc_list`` may be
    supplied to restrict the scan to a subset of locations.
    """
    locs = leaf_locations(ts, loc_list=loc_list)
    return _set_quantum_annotations(ts, [(loc, op) for loc in locs])


def annotate_leaf_state(ts, bitstring: str, *, loc_list=None) -> list[int]:
    """Set a simple product-state annotation on all leaf locations."""
    return annotate_leaf_operation(ts, quantum_state(bitstring), loc_list=loc_list)


def marker_locations(parse_result, marker: str) -> list[int]:
    """Return transition-system locations recorded for an inline QReach marker."""
    try:
        return list(parse_result.markers[marker])
    except KeyError as exc:
        raise KeyError(f"Unknown QReach mark: {marker}") from exc


def annotate_marker_operation(ts, parse_result, marker: str, op: pyqreach.QOperation) -> list[int]:
    """Set a quantum-operation annotation on all locations recorded for a marker."""
    locs = marker_locations(parse_result, marker)
    return _set_quantum_annotations(ts, [(loc, op) for loc in locs])


def annotate_marker_state(ts, parse_result, marker: str, bitstring: str) -> list[int]:
    """Set a simple product-state annotation on all locations recorded for a marker."""
    return annotate_marker_operation(ts, parse_result, marker, quantum_state(bitstring))


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


def _location_bound(ts, loc: int, bound: str = "lower") -> pyqreach.QOperation:
    if bound not in {"lower", "upper", "annotation"}:
        raise ValueError("bound must be one of: 'lower', 'upper', 'annotation'")
    if _is_symbolic_ts(ts):
        if bound == "annotation":
            return ts.getLocationAnnotation(loc)
        # TODO: SymTS currently exposes getLocationAnnotation rather than the
        # explicit lowerBound/upperBound pair.  Once symbolic pre/post workflow
        # is unified, decide how lower/upper snapshot names should map here.
        return ts.getLocationAnnotation(loc)
    if bound == "lower" or bound == "annotation":
        return ts.Locations[loc].lowerBound
    return ts.Locations[loc].upperBound


def _join_quantum_operations(ops: list[pyqreach.QOperation], qnum: int) -> pyqreach.QOperation:
    """Return the normalized span/disjunction of several QOperations."""
    if not ops:
        raise ValueError("Cannot join an empty QOperation list")
    return span_qops(ops)


def snapshot_operation(
    ts,
    locations: int | list[int],
    *,
    bound: str = "lower",
    merge: str = "single",
) -> pyqreach.QOperation | dict[int, pyqreach.QOperation]:
    """Build a quantum proposition from one or more transition-system locations.

    `merge` controls what happens when multiple locations are supplied:
    - 'single': require exactly one location.
    - 'first': use the first location.
    - 'join': build one QOperation spanning all selected location bounds.
    - 'per_location': return {loc: QOperation} without merging.
    """
    if isinstance(locations, int):
        locs = [locations]
    else:
        locs = list(locations)
    if not locs:
        raise ValueError("No snapshot locations were provided")

    if merge == "single":
        if len(locs) != 1:
            raise ValueError(f"Expected exactly one snapshot location, got {len(locs)}: {locs}")
        return _location_bound(ts, locs[0], bound)
    if merge == "first":
        return _location_bound(ts, locs[0], bound)
    if merge == "per_location":
        return {loc: _location_bound(ts, loc, bound) for loc in locs}
    if merge == "join":
        qnum = _infer_qnum(ts, locs[0])
        return _join_quantum_operations([_location_bound(ts, loc, bound) for loc in locs], qnum)
    raise ValueError("merge must be one of: 'single', 'first', 'join', 'per_location'")


def label_snapshot(
    ts,
    parse_result,
    marker: str,
    label: str | None = None,
    *,
    bound: str = "lower",
    merge: str = "single",
    locList: list | None = None,
) -> dict[str, list[int]]:
    """Label all TS locations satisfying the quantum proposition at a mark.

    This intentionally follows tsLabelling's semantics: after obtaining the
    snapshot QOperation from the marked location(s), it scans `locList` or all
    transition-system locations and labels every location satisfying that
    proposition.
    """
    if label is None:
        label = marker
    try:
        locations = parse_result.markers[marker]
    except KeyError as exc:
        raise KeyError(f"Unknown QReach mark: {marker}") from exc

    snapshot = snapshot_operation(ts, locations, bound=bound, merge=merge)
    labelled: dict[str, list[int]] = {}
    if isinstance(snapshot, dict):
        for loc, op in snapshot.items():
            loc_label = f"{label}_{loc}"
            before = {target: set(ts.getLabels(target)) for target in (range(ts.getLocationNum()) if locList is None else locList)}
            tsLabelling(ts, op, loc_label, locList=locList)
            labelled[loc_label] = [target for target, labels in before.items() if loc_label not in labels and loc_label in ts.getLabels(target)]
        return labelled

    iterList = list(range(ts.getLocationNum())) if locList is None else list(locList)
    before = {loc: set(ts.getLabels(loc)) for loc in iterList}
    tsLabelling(ts, snapshot, label, locList=iterList)
    labelled[label] = [loc for loc, labels in before.items() if label not in labels and label in ts.getLabels(loc)]
    return labelled

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
        targets = next_states if next_states else [s]
        smv += f"    state = {s} : {{ {', '.join(targets)} }};\n"
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


def _parse_counterexample_trace(cex_text: str) -> list[dict]:
    """Parse NuSMV counterexample text into a list of state dictionaries.

    Each dict contains ``state`` (location id as int) and any proposition
    assignments that appear on subsequent lines before the next state marker
    (e.g. ``leaf``, ``target``).
    """
    if not cex_text:
        return []

    # Split on state markers like "-> State: 1.17 <-"
    state_blocks = re.split(r"-> State: \d+\.\d+ <-", cex_text)
    trace: list[dict] = []

    for block in state_blocks:
        block = block.strip()
        if not block:
            continue
        entry: dict[str, object] = {}
        for line in block.splitlines():
            line = line.strip()
            m = re.match(r"state\s*=\s*(\d+)", line)
            if m:
                entry["state"] = int(m.group(1))
                continue
            m = re.match(r"(\w+)\s*=\s*(TRUE|FALSE)", line)
            if m:
                entry[m.group(1)] = m.group(2) == "TRUE"
        if "state" in entry:
            trace.append(entry)

    return trace


def _find_violating_step(trace: list[dict], antecedent: str, consequent: str) -> dict | None:
    """Return the first trace step where antecedent is TRUE and consequent is FALSE.

    In NuSMV counterexample output, propositions that are not listed on a state
    line are implicitly FALSE.
    """
    for step in trace:
        ant_val = step.get(antecedent, False)
        con_val = step.get(consequent, False)
        if ant_val is True and con_val is False:
            return step
    return None


def _trace_location_ids(trace: list[dict]) -> list[int]:
    """Extract the ordered list of QReach location ids from a parsed trace."""
    return [int(step["state"]) for step in trace]


def _identifier_for(ts, loc: int) -> str:
    """Return the parser-generated identifier for a location, or '' if none."""
    if hasattr(ts, "getIdentifier"):
        return ts.getIdentifier(loc) or ""
    try:
        return ts.Locations[loc].getIdentifier() or ""
    except Exception:
        return ""


def _labels_for(ts, loc: int) -> list[str]:
    """Return the labels set on a location."""
    try:
        return list(ts.getLabels(loc))
    except Exception:
        return []


def _analyse_counterexample(ts, cex_text: str, ctl_formula: str) -> dict | None:
    """Parse a NuSMV counterexample and map it back to QReach locations.

    Returns a dict with keys:
      - ``trace``: list of {location, identifier, labels} for each step
      - ``violating_location``: the first location where the implication fails
      - ``violating_identifier``: its parser-generated identifier
      - ``antecedent`` / ``consequent``: the two sides of the implication
      - ``edge_labels``: list of relation names along the trace path
    """
    trace = _parse_counterexample_trace(cex_text)
    if not trace:
        return None

    loc_ids = _trace_location_ids(trace)

    # Heuristic: extract antecedent/consequent from CTL formulas of the form
    #   AG (antecedent -> consequent)   or   AG (antecedent -> AX consequent)  etc.
    ant, con = None, None
    m = re.search(r"AG\s*\(\s*(\w+)\s*->\s*(?:\w+\s+)?(\w+)", ctl_formula)
    if m:
        ant, con = m.group(1), m.group(2)

    violating_step = None
    violating_loc = None
    if ant and con:
        violating_step = _find_violating_step(trace, ant, con)
        if violating_step:
            violating_loc = int(violating_step["state"])

    trace_info: list[dict] = []
    prev_loc: int | None = None
    edge_labels: list[str] = []

    for i, loc in enumerate(loc_ids):
        step_info: dict = {
            "step": i + 1,
            "location": loc,
            "identifier": _identifier_for(ts, loc),
            "labels": _labels_for(ts, loc),
        }
        # Capture proposition values from NuSMV trace when available
        if i < len(trace):
            for key, val in trace[i].items():
                if key not in ("state",):
                    step_info.setdefault("propositions", {})[key] = val
        trace_info.append(step_info)

        if prev_loc is not None:
            try:
                edge_labels.append(ts.getRelationName(prev_loc, loc))
            except Exception:
                edge_labels.append("")
        prev_loc = loc

    result: dict = {
        "trace": trace_info,
        "edge_labels": edge_labels,
    }

    if violating_loc is not None:
        result["violating_location"] = violating_loc
        result["violating_identifier"] = _identifier_for(ts, violating_loc)
        result["antecedent"] = ant
        result["consequent"] = con

    return result


def _format_counterexample_analysis(analysis: dict | None) -> str:
    """Render counterexample analysis as a human-readable string."""
    if analysis is None:
        return ""

    lines: list[str] = []
    viol = analysis.get("violating_location")
    if viol is not None:
        lines.append("=== Counterexample Analysis ===")
        lines.append(
            f"Violation: location {viol}"
            f" (identifier: {analysis.get('violating_identifier', '')!r})"
            f" satisfies {analysis.get('antecedent', '?')!r}"
            f" but NOT {analysis.get('consequent', '?')!r}."
        )
        lines.append("")

    lines.append("Trace (location → identifier → edge):")
    trace = analysis.get("trace", [])
    edges = analysis.get("edge_labels", [])
    for i, step in enumerate(trace):
        loc = step["location"]
        ident = step.get("identifier", "")
        edge = edges[i] if i < len(edges) else ""
        marker = " *** VIOLATION ***" if loc == viol else ""
        props = step.get("propositions", {})
        prop_str = " ".join(f"{k}={v}" for k, v in sorted(props.items())) if props else ""
        lines.append(
            f"  Step {step['step']:2d}: location {loc:4d}"
            f"  id={ident!r:30s}  --{edge}-->"
            + (f"  [{prop_str}]" if prop_str else "")
            + marker
        )

    return "\n".join(lines)


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
        return {'satisfied': None, 'counterexample': None, 'analysis': None, 'output': output}
    except FileNotFoundError:
        output = f'Error: {nusmv_cmd} not found. Please verify the path or ensure NuSMV is in PATH.'
        print(output)
        return {'satisfied': None, 'counterexample': None, 'analysis': None, 'output': output}
    except Exception as e:
        output = f'Error running NuSMV: {str(e)}'
        print(output)
        return {'satisfied': None, 'counterexample': None, 'analysis': None, 'output': output}
    finally:
        os.unlink(temp_file_name)

    # Parse the output to check if the specification is satisfied
    match = re.search(r"-- specification (.+) is (true|false)", output, re.MULTILINE | re.DOTALL)
    if match:
        satisfied = match.group(2) == 'true'
        counterexample = None
        analysis = None
        if not satisfied:
            # Extract counterexample if the specification is false
            cex_start = output.find("-- as demonstrated by the following execution sequence")
            if cex_start != -1:
                cex_end = output.find("********", cex_start)  # counterexample ends with a line of asterisks
                counterexample = output[cex_start:cex_end].strip() if cex_end != -1 else output[cex_start:].strip()
            analysis = _analyse_counterexample(ts, counterexample, ctl_formula)
        return {
            'satisfied': satisfied,
            'counterexample': counterexample,
            'analysis': analysis,
            'output': output,
        }
    else:
        print(f"Unexpected NuSMV output:\n{output}")
        return {'satisfied': None, 'counterexample': None, 'analysis': None, 'output': f'Unexpected output: {output}'}

