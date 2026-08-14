from qreach import symbolic_available

if not symbolic_available():
    print("SKIP: symbolic SymTS backend unavailable (LimTDD build)")
    raise SystemExit(0)

import pyqreach
from qiskit import QuantumCircuit, QuantumRegister, ClassicalRegister

from qreach.parse_qiskit import build_sym_ts_from_qiskit, parse_qiskit_cir
from qreach.qctl import tsLabellingClRegList, tsLabellingDefault


def build_demo_circuit() -> QuantumCircuit:
    q = QuantumRegister(2, "q")
    c = ClassicalRegister(1, "c")
    qc = QuantumCircuit(q, c)
    qc.h(q[0])
    qc.measure(q[0], c[0])
    return qc


def run_naive(qc: QuantumCircuit):
    pyqreach.initializeTransitionSystem()
    ts = pyqreach.TransitionSystem(False)
    end_locs = parse_qiskit_cir(qc, qc.num_qubits, ts)
    ts.setAnnotation([[0, pyqreach.QOperation(["0" * qc.num_qubits])]])
    ts.computingFixedPointPost()
    tsLabellingDefault(ts, "reachable")
    tsLabellingClRegList(ts, ["1"], "branch_one", locList=end_locs)
    return ts, end_locs


def run_symbolic(qc: QuantumCircuit):
    ts, end_locs = build_sym_ts_from_qiskit(qc, qnum=qc.num_qubits, max_locations=8)
    ts.setAnnotation(0, pyqreach.QOperation(["0" * qc.num_qubits]))
    ts.computingFixedPointPost()
    tsLabellingDefault(ts, "reachable")
    tsLabellingClRegList(ts, ["1"], "branch_one", locList=end_locs)
    return ts, end_locs


def summarize_naive(ts, end_locs):
    print("=== Naive TS ===")
    print("end_locs:", end_locs)
    for loc in end_locs:
        cp = ts.Locations[loc].cp.toString()
        identifier = ts.Locations[loc].getIdentifier()
        labels = ts.getLabels(loc)
        support = ts.Locations[loc].lowerBound.printFormal(False)
        dims = ts.printDims(loc)
        print(f"loc={loc} id={identifier} cp={cp} dims={dims} labels={labels}")
        print(f"  support={support}")


def summarize_symbolic(ts, end_locs):
    print("=== SymTS ===")
    print("end_locs:", end_locs)
    for loc in end_locs:
        cp = ts.getClassicalProposition(loc).toString()
        identifier = ts.getIdentifier(loc)
        labels = ts.getLabels(loc)
        support = ts.getLocationAnnotation(loc).printFormal(False)
        dims = ts.printDims(loc)
        print(f"loc={loc} id={identifier} cp={cp} dims={dims} labels={labels}")
        print(f"  support={support}")


def compare_results(naive_ts, naive_end_locs, sym_ts, sym_end_locs):
    print("=== Diff Check ===")
    naive_summary = {}
    for naive_loc in naive_end_locs:
        naive_cp = naive_ts.Locations[naive_loc].cp.toString()
        naive_summary[naive_cp] = {
            "loc": naive_loc,
            "identifier": naive_ts.Locations[naive_loc].getIdentifier(),
            "support": naive_ts.Locations[naive_loc].lowerBound.printFormal(False),
            "labels": sorted(naive_ts.getLabels(naive_loc)),
            "lower_dim": naive_ts.printDims(naive_loc)[1],
        }

    sym_summary = {}
    for sym_loc in sym_end_locs:
        sym_cp = sym_ts.getClassicalProposition(sym_loc).toString()
        sym_summary[sym_cp] = {
            "loc": sym_loc,
            "identifier": sym_ts.getIdentifier(sym_loc),
            "support": sym_ts.getLocationAnnotation(sym_loc).printFormal(False),
            "labels": sorted(sym_ts.getLabels(sym_loc)),
            "lower_dim": sym_ts.printDims(sym_loc)[1],
        }

    print("naive cp keys:", sorted(naive_summary.keys()))
    print("sym cp keys:", sorted(sym_summary.keys()))
    assert sorted(naive_summary.keys()) == sorted(sym_summary.keys()), "Classical proposition partitions differ between naive TS and SymTS"

    for cp_key in sorted(naive_summary.keys()):
        naive_info = naive_summary[cp_key]
        sym_info = sym_summary[cp_key]
        identifier_match = naive_info["identifier"] == sym_info["identifier"]
        support_match = naive_info["support"] == sym_info["support"]
        labels_match = naive_info["labels"] == sym_info["labels"]
        lower_dim_match = naive_info["lower_dim"] == sym_info["lower_dim"]
        print(
            f"cp={cp_key}: naive_loc={naive_info['loc']} sym_loc={sym_info['loc']} "
            f"id_match={identifier_match} support_match={support_match} "
            f"labels_match={labels_match} lower_dim_match={lower_dim_match}"
        )
        assert identifier_match, f"Identifier mismatch for cp={cp_key}"
        assert support_match, f"Support mismatch for cp={cp_key}"
        assert labels_match, f"Label mismatch for cp={cp_key}"
        assert lower_dim_match, f"Reachable dimension mismatch for cp={cp_key}"


def main():
    qc = build_demo_circuit()
    print(qc)

    naive_ts, naive_end_locs = run_naive(qc)
    sym_ts, sym_end_locs = run_symbolic(qc)

    summarize_naive(naive_ts, naive_end_locs)
    summarize_symbolic(sym_ts, sym_end_locs)
    compare_results(naive_ts, naive_end_locs, sym_ts, sym_end_locs)
    print("Minimal SymTS check passed.")


if __name__ == "__main__":
    main()