import pyqreach
from qiskit import QuantumCircuit, QuantumRegister, ClassicalRegister

from parse_qiskit import build_sym_ts_from_qiskit, parse_qiskit_cir
from qctl import tsLabellingClRegList, tsLabellingDefault


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
    if naive_end_locs != sym_end_locs:
        print("end location ids differ:", naive_end_locs, sym_end_locs)
    else:
        print("end location ids match:", naive_end_locs)

    for naive_loc, sym_loc in zip(naive_end_locs, sym_end_locs):
        naive_cp = naive_ts.Locations[naive_loc].cp.toString()
        sym_cp = sym_ts.getClassicalProposition(sym_loc).toString()
        naive_id = naive_ts.Locations[naive_loc].getIdentifier()
        sym_id = sym_ts.getIdentifier(sym_loc)
        naive_support = naive_ts.Locations[naive_loc].lowerBound.printFormal(False)
        sym_support = sym_ts.getLocationAnnotation(sym_loc).printFormal(False)
        print(f"loc {naive_loc}/{sym_loc}: cp_match={naive_cp == sym_cp}, id_match={naive_id == sym_id}, support_match={naive_support == sym_support}")


def main():
    qc = build_demo_circuit()
    print(qc)

    naive_ts, naive_end_locs = run_naive(qc)
    sym_ts, sym_end_locs = run_symbolic(qc)

    summarize_naive(naive_ts, naive_end_locs)
    summarize_symbolic(sym_ts, sym_end_locs)
    compare_results(naive_ts, naive_end_locs, sym_ts, sym_end_locs)


if __name__ == "__main__":
    main()