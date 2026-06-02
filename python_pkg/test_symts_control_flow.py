import pyqreach
import subprocess
import sys
from qiskit import QuantumCircuit, QuantumRegister, ClassicalRegister

from parse_qiskit import build_sym_ts_from_qiskit, parse_qiskit_cir


def build_rus_while_circuit() -> QuantumCircuit:
    qc = QuantumCircuit(2, 1)
    qc.x(1)
    qc.measure(1, 0)
    with qc.while_loop((0, 0b1)):
        qc.h(0)
        qc.measure(0, 0)
    return qc


def build_if_test_circuit() -> QuantumCircuit:
    q = QuantumRegister(2, "q")
    c = ClassicalRegister(1, "c")
    qc = QuantumCircuit(q, c)
    qc.h(q[0])
    qc.measure(q[0], c[0])
    with qc.if_test((c[0], 1)):
        qc.x(q[1])
    return qc


def run_naive(qc: QuantumCircuit):
    pyqreach.initializeTransitionSystem()
    ts = pyqreach.TransitionSystem(False)
    end_locs = parse_qiskit_cir(qc, qc.num_qubits, ts)
    ts.setAnnotation([[0, pyqreach.QOperation(["0" * qc.num_qubits])]])
    ts.computingFixedPointPost()
    return ts, end_locs


def run_symbolic(qc: QuantumCircuit, max_locations: int = 0):
    ts, end_locs = build_sym_ts_from_qiskit(qc, qnum=qc.num_qubits, max_locations=max_locations)
    ts.setAnnotation(0, pyqreach.QOperation(["0" * qc.num_qubits]))
    ts.computingFixedPointPost()
    return ts, end_locs


def collect_naive_end_summary(ts, end_locs):
    records = []
    for loc in end_locs:
        records.append(
            {
                "loc": loc,
                "cp": ts.Locations[loc].cp.toString(),
                "identifier": ts.Locations[loc].getIdentifier(),
                "support": ts.Locations[loc].lowerBound.printFormal(False),
                "lower_dim": ts.printDims(loc)[1],
            }
        )
    return sorted(records, key=lambda item: (item["loc"], item["cp"], item["identifier"], item["support"]))


def collect_symbolic_end_summary(ts, end_locs):
    records = []
    for loc in end_locs:
        records.append(
            {
                "loc": loc,
                "cp": ts.getClassicalProposition(loc).toString(),
                "identifier": ts.getIdentifier(loc),
                "support": ts.getLocationAnnotation(loc).printFormal(False),
                "lower_dim": ts.printDims(loc)[1],
            }
        )
    return sorted(records, key=lambda item: (item["loc"], item["cp"], item["identifier"], item["support"]))


def assert_match(name: str, naive_ts, naive_end_locs, sym_ts, sym_end_locs):
    naive_summary = collect_naive_end_summary(naive_ts, naive_end_locs)
    sym_summary = collect_symbolic_end_summary(sym_ts, sym_end_locs)

    print(f"=== {name} ===")
    print("naive locations:", naive_ts.getLocationNum(), "end_locs:", naive_end_locs)
    print("sym locations:", sym_ts.getLocationNum(), "end_locs:", sym_end_locs)
    print("naive summary:", naive_summary)
    print("sym summary:", sym_summary)

    assert naive_ts.getLocationNum() == sym_ts.getLocationNum(), f"{name}: location count differs"
    assert naive_end_locs == sym_end_locs, f"{name}: end location ids differ"
    assert naive_summary == sym_summary, f"{name}: end-location semantics differ"



def main():
    if len(sys.argv) > 1:
        case = sys.argv[1]
        if case == "while":
            while_qc = build_rus_while_circuit()
            naive_while, naive_while_end = run_naive(while_qc)
            sym_while, sym_while_end = run_symbolic(while_qc, max_locations=32)
            assert_match("while_loop", naive_while, naive_while_end, sym_while, sym_while_end)
            print("while_loop SymTS check passed.")
            return
        if case == "if":
            if_qc = build_if_test_circuit()
            naive_if, naive_if_end = run_naive(if_qc)
            sym_if, sym_if_end = run_symbolic(if_qc, max_locations=16)
            assert_match("if_test", naive_if, naive_if_end, sym_if, sym_if_end)
            print("if_test SymTS check passed.")
            return
        raise ValueError(f"Unsupported test case: {case}")

    for case in ("while", "if"):
        subprocess.run([sys.executable, __file__, case], check=True)

    print("Control-flow SymTS checks passed.")


if __name__ == "__main__":
    main()
