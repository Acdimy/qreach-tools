import pyqreach
import json
import subprocess
import sys
from time import time

from parse_qiskit import build_sym_ts_from_qiskit, parse_qiskit_cir
from circ_utils import prepare_steane_code
from qiskit import QuantumCircuit


def build_steane_if_test_circuit() -> QuantumCircuit:
    circ = QuantumCircuit(14, 14)
    prepare_steane_code(circ, list(range(7)))
    circ.h(13)
    prepare_steane_code(circ, list(range(7, 14)))
    for i in range(7):
        circ.cx(i, 7 + i)
    for i in range(7, 14):
        circ.measure(i, i)
    for i in range(7, 14):
        with circ.if_test((i, 1)):
            circ.x(i)
    return circ


def run_naive(qc: QuantumCircuit):
    pyqreach.initializeTransitionSystem()
    ts = pyqreach.TransitionSystem(False)
    cons_start = time()
    end_locs = parse_qiskit_cir(qc, qc.num_qubits, ts)
    cons_time = time() - cons_start

    ts.setAnnotation([[0, pyqreach.QOperation(["0" * qc.num_qubits])]])
    post_start = time()
    ts.computingFixedPointPost()
    post_time = time() - post_start
    return ts, end_locs, cons_time, post_time


def run_symbolic(qc: QuantumCircuit, max_locations: int = 0):
    cons_start = time()
    ts, end_locs = build_sym_ts_from_qiskit(qc, qnum=qc.num_qubits, max_locations=max_locations)
    cons_time = time() - cons_start

    ts.setAnnotation(0, pyqreach.QOperation(["0" * qc.num_qubits]))
    post_start = time()
    ts.computingFixedPointPost()
    post_time = time() - post_start
    return ts, end_locs, cons_time, post_time


def collect_end_summary_naive(ts, end_locs):
    records = []
    for loc in end_locs:
        records.append(
            {
                "cp": ts.Locations[loc].cp.toString(),
                "identifier": ts.Locations[loc].getIdentifier(),
                "lower_dim": ts.printDims(loc)[1],
            }
        )
    return sorted(records, key=lambda item: (item["cp"], item["identifier"]))


def collect_end_summary_symbolic(ts, end_locs):
    records = []
    for loc in end_locs:
        records.append(
            {
                "cp": ts.getClassicalProposition(loc).toString(),
                "identifier": ts.getIdentifier(loc),
                "lower_dim": ts.printDims(loc)[1],
            }
        )
    return sorted(records, key=lambda item: (item["cp"], item["identifier"]))


def main():
    if len(sys.argv) > 1:
        mode = sys.argv[1]
        qc = build_steane_if_test_circuit()

        if mode == "naive":
            ts, end_locs, cons_time, post_time = run_naive(qc)
            payload = {
                "kind": "naive",
                "construction_time": cons_time,
                "post_time": post_time,
                "location_count": ts.getLocationNum(),
                "end_count": len(end_locs),
                "summary": collect_end_summary_naive(ts, end_locs),
            }
            print("JSON_RESULT=" + json.dumps(payload, sort_keys=True))
            return

        if mode == "sym":
            ts, end_locs, cons_time, post_time = run_symbolic(qc, max_locations=0)
            payload = {
                "kind": "sym",
                "construction_time": cons_time,
                "post_time": post_time,
                "location_count": ts.getLocationNum(),
                "end_count": len(end_locs),
                "summary": collect_end_summary_symbolic(ts, end_locs),
            }
            print("JSON_RESULT=" + json.dumps(payload, sort_keys=True))
            return

        raise ValueError(f"Unsupported mode: {mode}")

    def run_case(mode: str):
        completed = subprocess.run(
            [sys.executable, __file__, mode],
            check=True,
            capture_output=True,
            text=True,
        )
        for line in completed.stdout.splitlines():
            if not line.startswith("JSON_RESULT="):
                print(line)
        for line in reversed(completed.stdout.splitlines()):
            if line.startswith("JSON_RESULT="):
                return json.loads(line[len("JSON_RESULT="):])
        raise RuntimeError(f"Missing JSON result for mode={mode}")

    naive_payload = run_case("naive")
    sym_payload = run_case("sym")

    print("=== Steane PostImage Comparison ===")
    print(f"naive construction time: {naive_payload['construction_time']:.2f}s")
    print(f"sym construction time: {sym_payload['construction_time']:.2f}s")
    print(f"naive post time: {naive_payload['post_time']:.2f}s")
    print(f"sym post time: {sym_payload['post_time']:.2f}s")
    print("naive location count:", naive_payload["location_count"])
    print("sym location count:", sym_payload["location_count"])
    print("naive end count:", naive_payload["end_count"])
    print("sym end count:", sym_payload["end_count"])

    assert naive_payload["location_count"] == sym_payload["location_count"], "Location count differs between naive TS and SymTS"
    assert naive_payload["end_count"] == sym_payload["end_count"], "End location count differs between naive TS and SymTS"
    assert naive_payload["summary"] == sym_payload["summary"], "End-location postImage results differ between naive TS and SymTS"

    print("Steane postImage SymTS check passed.")


if __name__ == "__main__":
    main()