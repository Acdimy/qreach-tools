from qreach import symbolic_available

if not symbolic_available():
    print("SKIP: symbolic SymTS backend unavailable (LimTDD build)")
    raise SystemExit(0)

import pyqreach
import json
import subprocess
import sys
from time import time

from qreach.parse_qiskit import build_sym_ts_from_qiskit, parse_qiskit_cir
from qreach.circ_utils import prepare_steane_code
from qiskit import QuantumCircuit


def build_steane_if_test_circuit(measure_count: int = 7) -> QuantumCircuit:
    assert 0 <= measure_count <= 7
    circ = QuantumCircuit(14, 14)
    prepare_steane_code(circ, list(range(7)))
    circ.h(13)
    prepare_steane_code(circ, list(range(7, 14)))
    for i in range(7):
        circ.cx(i, 7 + i)
    for i in range(7, 7 + measure_count):
        circ.measure(i, i)
    for i in range(7, 7 + measure_count):
        with circ.if_test((i, 1)):
            circ.x(i)
    return circ


def run_naive(qc: QuantumCircuit, run_post: bool = True):
    pyqreach.initializeTransitionSystem()
    ts = pyqreach.TransitionSystem(False)
    cons_start = time()
    end_locs = parse_qiskit_cir(qc, qc.num_qubits, ts)
    cons_time = time() - cons_start

    post_time = None
    if run_post:
        ts.setAnnotation([[0, pyqreach.QOperation(["0" * qc.num_qubits])]])
        post_start = time()
        ts.computingFixedPointPost()
        post_time = time() - post_start
    return ts, end_locs, cons_time, post_time


def run_symbolic(qc: QuantumCircuit, max_locations: int = 0, run_post: bool = True):
    cons_start = time()
    ts, end_locs = build_sym_ts_from_qiskit(qc, qnum=qc.num_qubits, max_locations=max_locations)
    cons_time = time() - cons_start

    post_time = None
    if run_post:
        ts.setAnnotation(0, pyqreach.QOperation(["0" * qc.num_qubits]))
        post_start = time()
        ts.computingFixedPointPost()
        post_time = time() - post_start
    return ts, end_locs, cons_time, post_time


def collect_ts_stats(ts, kind: str):
    if kind == "naive":
        return {
            "location_count": ts.getLocationNum(),
            "relation_count": len(ts.relations),
        }
    return {
        "location_count": ts.getLocationNum(),
        "relation_node_count": ts.getRelationNodeCount(),
        "annotation_node_count": ts.getAnnotationNodeCount(),
        "unique_node_count": ts.getTotalUniqueNodeCount(),
    }


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
        measure_count = int(sys.argv[2]) if len(sys.argv) > 2 else 7
        stage = sys.argv[3] if len(sys.argv) > 3 else "full"
        qc = build_steane_if_test_circuit(measure_count)
        run_post = stage == "full"

        if mode == "naive":
            ts, end_locs, cons_time, post_time = run_naive(qc, run_post=run_post)
            payload = {
                "kind": "naive",
                "measure_count": measure_count,
                "stage": stage,
                "construction_time": cons_time,
                "post_time": post_time,
                **collect_ts_stats(ts, "naive"),
                "end_count": len(end_locs),
                "summary": collect_end_summary_naive(ts, end_locs),
            }
            print("JSON_RESULT=" + json.dumps(payload, sort_keys=True))
            return

        if mode == "sym":
            ts, end_locs, cons_time, post_time = run_symbolic(qc, max_locations=0, run_post=run_post)
            payload = {
                "kind": "sym",
                "measure_count": measure_count,
                "stage": stage,
                "construction_time": cons_time,
                "post_time": post_time,
                **collect_ts_stats(ts, "sym"),
                "end_count": len(end_locs),
                "summary": collect_end_summary_symbolic(ts, end_locs),
            }
            print("JSON_RESULT=" + json.dumps(payload, sort_keys=True))
            return

        raise ValueError(f"Unsupported mode: {mode}")

    measure_count = int(sys.argv[1]) if len(sys.argv) > 1 else 7
    stage = sys.argv[2] if len(sys.argv) > 2 else "full"

    def run_case(mode: str):
        completed = subprocess.run(
            [sys.executable, __file__, mode, str(measure_count), stage],
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
    print(f"measure_count: {measure_count}")
    print(f"stage: {stage}")
    print(f"naive construction time: {naive_payload['construction_time']:.2f}s")
    print(f"sym construction time: {sym_payload['construction_time']:.2f}s")
    if naive_payload["post_time"] is not None:
        print(f"naive post time: {naive_payload['post_time']:.2f}s")
    if sym_payload["post_time"] is not None:
        print(f"sym post time: {sym_payload['post_time']:.2f}s")
    print("naive location count:", naive_payload["location_count"])
    print("sym location count:", sym_payload["location_count"])
    print("naive end count:", naive_payload["end_count"])
    print("sym end count:", sym_payload["end_count"])
    if "relation_count" in naive_payload:
        print("naive relation count:", naive_payload["relation_count"])
    if "relation_node_count" in sym_payload:
        print("sym relation nodes:", sym_payload["relation_node_count"])
        print("sym annotation nodes:", sym_payload["annotation_node_count"])
        print("sym unique nodes:", sym_payload["unique_node_count"])

    assert naive_payload["location_count"] == sym_payload["location_count"], "Location count differs between naive TS and SymTS"
    assert naive_payload["end_count"] == sym_payload["end_count"], "End location count differs between naive TS and SymTS"
    assert naive_payload["summary"] == sym_payload["summary"], "End-location postImage results differ between naive TS and SymTS"

    print("Steane postImage SymTS check passed.")


if __name__ == "__main__":
    main()