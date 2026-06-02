import json
import sys
from time import time

import pyqreach
from qiskit import QuantumCircuit

from circ_utils import prepare_steane_code
from parse_qiskit import build_sym_ts_from_qiskit


def build_steane_measure_only(measure_count: int) -> QuantumCircuit:
    circ = QuantumCircuit(14, 14)
    prepare_steane_code(circ, list(range(7)))
    circ.h(13)
    prepare_steane_code(circ, list(range(7, 14)))
    for index in range(7):
        circ.cx(index, 7 + index)
    for index in range(7, 7 + measure_count):
        circ.measure(index, index)
    return circ


def main():
    measure_counts = [int(sys.argv[1])] if len(sys.argv) > 1 else [0, 1, 3, 5, 7]
    for measure_count in measure_counts:
        qc = build_steane_measure_only(measure_count)
        cons_start = time()
        ts, end_locs = build_sym_ts_from_qiskit(qc, qnum=qc.num_qubits, max_locations=0)
        cons_time = time() - cons_start

        post_time = None
        error = None
        try:
            ts.setAnnotation(0, pyqreach.QOperation(["0" * qc.num_qubits]))
            post_start = time()
            ts.computingFixedPointPost()
            post_time = time() - post_start
        except Exception as exc:
            error = f"{type(exc).__name__}: {exc}"

        print(
            json.dumps(
                {
                    "measure_count": measure_count,
                    "construction_time": cons_time,
                    "post_time": post_time,
                    "location_count": ts.getLocationNum(),
                    "end_count": len(end_locs),
                    "relation_nodes": ts.getRelationNodeCount(),
                    "annotation_nodes": ts.getAnnotationNodeCount(),
                    "unique_nodes": ts.getTotalUniqueNodeCount(),
                    "error": error,
                },
                sort_keys=True,
            )
        )


if __name__ == "__main__":
    main()