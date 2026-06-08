import json
import sys
from pathlib import Path

import networkx as nx
import pyqreach

from qctl import nx2Graph_hierarchical
from test_symts_steane_post import build_steane_if_test_circuit, run_symbolic


def qadd_data_to_nx(data):
    graph = nx.DiGraph()
    for node in data["nodes"]:
        node_id = node["id"]
        kind = node["kind"]
        label = node["label"]
        highlight = kind == "terminal"
        graph.add_node(node_id, label=label, highlight=highlight, kind=kind)
    for edge in data["edges"]:
        graph.add_edge(edge["source"], edge["target"], label=str(edge["branch"]))
    return graph


def summarize_graph(data, name):
    internal = sum(1 for node in data["nodes"] if node["kind"] == "internal")
    terminal = sum(1 for node in data["nodes"] if node["kind"] == "terminal")
    max_var = max((node.get("var", -1) for node in data["nodes"]), default=-1)
    return {
        "name": name,
        "node_count": data["node_count"],
        "edge_count": data["edge_count"],
        "internal_nodes": internal,
        "terminal_nodes": terminal,
        "max_var": max_var,
    }


def main():
    measure_count = int(sys.argv[1]) if len(sys.argv) > 1 else 1
    stage = sys.argv[2] if len(sys.argv) > 2 else "full"
    output_dir = Path(sys.argv[3]) if len(sys.argv) > 3 else Path("output/steane_dd")
    output_dir.mkdir(parents=True, exist_ok=True)

    qc = build_steane_if_test_circuit(measure_count)
    ts, end_locs, cons_time, post_time = run_symbolic(qc, max_locations=0, run_post=(stage == "full"))

    relation_data = ts.exportRelationQADD()
    annotation_data = ts.exportAnnotationQADD()

    relation_graph = qadd_data_to_nx(relation_data)
    annotation_graph = qadd_data_to_nx(annotation_data)

    relation_path = output_dir / f"steane_relation_dd_m{measure_count}_{stage}"
    annotation_path = output_dir / f"steane_qoperation_dd_m{measure_count}_{stage}"
    nx2Graph_hierarchical(relation_graph, str(relation_path))
    nx2Graph_hierarchical(annotation_graph, str(annotation_path))

    summary = {
        "measure_count": measure_count,
        "stage": stage,
        "construction_time": cons_time,
        "post_time": post_time,
        "location_count": ts.getLocationNum(),
        "end_count": len(end_locs),
        "relation": summarize_graph(relation_data, "relation_dd"),
        "qoperation": summarize_graph(annotation_data, "qoperation_dd"),
        "relation_pdf": str(relation_path.with_suffix(".pdf")),
        "qoperation_pdf": str(annotation_path.with_suffix(".pdf")),
    }
    print(json.dumps(summary, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()