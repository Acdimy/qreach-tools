import json
import os
import subprocess
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parent
MEASURE_ONLY = ROOT / "debug_steane_measure_only.py"
STEANE_POST = ROOT / "test_symts_steane_post.py"


def run_json_command(argv, extra_env=None):
    env = os.environ.copy()
    if extra_env:
        env.update(extra_env)
    completed = subprocess.run(
        argv,
        check=True,
        capture_output=True,
        text=True,
        env=env,
    )
    for line in completed.stdout.splitlines():
        if not line.startswith("JSON_RESULT=") and not line.startswith("{"):
            print(line)
    for line in reversed(completed.stdout.splitlines()):
        if line.startswith("JSON_RESULT="):
            return json.loads(line[len("JSON_RESULT="):])
        if line.startswith("{"):
            return json.loads(line)
    raise RuntimeError(f"Missing JSON payload for command: {' '.join(argv)}")


def check_bounded_measure_only_regression():
    expected = {
        4: {"annotation_nodes": 98, "unique_nodes": 1862},
        8: {"annotation_nodes": 133, "unique_nodes": 1920},
    }
    payloads = []
    for limit in (4, 8):
        payload = run_json_command(
            [sys.executable, str(MEASURE_ONLY), "0"],
            extra_env={"TS_MAX_POST_ITER": str(limit)},
        )
        print(
            f"measure-only limit={limit}: post_time={payload['post_time']:.2f}s "
            f"annotation_nodes={payload['annotation_nodes']} unique_nodes={payload['unique_nodes']}"
        )
        assert payload["error"] is None, f"measure-only limit={limit} failed: {payload['error']}"
        assert payload["location_count"] == 37, f"Unexpected location count at limit={limit}"
        assert payload["end_count"] == 1, f"Unexpected end count at limit={limit}"
        assert payload["relation_nodes"] == 1206, f"Unexpected relation nodes at limit={limit}"
        assert payload["annotation_nodes"] == expected[limit]["annotation_nodes"], (
            f"Unexpected annotation nodes at limit={limit}: {payload['annotation_nodes']}"
        )
        assert payload["unique_nodes"] == expected[limit]["unique_nodes"], (
            f"Unexpected unique nodes at limit={limit}: {payload['unique_nodes']}"
        )
        payloads.append(payload)

    first, second = payloads
    assert first["annotation_nodes"] < second["annotation_nodes"], "Annotation nodes should grow with more iterations"
    assert first["unique_nodes"] < second["unique_nodes"], "Unique nodes should grow with more iterations"
    assert first["post_time"] < second["post_time"], "Bounded post time should grow with more iterations"


def check_if_construct_regression():
    naive = run_json_command([sys.executable, str(STEANE_POST), "naive", "1", "construct"])
    sym = run_json_command([sys.executable, str(STEANE_POST), "sym", "1", "construct"])
    print(
        "if-construct measure=1: "
        f"naive_locations={naive['location_count']} sym_locations={sym['location_count']}"
    )
    assert naive["location_count"] == sym["location_count"], "Construct-stage location count differs"
    assert naive["end_count"] == sym["end_count"], "Construct-stage end count differs"
    assert naive["summary"] == sym["summary"], "Construct-stage summaries differ"


def main():
    check_bounded_measure_only_regression()
    check_if_construct_regression()
    print("Steane bounded regression checks passed.")


if __name__ == "__main__":
    main()