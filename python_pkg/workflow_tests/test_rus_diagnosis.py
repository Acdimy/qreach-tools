from pathlib import Path
import sys
import subprocess

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import pyqreach
from qiskit import QuantumCircuit

from parse_qiskit import parse_qiskit_cir
from qctl import (
    amplitude_state,
    IRR_SQRT2,
    modelChecking,
    set_initial_state,
    _format_counterexample_analysis,
    _format_semantic_diagnosis,
)


def build_rus(*, corrected: bool) -> QuantumCircuit:
    qc = QuantumCircuit(3, 1)
    qc.x(2)
    qc.measure(2, 0)
    qc.reset(1)
    with qc.while_loop((0, 0b1)):
        if corrected:
            qc.reset(0)
        qc.h(0)
        qc.t(0)
        qc.cx(0, 1)
        qc.h(0)
        qc.cx(0, 1)
        qc.t(0)
        qc.h(0)
        qc.measure(0, 0)
    return qc


def run_rus_check(*, corrected: bool, with_target: bool):
    qc = build_rus(corrected=corrected)
    ts = pyqreach.TransitionSystem()
    parse_result = parse_qiskit_cir(qc, qc.num_qubits, ts, return_metadata=True)
    init_op = set_initial_state(ts, "000")
    ts.computingFixedPointPost()

    outloop_locs = [
        loc for loc in range(ts.getLocationNum())
        if "EW" in (ts.Locations[loc].getIdentifier() or "")
        and ts.printDims(loc)[1] > 0
    ]
    if not outloop_locs:
        raise RuntimeError("No outloop location with non-zero quantum state found")
    for loc in outloop_locs:
        ts.setLabel(loc, "outloop")

    diagnosis_config = None
    if with_target:
        target_s = amplitude_state({"001": 1.0 + 0j, "011": 1j * IRR_SQRT2})
        for loc in outloop_locs:
            ts.setAnnotation([[loc, target_s]])
        ts.computingFixedPointPre()
        ts.resetLocationBounds()
        set_initial_state(ts, "000")
        ts.computingFixedPointPost()
        from qctl import tsLabelling
        tsLabelling(ts, target_s, "s")
        diagnosis_config = {
            "initial_annotations": [(ts.getInitLocation(), init_op)],
            "target_annotations": [(loc, target_s) for loc in outloop_locs],
        }

    result = modelChecking(ts, "AG (outloop -> s)", semantic_diagnosis_config=diagnosis_config)
    return result


def _run_buggy_case() -> None:
    buggy = run_rus_check(corrected=False, with_target=True)
    assert buggy["satisfied"] is False
    assert buggy["diagnosis"] is not None
    assert buggy["diagnosis"]["status"] == "found"
    assert buggy["diagnosis"]["criterion"] == "incompatible_contribution_expansion"
    update = buggy["diagnosis"]["suspicious_update"]
    assert update["source_identifier"].startswith("S4.W")
    assert update["target_identifier"] == "S3" or update["target_identifier"].startswith("S4.W")
    assert update["sp_old_dim"] == 1
    assert update["sp_new_dim"] == 2
    assert update["contribution_compatible"] is False
    print(_format_counterexample_analysis(buggy["analysis"]))
    print(_format_semantic_diagnosis(buggy["diagnosis"]))


def _run_correct_case() -> None:
    correct = run_rus_check(corrected=True, with_target=True)
    assert correct["satisfied"] is True
    assert correct.get("diagnosis") is None


def _run_missing_target_case() -> None:
    missing_target = run_rus_check(corrected=False, with_target=False)
    assert missing_target["satisfied"] is False
    assert missing_target["diagnosis"] is None


def _run_unavailable_case() -> None:
    qc = build_rus(corrected=False)
    ts = pyqreach.TransitionSystem()
    parse_qiskit_cir(qc, qc.num_qubits, ts, return_metadata=True)
    set_initial_state(ts, "000")
    ts.computingFixedPointPost()
    for loc in range(ts.getLocationNum()):
        if "EW" in (ts.Locations[loc].getIdentifier() or "") and ts.printDims(loc)[1] > 0:
            ts.setLabel(loc, "outloop")
    direct_no_target = modelChecking(
        ts,
        "AG (outloop -> s)",
        semantic_diagnosis_config={"target_annotations": []},
    )
    assert direct_no_target["satisfied"] is False
    assert direct_no_target["diagnosis"]["status"] == "unavailable"


def _run_case_in_subprocess(case_name: str) -> None:
    cmd = [sys.executable, __file__, case_name]
    result = subprocess.run(cmd, check=False, text=True, capture_output=True)
    if result.stdout:
        print(result.stdout, end="")
    if result.stderr:
        print(result.stderr, end="", file=sys.stderr)
    assert result.returncode == 0, f"subprocess case {case_name!r} failed with {result.returncode}"


if __name__ == "__main__":
    case = sys.argv[1] if len(sys.argv) > 1 else "all"
    if case == "buggy":
        _run_buggy_case()
    elif case == "correct":
        _run_correct_case()
    elif case == "missing_target":
        _run_missing_target_case()
    elif case == "unavailable":
        _run_unavailable_case()
    else:
        for case_name in ["buggy", "correct", "missing_target", "unavailable"]:
            _run_case_in_subprocess(case_name)
        print("test_rus_diagnosis  OK")