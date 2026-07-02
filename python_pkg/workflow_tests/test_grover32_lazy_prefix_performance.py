from __future__ import annotations

from pathlib import Path

from debug_grover32_ccx_pathology import PYTHON_PKG, run_prefix_once


def _qasm(name: str) -> Path:
    return PYTHON_PKG / "benchmark" / "converted_qasm" / name


def test_grover32_plus_prefix_102_lazy_parse_under_budget() -> None:
    result = run_prefix_once(
        _qasm("single-it-grover32-plus.qasm"),
        max_instructions=102,
        initial_state=None,
        run_fixed_post=False,
    )

    assert result.status == "ok"
    assert result.num_locations == 103
    assert result.num_result_locations == 1
    assert result.time_parse != ""
    assert float(result.time_parse) < 5.0


def test_lazy_parse_controls_stay_fast_near_same_structure() -> None:
    linear = run_prefix_once(
        _qasm("single-it-grover32-plus-linear.qasm"),
        max_instructions=102,
        initial_state=None,
        run_fixed_post=False,
    )
    grover64 = run_prefix_once(
        _qasm("single-it-grover64-plus.qasm"),
        max_instructions=133,
        initial_state=None,
        run_fixed_post=False,
    )

    assert linear.status == "ok"
    assert linear.num_locations == 103
    assert linear.num_result_locations == 1
    assert linear.time_parse != ""
    assert float(linear.time_parse) < 1.0

    assert grover64.status == "ok"
    assert grover64.num_locations == 134
    assert grover64.num_result_locations == 1
    assert grover64.time_parse != ""
    assert float(grover64.time_parse) < 1.0
