from pathlib import Path
import sys

PYTHON_PKG = Path(__file__).resolve().parents[1]
if str(PYTHON_PKG) not in sys.path:
    sys.path.insert(0, str(PYTHON_PKG))

from workflow_tests.qasm_workflow_runner import (  # noqa: E402
    add_common_arguments,
    make_config_from_args,
    run_batch,
)


DEFAULT_INPUT_DIR = PYTHON_PKG / "benchmark"
DEFAULT_OUTPUT_DIR = PYTHON_PKG / "eval" / "scale_debug"
DEFAULT_OUTPUT_NAME = "qasm_benchmarks_lazy_results.csv"


def main() -> None:
    import argparse

    parser = argparse.ArgumentParser(
        description="Recursively run QisMC benchmark QASM checks. Use --no-lazy for ordinary parsing."
    )
    add_common_arguments(parser, default_error_injection=False)
    args = parser.parse_args()
    config = make_config_from_args(
        args,
        default_input_dir=DEFAULT_INPUT_DIR,
        default_output_dir=DEFAULT_OUTPUT_DIR,
        default_output_name=DEFAULT_OUTPUT_NAME,
    )
    run_batch(config)


if __name__ == "__main__":
    main()
