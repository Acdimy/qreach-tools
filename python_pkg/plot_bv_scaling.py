from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


PYTHON_PKG = Path(__file__).resolve().parent
DEFAULT_INPUT = PYTHON_PKG / "eval" / "scale_debug" / "bv_lazy_injected_results.csv"
DEFAULT_OUTPUT = PYTHON_PKG / "output" / "scale_bv_lazy_injected.pdf"
TITLE = "BV_DEBUG"


def generate_bv_scale_figure(
    filename: Path,
    *,
    output_file: Path,
    logscale_left: bool = True,
    logscale_right: bool = True,
    append_veri_time: bool = True,
) -> None:
    df = pd.read_csv(filename)
    df = df[df["status"] == "ok"].copy()
    if df.empty:
        raise ValueError(f"No successful rows found in {filename}")

    # Batch results are written in filename order, so sort numerically for scaling plots.
    df = df.sort_values("num_qubits").reset_index(drop=True)
    df = df[(df["num_qubits"] >= 5) & (df["num_qubits"] % 5 == 0)].reset_index(drop=True)
    if df.empty:
        raise ValueError("No BV rows found for qubit counts 5, 10, 15, ..., 100")

    x_vals = df["num_qubits"]
    x_idx = np.arange(len(x_vals))
    time_prepare = df["time_prepare_ts"].to_numpy(copy=True)
    time_post = df["time_fixed_point_post"].to_numpy(copy=True)
    if append_veri_time:
        time_post = time_post + df["verification_time"].to_numpy(copy=True)
    num_gates = df["num_gates"]
    num_locations = df["num_locations"]

    plt.style.use("default")
    plt.rcParams.update(
        {
            "font.size": 11,
            "axes.labelsize": 12,
            "axes.titlesize": 13,
            "xtick.direction": "in",
            "ytick.direction": "in",
            "axes.grid": True,
            "grid.linestyle": "--",
            "grid.alpha": 0.5,
        }
    )

    fig, ax1 = plt.subplots(figsize=(7, 4))

    width = 0.35
    ax1.bar(x_idx, time_prepare, width, label="Time (prepare_ts)", color="#7db5d8", alpha=0.8)
    ax1.bar(x_idx, time_post, width, bottom=time_prepare, label="Time (checking)", color="#f2a444", alpha=0.8)

    ax1.set_xlabel("Number of qubits")
    ax1.set_ylabel(f"Execution time (s{', log scale' if logscale_left else ''})")
    ax1.tick_params(axis="y")
    ax1.set_title(f"{TITLE} performance metrics")
    if logscale_left:
        ax1.set_yscale("log")
    ax1.grid(axis="y", linestyle="--", alpha=0.5)
    ax1.grid(axis="x", linestyle="--", alpha=0.5)

    ax2 = ax1.twinx()
    ax2.plot(x_idx, num_gates, color="#2b8cbe", marker="s", markersize=5, label="Num. of gates")
    ax2.plot(x_idx, num_locations, color="#de2d26", marker="^", markersize=5, label="Num. of locations")

    ax2.set_ylabel(f"Counts{' (log scale)' if logscale_right else ''}")
    ax2.tick_params(axis="y")
    if logscale_right:
        ax2.set_yscale("log")
    ax2.grid(False)

    lines_1, labels_1 = ax1.get_legend_handles_labels()
    lines_2, labels_2 = ax2.get_legend_handles_labels()
    ax1.legend(lines_1 + lines_2, labels_1 + labels_2, loc="upper left", frameon=False)

    ax1.set_xticks(x_idx)
    ax1.set_xticklabels([str(q) for q in x_vals])
    ax1.margins(x=0.05)
    fig.tight_layout()

    output_file.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_file, format=output_file.suffix.lstrip(".") or "pdf", bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Generate a BV scaling figure from bv_lazy_injected_results.csv."
    )
    parser.add_argument("--input-file", type=Path, default=DEFAULT_INPUT, help="Input CSV file")
    parser.add_argument("--output-file", type=Path, default=DEFAULT_OUTPUT, help="Output figure path")
    parser.add_argument("--no-logscale-left", dest="logscale_left", action="store_false", help="Disable log scale on the execution-time axis")
    parser.add_argument("--no-logscale-right", dest="logscale_right", action="store_false", help="Disable log scale on the counts axis")
    parser.add_argument("--no-verification-time", dest="append_veri_time", action="store_false", help="Exclude verification_time from the checking bar")
    parser.set_defaults(logscale_left=True, logscale_right=True, append_veri_time=True)
    args = parser.parse_args()

    generate_bv_scale_figure(
        args.input_file.expanduser().resolve(),
        output_file=args.output_file.expanduser().resolve(),
        logscale_left=args.logscale_left,
        logscale_right=args.logscale_right,
        append_veri_time=args.append_veri_time,
    )


if __name__ == "__main__":
    main()