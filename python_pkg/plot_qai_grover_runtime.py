from __future__ import annotations

import re
from pathlib import Path

import matplotlib
import matplotlib.pyplot as plt
import pandas as pd

CSV_PATH = Path("eval/qai_grover_case_qismc.csv")
OUTPUT_DIR = Path("output/qai_grover_runtime")
OUTPUT_BASENAME = "qai_grover_runtime_qismc"
TIMEOUT_S = 600.0

PATTERN = re.compile(
    r"single-it-grover(?P<size>\d+)-(?P<oracle>plus|zero)(?P<linear>-linear)?\.qasm$"
)

VARIANT_ORDER = ["plus", "zero", "plus-linear", "zero-linear"]
VARIANT_LABELS = {
    "plus": "Grover plus",
    "zero": "Grover zero",
    "plus-linear": "Grover plus linear",
    "zero-linear": "Grover zero linear",
}
VARIANT_COLORS = {
    "plus": "#2A9D8F",
    "zero": "#E76F51",
    "plus-linear": "#264653",
    "zero-linear": "#E9C46A",
}
VARIANT_MARKERS = {
    "plus": "o",
    "zero": "s",
    "plus-linear": "^",
    "zero-linear": "D",
}

matplotlib.rcParams.update(
    {
        "font.family": "serif",
        "font.size": 10,
        "axes.labelsize": 11,
        "legend.fontsize": 8,
        "figure.dpi": 150,
        "savefig.bbox": "tight",
        "savefig.pad_inches": 0.05,
    }
)


def parse_variant(name: str) -> pd.Series:
    match = PATTERN.match(name)
    if match is None:
        raise ValueError(f"Unexpected Grover filename: {name}")

    linear_suffix = "-linear" if match.group("linear") else ""
    variant = f"{match.group('oracle')}{linear_suffix}"
    problem_size = int(match.group("size"))
    return pd.Series({"variant": variant, "problem_size": problem_size})


def load_plot_data(csv_path: Path) -> pd.DataFrame:
    df = pd.read_csv(csv_path)
    df = df.drop_duplicates(subset=["relative_filename"], keep="last").copy()

    parsed = df["relative_filename"].apply(parse_variant)
    df = pd.concat([df, parsed], axis=1)

    df["time_total"] = pd.to_numeric(df["time_total"], errors="coerce")
    df["timeout_seconds"] = pd.to_numeric(df["timeout_seconds"], errors="coerce")
    df["plot_time_seconds"] = df["time_total"]

    timeout_mask = df["status"].eq("timeout") | (df["time_total"] >= TIMEOUT_S)
    df.loc[timeout_mask, "status"] = "timeout"
    df.loc[timeout_mask, "plot_time_seconds"] = TIMEOUT_S

    return df.sort_values(["variant", "problem_size"]).reset_index(drop=True)


def build_figure(df: pd.DataFrame) -> plt.Figure:
    plt.style.use("default")

    fig, ax = plt.subplots(figsize=(8.6, 4.8))

    for variant in VARIANT_ORDER:
        subset = df[df["variant"] == variant].sort_values("problem_size")
        if subset.empty:
            continue

        xs = subset["problem_size"].tolist()
        ys = subset["plot_time_seconds"].tolist()
        if xs:
            ax.plot(
                xs,
                ys,
                label=VARIANT_LABELS[variant],
                color=VARIANT_COLORS[variant],
                marker=VARIANT_MARKERS[variant],
                linewidth=1.2,
                markersize=4,
                zorder=3,
            )

        timeout_subset = subset[subset["status"] == "timeout"]
        if not timeout_subset.empty:
            tx = timeout_subset["problem_size"].tolist()
            ty = [TIMEOUT_S] * len(tx)
            ax.scatter(
                tx,
                ty,
                marker="x",
                s=36,
                linewidths=1.2,
                color=VARIANT_COLORS[variant],
                zorder=4,
            )

    ax.axhline(
        y=TIMEOUT_S,
        color="grey",
        linestyle="--",
        linewidth=0.7,
        alpha=0.6,
    )
    ax.text(
        0.98,
        TIMEOUT_S * 1.02,
        f"timeout ({int(TIMEOUT_S)}s)",
        color="grey",
        fontsize=7,
        ha="right",
        transform=ax.get_yaxis_transform(),
    )

    ax.set_xlabel("Grover problem size")
    ax.set_ylabel("Runtime (s)")
    ax.set_title("QisMC runtime on QAI Grover benchmarks", fontsize=11, pad=6)
    ax.set_yscale("log")
    ax.set_xticks(sorted(df["problem_size"].unique()))
    ax.grid(True, which="both", alpha=0.25, linewidth=0.5)
    ax.legend(loc="lower right", framealpha=0.9, edgecolor="grey", fontsize=7.5, ncol=2)

    fig.tight_layout()
    return fig


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    df = load_plot_data(CSV_PATH)
    fig = build_figure(df)

    pdf_path = OUTPUT_DIR / f"{OUTPUT_BASENAME}.pdf"
    png_path = OUTPUT_DIR / f"{OUTPUT_BASENAME}.png"
    fig.savefig(pdf_path, format="pdf", bbox_inches="tight")
    fig.savefig(png_path, dpi=300, bbox_inches="tight")

    print(df[["relative_filename", "variant", "problem_size", "status", "plot_time_seconds"]].to_string(index=False))
    print(f"Saved: {pdf_path}")
    print(f"Saved: {png_path}")


if __name__ == "__main__":
    main()
