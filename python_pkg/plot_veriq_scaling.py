#!/usr/bin/env python3
"""Generate scaling plots for VeriQ-bench experiments.

Reads all CSVs under ``eval/scale/`` and ``eval/scale_debug/`` and produces:

* ``output/veriq_scaling.pdf``  — total time vs n, all families, clean + injected
* ``output/veriq_scaling.png``  — same, PNG raster
* ``output/veriq_breakdown.pdf`` — stacked bar / subplot of parse vs fixed-point time

Dependencies: matplotlib (already in .venv).
"""
from __future__ import annotations

import csv
import re
from pathlib import Path

import matplotlib
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from matplotlib.lines import Line2D

# ---------------------------------------------------------------------------
# config
# ---------------------------------------------------------------------------
PYTHON_PKG = Path(__file__).resolve().parent
SCALE_DIR = PYTHON_PKG / "eval" / "scale"
SCALE_DEBUG_DIR = PYTHON_PKG / "eval" / "scale_debug"
OUTPUT_DIR = PYTHON_PKG / "output"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

TIMEOUT_S = 120.0  # mark rows at/above this as timed out

FAMILY_LABEL: dict[str, str] = {
    "dqc_pe": "DQC-PE",
    "dqc_qft": "DQC-QFT",
    "qft": "QFT",
    "pe": "PE",
    "grover": "Grover",
}

FAMILY_COLOR: dict[str, str] = {
    "dqc_pe": "#e41a1c",
    "dqc_qft": "#377eb8",
    "qft": "#4daf4a",
    "pe": "#984ea3",
    "grover": "#ff7f00",
}

FAMILY_MARKER: dict[str, str] = {
    "dqc_pe": "o",
    "dqc_qft": "s",
    "qft": "D",
    "pe": "^",
    "grover": "v",
}

matplotlib.rcParams.update({
    "font.family": "serif",
    "font.size": 10,
    "axes.labelsize": 11,
    "legend.fontsize": 8,
    "figure.dpi": 150,
    "savefig.bbox": "tight",
    "savefig.pad_inches": 0.05,
})


# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------
def _extract_n(filename: str) -> int:
    m = re.search(r"_(\d+)\.qasm$", filename)
    return int(m.group(1)) if m else 0


def _load_data(csv_dir: Path, mode_label: str) -> dict[str, list[dict]]:
    """Return {family: [sorted rows]} for all CSVs in *csv_dir*."""
    data: dict[str, list[dict]] = {}
    for csv_path in sorted(csv_dir.glob("*_lazy*results.csv")):
        fam = csv_path.stem.replace("_lazy_injected_results", "").replace("_lazy_results", "")
        rows: list[dict] = []
        with csv_path.open() as f:
            for r in csv.DictReader(f):
                rows.append(r)
        # sort by n
        rows.sort(key=lambda r: _extract_n(r.get("relative_filename", "")))
        data[fam] = rows
    return data


# ---------------------------------------------------------------------------
# Figure 1: scaling  (total time vs n)
# ---------------------------------------------------------------------------
def plot_scaling(
    clean_data: dict[str, list[dict]],
    injected_data: dict[str, list[dict]],
) -> Path:
    fig, (ax_clean, ax_inj) = plt.subplots(
        1, 2, figsize=(11, 5), sharey=True,
        gridspec_kw={"wspace": 0.06},
    )

    for ax, data, title in [
        (ax_clean, clean_data, "Clean (no error injection)"),
        (ax_inj, injected_data, "Injected (random Pauli / wrong init)"),
    ]:
        for fam, rows in data.items():
            if fam not in FAMILY_LABEL:
                continue
            color = FAMILY_COLOR[fam]
            marker = FAMILY_MARKER[fam]
            label = FAMILY_LABEL[fam]

            xs, ys = [], []
            for r in rows:
                if r.get("status") != "ok":
                    continue
                n = _extract_n(r.get("relative_filename", ""))
                t = float(r.get("time_total", 0))
                if t <= 0:
                    continue
                xs.append(n)
                ys.append(t)

            # timeout markers
            tx, ty = [], []
            for r in rows:
                if r.get("status") == "timeout":
                    n = _extract_n(r.get("relative_filename", ""))
                    tx.append(n)
                    ty.append(TIMEOUT_S)

            if xs:
                ax.plot(xs, ys, color=color, marker=marker, markersize=4,
                        linewidth=1.2, label=label, zorder=3)
            if tx:
                ax.scatter(tx, ty, marker="x", color=color, s=36,
                           linewidths=1.2, zorder=4)

        # timeout reference line
        ax.axhline(y=TIMEOUT_S, color="grey", linestyle="--", linewidth=0.7, alpha=0.6)
        ax.text(0.5, TIMEOUT_S * 1.02, "timeout (120s)", color="grey",
                fontsize=7, ha="right", transform=ax.get_yaxis_transform())

        ax.set_yscale("log")
        ax.set_title(title, fontsize=11, pad=6)
        ax.set_xlabel("n  (scale parameter)")
        ax.grid(True, which="both", alpha=0.25, linewidth=0.5)

    ax_clean.set_ylabel("Total time (s)")
    ax_clean.legend(loc="upper left", framealpha=0.9, edgecolor="grey",
                    fontsize=7.5)

    fig.suptitle("VeriQ-bench — Lazy QReach Scaling", fontsize=13, y=0.99)
    path = OUTPUT_DIR / "veriq_scaling.pdf"
    fig.savefig(path)
    fig.savefig(OUTPUT_DIR / "veriq_scaling.png")
    plt.close(fig)
    return path


# ---------------------------------------------------------------------------
# Figure 2: time breakdown  (parse vs fixed-point post)
# ---------------------------------------------------------------------------
def plot_breakdown(
    clean_data: dict[str, list[dict]],
    injected_data: dict[str, list[dict]],
) -> Path:
    families = [f for f in FAMILY_LABEL if f in clean_data]

    fig, axes = plt.subplots(
        2, len(families), figsize=(len(families) * 2.4, 8),
        sharex="col",
        gridspec_kw={"hspace": 0.35, "wspace": 0.25},
    )

    for col, fam in enumerate(families):
        for row_idx, (mode_label, data) in enumerate(
            [("clean", clean_data), ("injected", injected_data)]
        ):
            ax = axes[row_idx][col] if len(families) > 1 else axes[row_idx]
            rows = data.get(fam, [])
            if not rows:
                continue

            ns = []
            parse_times = []
            fixed_times = []
            for r in rows:
                if r.get("status") != "ok":
                    continue
                n = _extract_n(r.get("relative_filename", ""))
                prepare = float(r.get("time_prepare_ts", 0))
                fixed = float(r.get("time_fixed_point_post", 0))
                if prepare < 0 or fixed < 0:
                    continue
                ns.append(str(n))
                parse_times.append(prepare)
                fixed_times.append(fixed)

            if not ns:
                continue

            x = range(len(ns))
            w = 0.35
            color = FAMILY_COLOR[fam]
            ax.bar([i - w / 2 for i in x], parse_times, w,
                   color=color, alpha=0.75, label="parse / build", linewidth=0)
            ax.bar([i + w / 2 for i in x], fixed_times, w,
                   color=color, alpha=0.35, label="fixed-point post", linewidth=0)

            ax.set_xticks(list(x))
            ax.set_xticklabels(ns, rotation=45, ha="right", fontsize=7)
            ax.set_yscale("log")
            ax.grid(axis="y", alpha=0.25, linewidth=0.5)

            if col == 0:
                ax.set_ylabel(f"{mode_label}\nTime (s)", fontsize=9)
            if row_idx == 0:
                ax.set_title(FAMILY_LABEL[fam], fontsize=10, pad=4)

    # single shared legend
    handles = [
        Line2D([0], [0], color="grey", alpha=0.75, lw=6, label="parse / build"),
        Line2D([0], [0], color="grey", alpha=0.35, lw=6, label="fixed-point post"),
    ]
    fig.legend(handles=handles, loc="lower center", ncol=2, fontsize=8,
               framealpha=0.9, edgecolor="grey")

    fig.suptitle("VeriQ-bench — Time Breakdown", fontsize=13, y=1.00)
    path = OUTPUT_DIR / "veriq_breakdown.pdf"
    fig.savefig(path)
    fig.savefig(OUTPUT_DIR / "veriq_breakdown.png")
    plt.close(fig)
    return path


# ---------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------
def main() -> None:
    clean_data = _load_data(SCALE_DIR, "clean")
    injected_data = _load_data(SCALE_DEBUG_DIR, "injected")

    # report what was found
    for label, data in [("clean", clean_data), ("injected", injected_data)]:
        for fam, rows in sorted(data.items()):
            ok = sum(1 for r in rows if r.get("status") == "ok")
            to = sum(1 for r in rows if r.get("status") == "timeout")
            print(f"  {label:8s} {fam:8s}  {ok:3d} ok  {to:3d} timeout")

    p1 = plot_scaling(clean_data, injected_data)
    print(f"\nscaling plot  → {p1}")
    p2 = plot_breakdown(clean_data, injected_data)
    print(f"breakdown plot → {p2}")


if __name__ == "__main__":
    main()
