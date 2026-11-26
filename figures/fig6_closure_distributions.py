"""Attachment Figure 3: annual closure-hour distributions with caps (experiment E4)."""

from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np

from common import CADENCE_COLORS, OI, read_csv, read_parquet, save

DIST_LABELS = {
    "record_flights_1_to_4": "record, Flights 1--4 (12.9 h)",
    "expected_upper_20h": "hypothetical upper case (20 h)",
    "recent_flights_2_to_4": "record, Flights 2--4 (9 h)",
}
DIST_COLORS = {
    "record_flights_1_to_4": OI["blue"],
    "expected_upper_20h": OI["vermillion"],
    "recent_flights_2_to_4": OI["green"],
}
CAPS = ((500, ":"), (800, "--"), (2000, "-"))
RECORD = (
    (24.5, "Flight 1 (24.5 h)", OI["vermillion"]),
    (10.5, "Flight 2 (10.5 h)", OI["orange"]),
    (8.75, "Flight 3 (8.75 h)", OI["purple"]),
    (7.75, "Flight 4 (7.75 h)", OI["green"]),
)


def main() -> None:
    raw = read_parquet("e4_raw.parquet", columns=["distribution", "cadence", "closure_anom_h"])
    tab = read_csv("e4_cap_probabilities.csv")
    dists = [d for d in DIST_LABELS if d in set(raw["distribution"])]
    fig = plt.figure(figsize=(6.5, 5.0))
    gs = fig.add_gridspec(2, 3, height_ratios=[1.0, 0.7], hspace=0.55, wspace=0.6)
    axes = [fig.add_subplot(gs[0, i]) for i in range(3)]
    legend_axes = [fig.add_subplot(gs[1, i]) for i in range(3)]
    for lax in legend_axes:
        lax.axis("off")

    ax = axes[0]
    rec = raw[raw["distribution"] == "record_flights_1_to_4"]
    cads = sorted(rec["cadence"].unique())
    bins = np.linspace(0, max(rec["closure_anom_h"].max(), 2100), 80)
    for c in cads:
        h = rec[rec["cadence"] == c]["closure_anom_h"]
        ax.hist(
            h,
            bins=bins,
            histtype="step",
            color=CADENCE_COLORS[int(c)],
            label=f"{int(c)} launches per year",
            density=True,
        )
    for cap, ls in CAPS:
        ax.axvline(cap, color="gray", linestyle=ls, linewidth=0.8, label=f"{cap:,} h cap")
    ax.set_xlabel("annual hours (operations\nand anomalies)")
    ax.set_xticks([0, 3000, 6000])
    ax.set_ylabel("probability density")
    ax.set_title(r"\textbf{(a)} Annual totals, record" + "\nper-launch window (12.9 h)", fontsize=8, loc="left")
    ax.set_yticks([])
    legend_axes[0].legend(*ax.get_legend_handles_labels(), fontsize=8, loc="upper left", ncol=1)

    ax = axes[1]
    for d in dists:
        t = tab[tab["distribution"] == d].sort_values("cadence")
        ax.plot(
            t["cadence"], t["closure_anom_h_mean"], marker="o", markersize=3, color=DIST_COLORS[d], label=DIST_LABELS[d]
        )
    for cap, ls in CAPS:
        ax.axhline(cap, color="gray", linestyle=ls, linewidth=0.8, label=f"{cap:,} h cap")
    ax.set_xlabel("launches per year")
    ax.set_ylabel("expected annual closure hours")
    ax.set_yticks([0, 2000, 4000, 6000])
    ax.set_title(r"\textbf{(b)} Expected annual hours" + "\nagainst cadence", fontsize=8, loc="left")
    ax.set_xticks([25, 100, 200, 456])
    legend_axes[1].legend(*ax.get_legend_handles_labels(), fontsize=8, loc="upper left", ncol=1)

    ax = axes[2]
    cad = np.linspace(25, 456, 200)
    for cap, ls in CAPS:
        ax.plot(cad, cap / cad, color=OI["blue"], linestyle=ls, label=f"{cap:,} h cap / cadence")
    for hrs, lab, col in RECORD:
        ax.axhline(hrs, color=col, linewidth=0.9, label=lab)
    ax.set_xlabel("launches per year")
    ax.set_ylabel("hours per launch")
    ax.set_ylim(0, 30)
    ax.set_xticks([25, 100, 200, 456])
    ax.set_title(
        r"\textbf{(c)} Per-launch budget implied" + "\nby each cap, and recorded hours", fontsize=8, loc="left"
    )
    legend_axes[2].legend(*ax.get_legend_handles_labels(), fontsize=8, loc="upper left", ncol=1)
    save(fig, "fig6_closure_distributions")


if __name__ == "__main__":
    main()
