"""Figure 6: annual closure-hour distributions with caps (experiment E4)."""

from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np

from common import OI, SERIES, read_csv, read_parquet, save


def main() -> None:
    raw = read_parquet("e4_raw.parquet")
    tab = read_csv("e4_cap_probabilities.csv")
    dists = [d for d in raw["distribution"].unique() if not d.startswith("record_weekday")]
    fig, axes = plt.subplots(1, 3, figsize=(7.2, 2.4))

    ax = axes[0]
    rec = raw[raw["distribution"] == dists[0]]
    cads = sorted(rec["cadence"].unique())
    bins = np.linspace(0, max(rec["closure_anom_h"].max(), 2100), 80)
    for c, col in zip(cads, SERIES, strict=False):
        h = rec[rec["cadence"] == c]["closure_anom_h"]
        ax.hist(h, bins=bins, histtype="step", color=col, label=f"{int(c)}/yr", density=True)
    for cap, ls in ((500, ":"), (800, "--"), (2000, "-")):
        ax.axvline(cap, color="gray", linestyle=ls, linewidth=0.8)
    ax.set_xlabel("annual access-restriction hours (operations + anomalies)")
    ax.set_ylabel("density")
    ax.set_title(f"Per-launch window: {dists[0].replace('_', ' ')}", fontsize=6.5)
    ax.legend(fontsize=5.5)

    ax = axes[1]
    for d, col in zip(dists, SERIES, strict=False):
        t = tab[tab["distribution"] == d].sort_values("cadence")
        ax.plot(
            t["cadence"],
            t["closure_anom_h_mean"],
            marker="o",
            markersize=3,
            color=col,
            label=f"{d.replace('_', ' ')} (mean {t['mean_h'].iloc[0]:g} h/launch)",
        )
    for cap, ls in ((500, ":"), (800, "--"), (2000, "-")):
        ax.axhline(cap, color="gray", linestyle=ls, linewidth=0.8)
        ax.annotate(f"{cap} h", (30, cap * 1.05), fontsize=5.5, color="gray")
    ax.set_xlabel("launches per year")
    ax.set_ylabel("mean annual closure hours")
    ax.set_title("Expected closure hours vs cadence", fontsize=6.5)
    ax.legend(fontsize=5)

    ax = axes[2]
    cad = np.linspace(25, 456, 100)
    for cap, ls in ((500, ":"), (800, "--"), (2000, "-")):
        ax.plot(cad, cap / cad, color=OI["blue"], linestyle=ls, label=f"{cap} h cap")
    for hrs, lab, col in (
        (24.5, "Flight 1 (24.5 h)", OI["vermillion"]),
        (10.5, "Flight 2 (10.5 h)", OI["orange"]),
        (7.75, "Flight 4 (7.75 h)", OI["green"]),
    ):
        ax.axhline(hrs, color=col, linewidth=0.8, label=lab)
    ax.set_yscale("log")
    ax.set_xlabel("launches per year")
    ax.set_ylabel("closure budget per launch (h)")
    ax.set_title("Cap divided by cadence vs recorded per-launch hours", fontsize=6.5)
    ax.legend(fontsize=5)
    fig.suptitle(
        "Public-access restriction hours: Monte Carlo distributions "
        "(hypothetical per-event closure parameters bracketing the FAA record)",
        fontsize=7,
    )
    fig.tight_layout()
    save(fig, "fig6_closure_distributions")


if __name__ == "__main__":
    main()
