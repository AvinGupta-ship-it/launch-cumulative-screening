"""Attachment Figure 1: metrics against cadence with extrapolation overlays (experiments E1 and E2)."""

from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np

from common import (
    COMMUNITY_EDGE,
    COMMUNITY_SPI,
    EXTRAP_COLOR,
    HABITAT,
    MODEL_COLORS,
    MODEL_LABELS,
    read_parquet,
    save,
)

BUSINESS_DAYS_2030 = 250
PANELS = [
    ("deluge_gal", "Deluge water\n(million gal per yr)", 1e-6, "linear"),
    (f"ldn|{COMMUNITY_EDGE}", "DNL at Port Isabel edge\n(dB)", 1.0, "energy"),
    (f"na100|{COMMUNITY_EDGE}", "Events $>$ 100 dBA $L_{A\\max}$,\nPort Isabel edge (per yr)", 1.0, "linear"),
    (f"days90|{COMMUNITY_SPI}", "Days with $>$ 90 dBA,\nSouth Padre Island edge (per yr)", 1.0, "linear"),
    (f"boom2|{COMMUNITY_EDGE}", "Booms $>$ 2 psf,\nPort Isabel edge (per yr)", 1.0, "linear"),
    (f"disturbed|{HABITAT}", "Disturbed fraction of season,\nhabitat (1 h recovery)", 1.0, "linear"),
    (f"longest_gap|{HABITAT}", "Longest undisturbed gap\nin season, habitat (h)", 24.0, "none"),
    ("closure_h", "Access-restriction hours\nper year (operations)", 1.0, "linear"),
    ("weekend_avail", "Weekend and holiday daylight\navailable (fraction)", 1.0, "none"),
]


def main() -> None:
    raw = read_parquet("e1_raw.parquet")
    sub = raw[raw["external_fraction"] == 0.0]
    g = sub.groupby("cadence")
    cad = np.array(sorted(sub["cadence"].unique()))
    ref = cad.min()
    e2 = read_parquet("e2_raw.parquet")
    e2g = e2.groupby(["model", "cadence"])
    fig, axes = plt.subplots(3, 3, figsize=(6.5, 6.6))
    letters = "abcdefghi"
    for ax, (col, label, scale, mode), letter in zip(axes.ravel(), PANELS, letters, strict=True):
        mean = g[col].mean().reindex(cad).to_numpy() * scale
        lo = g[col].quantile(0.025).reindex(cad).to_numpy() * scale
        hi = g[col].quantile(0.975).reindex(cad).to_numpy() * scale
        ax.fill_between(cad, lo, hi, color=MODEL_COLORS["poisson"], alpha=0.15, linewidth=0)
        ax.plot(cad, mean, marker="o", color=MODEL_COLORS["poisson"], markersize=3.5, label="Poisson, all days")
        if mode == "linear":
            ax.plot(cad, mean[0] * cad / ref, color=EXTRAP_COLOR, linestyle="--", linewidth=1.0)
        elif mode == "energy":
            ax.plot(cad, mean[0] + 10 * np.log10(cad / ref), color=EXTRAP_COLOR, linestyle="--", linewidth=1.0)
        for model, marker in (("weekday_renewal", "s"), ("campaign", "^")):
            xs, ys = [], []
            for c in sorted(e2["cadence"].unique()):
                xs.append(c)
                ys.append(float(e2g[col].mean().loc[(model, c)]) * scale)
            ax.plot(
                xs,
                ys,
                marker=marker,
                linestyle="none",
                markersize=4.5,
                markerfacecolor="white",
                markeredgewidth=1.2,
                color=MODEL_COLORS[model],
                label=MODEL_LABELS[model],
                zorder=4,
            )
        if col.startswith("days90"):
            ax.axhline(365, color="gray", linewidth=0.7, linestyle=":")
            ax.axhline(BUSINESS_DAYS_2030, color="gray", linewidth=0.7, linestyle="-.")
            # Labels sit where no data lie: the extrapolation line and the band are below 365 at the left and
            # above 250 at the right.
            ax.text(30, 365 + 12, "365 days", fontsize=8, color="gray", va="bottom", ha="left")
            ax.text(
                420,
                BUSINESS_DAYS_2030 - 14,
                f"{BUSINESS_DAYS_2030} business\ndays",
                fontsize=8,
                color="gray",
                va="top",
                ha="right",
                linespacing=1.0,
            )
        ax.set_title(
            rf"\textbf{{({letter})}} " + label.split("\n")[0] + "\n" + label.split("\n")[1], fontsize=8, loc="left"
        )
        ax.set_xlabel("launches per year")
        ax.set_xticks([25, 100, 200, 456])
    handles, labels = axes[0, 0].get_legend_handles_labels()
    handles.append(plt.Line2D([], [], color=EXTRAP_COLOR, linestyle="--"))
    labels.append("linear or energy-average extrapolation from 25/yr")
    handles.insert(1, plt.Rectangle((0, 0), 1, 1, color=MODEL_COLORS["poisson"], alpha=0.15))
    labels.insert(1, r"95\% Monte Carlo interval (Poisson)")
    fig.legend(handles, labels, loc="lower center", ncol=3, bbox_to_anchor=(0.5, -0.045), fontsize=8)
    fig.tight_layout(h_pad=1.0, w_pad=1.0)
    save(fig, "fig3_metric_vs_cadence")


if __name__ == "__main__":
    main()
