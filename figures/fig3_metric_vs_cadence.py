"""Figure 3: metrics versus cadence with linear-extrapolation overlays (experiment E1)."""

from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np

from common import COMMUNITY_EDGE, COMMUNITY_SPI, HABITAT, OI, read_parquet, save


def main() -> None:
    raw = read_parquet("e1_raw.parquet")
    sub = raw[raw["external_fraction"] == 0.0]
    g = sub.groupby("cadence")
    cad = np.array(sorted(sub["cadence"].unique()))
    ref = cad.min()
    panels = [
        ("deluge_gal", "Deluge water (million gal/yr)", 1e-6, "linear"),
        (f"ldn|{COMMUNITY_EDGE}", "$L_{dn}$ at Port Isabel edge (dB)", 1.0, "energy"),
        (f"na100|{COMMUNITY_EDGE}", "Events > 100 dBA $L_{Amax}$, Port Isabel edge", 1.0, "linear"),
        (f"days90|{COMMUNITY_SPI}", "Days with > 90 dBA, South Padre Island edge", 1.0, "linear"),
        (f"boom2|{COMMUNITY_EDGE}", "Booms > 2 psf, Port Isabel edge", 1.0, "linear"),
        (f"disturbed|{HABITAT}", "Disturbed fraction, habitat receptor (1 h recovery)", 1.0, "linear"),
        (f"longest_gap|{HABITAT}", "Longest undisturbed gap in season (days)", 1.0, "none"),
        ("closure_h", "Access-restriction hours per year", 1.0, "linear"),
        ("weekend_avail", "Weekend/holiday daylight available (fraction)", 1.0, "none"),
    ]
    fig, axes = plt.subplots(3, 3, figsize=(7.2, 6.2))
    for ax, (col, label, scale, mode) in zip(axes.ravel(), panels, strict=True):
        mean = g[col].mean().reindex(cad).to_numpy() * scale
        lo = g[col].quantile(0.025).reindex(cad).to_numpy() * scale
        hi = g[col].quantile(0.975).reindex(cad).to_numpy() * scale
        ax.fill_between(cad, lo, hi, color=OI["blue"], alpha=0.15, linewidth=0)
        ax.plot(cad, mean, marker="o", color=OI["blue"], markersize=3, label="simulated mean (95% interval)")
        if mode == "linear":
            ax.plot(
                cad,
                mean[0] * cad / ref,
                color=OI["vermillion"],
                linestyle="--",
                linewidth=1.0,
                label="linear extrapolation from 25/yr",
            )
        elif mode == "energy":
            ax.plot(
                cad,
                mean[0] + 10 * np.log10(cad / ref),
                color=OI["vermillion"],
                linestyle="--",
                linewidth=1.0,
                label="energy-average extrapolation",
            )
        ax.set_title(label, fontsize=7)
        ax.set_xlabel("launches per year")
        ax.set_xticks([25, 100, 200, 456])
    axes[0, 0].legend(fontsize=5.5, loc="upper left")
    axes[0, 1].legend(fontsize=5.5, loc="lower right")
    fig.suptitle(
        "Scaling with cadence (Poisson timing, 50% night share, no external landings); hypothetical closure parameters",
        fontsize=7.5,
    )
    fig.tight_layout()
    save(fig, "fig3_metric_vs_cadence")


if __name__ == "__main__":
    main()
