"""Comment Figure 1: exposure days against cadence under three timing models (experiments E1 and E2)."""

from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np

from common import COMMUNITY_SPI, EXTRAP_COLOR, MODEL_COLORS, MODEL_LABELS, read_parquet, save

BUSINESS_DAYS_2030 = 250


def main() -> None:
    col = f"days90|{COMMUNITY_SPI}"
    raw = read_parquet("e1_raw.parquet", columns=["cadence", "external_fraction", col])
    sub = raw[raw["external_fraction"] == 0.0]
    g = sub.groupby("cadence")[col]
    cad = np.array(sorted(sub["cadence"].unique()))
    mean = g.mean().reindex(cad).to_numpy()
    lo = g.quantile(0.025).reindex(cad).to_numpy()
    hi = g.quantile(0.975).reindex(cad).to_numpy()
    e2 = read_parquet("e2_raw.parquet", columns=["cadence", "model", col])
    fig, ax = plt.subplots(figsize=(4.8, 3.3))
    ax.plot(cad, mean[0] * cad / cad.min(), color=EXTRAP_COLOR, linestyle="--", label="linear extrapolation from 25/yr")
    ax.fill_between(cad, lo, hi, color=MODEL_COLORS["poisson"], alpha=0.15, linewidth=0)
    ax.plot(
        cad, mean, marker="o", markersize=4, color=MODEL_COLORS["poisson"], label="simulated, Poisson timing, all days"
    )
    for model, marker in (("weekday_renewal", "s"), ("campaign", "^")):
        m = e2[e2["model"] == model].groupby("cadence")[col].mean()
        ax.plot(
            m.index,
            m.to_numpy(),
            marker=marker,
            linestyle="none",
            markersize=5.5,
            markerfacecolor="white",
            markeredgewidth=1.3,
            color=MODEL_COLORS[model],
            label=f"simulated, {MODEL_LABELS[model]}",
            zorder=4,
        )
    ax.axhline(365, color="gray", linewidth=0.7, linestyle=":")
    ax.axhline(BUSINESS_DAYS_2030, color="gray", linewidth=0.7, linestyle="-.")
    ax.text(30, 373, "365 days in the year", fontsize=8, color="gray", va="bottom")
    ax.text(30, BUSINESS_DAYS_2030 - 22, f"{BUSINESS_DAYS_2030} business days", fontsize=8, color="gray", va="top")
    ax.set_xlabel("launches per year")
    ax.set_ylabel(r"days per year with an event $>$ 90 dBA" + "\n" + "at the South Padre Island edge")
    ax.set_xticks([25, 100, 200, 456])
    ax.set_ylim(0, 850)
    ax.legend(loc="upper left", fontsize=8, bbox_to_anchor=(0.0, 0.98))
    fig.tight_layout()
    save(fig, "figc1_exposure_days")


if __name__ == "__main__":
    main()
