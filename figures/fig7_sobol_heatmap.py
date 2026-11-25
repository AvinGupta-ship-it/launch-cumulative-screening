"""Figure 7: Sobol total-order indices (metric x input x cadence), experiment E6."""

from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np

from common import read_csv, save

OUTPUT_LABELS = {
    "ldn": "$L_{dn}$ (community)",
    "cdnl": "CDNL (community)",
    "na100": "N > 100 dBA (community)",
    "days90": "days > 90 dBA (community)",
    "boom2": "booms > 2 psf (community)",
    "disturbed": "disturbed fraction (habitat)",
    "longest_gap": "longest gap (habitat)",
    "closure_h": "closure hours",
    "weekend_avail": "weekend availability",
}


def main() -> None:
    si = read_csv("e6_sobol_indices.csv")
    cads = sorted(si["cadence"].unique())
    factors = list(dict.fromkeys(si["factor"]))
    outputs = [o for o in OUTPUT_LABELS if o in set(si["output"])]
    fig, axes = plt.subplots(1, len(cads), figsize=(7.2, 3.6), sharey=True)
    for ax, c in zip(np.atleast_1d(axes), cads, strict=True):
        sub = si[si["cadence"] == c]
        z = np.array(
            [[float(sub[(sub["output"] == o) & (sub["factor"] == f)]["ST"].iloc[0]) for f in factors] for o in outputs]
        )
        im = ax.imshow(np.clip(z, 0, 1), cmap="viridis", vmin=0, vmax=1, aspect="auto")
        ax.set_xticks(range(len(factors)))
        ax.set_xticklabels([f.replace("_", " ") for f in factors], rotation=90, fontsize=5.5)
        ax.set_yticks(range(len(outputs)))
        ax.set_yticklabels([OUTPUT_LABELS[o] for o in outputs], fontsize=6)
        ax.set_title(f"{int(c)} launches/yr", fontsize=7)
        for i in range(len(outputs)):
            for j in range(len(factors)):
                if z[i, j] >= 0.3:
                    ax.text(
                        j,
                        i,
                        f"{min(z[i, j], 1.0):.1f}",
                        ha="center",
                        va="center",
                        fontsize=4.5,
                        color="white" if z[i, j] < 0.6 else "black",
                    )
    fig.colorbar(im, ax=axes, shrink=0.7, pad=0.02, label="total-order Sobol index $S_T$")
    fig.suptitle(
        "Dominant uncertainties by metric and cadence "
        "(campaign timing on business days; inputs sampled over documented ranges)",
        fontsize=7,
    )
    save(fig, "fig7_sobol_heatmap")


if __name__ == "__main__":
    main()
