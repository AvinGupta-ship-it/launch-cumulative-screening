"""Attachment Figure 5: Sobol total-order indices (metric x input x cadence), experiment E6."""

from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np

from common import FACTOR_LABELS, read_csv, save

OUTPUT_LABELS = {
    "ldn": "DNL (community)",
    "cdnl": "CDNL (community)",
    "na100": r"events $>$ 100 dBA (community)",
    "days90": r"days $>$ 90 dBA (community)",
    "boom2": r"booms $>$ 2 psf (community)",
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
    fig, axes = plt.subplots(len(cads), 1, figsize=(6.5, 8.2), sharex=True)
    im = None
    for ax, c, letter in zip(np.atleast_1d(axes), cads, "abc", strict=True):
        sub = si[si["cadence"] == c]
        z = np.array(
            [[float(sub[(sub["output"] == o) & (sub["factor"] == f)]["ST"].iloc[0]) for f in factors] for o in outputs]
        )
        im = ax.imshow(np.clip(z, 0, 1), cmap="viridis", vmin=0, vmax=1, aspect="auto")
        ax.set_xticks(range(len(factors)))
        ax.set_xticklabels(
            [FACTOR_LABELS.get(f, f.replace("_", " ")) for f in factors], rotation=60, ha="right", fontsize=8
        )
        ax.set_yticks(range(len(outputs)))
        ax.set_yticklabels([OUTPUT_LABELS[o] for o in outputs], fontsize=8)
        ax.set_title(rf"\textbf{{({letter})}} {int(c)} launches per year", fontsize=8, loc="left")
        for i in range(len(outputs)):
            for j in range(len(factors)):
                if z[i, j] >= 0.2:
                    ax.text(
                        j,
                        i,
                        f"{min(z[i, j], 1.0):.1f}",
                        ha="center",
                        va="center",
                        fontsize=8,
                        color="white" if z[i, j] < 0.6 else "black",
                    )
        ax.tick_params(axis="x", labelbottom=(c == cads[-1]))
    cb = fig.colorbar(im, ax=axes, shrink=0.5, pad=0.02, aspect=30)
    cb.set_label(r"total-order Sobol index $S_T$")
    save(fig, "fig7_sobol_heatmap")


if __name__ == "__main__":
    main()
