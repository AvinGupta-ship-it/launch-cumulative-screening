"""Figure 4: L_dn contour maps at 25, 100 and 456 launches per year, with night-share comparison."""

from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import xarray as xr

from common import OI, RESULTS, save


def main() -> None:
    ds = xr.open_dataset(RESULTS / "ldn_fields.nc")
    rec = pd.read_csv(RESULTS / "receptors_used.csv")
    pad_x, pad_y = float(ds.attrs["pad_x"]), float(ds.attrs["pad_y"])
    x = (ds["x"].values - pad_x) / 1000.0
    y = (ds["y"].values - pad_y) / 1000.0
    cases = [
        ("existing_ldn", "25/yr, 12% night (existing)"),
        ("phase1_ldn", "100/yr, 50% night (Phase 1)"),
        ("phase2_ldn", "456/yr, 50% night (Phase 2)"),
        ("phase2_fn012_ldn", "456/yr, 12% night"),
    ]
    levels = [55, 60, 65, 70, 75]
    fig, axes = plt.subplots(1, 4, figsize=(7.2, 2.3), sharex=True, sharey=True)
    for ax, (key, title) in zip(axes, cases, strict=True):
        f = ds[key].values
        cf_ = ax.contourf(x, y, f, levels=[50, 55, 60, 65, 70, 75, 90], cmap="viridis", alpha=0.85)
        cs = ax.contour(x, y, f, levels=levels, colors="white", linewidths=0.5)
        ax.clabel(cs, fmt="%d", fontsize=5, inline=True)
        ax.plot(0, 0, marker="X", color=OI["vermillion"], markersize=6)
        comm = rec[rec["kind"] == "community"]
        ax.plot(
            (comm["x"] - pad_x) / 1000,
            (comm["y"] - pad_y) / 1000,
            "o",
            color="white",
            markersize=2.5,
            markeredgecolor="black",
            markeredgewidth=0.4,
        )
        ax.set_title(title, fontsize=6.5)
        ax.set_xlim(-30, 30)
        ax.set_ylim(-30, 30)
        ax.set_aspect("equal")
        ax.set_xlabel("km east of pad")
    axes[0].set_ylabel("km north of pad")
    cb = fig.colorbar(cf_, ax=axes, shrink=0.8, pad=0.02)
    cb.set_label("$L_{dn}$ (dB), expected-count field")
    fig.suptitle(
        "Day-night average level from fitted FAA single-event SEL footprints; white dots are community receptors",
        fontsize=7.5,
    )
    save(fig, "fig4_ldn_maps")
    _ = np.nan


if __name__ == "__main__":
    main()
