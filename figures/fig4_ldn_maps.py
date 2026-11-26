"""Attachment Figure 2: DNL fields at 25, 100 and 456 launches per year, with night-share comparison."""

from __future__ import annotations

import geopandas as gpd
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import xarray as xr
from matplotlib.patches import Circle

from common import OI, RECEPTOR_LETTERS, RESULTS, ROOT, panel_label, save

GPKG = ROOT / "data" / "processed" / "site_layers.gpkg"
CASES = [
    ("existing_ldn", r"25/yr, 12\% night (existing)"),
    ("phase1_ldn", r"100/yr, 50\% night (Phase 1)"),
    ("phase2_ldn", r"456/yr, 50\% night (Phase 2)"),
    ("phase2_fn012_ldn", r"456/yr, 12\% night"),
]
RADIUS_KM = 30.0
LETTER_OFFSETS = {
    "A": (-2.8, 1.2),
    "B": (1.2, 0.8),
    "C": (-2.6, 0.6),
    "D": (-0.2, 1.6),
    "E": (-2.6, -0.6),
    "F": (-7.0, -9.0),
    "G": (1.2, 0.6),
    "H": (-1.0, -11.0),
    "I": (-11.0, -4.0),
}
LEADER = {"F", "H", "I"}  # receptors within 4 km of the pad are labelled with a leader line


def main() -> None:
    ds = xr.open_dataset(RESULTS / "ldn_fields.nc")
    rec = pd.read_csv(RESULTS / "receptors_used.csv").set_index("name")
    pad_x, pad_y = float(ds.attrs["pad_x"]), float(ds.attrs["pad_y"])
    x = (ds["x"].values - pad_x) / 1000.0
    y = (ds["y"].values - pad_y) / 1000.0
    water = gpd.read_file(GPKG, layer="water")
    # One dissolved water body so that only shorelines are drawn; the county's offshore limit (a straight
    # north-south line) is dropped by masking vertices east of it.
    # One dissolved water body so that only shorelines are drawn (bays, lagoons and the Gulf; not ponds). The
    # county's offshore limit, a straight north-south line, is not a shoreline and is masked out.
    gulf_east = max(
        g.bounds[2] for g, n in zip(water.geometry, water["FULLNAME"], strict=True) if n == "Gulf of America"
    )
    dissolved = water[water["AWATER"] > 3.0e6].geometry.union_all()
    water_geoms = list(getattr(dissolved, "geoms", [dissolved]))
    fig, axes = plt.subplots(2, 2, figsize=(6.5, 6.6), sharex=True, sharey=True)
    cf_ = None
    for ax, (key, title), letter in zip(axes.ravel(), CASES, "abcd", strict=True):
        f = ds[key].values
        clip = Circle((0, 0), RADIUS_KM, transform=ax.transData)
        ax.add_patch(Circle((0, 0), RADIUS_KM, facecolor="white", edgecolor="none", zorder=0))
        cf_ = ax.contourf(x, y, f, levels=[50, 55, 60, 65, 70, 75, 90], cmap="viridis", zorder=2)
        cf_.set_clip_path(clip)
        cs = ax.contour(x, y, f, levels=[55, 65, 75], colors="white", linewidths=0.7, zorder=3)
        cs.set_clip_path(clip)
        ax.clabel(cs, fmt="%d", fontsize=8, inline=True, inline_spacing=2)
        for label in cs.labelTexts:  # drop labels that the circular boundary would cut
            lx, ly = label.get_position()
            if np.hypot(lx, ly) > RADIUS_KM - 3.0:
                label.set_visible(False)
        for geom in water_geoms:
            for ring in [geom.exterior, *geom.interiors]:
                gx, gy = (np.asarray(v) for v in ring.xy)
                gx = np.where(gx > gulf_east - 1500.0, np.nan, gx)
                line = ax.plot((gx - pad_x) / 1000, (gy - pad_y) / 1000, color="black", lw=0.35, zorder=4)[0]
                line.set_clip_path(clip)
        ax.add_patch(Circle((0, 0), RADIUS_KM, facecolor="none", edgecolor="black", lw=0.5, zorder=5))
        ax.plot(0, 0, marker="X", color=OI["vermillion"], markersize=7, zorder=6)
        for let, name in RECEPTOR_LETTERS.items():
            rx = (rec.loc[name, "x"] - pad_x) / 1000
            ry = (rec.loc[name, "y"] - pad_y) / 1000
            ax.plot(rx, ry, "o", color="white", markersize=3.5, markeredgecolor="black", markeredgewidth=0.5, zorder=6)
            dx, dy = LETTER_OFFSETS[let]
            box = {"boxstyle": "square,pad=0.08", "facecolor": "white", "edgecolor": "none", "alpha": 0.8}
            if let in LEADER:
                ax.annotate(
                    let,
                    (rx, ry),
                    (rx + dx, ry + dy),
                    fontsize=8,
                    ha="center",
                    va="center",
                    zorder=7,
                    arrowprops={"arrowstyle": "-", "lw": 0.5, "color": "black", "shrinkB": 2},
                    bbox=box,
                )
            else:
                ax.text(rx + dx, ry + dy, let, fontsize=8, ha="left", va="center", zorder=7, bbox=box)
        ax.set_title(title, fontsize=8)
        ax.set_xlim(-31, 31)
        ax.set_ylim(-31, 31)
        ax.set_aspect("equal")
        ax.set_xticks([-20, 0, 20])
        ax.set_yticks([-20, 0, 20])
        panel_label(ax, letter, x=-0.12, y=1.02)
    for ax in axes[1]:
        ax.set_xlabel("km east of pad")
    for ax in axes[:, 0]:
        ax.set_ylabel("km north of pad")
    cb = fig.colorbar(cf_, ax=axes, shrink=0.75, pad=0.03)
    cb.set_label("DNL (dB), expected-count field")
    save(fig, "fig4_ldn_maps")


if __name__ == "__main__":
    main()
