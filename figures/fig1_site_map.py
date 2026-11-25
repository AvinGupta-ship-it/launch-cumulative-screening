"""Figure 1: site map with receptors, protected areas, habitat and footprints."""

from __future__ import annotations

import geopandas as gpd
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import xarray as xr

from common import OI, RESULTS, ROOT, save

GPKG = ROOT / "data" / "processed" / "site_layers.gpkg"


def main() -> None:
    rec = pd.read_csv(RESULTS / "receptors_used.csv")
    blocks = gpd.read_file(GPKG, layer="blocks")
    refuges = gpd.read_file(GPKG, layer="refuges")
    ch = gpd.read_file(GPKG, layer="critical_habitat")
    parks = gpd.read_file(GPKG, layer="parks")
    ds = xr.open_dataset(RESULTS / "ldn_fields.nc")
    pad_x, pad_y = float(ds.attrs["pad_x"]), float(ds.attrs["pad_y"])

    fig, ax = plt.subplots(figsize=(6.5, 6.0))
    land = blocks[blocks["AREALAND"] > 0]
    water = blocks[blocks["AREALAND"] == 0]
    water.plot(ax=ax, color="#dbe9f4", edgecolor="none")
    land.plot(ax=ax, color="#f2efe9", edgecolor="none")
    refuges.plot(ax=ax, color="#cfe8cf", edgecolor="#4f8f4f", linewidth=0.3, alpha=0.9)
    ch[ch["comname"] == "Piping Plover"].plot(
        ax=ax, facecolor="none", edgecolor=OI["purple"], linewidth=0.8, hatch="///"
    )
    ch[ch["comname"] != "Piping Plover"].plot(
        ax=ax, facecolor="none", edgecolor=OI["orange"], linewidth=0.8, linestyle="--"
    )
    parks.plot(ax=ax, facecolor="#e8d8a8", edgecolor="#8a6d1a", linewidth=0.5)

    x = ds["x"].values
    y = ds["y"].values
    for key, color in (("existing_ldn", OI["blue"]), ("phase2_ldn", OI["vermillion"])):
        ax.contour(x, y, ds[key].values, levels=[65.0], colors=[color], linewidths=1.3)
    ax.plot([], [], color=OI["blue"], label="DNL 65 dB, existing (25/yr)")
    ax.plot([], [], color=OI["vermillion"], label="DNL 65 dB, Phase 2 (456/yr, 50% night)")
    cs2 = ax.contour(
        x, y, ds["phase2_cdnl"].values, levels=[60.0], colors=[OI["green"]], linewidths=1.0, linestyles="dotted"
    )
    ax.plot([], [], color=OI["green"], linestyle="dotted", label="CDNL 60 dB, Phase 2")

    kinds = {
        "community": ("o", OI["black"]),
        "community_edge": ("s", OI["black"]),
        "park": ("^", "#8a6d1a"),
        "refuge": ("D", "#2f6f2f"),
        "critical_habitat": ("*", OI["purple"]),
    }
    for _, r in rec.iterrows():
        m, c = kinds.get(r["kind"], ("o", "k"))
        ax.plot(
            r["x"],
            r["y"],
            marker=m,
            color=c,
            markersize=5 if m != "*" else 8,
            linestyle="none",
            markeredgecolor="white",
            markeredgewidth=0.4,
        )
        if r["kind"] in ("community", "park") and r["r_km"] < 32:
            ax.annotate(
                r["name"].replace(" (near pad)", "").replace(" (nearest)", ""),
                (r["x"], r["y"]),
                xytext=(3, 3),
                textcoords="offset points",
                fontsize=5.5,
            )
    ax.plot(
        pad_x,
        pad_y,
        marker="X",
        color=OI["vermillion"],
        markersize=9,
        linestyle="none",
        label="Launch pad (reference point)",
    )
    for k_km in (10, 20):
        th = np.linspace(0, 2 * np.pi, 200)
        ax.plot(
            pad_x + 1000 * k_km * np.cos(th),
            pad_y + 1000 * k_km * np.sin(th),
            color="gray",
            linewidth=0.4,
            linestyle=":",
        )
        ax.annotate(f"{k_km} km", (pad_x + 1000 * k_km * 0.72, pad_y + 1000 * k_km * 0.72), fontsize=5.5, color="gray")
    ax.set_xlim(pad_x - 30_000, pad_x + 30_000)
    ax.set_ylim(pad_y - 30_000, pad_y + 30_000)
    ax.set_aspect("equal")
    ax.set_xlabel("UTM 14N easting (m)")
    ax.set_ylabel("UTM 14N northing (m)")
    ax.legend(loc="lower left", fontsize=6)
    ax.set_title("Study domain, receptors and cumulative-noise contours (expected-count fields)", fontsize=8)
    save(fig, "fig1_site_map")


if __name__ == "__main__":
    main()
