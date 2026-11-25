"""Shared figure settings: paths, palette and helpers.

Figures read saved results only; no computation of scientific content
happens here. The results directory is ``results/paper`` unless the
environment variable ``LAUNCHCUM_RESULTS`` points elsewhere.
"""

from __future__ import annotations

import json
import os
from pathlib import Path

import matplotlib as mpl
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
RESULTS = Path(os.environ.get("LAUNCHCUM_RESULTS", ROOT / "results" / "paper"))
FIGDIR = ROOT / "paper" / "figures"
FIGDIR.mkdir(parents=True, exist_ok=True)

mpl.rc_file(ROOT / "figures" / "matplotlibrc")

# Okabe-Ito palette (colorblind-safe)
OI = {
    "black": "#000000",
    "orange": "#E69F00",
    "skyblue": "#56B4E9",
    "green": "#009E73",
    "yellow": "#F0E442",
    "blue": "#0072B2",
    "vermillion": "#D55E00",
    "purple": "#CC79A7",
}
SERIES = [OI["blue"], OI["vermillion"], OI["green"], OI["orange"], OI["purple"], OI["skyblue"], OI["black"]]

COMMUNITY_EDGE = "Port Isabel (nearest edge)"
COMMUNITY_SPI = "South Padre Island (nearest edge)"
HABITAT = "Piping plover critical habitat (near pad)"
REFUGE = "Lower Rio Grande Valley NWR (near pad)"
PARK = "Boca Chica State Park"


def read_csv(name: str) -> pd.DataFrame:
    return pd.read_csv(RESULTS / name)


def read_parquet(name: str) -> pd.DataFrame:
    return pd.read_parquet(RESULTS / name)


def summary() -> dict:
    return json.loads((RESULTS / "summary.json").read_text())["data"]


def save(fig, name: str) -> None:
    fig.savefig(FIGDIR / f"{name}.pdf")
    fig.savefig(FIGDIR / f"{name}.png", dpi=200)
    print("wrote", FIGDIR / f"{name}.pdf")


def stage_label(df_label: str) -> str:
    return {
        "scenario_existing": "Existing (25/yr)",
        "scenario_phase1": "Phase 1 (100/yr)",
        "scenario_phase2": "Phase 2 (456/yr, weekday)",
        "scenario_phase2_clustered": "Phase 2, clustered",
        "scenario_phase2_hypothetical_upscaled_vehicle": "Phase 2, hypothetical upscaled vehicle",
    }.get(df_label, df_label)
