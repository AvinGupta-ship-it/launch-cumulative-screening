"""Shared figure settings: paths, palette and helpers.

Figures read saved results only; no computation of scientific content
happens here. The results directory is ``results/paper`` unless the
environment variable ``LAUNCHCUM_RESULTS`` points elsewhere. Text is set
with LaTeX so that figure fonts match the body font of the documents.
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

# One color per cadence level and one per timing model, used in every figure.
CADENCE_COLORS = {25: OI["black"], 50: OI["skyblue"], 100: OI["orange"], 200: OI["vermillion"], 456: OI["purple"]}
MODEL_COLORS = {"poisson": OI["blue"], "weekday_renewal": OI["green"], "campaign": OI["purple"]}
MODEL_LABELS = {"poisson": "Poisson, all days", "weekday_renewal": "weekday-only renewal", "campaign": "campaign"}
EXTRAP_COLOR = OI["vermillion"]

COMMUNITY_EDGE = "Port Isabel (nearest edge)"
COMMUNITY_SPI = "South Padre Island (nearest edge)"
HABITAT = "Piping plover critical habitat (near pad)"
REFUGE = "Lower Rio Grande Valley NWR (near pad)"
PARK = "Boca Chica State Park"

# Named receptors reported in the documents, keyed by the letter used in the maps.
RECEPTOR_LETTERS = {
    "A": COMMUNITY_EDGE,
    "B": COMMUNITY_SPI,
    "C": "Laguna Vista (nearest edge)",
    "D": "Laguna Heights (nearest edge)",
    "E": "Brownsville (nearest edge)",
    "F": "Boca Chica Village",
    "G": "Isla Blanca Park",
    "H": PARK,
    "I": HABITAT,
}

# Readable names of the sensitivity-analysis inputs (factor table in the configuration).
FACTOR_LABELS = {
    "night_fraction": "night share",
    "scrub_probability": "scrub probability",
    "turnaround_min_days": "minimum pad turnaround",
    "n_pads": "number of pads",
    "cluster_rate_per_year": "campaign rate",
    "cluster_spread_days": "campaign spread",
    "level_offset_db": "vehicle level offset",
    "boom_scale": "boom scale factor",
    "closure_launch_mean_h": "launch closure mean",
    "closure_launch_cv": "launch closure CV",
    "anomaly_probability": "anomaly probability",
    "static_fires_per_launch": "static fires per launch",
    "external_fraction": "external landing share",
    "p_ss_landing_local": "local Starship landing share",
    "recovery_days": "recovery time",
    "window_days": "required window",
    "random_seed": "random seed (Monte Carlo)",
}


def read_csv(name: str) -> pd.DataFrame:
    return pd.read_csv(RESULTS / name)


def read_parquet(name: str, columns: list[str] | None = None) -> pd.DataFrame:
    return pd.read_parquet(RESULTS / name, columns=columns)


def summary() -> dict:
    return json.loads((RESULTS / "summary.json").read_text())["data"]


def save(fig, name: str) -> None:
    fig.savefig(FIGDIR / f"{name}.pdf")
    fig.savefig(FIGDIR / f"{name}.png", dpi=200)
    print("wrote", FIGDIR / f"{name}.pdf")


def panel_label(ax, letter: str, x: float = -0.18, y: float = 1.04) -> None:
    """Bold panel letter outside the top-left corner of an axes."""
    ax.text(x, y, rf"\textbf{{({letter})}}", transform=ax.transAxes, ha="left", va="bottom", fontsize=9)


def tex(s: str) -> str:
    """Escape a data string for LaTeX text."""
    return s.replace("_", " ").replace("%", r"\%").replace("&", r"\&")


def stage_label(df_label: str) -> str:
    return {
        "scenario_existing": "Existing (25/yr)",
        "scenario_phase1": "Phase 1 (100/yr)",
        "scenario_phase2": "Phase 2 (456/yr, weekday)",
        "scenario_phase2_clustered": "Phase 2, clustered",
        "scenario_phase2_hypothetical_upscaled_vehicle": "Phase 2, hypothetical upscaled vehicle",
    }.get(df_label, df_label)
