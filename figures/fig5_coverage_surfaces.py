"""Figure 5: coverage and longest-gap response surfaces over rate and recovery time (experiment E5)."""

from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np
from common import HABITAT, OI, SERIES, read_csv, save


def main() -> None:
    s = read_csv("e5_surfaces.csv")
    s = s[s["receptor"] == HABITAT]
    cads = sorted(s["cadence"].unique())
    fig, axes = plt.subplots(1, 3, figsize=(7.2, 2.4))

    ax = axes[0]
    for c, col in zip(cads, SERIES, strict=False):
        p = s[(s["model"] == "poisson") & (s["cadence"] == c)].sort_values("recovery_days")
        w = s[(s["model"] == "weekday_renewal") & (s["cadence"] == c)].sort_values("recovery_days")
        ax.plot(p["recovery_days"] * 24, p["disturbed"], color=col, label=f"{int(c)}/yr Poisson")
        ax.plot(w["recovery_days"] * 24, w["disturbed"], color=col, linestyle="--", linewidth=0.9)
        ax.plot(p["recovery_days"] * 24, p["coverage_closed_form"], color=col, linestyle=":", linewidth=0.7)
    ax.set_xscale("log")
    ax.set_xlabel("recovery time $\\tau$ (h)")
    ax.set_ylabel("disturbed fraction of season")
    ax.set_ylim(0, 1)
    ax.set_title("Coverage (solid Poisson, dashed weekday; dotted closed form)", fontsize=6.5)
    ax.legend(fontsize=5, loc="upper left")

    ax = axes[1]
    for c, col in zip(cads, SERIES, strict=False):
        p = s[(s["model"] == "poisson") & (s["cadence"] == c)].sort_values("recovery_days")
        over = (p["coverage_linear"] - p["disturbed"]) / p["disturbed"]
        ax.plot(p["lambda_tau"], over, color=col, label=f"{int(c)}/yr")
    ax.axhline(0.1, color="gray", linewidth=0.6, linestyle=":")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("$\\lambda E[\\tau]$")
    ax.set_ylabel("linear overstatement of coverage")
    ax.set_title("Linear-extrapolation error of coverage", fontsize=6.5)

    ax = axes[2]
    wcols = [c for c in s.columns if c.startswith("has_window|")]
    ws = [float(c.split("|")[1]) for c in wcols]
    top = s[(s["model"] == "weekday_renewal") & (s["cadence"] == max(cads))].sort_values("recovery_days")
    z = np.array([[1 - top[col].to_numpy()[i] for col in wcols] for i in range(len(top))])
    im = ax.imshow(z, aspect="auto", origin="lower", cmap="magma", vmin=0, vmax=1,
                   extent=[-0.5, len(ws) - 0.5, -0.5, len(top) - 0.5])
    ax.set_xticks(range(len(ws)))
    ax.set_xticklabels([f"{w*24:g}h" if w < 1 else f"{w:g}d" for w in ws], fontsize=5.5)
    ax.set_yticks(range(len(top)))
    ax.set_yticklabels([f"{t*24:g}h" if t < 1 else f"{t:g}d" for t in top["recovery_days"]], fontsize=5.5)
    ax.set_xlabel("required undisturbed window $w$")
    ax.set_ylabel("recovery time $\\tau$")
    ax.set_title(f"P(no window $\\geq w$ in season), {int(max(cads))}/yr weekday", fontsize=6.5)
    fig.colorbar(im, ax=ax, shrink=0.8)
    fig.suptitle("Disturbance-coverage regimes at the piping plover critical-habitat receptor (events > 90 dBA or > 1 psf); "
                 "$\\tau$ and $w$ swept, not species-specific", fontsize=7)
    fig.tight_layout()
    save(fig, "fig5_coverage_surfaces")
    _ = OI


if __name__ == "__main__":
    main()
