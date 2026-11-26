"""Attachment Figure 4: coverage and interval metrics over recovery time and window (experiment E5)."""

from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np

from common import CADENCE_COLORS, HABITAT, read_csv, save


def duration_label(days: float) -> str:
    minutes = days * 1440.0
    if minutes < 59:
        return f"{minutes:.0f} min"
    hours = days * 24.0
    if hours < 23.5:
        return f"{hours:.0f} h"
    return f"{days:.0f} d"


def main() -> None:
    s = read_csv("e5_surfaces.csv")
    s = s[s["receptor"] == HABITAT]
    cads = sorted(s["cadence"].unique())
    fig, axes = plt.subplots(2, 2, figsize=(6.5, 6.0))

    ax = axes[0, 0]
    for c in cads:
        col = CADENCE_COLORS[int(c)]
        p = s[(s["model"] == "poisson") & (s["cadence"] == c)].sort_values("recovery_days")
        w = s[(s["model"] == "weekday_renewal") & (s["cadence"] == c)].sort_values("recovery_days")
        ax.plot(p["recovery_days"] * 24, p["disturbed"], color=col, label=f"{int(c)} launches per year")
        ax.plot(w["recovery_days"] * 24, w["disturbed"], color=col, linestyle="--", linewidth=0.9)
        ax.plot(p["recovery_days"] * 24, p["coverage_closed_form"], color=col, linestyle=":", linewidth=0.8)
    ax.set_xscale("log")
    ax.set_xlabel(r"recovery time $\tau$ (h)")
    ax.set_ylabel("disturbed fraction of season")
    ax.set_ylim(0, 1)
    ax.set_title(r"\textbf{(a)} Coverage against recovery time", fontsize=8, loc="left")
    handles, labels = ax.get_legend_handles_labels()

    ax = axes[0, 1]
    for c in cads:
        p = s[(s["model"] == "poisson") & (s["cadence"] == c)].sort_values("recovery_days")
        over = ((p["coverage_linear"] - p["disturbed"]) / p["disturbed"]).to_numpy()
        lt = p["lambda_tau"].to_numpy()
        keep = over > 0
        ax.plot(lt[keep], over[keep], color=CADENCE_COLORS[int(c)])
    ax.axhline(0.1, color="gray", linewidth=0.7, linestyle=":")
    ax.text(0.6, 0.085, r"10\% overstatement", fontsize=8, color="gray", va="top", ha="left")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel(r"$\lambda\,E[\tau]$ (Poisson timing)")
    ax.set_ylabel(r"overstatement of coverage by $\lambda\tau$")
    ax.set_title(r"\textbf{(b)} Error of the per-event sum $\lambda\tau$", fontsize=8, loc="left")

    ax = axes[1, 0]
    wcols = [c for c in s.columns if c.startswith("has_window|")]
    ws = [float(c.split("|")[1]) for c in wcols]
    top = s[(s["model"] == "weekday_renewal") & (s["cadence"] == max(cads))].sort_values("recovery_days")
    z = np.array([[1 - top[col].to_numpy()[i] for col in wcols] for i in range(len(top))])
    im = ax.imshow(
        z,
        aspect="auto",
        origin="lower",
        cmap="magma",
        vmin=0,
        vmax=1,
        extent=[-0.5, len(ws) - 0.5, -0.5, len(top) - 0.5],
    )
    ax.set_xticks(range(len(ws)))
    ax.set_xticklabels([duration_label(w) for w in ws], fontsize=8)
    ax.set_yticks(range(len(top)))
    ax.set_yticklabels([duration_label(t) for t in top["recovery_days"]], fontsize=8)
    ax.set_xlabel(r"required undisturbed window $w$")
    ax.set_ylabel(r"recovery time $\tau$")
    ax.set_title(
        f"\\textbf{{(c)}} P(no window $\\geq w$ in season), {int(max(cads))}/yr, weekday-only", fontsize=8, loc="left"
    )
    cb = fig.colorbar(im, ax=ax, shrink=0.9, pad=0.04)
    cb.ax.tick_params(labelsize=8)

    ax = axes[1, 1]
    ax.axis("off")
    handles.append(plt.Line2D([], [], color="gray", linestyle="-"))
    labels.append("Poisson timing, all days (a)")
    handles.append(plt.Line2D([], [], color="gray", linestyle="--"))
    labels.append("weekday-only renewal (a)")
    handles.append(plt.Line2D([], [], color="gray", linestyle=":"))
    labels.append(r"closed form $1-e^{-\lambda\tau}$ (a)")
    ax.legend(handles, labels, loc="center left", fontsize=8, frameon=False, bbox_to_anchor=(0.05, 0.5))
    fig.tight_layout(w_pad=1.5, h_pad=1.5)
    save(fig, "fig5_coverage_surfaces")


if __name__ == "__main__":
    main()
