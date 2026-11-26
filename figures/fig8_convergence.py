"""Supporting-information figure: Monte Carlo convergence of the headline quantities."""

from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np

from common import OI, read_csv, read_parquet, save


def main() -> None:
    se = read_csv("monte_carlo_standard_errors.csv")
    raw = read_parquet("scenarios_raw.parquet")
    n_pts = np.unique(np.geomspace(50, raw.groupby("scenario").size().max(), 40).astype(int))
    fig, ax = plt.subplots(figsize=(5.0, 3.4))
    worst = None
    for _, r in se.iterrows():
        v = raw[raw["scenario"] == r["scenario"]][r["metric"]].to_numpy(float)
        if v.std(ddof=1) == 0.0:
            continue  # a deterministic count (every realization identical) has no Monte Carlo error
        ns = n_pts[n_pts <= v.size]
        frac = np.array([v[:n].std(ddof=1) / np.sqrt(n) / abs(v[:n].mean()) for n in ns])
        ax.plot(ns, 100 * frac, color=OI["blue"], alpha=0.35, linewidth=0.8)
        if worst is None or frac[-1] > worst[1][-1]:
            worst = (ns, frac, f"{r['scenario']}: {r['metric']}")
    if worst is not None:
        ax.plot(
            worst[0],
            100 * worst[1],
            color=OI["vermillion"],
            linewidth=1.3,
            label=f"largest at the full sample ({100 * worst[1][-1]:.2f}\\%)",
        )
    ax.axhline(2.0, color="gray", linestyle="--", linewidth=0.8, label=r"2\% criterion")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("number of simulated years")
    ax.set_ylabel(r"standard error of the mean (\% of the mean)")
    ax.set_ylim(1e-3, 20)
    ax.legend(fontsize=8, loc="upper right")
    fig.tight_layout()
    save(fig, "fig8_convergence")


if __name__ == "__main__":
    main()
