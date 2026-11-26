"""Comment Figure 2: per-launch closure budget implied by annual caps against the recorded per-launch hours."""

from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np

from common import OI, save

RECORD = (
    (24.5, "Flight 1 (24.5 h)", OI["vermillion"]),
    (10.5, "Flight 2 (10.5 h)", OI["orange"]),
    (8.75, "Flight 3 (8.75 h)", OI["purple"]),
    (7.75, "Flight 4 (7.75 h)", OI["green"]),
)
CAPS = ((500, ":"), (800, "--"), (2000, "-"))


def main() -> None:
    cad = np.linspace(25, 456, 200)
    fig, ax = plt.subplots(figsize=(4.8, 3.3))
    for cap, ls in CAPS:
        ax.plot(cad, cap / cad, color=OI["blue"], linestyle=ls, label=f"{cap:,} h per year cap divided by cadence")
    for hrs, lab, col in RECORD:
        ax.axhline(hrs, color=col, linewidth=0.9, label=lab)
    ax.set_xlabel("launches per year")
    ax.set_ylabel("access-restriction hours per launch")
    ax.set_xticks([25, 100, 200, 456])
    ax.set_ylim(0, 30)
    # The legend sits between the Flight 1 line (24.5 h) and the Flight 2 line (10.5 h) at high cadence,
    # where no curve passes.
    ax.legend(loc="center right", fontsize=8, ncol=1, bbox_to_anchor=(0.99, 0.58))
    fig.tight_layout()
    save(fig, "figc2_closure_budget")


if __name__ == "__main__":
    main()
