"""Figure 2: schematic of the four aggregation classes (closed forms, illustrative parameters)."""

from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np

from common import OI, save
from launchcum import closed_form as cf


def main() -> None:
    n = np.linspace(1, 500, 400)
    n0 = 25.0
    fig, axes = plt.subplots(1, 4, figsize=(7.2, 2.1))

    ax = axes[0]
    ax.plot(n, n / n0, color=OI["blue"])
    ax.set_title("(i) Additive", fontsize=6.5)
    ax.set_ylabel("metric / metric at 25")

    ax = axes[1]
    ax.plot(n, 10 * np.log10(n / n0), color=OI["blue"], label="count term")
    ax.plot(
        n,
        10 * np.log10(n / n0) + cf.ldn_night_term(0.12, 0.5),
        color=OI["vermillion"],
        label="with night share 0.12 to 0.5",
    )
    ax.set_title("(ii) Energy-averaged ($L_{dn}$, dB)", fontsize=6.5)
    ax.legend(fontsize=5.5, loc="lower right")

    ax = axes[2]
    for tau_h, c in ((0.5, OI["green"]), (6.0, OI["blue"]), (24.0, OI["vermillion"])):
        lam = n / 365.0
        ax.plot(n, cf.poisson_coverage(lam, tau_h / 24.0), color=c, label=f"recovery {tau_h:g} h")
        ax.plot(n, np.minimum(cf.linear_coverage(lam, tau_h / 24.0), 1.0), color=c, linestyle=":", linewidth=0.8)
    ax.set_ylim(0, 1)
    ax.set_title("(iii) Coverage", fontsize=6.5)
    ax.legend(fontsize=5.5, loc="upper left")

    ax = axes[3]
    ax.plot(n, [cf.expected_exceedance_days(v, 365) for v in n], color=OI["blue"], label="any day")
    ax.plot(n, [cf.expected_exceedance_days(v, 250) for v in n], color=OI["vermillion"], label="business days only")
    ax.plot(n, n, color="gray", linestyle=":", linewidth=0.8, label="linear")
    ax.set_ylim(0, 400)
    ax.set_title("(iv) Exceedance days", fontsize=6.5)
    ax.legend(fontsize=5.5, loc="upper left")

    for ax in axes:
        ax.set_xlabel("events per year")
    fig.suptitle(
        "Scaling of four aggregation classes under Poisson operations (schematic; dotted = linear extrapolation)",
        fontsize=7.5,
        y=1.06,
    )
    fig.subplots_adjust(wspace=0.35)
    save(fig, "fig2_aggregation_classes")


if __name__ == "__main__":
    main()
