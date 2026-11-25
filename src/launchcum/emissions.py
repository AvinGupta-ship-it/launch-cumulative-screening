"""Fleet stratospheric emission inventory (bounded sub-analysis).

Annual emissions above the tropopause from methane-oxygen combustion are
estimated as

    E_s = N * M_prop * f_strat * EI_s

for each species ``s``, where ``N`` is the annual launch count, ``M_prop``
the propellant mass per launch, ``f_strat`` the fraction of propellant
consumed above the tropopause and ``EI_s`` the emission index. Ranges are
propagated by Monte Carlo sampling of uniform distributions over the
documented bounds. No atmospheric chemistry or transport is modelled.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd
from numpy.typing import NDArray

FloatArray = NDArray[np.float64]


@dataclass(frozen=True)
class Range:
    low: float
    high: float
    source: str

    def sample(self, rng: np.random.Generator, n: int) -> FloatArray:
        return np.asarray(rng.uniform(self.low, self.high, size=n))

    @property
    def mid(self) -> float:
        return 0.5 * (self.low + self.high)


@dataclass(frozen=True)
class SiteCadence:
    name: str
    launches_per_year: float
    source: str


@dataclass(frozen=True)
class EmissionInputs:
    propellant_per_launch_t: Range
    f_stratosphere: Range
    ei_g_per_kg: dict[str, Range]


def inventory(
    sites: list[SiteCadence],
    inputs: EmissionInputs,
    rng: np.random.Generator,
    n_samples: int = 20000,
) -> pd.DataFrame:
    """Monte Carlo inventory per site and for the fleet total.

    Returns a table of species by site with median and 5-95 percentile
    annual emissions in tonnes above the tropopause.
    """
    m = inputs.propellant_per_launch_t.sample(rng, n_samples)
    f = inputs.f_stratosphere.sample(rng, n_samples)
    rows = []
    total: dict[str, FloatArray] = {s: np.zeros(n_samples) for s in inputs.ei_g_per_kg}
    for site in sites:
        for species, ei in inputs.ei_g_per_kg.items():
            e = site.launches_per_year * m * f * ei.sample(rng, n_samples) * 1e-3  # t
            total[species] += e
            rows.append(_row(site.name, species, site.launches_per_year, e))
    n_fleet = sum(s.launches_per_year for s in sites)
    for species, e in total.items():
        rows.append(_row("fleet_total", species, n_fleet, e))
    return pd.DataFrame(rows)


def _row(site: str, species: str, n: float, e: FloatArray) -> dict[str, float | str]:
    q = np.percentile(e, [5, 50, 95])
    return {
        "site": site,
        "species": species,
        "launches_per_year": n,
        "p05_t": q[0],
        "median_t": q[1],
        "p95_t": q[2],
        "mean_t": float(e.mean()),
    }
