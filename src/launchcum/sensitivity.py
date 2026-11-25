"""Global sensitivity analysis wrappers around SALib."""

from __future__ import annotations

from collections.abc import Callable
from dataclasses import dataclass

import numpy as np
import pandas as pd
from numpy.typing import NDArray
from SALib.analyze import morris as morris_analyze
from SALib.analyze import sobol as sobol_analyze
from SALib.sample import morris as morris_sample
from SALib.sample import sobol as sobol_sample

FloatArray = NDArray[np.float64]


@dataclass(frozen=True)
class Factor:
    name: str
    low: float
    high: float
    source_rank: int
    note: str = ""


def problem_from_factors(factors: list[Factor]) -> dict[str, object]:
    return {
        "num_vars": len(factors),
        "names": [f.name for f in factors],
        "bounds": [[f.low, f.high] for f in factors],
    }


def saltelli_design(factors: list[Factor], n_base: int, seed: int) -> FloatArray:
    """Saltelli sample of size ``n_base (k + 2)`` with first/total indices only."""
    problem = problem_from_factors(factors)
    return np.asarray(sobol_sample.sample(problem, n_base, calc_second_order=False, seed=seed))


def sobol_indices(factors: list[Factor], y: FloatArray, n_boot: int = 200, seed: int = 0) -> pd.DataFrame:
    """First- and total-order Sobol indices with bootstrap confidence intervals."""
    problem = problem_from_factors(factors)
    y = np.asarray(y, dtype=float)
    if not np.all(np.isfinite(y)):
        y = np.where(np.isfinite(y), y, np.nanmedian(y[np.isfinite(y)]) if np.isfinite(y).any() else 0.0)
    if np.std(y) == 0:
        return pd.DataFrame({"factor": problem["names"], "S1": 0.0, "S1_conf": 0.0, "ST": 0.0, "ST_conf": 0.0})
    si = sobol_analyze.analyze(problem, y, calc_second_order=False, num_resamples=n_boot, conf_level=0.95, seed=seed)
    return pd.DataFrame(
        {"factor": problem["names"], "S1": si["S1"], "S1_conf": si["S1_conf"], "ST": si["ST"], "ST_conf": si["ST_conf"]}
    )


def morris_design(factors: list[Factor], n_traj: int, seed: int) -> FloatArray:
    problem = problem_from_factors(factors)
    return np.asarray(morris_sample.sample(problem, n_traj, num_levels=4, seed=seed))


def morris_effects(factors: list[Factor], x: FloatArray, y: FloatArray) -> pd.DataFrame:
    problem = problem_from_factors(factors)
    si = morris_analyze.analyze(problem, x, np.asarray(y, float), conf_level=0.95, num_resamples=200)
    return pd.DataFrame(
        {"factor": problem["names"], "mu_star": si["mu_star"], "sigma": si["sigma"], "mu_star_conf": si["mu_star_conf"]}
    )


def ishigami(x: FloatArray, a: float = 7.0, b: float = 0.1) -> FloatArray:
    """Ishigami test function on ``[-pi, pi]^3``."""
    return np.asarray(np.sin(x[:, 0]) + a * np.sin(x[:, 1]) ** 2 + b * x[:, 2] ** 4 * np.sin(x[:, 0]))


def ishigami_analytic(a: float = 7.0, b: float = 0.1) -> dict[str, FloatArray]:
    """Analytical first- and total-order indices of the Ishigami function."""
    pi = np.pi
    v1 = 0.5 * (1 + b * pi**4 / 5) ** 2
    v2 = a**2 / 8
    v13 = b**2 * pi**8 * 8 / 225
    v = v1 + v2 + v13
    s1 = np.asarray([v1 / v, v2 / v, 0.0])
    st = np.asarray([(v1 + v13) / v, v2 / v, v13 / v])
    return {"S1": s1, "ST": st}


def evaluate_rows(fn: Callable[[FloatArray], float], x: FloatArray) -> FloatArray:
    return np.asarray([fn(row) for row in x], dtype=float)
