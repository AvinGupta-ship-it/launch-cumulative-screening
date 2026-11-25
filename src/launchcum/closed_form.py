"""Closed-form results for Poisson operations.

These expressions serve as test oracles for the simulation code and are
quoted in the manuscript. All rates are per day unless stated otherwise.
"""

from __future__ import annotations

import numpy as np
from numpy.typing import ArrayLike, NDArray

FloatArray = NDArray[np.float64]

SECONDS_PER_YEAR = 86400.0 * 365.0
NIGHT_PENALTY_DB = 10.0


def ldn_from_events(sel_db: ArrayLike, night: ArrayLike, seconds: float = SECONDS_PER_YEAR) -> float:
    """Day-night average level from per-event sound exposure levels.

    ``sel_db`` are per-event SEL values at one receptor; ``night`` flags the
    events receiving the 10 dB penalty. ``seconds`` is the averaging period.
    """
    sel = np.asarray(sel_db, dtype=float)
    nt = np.asarray(night, dtype=bool)
    if sel.size == 0:
        return -np.inf
    e = np.power(10.0, (sel + NIGHT_PENALTY_DB * nt) / 10.0)
    return float(10.0 * np.log10(e.sum() / seconds))


def ldn_increment(n1: float, n2: float, fn1: float, fn2: float) -> float:
    """Change in L_dn when counts go ``n1`` to ``n2`` and the night fraction
    ``fn1`` to ``fn2`` with fixed per-event levels.

    The decomposition is ``10 log10(n2/n1) + 10 log10((1+9 fn2)/(1+9 fn1))``.
    """
    return float(10.0 * np.log10(n2 / n1) + 10.0 * np.log10((1 + 9 * fn2) / (1 + 9 * fn1)))


def ldn_count_term(n1: float, n2: float) -> float:
    return float(10.0 * np.log10(n2 / n1))


def ldn_night_term(fn1: float, fn2: float) -> float:
    return float(10.0 * np.log10((1 + 9 * fn2) / (1 + 9 * fn1)))


def poisson_coverage(rate_per_day: ArrayLike, mean_recovery_days: ArrayLike) -> FloatArray:
    """Stationary disturbed fraction for a Poisson process with recovery time.

    Each event starts a disturbed interval of random length with mean
    ``mean_recovery_days``; the union of intervals is a busy period of an
    M/G/infinity system, so the fraction of time in the disturbed state is
    ``1 - exp(-lambda E[tau])``.
    """
    lam = np.asarray(rate_per_day, dtype=float)
    tau = np.asarray(mean_recovery_days, dtype=float)
    return np.asarray(1.0 - np.exp(-lam * tau), dtype=float)


def linear_coverage(rate_per_day: ArrayLike, mean_recovery_days: ArrayLike) -> FloatArray:
    """Linear (non-overlapping) extrapolation ``lambda E[tau]`` of the disturbed
    fraction, the quantity implied by summing per-event footprints."""
    return np.asarray(np.asarray(rate_per_day, float) * np.asarray(mean_recovery_days, float))


def expected_exceedance_days(n_events: float, n_days: int) -> float:
    """Expected number of days with at least one event when ``n_events`` are
    placed as a Poisson process over ``n_days`` eligible days."""
    lam_d = n_events / n_days
    return float(n_days * (1.0 - np.exp(-lam_d)))


def poisson_gap_survival(rate_per_day: float, gap_days: ArrayLike) -> FloatArray:
    """Probability that an inter-event interval exceeds ``gap_days``."""
    g = np.asarray(gap_days, dtype=float)
    return np.asarray(np.exp(-rate_per_day * g))


def prob_no_window(rate_per_day: float, season_days: float, window_days: float, recovery_days: float) -> float:
    """Approximate probability that a season contains no undisturbed window.

    An undisturbed window of length ``window_days`` requires an inter-event
    gap of at least ``window_days + recovery_days`` (fixed recovery time).
    Gaps in a Poisson process are exponential; with roughly ``lambda S``
    gaps in a season of length ``S`` the probability that none exceeds the
    required length is approximately ``exp(-lambda S exp(-lambda (w+tau)))``.
    The expression is asymptotic and is checked against simulation.
    """
    lam = rate_per_day
    need = window_days + recovery_days
    return float(np.exp(-lam * season_days * np.exp(-lam * need)))


def expected_longest_gap(rate_per_day: float, season_days: float) -> float:
    """Leading-order expectation of the longest gap among ``lambda S`` exponential
    gaps: ``(ln(lambda S) + gamma)/lambda`` with Euler's constant ``gamma``."""
    n = rate_per_day * season_days
    if n <= 1:
        return season_days
    return float((np.log(n) + np.euler_gamma) / rate_per_day)


def linear_extrapolation_error(y_ref: float, n_ref: float, y_n: ArrayLike, n: ArrayLike) -> FloatArray:
    """Relative error of the linear extrapolation ``(n/n_ref) y_ref`` against ``y_n``."""
    yn = np.asarray(y_n, dtype=float)
    nn = np.asarray(n, dtype=float)
    lin = nn / n_ref * y_ref
    with np.errstate(divide="ignore", invalid="ignore"):
        return np.asarray((lin - yn) / yn)


def energy_extrapolation_error(
    ldn_ref: float, n_ref: float, ldn_n: ArrayLike, n: ArrayLike, fn_ref: float, fn_n: ArrayLike
) -> FloatArray:
    """Error in dB of the energy-average extrapolation of L_dn from a reference
    cadence, relative to the simulated L_dn at cadence ``n``."""
    nn = np.asarray(n, dtype=float)
    fn = np.asarray(fn_n, dtype=float)
    pred = ldn_ref + 10 * np.log10(nn / n_ref) + 10 * np.log10((1 + 9 * fn) / (1 + 9 * fn_ref))
    return np.asarray(pred - np.asarray(ldn_n, dtype=float))


def cap_probability_normal(mean_hours: float, sd_hours: float, cap_hours: float) -> float:
    """Normal approximation to ``P(H > cap)`` for annual closure hours ``H``."""
    from scipy.stats import norm

    return float(norm.sf(cap_hours, loc=mean_hours, scale=sd_hours))
