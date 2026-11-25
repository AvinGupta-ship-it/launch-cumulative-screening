"""Stochastic operations scheduler.

Operations are generated as a marked point process on the calendar of one
scenario year. Three temporal models are implemented:

``poisson``
    Homogeneous Poisson process (annual count Poisson-distributed about the
    configured mean). The reference case for the closed forms.
``weekday_renewal``
    Launches on business days only; same-pad intervals follow a shifted
    gamma distribution with a minimum turnaround; exactly the configured
    number of launches per year.
``campaign``
    Neyman-Scott cluster process matched to the same mean annual count.

Each launch spawns a Super Heavy landing shortly after lift-off and, with a
configurable probability, a Starship landing hours to days later. Landings
of vehicles launched elsewhere arrive as an independent process. Scrubs and
static fires are attached to launches. Every event carries a closure window.

The generator works on plain NumPy arrays for speed; :func:`events_to_frame`
converts a realization to a typed ``pandas.DataFrame``.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Literal

import numpy as np
import pandas as pd
from numpy.typing import NDArray

from .opcalendar import (
    OpCalendar,
    build_calendar,
    business_time_to_calendar,
    draw_time_of_day,
    is_night,
)

FloatArray = NDArray[np.float64]
IntArray = NDArray[np.int64]
BoolArray = NDArray[np.bool_]

EVENT_TYPES = ("launch", "sh_landing", "ss_landing", "static_fire", "scrub")
TYPE_CODE = {name: i for i, name in enumerate(EVENT_TYPES)}
TemporalModel = Literal["poisson", "weekday_renewal", "campaign"]


@dataclass(frozen=True)
class ClosureModel:
    """Per-event closure window parameters (hours).

    Durations are lognormal with the given mean and coefficient of
    variation. ``lead_fraction`` is the share of the window preceding the
    event. ``anomaly_probability`` is the per-launch probability of an
    anomaly adding ``anomaly_hours`` of access restriction.
    """

    launch_mean_h: float = 12.9
    launch_cv: float = 0.55
    ss_landing_mean_h: float = 3.0
    ss_landing_cv: float = 0.5
    static_fire_mean_h: float = 4.0
    static_fire_cv: float = 0.5
    scrub_fraction: float = 0.7
    lead_fraction: float = 0.75
    anomaly_probability: float = 0.02
    anomaly_hours: float = 48.0
    external_landing_mean_h: float = 3.0


@dataclass(frozen=True)
class ScheduleParams:
    """Parameters of one operations scenario."""

    n_launches: float
    model: TemporalModel = "poisson"
    year: int = 2030
    weekday_only: bool = False
    night_fraction: float = 0.5
    n_pads: int = 1
    turnaround_min_days: float = 0.5
    gamma_shape: float = 2.0
    renewal_slack: float = 0.1
    cluster_rate_per_year: float = 60.0
    cluster_spread_days: float = 2.0
    p_sh_landing_local: float = 1.0
    p_ss_landing_local: float = 1.0
    sh_landing_delay_min: float = 7.0
    ss_landing_delay_days_min: float = 0.04
    ss_landing_delay_days_max: float = 3.0
    n_external_landings: float = 0.0
    external_weekday_only: bool = False
    static_fires_per_launch: float = 0.2
    static_fire_lead_days: tuple[float, float] = (1.0, 14.0)
    scrub_probability: float = 0.15
    scrub_recycle_days: tuple[float, float] = (1.0, 3.0)
    closures: ClosureModel = field(default_factory=ClosureModel)


@dataclass
class Realization:
    """One annual realization of the operations process (struct of arrays)."""

    t: FloatArray
    type_code: IntArray
    night: BoolArray
    external: BoolArray
    pad: IntArray
    closure_start: FloatArray
    closure_hours: FloatArray
    anomaly_hours: float
    feasible: bool
    renewal_rejections: int

    def select(self, *types: str) -> BoolArray:
        codes = [TYPE_CODE[t] for t in types]
        return np.isin(self.type_code, codes)

    @property
    def n_events(self) -> int:
        return int(self.t.shape[0])

    def count(self, event_type: str) -> int:
        return int((self.type_code == TYPE_CODE[event_type]).sum())


def _lognormal(rng: np.random.Generator, mean: float, cv: float, size: int) -> FloatArray:
    if size == 0:
        return np.zeros(0)
    if mean <= 0:
        return np.zeros(size)
    sigma2 = np.log(1.0 + cv**2)
    mu = np.log(mean) - 0.5 * sigma2
    return np.asarray(rng.lognormal(mu, np.sqrt(sigma2), size))


def _launch_days_poisson(rng: np.random.Generator, p: ScheduleParams, cal: OpCalendar) -> IntArray:
    n = int(rng.poisson(p.n_launches))
    if p.weekday_only:
        b = rng.random(n) * cal.n_business
        return business_time_to_calendar(cal, b)
    return np.asarray(np.floor(rng.random(n) * cal.n_days), dtype=np.int64)


def _launch_days_renewal(
    rng: np.random.Generator, p: ScheduleParams, cal: OpCalendar
) -> tuple[IntArray, IntArray, bool, int]:
    """Exactly ``n_launches`` launches on business days with per-pad turnaround."""
    n = int(round(p.n_launches))
    axis_len = float(cal.n_business if p.weekday_only else cal.n_days)
    pads = rng.integers(0, p.n_pads, size=n)
    days = np.empty(n, dtype=np.int64)
    rejections = 0
    feasible = True
    for pad in range(p.n_pads):
        idx = np.flatnonzero(pads == pad)
        m = idx.shape[0]
        if m == 0:
            continue
        mean_interval = axis_len * (1.0 - p.renewal_slack) / m
        excess = mean_interval - p.turnaround_min_days
        if excess <= 0:
            # The cadence cannot be met with this turnaround: fall back to an
            # evenly spaced schedule at the minimum turnaround and flag it.
            feasible = False
            pos = p.turnaround_min_days * np.arange(m) + rng.random() * max(
                axis_len - p.turnaround_min_days * (m - 1), 0.0
            )
            pos = np.minimum(pos, axis_len - 1e-9)
        else:
            theta = excess / p.gamma_shape
            pos = None
            for _ in range(50):
                gaps = p.turnaround_min_days + rng.gamma(p.gamma_shape, theta, size=m)
                total = gaps[1:].sum()
                if total <= axis_len:
                    start = rng.random() * (axis_len - total)
                    pos = start + np.concatenate([[0.0], np.cumsum(gaps[1:])])
                    break
                rejections += 1
            if pos is None:
                gaps = p.turnaround_min_days + rng.gamma(p.gamma_shape, theta, size=m)
                gaps = gaps[1:] * (axis_len * 0.98 / gaps[1:].sum())
                pos = np.concatenate([[0.0], np.cumsum(gaps)]) + 0.01 * axis_len
                feasible = False
        if p.weekday_only:
            days[idx] = business_time_to_calendar(cal, np.asarray(pos))
        else:
            days[idx] = np.clip(np.floor(pos).astype(np.int64), 0, cal.n_days - 1)
    return days, pads.astype(np.int64), feasible, rejections


def _launch_days_campaign(rng: np.random.Generator, p: ScheduleParams, cal: OpCalendar) -> IntArray:
    axis_len = float(cal.n_business if p.weekday_only else cal.n_days)
    n_clusters = int(rng.poisson(p.cluster_rate_per_year))
    if n_clusters == 0:
        return np.zeros(0, dtype=np.int64)
    mean_size = p.n_launches / p.cluster_rate_per_year
    centres = rng.random(n_clusters) * axis_len
    sizes = rng.poisson(mean_size, size=n_clusters)
    parents = np.repeat(centres, sizes)
    pos = parents + rng.normal(0.0, p.cluster_spread_days, size=parents.shape[0])
    pos = np.mod(pos, axis_len)
    if p.weekday_only:
        return business_time_to_calendar(cal, pos)
    return np.asarray(np.floor(pos).astype(np.int64))


def generate(params: ScheduleParams, rng: np.random.Generator, cal: OpCalendar | None = None) -> Realization:
    """Generate one annual realization of operations."""
    p = params
    cal = cal or build_calendar(p.year)
    rejections = 0
    feasible = True
    launch_day: IntArray
    pads: IntArray
    if p.model == "poisson":
        launch_day = _launch_days_poisson(rng, p, cal)
        pads = np.asarray(rng.integers(0, p.n_pads, size=launch_day.shape[0]), dtype=np.int64)
    elif p.model == "weekday_renewal":
        launch_day, pads, feasible, rejections = _launch_days_renewal(rng, p, cal)
    elif p.model == "campaign":
        launch_day = _launch_days_campaign(rng, p, cal)
        pads = np.asarray(rng.integers(0, p.n_pads, size=launch_day.shape[0]), dtype=np.int64)
    else:  # pragma: no cover - guarded by pydantic
        raise ValueError(f"unknown temporal model {p.model!r}")

    n_l = launch_day.shape[0]
    launch_night = np.asarray(rng.random(n_l) < p.night_fraction, dtype=bool)
    t_launch = launch_day + draw_time_of_day(rng, launch_night) / 24.0
    t_launch.sort()

    # Super Heavy landings a few minutes after lift-off.
    sh_mask = np.asarray(rng.random(n_l) < p.p_sh_landing_local, dtype=bool)
    t_sh = t_launch[sh_mask] + p.sh_landing_delay_min / 1440.0
    # Starship landings hours to days later with their own day/night draw.
    ss_mask = np.asarray(rng.random(n_l) < p.p_ss_landing_local, dtype=bool)
    n_ss = int(ss_mask.sum())
    delay = rng.uniform(p.ss_landing_delay_days_min, p.ss_landing_delay_days_max, size=n_ss)
    ss_day = np.floor(t_launch[ss_mask] + delay)
    ss_night = np.asarray(rng.random(n_ss) < p.night_fraction, dtype=bool)
    t_ss = ss_day + draw_time_of_day(rng, ss_night) / 24.0
    # Landings of vehicles launched from other sites.
    n_ext = int(rng.poisson(p.n_external_landings)) if p.n_external_landings > 0 else 0
    if p.external_weekday_only:
        ext_day = business_time_to_calendar(cal, rng.random(n_ext) * cal.n_business)
    else:
        ext_day = np.floor(rng.random(n_ext) * cal.n_days)
    ext_night = np.asarray(rng.random(n_ext) < p.night_fraction, dtype=bool)
    t_ext = ext_day + draw_time_of_day(rng, ext_night) / 24.0
    # Static fires before launch (daytime, per the current licence terms).
    n_sf_per = rng.poisson(p.static_fires_per_launch, size=n_l)
    sf_parent = np.repeat(np.arange(n_l), n_sf_per)
    lo, hi = p.static_fire_lead_days
    sf_day = np.floor(t_launch[sf_parent] - rng.uniform(lo, hi, size=sf_parent.shape[0]))
    t_sf = sf_day + draw_time_of_day(rng, np.zeros(sf_parent.shape[0], dtype=bool)) / 24.0
    # Scrubs: geometric number of scrubbed attempts before each launch.
    n_scrub_per = (
        rng.geometric(1.0 - p.scrub_probability, size=n_l) - 1
        if p.scrub_probability > 0
        else np.zeros(n_l, dtype=np.int64)
    )
    n_scrub_per = np.asarray(n_scrub_per, dtype=np.int64)
    sc_parent = np.repeat(np.arange(n_l), n_scrub_per)
    rlo, rhi = p.scrub_recycle_days
    attempt_index = np.concatenate([np.arange(1, k + 1) for k in n_scrub_per]) if sc_parent.shape[0] else np.zeros(0)
    t_sc = t_launch[sc_parent] - attempt_index * rng.uniform(rlo, rhi, size=sc_parent.shape[0])

    t_all = np.concatenate([t_launch, t_sh, t_ss, t_ext, t_sf, t_sc])
    codes = np.concatenate(
        [
            np.full(n_l, TYPE_CODE["launch"]),
            np.full(t_sh.shape[0], TYPE_CODE["sh_landing"]),
            np.full(n_ss, TYPE_CODE["ss_landing"]),
            np.full(n_ext, TYPE_CODE["ss_landing"]),
            np.full(t_sf.shape[0], TYPE_CODE["static_fire"]),
            np.full(t_sc.shape[0], TYPE_CODE["scrub"]),
        ]
    ).astype(np.int64)
    external = np.concatenate(
        [
            np.zeros(n_l + t_sh.shape[0] + n_ss, dtype=bool),
            np.ones(n_ext, dtype=bool),
            np.zeros(t_sf.shape[0] + t_sc.shape[0], dtype=bool),
        ]
    )
    pad_all = np.concatenate(
        [
            pads,
            pads[sh_mask],
            pads[ss_mask],
            np.zeros(n_ext, dtype=np.int64),
            pads[sf_parent],
            pads[sc_parent],
        ]
    ).astype(np.int64)

    # Keep events inside the scenario year (static fires and scrubs attached
    # to early-January launches may fall into the previous year).
    keep = (t_all >= 0.0) & (t_all < cal.n_days)
    t_all, codes, external, pad_all = t_all[keep], codes[keep], external[keep], pad_all[keep]
    order = np.argsort(t_all, kind="stable")
    t_all, codes, external, pad_all = t_all[order], codes[order], external[order], pad_all[order]
    night = is_night(t_all)

    # Closure windows. The launch window covers the launch and the Super
    # Heavy landing that follows it, so Super Heavy landings carry no window.
    c = p.closures
    hours = np.zeros(t_all.shape[0])
    for name, mean, cv in (
        ("launch", c.launch_mean_h, c.launch_cv),
        ("ss_landing", c.ss_landing_mean_h, c.ss_landing_cv),
        ("static_fire", c.static_fire_mean_h, c.static_fire_cv),
        ("scrub", c.launch_mean_h * c.scrub_fraction, c.launch_cv),
    ):
        m = codes == TYPE_CODE[name]
        hours[m] = _lognormal(rng, mean, cv, int(m.sum()))
    ext_m = (codes == TYPE_CODE["ss_landing"]) & external
    hours[ext_m] = _lognormal(rng, c.external_landing_mean_h, c.ss_landing_cv, int(ext_m.sum()))
    start = t_all - c.lead_fraction * hours / 24.0
    n_anom = int(rng.binomial(n_l, c.anomaly_probability)) if n_l else 0
    anomaly_hours = float(n_anom * c.anomaly_hours)

    return Realization(
        t=t_all,
        type_code=codes,
        night=night,
        external=external,
        pad=pad_all,
        closure_start=start,
        closure_hours=hours,
        anomaly_hours=anomaly_hours,
        feasible=feasible,
        renewal_rejections=rejections,
    )


def generate_many(params: ScheduleParams, seed: int | np.random.SeedSequence, n: int) -> list[Realization]:
    """Generate ``n`` independent realizations from spawned seed sequences."""
    ss = seed if isinstance(seed, np.random.SeedSequence) else np.random.SeedSequence(seed)
    cal = build_calendar(params.year)
    return [generate(params, np.random.default_rng(child), cal) for child in ss.spawn(n)]


def events_to_frame(real: Realization, realization_id: int = 0) -> pd.DataFrame:
    """Convert a realization to a typed ``DataFrame``."""
    return pd.DataFrame(
        {
            "realization": np.full(real.n_events, realization_id, dtype=np.int64),
            "t_days": real.t,
            "event_type": pd.Categorical([EVENT_TYPES[int(c)] for c in real.type_code], categories=list(EVENT_TYPES)),
            "night": real.night,
            "external_origin": real.external,
            "pad": real.pad,
            "closure_start_days": real.closure_start,
            "closure_hours": real.closure_hours,
        }
    )
