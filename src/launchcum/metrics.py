"""Metric algebra applied to realizations of the operations process.

Four classes of cumulative-effect metric are implemented:

(i)   additive totals (water, emissions);
(ii)  energy-averaged noise (``L_dn`` from SEL, ``CDNL`` from boom overpressure);
(iii) number-above, exceedance-day and boom-count metrics;
(iv)  coverage and interval metrics for disturbance-sensitive receptors;

plus public-access availability metrics derived from closure windows.
All functions operate on one :class:`~launchcum.schedule.Realization`
and a :class:`ReceptorLevels` table of per-event levels at receptors.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray

from .closed_form import NIGHT_PENALTY_DB, SECONDS_PER_YEAR
from .footprints import FootprintLibrary, csel_from_psf
from .opcalendar import NIGHT_END_H, NIGHT_START_H, OpCalendar
from .schedule import EVENT_TYPES, TYPE_CODE, Realization

FloatArray = NDArray[np.float64]
IntArray = NDArray[np.int64]
BoolArray = NDArray[np.bool_]

NOISE_TYPES = ("launch", "sh_landing", "ss_landing", "static_fire")
BOOM_TYPES = ("sh_landing", "ss_landing")


@dataclass
class ReceptorLevels:
    """Per-event-type levels at a set of receptors (rows: event types)."""

    names: list[str]
    r_km: FloatArray
    azimuth_deg: FloatArray
    lamax: FloatArray
    sel: FloatArray
    psf: FloatArray

    @property
    def csel(self) -> FloatArray:
        out = np.full_like(self.psf, -np.inf)
        m = self.psf > 0
        out[m] = csel_from_psf(self.psf[m])
        return out

    @property
    def n_receptors(self) -> int:
        return len(self.names)


def receptor_levels(
    lib: FootprintLibrary,
    names: list[str],
    r_km: FloatArray,
    azimuth_deg: FloatArray,
    level_offset_db: float = 0.0,
    boom_scale: float = 1.0,
) -> ReceptorLevels:
    """Evaluate fitted footprints at receptors for every event type."""
    n_t = len(EVENT_TYPES)
    k = len(names)
    lamax = np.full((n_t, k), -np.inf)
    sel = np.full((n_t, k), -np.inf)
    psf = np.zeros((n_t, k))
    for et in NOISE_TYPES:
        i = TYPE_CODE[et]
        lamax[i] = lib.level(et, "lamax", r_km, level_offset_db)
        sel[i] = lib.level(et, "sel", r_km, level_offset_db)
    for et in BOOM_TYPES:
        psf[TYPE_CODE[et]] = lib.boom(et, r_km, azimuth_deg, boom_scale)
    return ReceptorLevels(list(names), np.asarray(r_km, float), np.asarray(azimuth_deg, float), lamax, sel, psf)


# --------------------------------------------------------------------------
# (i) additive totals
# --------------------------------------------------------------------------


def additive_totals(real: Realization, per_operation: dict[str, dict[str, float]]) -> dict[str, float]:
    """Sum per-operation quantities over events.

    ``per_operation`` maps quantity name to a mapping of event type to the
    per-event value, e.g. ``{"deluge_gal": {"launch": 422000, ...}}``.
    """
    out: dict[str, float] = {}
    for qname, by_type in per_operation.items():
        total = 0.0
        for et, q in by_type.items():
            total += q * real.count(et)
        out[qname] = total
    return out


# --------------------------------------------------------------------------
# (ii) energy-averaged noise
# --------------------------------------------------------------------------


def ldn(real: Realization, lv: ReceptorLevels, seconds: float = SECONDS_PER_YEAR) -> FloatArray:
    """Day-night average level at each receptor from per-event SEL."""
    sel = lv.sel[real.type_code]  # (n_events, k)
    pen = NIGHT_PENALTY_DB * real.night[:, None]
    with np.errstate(over="ignore"):
        e = np.power(10.0, (sel + pen) / 10.0)
    e[~np.isfinite(sel)] = 0.0
    tot = e.sum(axis=0)
    with np.errstate(divide="ignore"):
        return np.asarray(10.0 * np.log10(tot / seconds))


def cdnl(real: Realization, lv: ReceptorLevels, seconds: float = SECONDS_PER_YEAR) -> FloatArray:
    """C-weighted day-night level from sonic booms (FAA Desk Reference form).

    ``CDNL = 10 log10( sum_i 10^((CSEL_i + 10*night_i)/10) / T )`` with ``T``
    the number of seconds in the year; for ``N`` identical daytime booms this
    reduces to ``CSEL + 10 log10(N/365) - 49.4``.
    """
    csel = lv.csel[real.type_code]
    pen = NIGHT_PENALTY_DB * real.night[:, None]
    with np.errstate(over="ignore"):
        e = np.power(10.0, (csel + pen) / 10.0)
    e[~np.isfinite(csel)] = 0.0
    tot = e.sum(axis=0)
    with np.errstate(divide="ignore"):
        return np.asarray(10.0 * np.log10(tot / seconds))


def ldn_from_counts(
    sel_by_type: FloatArray, n_day: FloatArray, n_night: FloatArray, seconds: float = SECONDS_PER_YEAR
) -> FloatArray:
    """L_dn field from expected day/night counts per type and SEL fields.

    ``sel_by_type`` has shape ``(n_types, ...)``; ``n_day`` and ``n_night``
    have shape ``(n_types,)``. Non-finite SEL values contribute nothing.
    """
    with np.errstate(over="ignore"):
        e = np.power(10.0, sel_by_type / 10.0)
    e[~np.isfinite(sel_by_type)] = 0.0
    w = (n_day + 10.0 * n_night).reshape((-1,) + (1,) * (sel_by_type.ndim - 1))
    tot = (e * w).sum(axis=0)
    with np.errstate(divide="ignore"):
        return np.asarray(10.0 * np.log10(tot / seconds))


# --------------------------------------------------------------------------
# (iii) number-above, exceedance days, boom counts
# --------------------------------------------------------------------------


def number_above(real: Realization, lv: ReceptorLevels, threshold_db: float, night_only: bool = False) -> IntArray:
    lam = lv.lamax[real.type_code]
    above = lam > threshold_db
    if night_only:
        above &= real.night[:, None]
    return np.asarray(above.sum(axis=0), dtype=np.int64)


def exceedance_days(real: Realization, lv: ReceptorLevels, threshold_db: float) -> IntArray:
    """Number of calendar days with at least one event above ``threshold_db``."""
    lam = lv.lamax[real.type_code]
    above = lam > threshold_db  # (n_events, k)
    day = np.floor(real.t).astype(np.int64)
    out = np.zeros(lv.n_receptors, dtype=np.int64)
    for k in range(lv.n_receptors):
        out[k] = np.unique(day[above[:, k]]).shape[0]
    return out


def boom_counts(real: Realization, lv: ReceptorLevels, psf_threshold: float) -> IntArray:
    p = lv.psf[real.type_code]
    return np.asarray((p > psf_threshold).sum(axis=0), dtype=np.int64)


# --------------------------------------------------------------------------
# (iv) coverage and interval metrics
# --------------------------------------------------------------------------


def merge_intervals(start: FloatArray, end: FloatArray) -> tuple[FloatArray, FloatArray]:
    """Union of intervals as sorted disjoint ``(start, end)`` arrays."""
    if start.shape[0] == 0:
        return np.zeros(0), np.zeros(0)
    order = np.argsort(start, kind="stable")
    s = start[order]
    e = end[order]
    cummax = np.maximum.accumulate(e)
    new_group = np.concatenate([[True], s[1:] > cummax[:-1]])
    gid = np.cumsum(new_group) - 1
    ms = s[new_group]
    me = np.zeros(ms.shape[0])
    np.maximum.at(me, gid, e)
    return ms, me


def union_length(start: FloatArray, end: FloatArray, lo: float, hi: float) -> float:
    """Length of the union of intervals clipped to ``[lo, hi]``."""
    ms, me = merge_intervals(np.clip(start, lo, hi), np.clip(end, lo, hi))
    return float(np.sum(me - ms))


def gaps_within(start: FloatArray, end: FloatArray, lo: float, hi: float) -> FloatArray:
    """Lengths of undisturbed gaps inside ``[lo, hi]`` (including the ends)."""
    ms, me = merge_intervals(np.clip(start, lo, hi), np.clip(end, lo, hi))
    if ms.shape[0] == 0:
        return np.asarray([hi - lo])
    edges_s = np.concatenate([[lo], me])
    edges_e = np.concatenate([ms, [hi]])
    g = edges_e - edges_s
    return np.asarray(g[g > 0])


@dataclass
class CoverageResult:
    disturbed_fraction: float
    longest_gap_days: float
    n_gaps: int
    has_window: dict[float, bool]
    n_windows: dict[float, int]


def coverage(
    t_events: FloatArray,
    recovery_days: FloatArray,
    season: tuple[float, float],
    windows_days: tuple[float, ...] = (),
) -> CoverageResult:
    """Coverage statistics of the union of ``[t_i, t_i + tau_i]`` within a season."""
    lo, hi = season
    start = t_events
    end = t_events + recovery_days
    inside = (end > lo) & (start < hi)
    start, end = start[inside], end[inside]
    length = hi - lo
    covered = union_length(start, end, lo, hi) if start.shape[0] else 0.0
    g = gaps_within(start, end, lo, hi)
    longest = float(np.max(g)) if g.size else 0.0
    return CoverageResult(
        disturbed_fraction=covered / length,
        longest_gap_days=longest,
        n_gaps=int(g.size),
        has_window={w: bool(longest >= w) for w in windows_days},
        n_windows={w: int(np.sum(np.floor(g / w))) for w in windows_days},
    )


def disturbance_times(
    real: Realization,
    lv: ReceptorLevels,
    receptor: int,
    threshold_db: float | None,
    psf_threshold: float | None = None,
    include_scrubs: bool = False,
) -> FloatArray:
    """Event times that disturb ``receptor``.

    An event disturbs the receptor if its LAmax exceeds ``threshold_db`` or
    its boom exceeds ``psf_threshold``. With both thresholds ``None`` every
    noise-producing event counts.
    """
    types = real.type_code
    lam = lv.lamax[types, receptor]
    psf = lv.psf[types, receptor]
    if threshold_db is None and psf_threshold is None:
        m = np.isin(types, [TYPE_CODE[t] for t in NOISE_TYPES])
    else:
        m = np.zeros(types.shape[0], dtype=bool)
        if threshold_db is not None:
            m |= lam > threshold_db
        if psf_threshold is not None:
            m |= psf > psf_threshold
    if include_scrubs:
        m |= types == TYPE_CODE["scrub"]
    return np.asarray(real.t[m])


# --------------------------------------------------------------------------
# (v) public-access availability
# --------------------------------------------------------------------------


@dataclass
class AccessResult:
    closure_hours: float
    closure_hours_with_anomaly: float
    n_closure_windows: int
    weekend_daylight_available: float
    summer_daylight_available: float
    protected_holiday_hours_closed: float


def _hourly_closed_mask(start_days: FloatArray, hours: FloatArray, n_days: int) -> BoolArray:
    n_h = n_days * 24
    closed = np.zeros(n_h, dtype=bool)
    s = np.floor(start_days * 24.0).astype(np.int64)
    e = np.ceil((start_days + hours / 24.0) * 24.0).astype(np.int64)
    for a, b in zip(s, e, strict=True):
        a = max(a, 0)
        b = min(b, n_h)
        if b > a:
            closed[a:b] = True
    return closed


def access_metrics(real: Realization, cal: OpCalendar, apply_pea_measures: bool = False) -> AccessResult:
    """Closure totals and availability of weekend, holiday and summer daylight.

    With ``apply_pea_measures`` the closure windows that fall on protected
    holidays or on summer weekends are removed before totals are computed,
    mirroring the 2022 PEA mitigation measures; the removed events are
    assumed to be rescheduled and are not re-placed.
    """
    m = real.closure_hours > 0
    start = real.closure_start[m]
    hours = real.closure_hours[m]
    if apply_pea_measures and start.shape[0]:
        no_go = cal.protected_holiday | (cal.summer & cal.weekend)
        d0 = np.clip(np.floor(start).astype(np.int64), 0, cal.n_days - 1)
        d1 = np.clip(np.floor(start + hours / 24.0).astype(np.int64), 0, cal.n_days - 1)
        blocked = np.asarray([bool(no_go[a : b + 1].any()) for a, b in zip(d0, d1, strict=True)], dtype=bool)
        start, hours = start[~blocked], hours[~blocked]
    total = union_length(start, start + hours / 24.0, 0.0, float(cal.n_days)) * 24.0
    closed = _hourly_closed_mask(start, hours, cal.n_days)
    hour_of_day = np.arange(cal.n_days * 24) % 24
    day_idx = np.arange(cal.n_days * 24) // 24
    daylight = (hour_of_day >= NIGHT_END_H) & (hour_of_day < NIGHT_START_H)
    weekend_or_holiday = cal.weekend[day_idx] | cal.holiday[day_idx]
    elig_w = daylight & weekend_or_holiday
    elig_s = daylight & cal.summer[day_idx]
    prot = daylight & cal.protected_holiday[day_idx]
    return AccessResult(
        closure_hours=total,
        closure_hours_with_anomaly=total + real.anomaly_hours,
        n_closure_windows=int(start.shape[0]),
        weekend_daylight_available=1.0 - closed[elig_w].mean() if elig_w.any() else 1.0,
        summer_daylight_available=1.0 - closed[elig_s].mean() if elig_s.any() else 1.0,
        protected_holiday_hours_closed=float(closed[prot].sum()),
    )
