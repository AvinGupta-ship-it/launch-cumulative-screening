import numpy as np
import pytest
from hypothesis import given, settings
from hypothesis import strategies as st

from launchcum import closed_form as cf
from launchcum.metrics import (
    ReceptorLevels,
    access_metrics,
    cdnl,
    coverage,
    exceedance_days,
    gaps_within,
    ldn,
    merge_intervals,
    number_above,
    union_length,
)
from launchcum.opcalendar import build_calendar
from launchcum.schedule import EVENT_TYPES, TYPE_CODE, ScheduleParams, generate


def _levels(k=2):
    n_t = len(EVENT_TYPES)
    lamax = np.full((n_t, k), -np.inf)
    sel = np.full((n_t, k), -np.inf)
    psf = np.zeros((n_t, k))
    lamax[TYPE_CODE["launch"]] = [105.0, 95.0]
    sel[TYPE_CODE["launch"]] = [110.0, 100.0]
    lamax[TYPE_CODE["sh_landing"]] = [95.0, 85.0]
    sel[TYPE_CODE["sh_landing"]] = [100.0, 90.0]
    psf[TYPE_CODE["sh_landing"]] = [6.0, 1.5]
    return ReceptorLevels(["a", "b"], np.array([5.0, 10.0]), np.zeros(k), lamax, sel, psf)


def test_merge_and_union():
    s = np.array([0.0, 1.0, 5.0, 2.5])
    e = np.array([2.0, 3.0, 6.0, 2.8])
    ms, me = merge_intervals(s, e)
    assert ms.tolist() == [0.0, 5.0] and me.tolist() == [3.0, 6.0]
    assert union_length(s, e, 0.0, 10.0) == pytest.approx(4.0)
    gaps = gaps_within(s, e, 0.0, 10.0)
    assert sorted(gaps.tolist()) == pytest.approx([2.0, 4.0])


def test_coverage_recovers_poisson_closed_form():
    rng = np.random.default_rng(0)
    lam, tau, season = 2.0, 0.5, (0.0, 2000.0)
    phi = []
    for _ in range(20):
        n = rng.poisson(lam * season[1])
        t = np.sort(rng.random(n) * season[1])
        cov = coverage(t, rng.exponential(tau, n), season)
        phi.append(cov.disturbed_fraction)
    assert np.mean(phi) == pytest.approx(float(cf.poisson_coverage(lam, tau)), abs=0.01)


@settings(max_examples=30, deadline=None)
@given(st.floats(0.01, 5.0), st.floats(0.001, 2.0))
def test_coverage_monotone_in_recovery_time(lam, tau):
    rng = np.random.default_rng(1)
    n = rng.poisson(lam * 365)
    t = np.sort(rng.random(n) * 365)
    u = rng.random(n)
    c1 = coverage(t, -tau * np.log(1 - u), (0.0, 365.0)).disturbed_fraction
    c2 = coverage(t, -2 * tau * np.log(1 - u), (0.0, 365.0)).disturbed_fraction
    assert 0.0 <= c1 <= c2 <= 1.0


def test_ldn_and_counts_against_direct_formula():
    cal = build_calendar(2030)
    p = ScheduleParams(
        n_launches=200,
        model="poisson",
        night_fraction=0.5,
        static_fires_per_launch=0,
        scrub_probability=0,
        p_ss_landing_local=0.0,
    )
    real = generate(p, np.random.default_rng(2), cal)
    lv = _levels()
    ldn_vals = ldn(real, lv)
    sel = lv.sel[real.type_code]
    pen = 10 * real.night[:, None]
    e = np.where(np.isfinite(sel), 10 ** ((sel + pen) / 10), 0.0).sum(axis=0)
    assert np.allclose(ldn_vals, 10 * np.log10(e / cf.SECONDS_PER_YEAR))
    na = number_above(real, lv, 100.0)
    assert na[0] == real.count("launch")  # only launches exceed 100 at receptor a
    assert na[1] == 0
    d = exceedance_days(real, lv, 90.0)
    assert d[0] <= 365 and d[0] <= real.count("launch") + real.count("sh_landing")
    c = cdnl(real, lv)
    assert c[0] > c[1]


def test_access_metrics_respect_caps_and_availability():
    cal = build_calendar(2030)
    p = ScheduleParams(n_launches=456, model="weekday_renewal", weekday_only=True, n_pads=2)
    real = generate(p, np.random.default_rng(4), cal)
    acc = access_metrics(real, cal)
    assert acc.closure_hours > 0
    assert acc.closure_hours <= real.closure_hours.sum() + 1e-9  # union never exceeds the sum
    assert 0.0 <= acc.weekend_daylight_available <= 1.0
    acc2 = access_metrics(real, cal, apply_pea_measures=True)
    assert acc2.protected_holiday_hours_closed == 0.0
    assert acc2.closure_hours <= acc.closure_hours + 1e-9
