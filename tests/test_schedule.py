import numpy as np
import pytest
from hypothesis import given, settings
from hypothesis import strategies as st

from launchcum.opcalendar import build_calendar, is_night
from launchcum.schedule import ScheduleParams, events_to_frame, generate, generate_many


@pytest.fixture(scope="module")
def cal():
    return build_calendar(2030)


def test_calendar_counts(cal):
    assert cal.n_days == 365
    assert 248 <= cal.n_business <= 254
    assert cal.weekend.sum() == 104
    assert cal.summer.sum() > 90
    assert cal.protected_holiday.sum() >= 6


@settings(max_examples=20, deadline=None)
@given(st.integers(0, 10_000))
def test_weekday_generator_never_places_launches_on_weekends_or_holidays(seed):
    cal = build_calendar(2030)
    p = ScheduleParams(n_launches=456, model="weekday_renewal", weekday_only=True, n_pads=2, turnaround_min_days=0.5)
    real = generate(p, np.random.default_rng(seed), cal)
    launches = real.type_code == 0
    days = np.floor(real.t[launches]).astype(int)
    assert not cal.weekend[days].any()
    assert not cal.holiday[days].any()
    assert real.count("launch") == 456
    assert real.feasible


def test_campaign_generator_respects_business_days():
    cal = build_calendar(2030)
    p = ScheduleParams(
        n_launches=456, model="campaign", weekday_only=True, cluster_rate_per_year=60, cluster_spread_days=2.0
    )
    real = generate(p, np.random.default_rng(3), cal)
    days = np.floor(real.t[real.type_code == 0]).astype(int)
    assert not cal.weekend[days].any()


def test_poisson_count_in_distribution():
    reals = generate_many(ScheduleParams(n_launches=100, model="poisson"), 7, 400)
    counts = np.array([r.count("launch") for r in reals])
    # mean within 3 standard errors of 100, variance consistent with Poisson
    assert abs(counts.mean() - 100) < 3 * np.sqrt(100 / 400)
    assert 60 < counts.var() < 150


def test_landings_follow_launches_and_night_fraction():
    cal = build_calendar(2030)
    p = ScheduleParams(n_launches=400, model="poisson", night_fraction=0.5, sh_landing_delay_min=7.0)
    real = generate(p, np.random.default_rng(11), cal)
    assert real.count("sh_landing") == real.count("launch")
    launches = real.t[real.type_code == 0]
    sh = real.t[real.type_code == 1]
    assert np.allclose(np.sort(sh) - np.sort(launches), 7.0 / 1440.0)
    night_frac = real.night[real.type_code == 0].mean()
    assert abs(night_frac - 0.5) < 0.1
    assert np.array_equal(real.night, is_night(real.t))


def test_night_definition_matches_noi():
    t = np.array([0.0, 6.9 / 24, 7.0 / 24, 21.9 / 24, 22.0 / 24, 23.5 / 24])
    assert is_night(t).tolist() == [True, True, False, False, True, True]


def test_external_landings_and_frame():
    cal = build_calendar(2030)
    p = ScheduleParams(n_launches=50, model="poisson", n_external_landings=200)
    real = generate(p, np.random.default_rng(5), cal)
    assert real.external.sum() > 150
    df = events_to_frame(real)
    assert set(df.columns) >= {"t_days", "event_type", "night", "closure_hours"}
    assert (df["closure_hours"] >= 0).all()
    assert df["t_days"].is_monotonic_increasing


def test_seed_reproducibility():
    p = ScheduleParams(n_launches=120, model="weekday_renewal", weekday_only=True)
    a = generate(p, np.random.default_rng(42))
    b = generate(p, np.random.default_rng(42))
    assert np.array_equal(a.t, b.t) and np.array_equal(a.type_code, b.type_code)


def test_infeasible_turnaround_is_flagged():
    cal = build_calendar(2030)
    p = ScheduleParams(n_launches=456, model="weekday_renewal", weekday_only=True, n_pads=1, turnaround_min_days=2.0)
    real = generate(p, np.random.default_rng(1), cal)
    assert not real.feasible


def test_generate_without_scrubs_or_static_fires() -> None:
    from dataclasses import replace

    from launchcum.opcalendar import build_calendar
    from launchcum.schedule import TYPE_CODE, ScheduleParams, generate

    p = replace(ScheduleParams(n_launches=30.0), scrub_probability=0.0, static_fires_per_launch=0.0)
    real = generate(p, np.random.default_rng(3), build_calendar(2030))
    assert not np.any(real.type_code == TYPE_CODE["scrub"])
    assert not np.any(real.type_code == TYPE_CODE["static_fire"])
    assert np.all(np.diff(real.t) >= 0)


def test_scrub_attempts_recede_in_order_and_landings_follow_launches() -> None:
    from dataclasses import replace

    from launchcum.opcalendar import build_calendar
    from launchcum.schedule import TYPE_CODE, ScheduleParams, generate

    p = replace(ScheduleParams(n_launches=200.0), scrub_probability=0.5)
    real = generate(p, np.random.default_rng(5), build_calendar(2030))
    t_launch = real.t[real.type_code == TYPE_CODE["launch"]]
    t_ss = real.t[(real.type_code == TYPE_CODE["ss_landing"]) & ~real.external]
    # every local Starship landing follows some launch by at least the minimum delay
    assert np.all(t_ss >= t_launch.min() + p.ss_landing_delay_days_min)
    t_sc = real.t[real.type_code == TYPE_CODE["scrub"]]
    assert np.all(t_sc < t_launch.max())
