import numpy as np
import pytest
from hypothesis import given, settings
from hypothesis import strategies as st

from launchcum import closed_form as cf
from launchcum.footprints import csel_from_psf, psf_from_csel


def test_ldn_increment_decomposition():
    # 25 -> 456 launches with the night share moving from 0.12 to 0.5
    count = cf.ldn_count_term(25, 456)
    night = cf.ldn_night_term(0.12, 0.5)
    assert count == pytest.approx(12.61, abs=0.01)
    assert 3.0 < night < 4.5
    assert cf.ldn_increment(25, 456, 0.12, 0.5) == pytest.approx(count + night)


def test_ldn_from_events_matches_closed_form():
    rng = np.random.default_rng(0)
    sel = np.full(100, 110.0)
    night = rng.random(100) < 0.3
    base = cf.ldn_from_events(sel, night)
    # doubling the count at the same night share adds 10 log10 2
    doubled = cf.ldn_from_events(np.tile(sel, 2), np.tile(night, 2))
    assert doubled - base == pytest.approx(10 * np.log10(2), abs=1e-9)


@given(st.floats(1e-4, 5.0), st.floats(1e-4, 5.0))
def test_poisson_coverage_bounds(lam, tau):
    phi = float(cf.poisson_coverage(lam, tau))
    lin = float(cf.linear_coverage(lam, tau))
    assert 0.0 <= phi < 1.0
    assert phi <= lin + 1e-12  # overlap only reduces coverage


@settings(max_examples=50)
@given(st.floats(1.0, 2000.0), st.integers(100, 366))
def test_exceedance_days_bounded_by_days(n_events, n_days):
    d = cf.expected_exceedance_days(n_events, n_days)
    assert 0.0 < d < n_days
    assert d <= n_events


def test_linear_extrapolation_error_sign():
    # A saturating metric extrapolated linearly is overstated.
    lam = np.array([25, 100, 456]) / 365.0
    phi = cf.poisson_coverage(lam, 1.0)
    err = cf.linear_extrapolation_error(phi[0], 25, phi, np.array([25, 100, 456]))
    assert err[0] == pytest.approx(0.0, abs=1e-12)
    assert np.all(err[1:] > 0)


def test_csel_conversion_roundtrip_and_values():
    # FAA Desk Reference: Lpk = 127.6 + 20 log10(psf); Lpk - CSEL = 26 dB
    assert float(csel_from_psf(1.0)) == pytest.approx(101.6)
    assert float(csel_from_psf(10.0)) == pytest.approx(121.6)
    p = np.array([0.5, 2.0, 21.0])
    assert np.allclose(psf_from_csel(csel_from_psf(p)), p)


def test_cdnl_desk_reference_form():
    # 22 identical daytime booms per year at 10 psf give CDNL = CSEL + 10 log10(22/365) - 49.4
    csel = float(csel_from_psf(10.0))
    expected = csel + 10 * np.log10(22 / 365) - 49.4
    sim = 10 * np.log10(22 * 10 ** (csel / 10) / cf.SECONDS_PER_YEAR)
    assert sim == pytest.approx(expected, abs=0.05)  # the Desk Reference rounds 10 log10(86400) to 49.4
