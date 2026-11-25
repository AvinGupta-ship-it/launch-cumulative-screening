from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from launchcum.data.transcribed import TABLES, TranscriptionMismatch, reconcile, reconcile_all
from launchcum.footprints import MILES_TO_KM, build_library
from launchcum.spatial import Grid, ReceptorSet, pad_utm

ROOT = Path(__file__).resolve().parents[1]


@pytest.fixture(scope="module")
def tables():
    return reconcile_all()


@pytest.fixture(scope="module")
def receptors():
    px, py = pad_utm(-97.15456, 25.9963)
    return ReceptorSet.from_lonlat(pd.read_csv(ROOT / "data/processed/receptors.csv"), px, py)


def test_every_transcribed_table_has_two_entries_that_reconcile(tables):
    base = ROOT / "data" / "transcribed"
    for name in TABLES:
        assert (base / f"{name}_entry_a.csv").exists()
        assert (base / f"{name}_entry_b.csv").exists()
        assert len(tables[name]) > 0


def test_reconcile_detects_mismatch(tmp_path):
    base = ROOT / "data" / "transcribed"
    a = (base / "ghg_table4_entry_a.csv").read_text()
    (tmp_path / "ghg_table4_entry_a.csv").write_text(a)
    (tmp_path / "ghg_table4_entry_b.csv").write_text(a.replace("97342", "97324"))
    with pytest.raises(TranscriptionMismatch):
        reconcile("ghg_table4", tmp_path)


def test_transcribed_rows_cite_pages(tables):
    for name, df in tables.items():
        assert "page" in df.columns, name
        assert df["page"].notna().all(), name


def test_footprint_interpolants_reproduce_points_within_tolerance(tables, receptors):
    from launchcum.pipeline import footprint_points_with_digitized

    pts = footprint_points_with_digitized(tables["footprint_points"])
    lib = build_library(pts, receptors.distance_map())
    res = lib.residual_table()
    dig_ids = set(pts[pts["source"].astype(str).str.startswith("digitized")]["point_id"])
    for (et, m), fit in lib.noise.items():
        g = pts[(pts["event_type"] == et) & (pts["metric"] == m)]
        if fit.interp_r_km.size >= 2:
            # digitized contour radii are reproduced to within 1 dB (monotone interpolation)
            dig = g[g["point_id"].isin(dig_ids)]
            r = dig["distance_miles"].to_numpy(float) * MILES_TO_KM
            lvl = dig["value"].to_numpy(float)
            assert np.max(np.abs(fit(r) - lvl)) < 1.0 + 1e-9, (et, m)
            # and the level is monotone non-increasing in distance
            rr = np.linspace(0.3, 60, 300)
            assert np.all(np.diff(fit(rr)) <= 1e-9), (et, m)
    # residuals of text-quoted extents are reported; the largest disagreement is bounded
    assert res["residual_db"].abs().max() < 15.0
    for et, bf in lib.booms.items():
        assert np.all(np.diff(bf(np.linspace(1, 60, 50))) <= 1e-12), et
        assert bf(0.5) <= bf.psf_cap * 1.0000001


def test_boom_cap_and_extent(tables, receptors):
    lib = build_library(tables["footprint_points"], receptors.distance_map())
    sh = lib.booms["sh_landing"]
    assert sh(0.5) == pytest.approx(21.0)
    # transcribed 1 psf extent is 28 miles; the fitted law must put 1 psf within a factor 2 of it
    d1 = sh.distance_for_psf(1.0) / MILES_TO_KM
    assert 14.0 < d1 < 56.0


def test_grid_geometry():
    g = Grid(0.0, 0.0, 30_000.0, 250.0)
    assert g.n == 240
    assert g.r_km.min() < 0.2
    assert abs(g.r_km.max() - np.hypot(29.875, 29.875)) < 0.01
    field = (g.r_km < 5.0).astype(float)
    assert g.area_above(field, 0.5) == pytest.approx(np.pi * 25.0, rel=0.02)
    assert g.max_extent_km(field, 0.5) < 5.0


def test_receptor_distances(receptors):
    names = receptors.names
    assert "Port Isabel (nearest edge)" in names and "Piping plover critical habitat (near pad)" in names
    r = receptors.distance_map()
    assert 5.0 < r["Port Isabel (nearest edge)"] < 12.0
    assert r["Brownsville"] > r["Brownsville (nearest edge)"]
