"""Per-event footprints digitized from the FAA environmental record.

Footprints are not computed from a new acoustic model. Contour extents and
receptor-level bands stated in the 2025 Tiered Environmental Assessment and
its noise appendix are transcribed to ``data/transcribed`` with page
citations; this module fits simple propagation forms to those points and
reports residuals.

Noise levels follow ``L(r) = L0 - 20 log10(r / r0) - alpha (r - r0)`` with
``r`` in kilometres, ``r0 = 1`` km. Sonic boom overpressure follows a
piecewise log-linear interpolant in ``(log r, log p)`` capped at the maximum
overpressure stated in the record. Hypothetical vehicle variants enter
through additive level offsets ``s_L`` (dB) and multiplicative boom factors
``s_p``.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import pandas as pd
from numpy.typing import ArrayLike, NDArray

FloatArray = NDArray[np.float64]

MILES_TO_KM = 1.609344
R0_KM = 1.0

# FAA Order 1050.1F Desk Reference, chapter 11 (Noise), "Calculation of CDNL":
# Lpk = 127.6 + 20 log10(psf); Lpk - CSEL = 26 dB.
LPK_OFFSET_DB = 127.6
LPK_MINUS_CSEL_DB = 26.0


@dataclass
class PropagationFit:
    """Fitted level-distance relation for one event type and metric."""

    event_type: str
    metric: str
    l0: float
    alpha_per_km: float
    r_km: FloatArray
    level_db: FloatArray
    weight: FloatArray
    residual_db: FloatArray
    shared_alpha: bool
    r_min_km: float = 0.3
    level_cap_db: float | None = None
    interp_r_km: FloatArray = field(default_factory=lambda: np.zeros(0))
    interp_level_db: FloatArray = field(default_factory=lambda: np.zeros(0))

    def law(self, r_km: ArrayLike, offset_db: float = 0.0) -> FloatArray:
        """The fitted propagation law alone."""
        r = np.maximum(np.asarray(r_km, dtype=float), self.r_min_km)
        lvl = self.l0 - 20.0 * np.log10(r / R0_KM) - self.alpha_per_km * (r - R0_KM) + offset_db
        if self.level_cap_db is not None:
            lvl = np.minimum(lvl, self.level_cap_db + offset_db)
        return np.asarray(lvl)

    def __call__(self, r_km: ArrayLike, offset_db: float = 0.0) -> FloatArray:
        """Level at distance ``r_km``.

        Inside the range of digitized contours the level is interpolated
        monotonically in ``log r`` through the contour radii; beyond the
        outermost contour the fitted law, shifted to match the outermost
        contour, is used. Inside the innermost contour the level is held at
        that contour's value.
        """
        r = np.asarray(r_km, dtype=float)
        if self.interp_r_km.size < 2:
            return self.law(r, offset_db)
        lr = np.log10(np.maximum(r, 1e-3))
        xs = np.log10(self.interp_r_km)
        ys = self.interp_level_db
        inner = ys[0]
        out = np.interp(lr, xs, ys)
        beyond = r > self.interp_r_km[-1]
        if np.any(beyond):
            shift = ys[-1] - self.law(self.interp_r_km[-1])[()]
            out = np.where(beyond, self.law(r) + shift, out)
        out = np.where(r < self.interp_r_km[0], inner, out)
        return np.asarray(out + offset_db)

    def distance_for_level(self, level_db: float, offset_db: float = 0.0) -> float:
        """Distance (km) at which the fitted level equals ``level_db`` (bisection)."""
        target = level_db
        lo, hi = self.r_min_km, 1000.0
        if self(lo, offset_db) < target:
            return float(lo)
        if self(hi, offset_db) > target:
            return float(hi)
        for _ in range(100):
            mid = 0.5 * (lo + hi)
            if self(mid, offset_db) > target:
                lo = mid
            else:
                hi = mid
        return float(0.5 * (lo + hi))

    @property
    def rms_residual_db(self) -> float:
        return float(np.sqrt(np.mean(self.residual_db**2))) if self.residual_db.size else float("nan")


@dataclass
class BoomFit:
    """Power-law overpressure fit ``log10 p = a + b log10 r`` with a near-field cap.

    A single log-log line is fitted by weighted least squares to the
    transcribed contour extents and receptor values; it is monotone in
    distance and extrapolates inward until the stated maximum overpressure.
    """

    event_type: str
    r_km: FloatArray
    psf: FloatArray
    psf_cap: float
    weight: FloatArray = field(default_factory=lambda: np.zeros(0))
    a: float = 0.0
    b: float = -1.0
    residual_log10: FloatArray = field(default_factory=lambda: np.zeros(0))
    elongation: float = 0.0
    axis_deg: float = 0.0

    def __post_init__(self) -> None:
        if self.weight.shape[0] != self.r_km.shape[0]:
            self.weight = np.ones_like(self.r_km)
        x = np.log10(self.r_km)
        y = np.log10(self.psf)
        if x.size >= 2 and np.ptp(x) > 0:
            w = np.sqrt(self.weight)
            a_mat = np.column_stack([np.ones_like(x), x]) * w[:, None]
            sol, *_ = np.linalg.lstsq(a_mat, y * w, rcond=None)
            self.a, self.b = float(sol[0]), float(min(sol[1], -0.05))
            if sol[1] > -0.05:
                # Enforce a decreasing law: refit the intercept with the slope pinned.
                self.a = float(np.average(y - self.b * x, weights=self.weight))
        elif x.size == 1:
            self.a, self.b = float(y[0] + x[0]), -1.0
        self.residual_log10 = (self.a + self.b * x) - y

    def __call__(self, r_km: ArrayLike, azimuth_deg: ArrayLike | None = None, scale: float = 1.0) -> FloatArray:
        r = np.asarray(r_km, dtype=float)
        if azimuth_deg is not None and self.elongation > 0:
            # Hypothetical directivity: the footprint is stretched along the
            # approach axis by (1 + e) and compressed across it by (1 - e).
            th = np.deg2rad(np.asarray(azimuth_deg, dtype=float) - self.axis_deg)
            g = np.sqrt((np.cos(th) / (1 + self.elongation)) ** 2 + (np.sin(th) / (1 - self.elongation)) ** 2)
            r = r * g
        lr = np.log10(np.maximum(r, 1e-3))
        p = np.minimum(np.power(10.0, self.a + self.b * lr), self.psf_cap) * scale
        return np.asarray(p)

    def distance_for_psf(self, psf: float, scale: float = 1.0) -> float:
        lo, hi = 0.05, 2000.0
        if self(lo, scale=scale) < psf:
            return float(lo)
        if self(hi, scale=scale) > psf:
            return float(hi)
        for _ in range(100):
            mid = 0.5 * (lo + hi)
            if self(mid, scale=scale) > psf:
                lo = mid
            else:
                hi = mid
        return float(0.5 * (lo + hi))


def csel_from_psf(psf: ArrayLike) -> FloatArray:
    """C-weighted sound exposure level from peak overpressure (FAA Desk Reference)."""
    p = np.maximum(np.asarray(psf, dtype=float), 1e-6)
    return np.asarray(LPK_OFFSET_DB + 20.0 * np.log10(p) - LPK_MINUS_CSEL_DB)


def psf_from_csel(csel_db: ArrayLike) -> FloatArray:
    return np.asarray(np.power(10.0, (np.asarray(csel_db, float) + LPK_MINUS_CSEL_DB - LPK_OFFSET_DB) / 20.0))


def _fit_shared_alpha(groups: list[tuple[FloatArray, FloatArray, FloatArray]]) -> tuple[list[float], float]:
    """Least-squares fit of one alpha and a per-group L0.

    Each group supplies ``(r_km, level_db, weight)``. The model is linear in
    the unknowns ``(L0_g, alpha)`` once spherical spreading is removed.
    """
    rows = []
    rhs = []
    wts = []
    n_g = len(groups)
    for g, (r, lvl, w) in enumerate(groups):
        for ri, li, wi in zip(r, lvl, w, strict=True):
            row = np.zeros(n_g + 1)
            row[g] = 1.0
            row[-1] = -(ri - R0_KM)
            rows.append(row)
            rhs.append(li + 20.0 * np.log10(ri / R0_KM))
            wts.append(wi)
    a = np.asarray(rows) * np.sqrt(np.asarray(wts))[:, None]
    b = np.asarray(rhs) * np.sqrt(np.asarray(wts))
    sol, *_ = np.linalg.lstsq(a, b, rcond=None)
    alpha = float(max(sol[-1], 0.0))
    if alpha != sol[-1]:
        # Re-fit L0 with alpha clipped at zero.
        sol = [float(np.average(lvl + 20 * np.log10(r / R0_KM), weights=w)) for r, lvl, w in groups] + [0.0]
    return [float(s) for s in sol[:-1]], alpha


@dataclass
class FootprintLibrary:
    """Collection of fitted footprints indexed by ``(event_type, metric)``."""

    noise: dict[tuple[str, str], PropagationFit]
    booms: dict[str, BoomFit]
    points: pd.DataFrame

    def level(self, event_type: str, metric: str, r_km: ArrayLike, offset_db: float = 0.0) -> FloatArray:
        key = (event_type, metric)
        if key not in self.noise:
            return np.full(np.shape(r_km), -np.inf)
        return self.noise[key](r_km, offset_db)

    def boom(
        self, event_type: str, r_km: ArrayLike, azimuth_deg: ArrayLike | None = None, scale: float = 1.0
    ) -> FloatArray:
        if event_type not in self.booms:
            return np.zeros(np.shape(r_km))
        return self.booms[event_type](r_km, azimuth_deg, scale)

    def residual_table(self) -> pd.DataFrame:
        recs = []
        for (et, met), fit in self.noise.items():
            law_res = fit.law(fit.r_km) - fit.level_db
            for r, lv, res, lres in zip(fit.r_km, fit.level_db, fit.residual_db, law_res, strict=True):
                recs.append(
                    {
                        "event_type": et,
                        "metric": met,
                        "r_km": r,
                        "transcribed": lv,
                        "residual_db": res,
                        "law_residual_db": lres,
                        "l0": fit.l0,
                        "alpha_per_km": fit.alpha_per_km,
                    }
                )
        for et, bf in self.booms.items():
            for r, p, res in zip(bf.r_km, bf.psf, bf.residual_log10, strict=True):
                recs.append(
                    {
                        "event_type": et,
                        "metric": "boom_psf",
                        "r_km": r,
                        "transcribed": p,
                        "residual_db": 20 * res,
                        "law_residual_db": 20 * res,
                        "l0": np.nan,
                        "alpha_per_km": np.nan,
                    }
                )
        return pd.DataFrame(recs)


def build_library(points: pd.DataFrame, receptor_distance_km: dict[str, float] | None = None) -> FootprintLibrary:
    """Fit footprints from the reconciled transcription table.

    ``points`` must contain columns ``event_type, metric, value, value_hi,
    distance_miles, location, weight``. Rows with a ``location`` instead of
    a distance are resolved through ``receptor_distance_km``; rows whose
    location is unknown are skipped. Band rows (``value_hi`` present) use
    the band midpoint with weight reduced by the band half-width.
    """
    df = points.copy()
    dist = np.where(df["distance_miles"].notna(), df["distance_miles"].astype(float) * MILES_TO_KM, np.nan)
    if receptor_distance_km:
        for i, loc in enumerate(df["location"].fillna("")):
            if np.isnan(dist[i]) and loc in receptor_distance_km:
                dist[i] = receptor_distance_km[loc]
    df["r_km"] = dist
    df = df[df["r_km"].notna() & (df["r_km"] > 0)].copy()
    hi = df["value_hi"].astype(float)
    mid = np.where(hi.notna(), 0.5 * (df["value"].astype(float) + hi.fillna(0)), df["value"].astype(float))
    df["level"] = mid
    halfwidth = np.where(hi.notna(), 0.5 * (hi.fillna(0) - df["value"].astype(float)), 0.0)
    base_w = df["weight"].fillna(1.0).astype(float).to_numpy()
    df["w"] = base_w / (1.0 + halfwidth**2 / 4.0)

    noise: dict[tuple[str, str], PropagationFit] = {}
    for metric in sorted(set(df["metric"]) - {"boom_psf"}):
        sub = df[df["metric"] == metric]
        groups = []
        keys = []
        for et, g in sub.groupby("event_type", sort=True):
            groups.append((g["r_km"].to_numpy(float), g["level"].to_numpy(float), g["w"].to_numpy(float)))
            keys.append(str(et))
        l0s, alpha = _fit_shared_alpha(groups)
        for et, (r, lvl, w), l0 in zip(keys, groups, l0s, strict=True):
            cap = _metric_cap(metric)
            fit = PropagationFit(
                et, metric, l0, alpha, r, lvl, w, np.zeros_like(lvl), shared_alpha=True, level_cap_db=cap
            )
            g = sub[sub["event_type"] == et]
            dig = g[g["source"].astype(str).str.startswith("digitized")]
            if len(dig) >= 2:
                order = np.argsort(dig["r_km"].to_numpy(float))
                rr = dig["r_km"].to_numpy(float)[order]
                ll = dig["level"].to_numpy(float)[order]
                # enforce monotone decrease with distance
                ll = np.minimum.accumulate(ll)
                keep = np.concatenate([[True], np.diff(rr) > 0])
                fit.interp_r_km = rr[keep]
                fit.interp_level_db = ll[keep]
            fit.residual_db = fit(r) - lvl
            noise[(et, metric)] = fit

    booms: dict[str, BoomFit] = {}
    bsub = df[df["metric"] == "boom_psf"]
    for et, g in bsub.groupby("event_type", sort=True):
        r = g["r_km"].to_numpy(float)
        p = g["level"].to_numpy(float)
        w = g["w"].to_numpy(float)
        cap_rows = points[(points["event_type"] == et) & (points["metric"] == "boom_psf_cap")]
        cap = float(cap_rows["value"].max()) if len(cap_rows) else float(p.max())
        order = np.argsort(r)
        booms[str(et)] = BoomFit(str(et), r[order], p[order], cap, weight=w[order])
    return FootprintLibrary(noise=noise, booms=booms, points=df)


def _metric_cap(metric: str) -> float | None:
    # Innermost contour levels stated in the record bound the fitted forms
    # near the pad, where the point-source law is not meaningful.
    return {"lamax": 150.0, "sel": 160.0, "lmax": 160.0}.get(metric)


def load_points(path: Path) -> pd.DataFrame:
    return pd.read_csv(path)
