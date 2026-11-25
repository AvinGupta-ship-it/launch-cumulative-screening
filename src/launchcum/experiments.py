"""Computational experiments E1-E8 and the model-verification case.

Each function returns tidy ``pandas`` tables; orchestration, provenance and
file output live in :mod:`launchcum.pipeline`.
"""

from __future__ import annotations

import dataclasses
from collections.abc import Mapping, Sequence
from concurrent.futures import ProcessPoolExecutor
from dataclasses import dataclass, replace

import numpy as np
import pandas as pd
from numpy.typing import NDArray

from . import closed_form as cf
from .emissions import EmissionInputs, Range, SiteCadence, inventory
from .footprints import FootprintLibrary
from .io import PipelineConfig, ScenarioConfig
from .metrics import (
    ReceptorLevels,
    access_metrics,
    additive_totals,
    boom_counts,
    cdnl,
    coverage,
    disturbance_times,
    exceedance_days,
    ldn,
    ldn_from_counts,
    number_above,
    receptor_levels,
)
from .opcalendar import OpCalendar, build_calendar
from .schedule import EVENT_TYPES, TYPE_CODE, ClosureModel, Realization, ScheduleParams, generate
from .sensitivity import Factor, saltelli_design, sobol_indices
from .spatial import Grid, ReceptorSet

FloatArray = NDArray[np.float64]

DEFAULT_RECOVERY_DAYS = 1.0 / 24.0  # one hour
DEFAULT_WINDOW_DAYS = 4.0 / 24.0  # four hours


@dataclass(frozen=True)
class MetricSpec:
    """What to compute per realization."""

    lamax_thresholds: tuple[float, ...] = (90.0, 100.0, 110.0, 115.0)
    psf_thresholds: tuple[float, ...] = (1.0, 2.0, 4.0, 6.0)
    disturbance_lamax_db: float = 90.0
    recovery_days: float = DEFAULT_RECOVERY_DAYS
    window_days: float = DEFAULT_WINDOW_DAYS
    season: tuple[float, float] = (59.0, 243.0)
    coverage_receptors: tuple[int, ...] = ()
    per_operation: dict[str, dict[str, float]] = dataclasses.field(default_factory=dict)
    apply_pea_measures: bool = False


def realization_metrics(
    real: Realization, lv: ReceptorLevels, cal: OpCalendar, spec: MetricSpec, rng: np.random.Generator
) -> dict[str, float]:
    """All metrics of one realization as a flat mapping."""
    out: dict[str, float] = {
        "n_launch": real.count("launch"),
        "n_sh_landing": real.count("sh_landing"),
        "n_ss_landing": real.count("ss_landing"),
        "n_static_fire": real.count("static_fire"),
        "n_scrub": real.count("scrub"),
        "n_night_launch": int(np.sum(real.night & (real.type_code == TYPE_CODE["launch"]))),
        "n_external": int(real.external.sum()),
        "feasible": float(real.feasible),
        "renewal_rejections": real.renewal_rejections,
    }
    for q, v in additive_totals(real, spec.per_operation).items():
        out[q] = v
    names = lv.names
    ldn_vals = ldn(real, lv)
    c = cdnl(real, lv)
    for k, nm in enumerate(names):
        out[f"ldn|{nm}"] = ldn_vals[k]
        out[f"cdnl|{nm}"] = c[k]
    for thr in spec.lamax_thresholds:
        na = number_above(real, lv, thr)
        nan = number_above(real, lv, thr, night_only=True)
        for k, nm in enumerate(names):
            out[f"na{int(thr)}|{nm}"] = na[k]
            out[f"na{int(thr)}night|{nm}"] = nan[k]
    for thr in (spec.lamax_thresholds[0], spec.lamax_thresholds[1]):
        d = exceedance_days(real, lv, thr)
        for k, nm in enumerate(names):
            out[f"days{int(thr)}|{nm}"] = d[k]
    for p in spec.psf_thresholds:
        b = boom_counts(real, lv, p)
        for k, nm in enumerate(names):
            out[f"boom{p:g}|{nm}"] = b[k]
    for k in spec.coverage_receptors:
        nm = names[k]
        t = disturbance_times(real, lv, k, spec.disturbance_lamax_db, psf_threshold=1.0)
        tau = np.asarray(rng.exponential(spec.recovery_days, size=t.shape[0]), dtype=float)
        cov = coverage(t, tau, spec.season, (spec.window_days,))
        out[f"disturbed|{nm}"] = cov.disturbed_fraction
        out[f"longest_gap|{nm}"] = cov.longest_gap_days
        out[f"has_window|{nm}"] = float(cov.has_window[spec.window_days])
        out[f"n_dist|{nm}"] = float(np.sum((t >= spec.season[0]) & (t < spec.season[1])))
    acc = access_metrics(real, cal, spec.apply_pea_measures)
    out["closure_h"] = acc.closure_hours
    out["closure_anom_h"] = acc.closure_hours_with_anomaly
    out["n_closure_windows"] = acc.n_closure_windows
    out["weekend_avail"] = acc.weekend_daylight_available
    out["summer_avail"] = acc.summer_daylight_available
    out["protected_holiday_closed_h"] = acc.protected_holiday_hours_closed
    return out


def _chunk(
    args: tuple[ScheduleParams, ReceptorLevels, MetricSpec, list[np.random.SeedSequence], int],
) -> list[dict[str, float]]:
    params, lv, spec, seeds, offset = args
    cal = build_calendar(params.year)
    rows = []
    for i, s in enumerate(seeds):
        rng = np.random.default_rng(s)
        real = generate(params, rng, cal)
        row = realization_metrics(real, lv, cal, spec, rng)
        row["realization"] = offset + i
        rows.append(row)
    return rows


def simulate(
    params: ScheduleParams,
    lv: ReceptorLevels,
    spec: MetricSpec,
    n: int,
    seed: np.random.SeedSequence,
    n_workers: int = 1,
) -> pd.DataFrame:
    """Monte Carlo table of metrics (one row per realization)."""
    children = seed.spawn(n)
    if n_workers <= 1 or n < 4 * n_workers:
        rows = _chunk((params, lv, spec, children, 0))
    else:
        bounds = np.linspace(0, n, n_workers + 1).astype(int)
        tasks = [(params, lv, spec, children[a:b], int(a)) for a, b in zip(bounds[:-1], bounds[1:], strict=True)]
        rows = []
        with ProcessPoolExecutor(max_workers=n_workers) as ex:
            for part in ex.map(_chunk, tasks):
                rows.extend(part)
    return pd.DataFrame(rows)


def summarize(df: pd.DataFrame, group: list[str] | None = None) -> pd.DataFrame:
    """Mean, standard error, standard deviation and 2.5/97.5 percentiles per column."""
    num = df.select_dtypes("number")
    if group:
        g = num.groupby([df[c] for c in group])
        mean = g.mean()
        sd = g.std(ddof=1)
        n = g.size()
        se = sd.div(np.sqrt(n), axis=0)
        q025 = g.quantile(0.025)
        q975 = g.quantile(0.975)
        out = pd.concat({"mean": mean, "se": se, "sd": sd, "q025": q025, "q975": q975}, axis=1)
        return out
    mean_s = num.mean()
    sd_s = num.std(ddof=1)
    se_s = sd_s / np.sqrt(len(num))
    return pd.DataFrame(
        {"mean": mean_s, "se": se_s, "sd": sd_s, "q025": num.quantile(0.025), "q975": num.quantile(0.975)}
    )


# --------------------------------------------------------------------------
# verification against the 2025 Tiered EA (model verification level)
# --------------------------------------------------------------------------


def verification_fields(lib: FootprintLibrary, grid: Grid, static_fire_ss_available: bool) -> dict[str, FloatArray]:
    """Deterministic DNL and CDNL fields for the 2025 EA cumulative-noise scenario.

    Scenario (Appendix B section 6.1; sonic boom analysis p.57): 25 launches
    (22 day / 3 night), 22 daytime Super Heavy and 22 daytime Starship
    landings at the site (the 3 night landings of each are offshore), 5
    booster and 6 ship static fires (daytime); CDNL from 44 daytime booms.
    """
    sel = np.stack([lib.level(et, "sel", grid.r_km) for et in EVENT_TYPES])
    n_day = np.zeros(len(EVENT_TYPES))
    n_night = np.zeros(len(EVENT_TYPES))
    n_day[TYPE_CODE["launch"]] = 22
    n_night[TYPE_CODE["launch"]] = 3
    n_day[TYPE_CODE["sh_landing"]] = 22
    n_day[TYPE_CODE["ss_landing"]] = 22
    n_day[TYPE_CODE["static_fire"]] = 5 if static_fire_ss_available else 11
    dnl = ldn_from_counts(sel, n_day, n_night)
    if static_fire_ss_available:
        sel_ss = lib.level("static_fire_ss", "sel", grid.r_km)[None]
        e_ss = 6 * np.power(10.0, sel_ss / 10.0)
        e_ss[~np.isfinite(sel_ss)] = 0.0
        dnl = 10 * np.log10(np.power(10.0, dnl / 10.0) + e_ss[0] / cf.SECONDS_PER_YEAR)
    from .footprints import csel_from_psf

    csel_sh = csel_from_psf(lib.boom("sh_landing", grid.r_km, grid.azimuth_deg))
    csel_ss = csel_from_psf(lib.boom("ss_landing", grid.r_km, grid.azimuth_deg))
    e = 22 * np.power(10.0, csel_sh / 10.0) + 22 * np.power(10.0, csel_ss / 10.0)
    cd = 10 * np.log10(e / cf.SECONDS_PER_YEAR)
    return {"dnl": dnl, "cdnl": cd}


def verification_table(
    lib: FootprintLibrary,
    grid: Grid,
    targets: pd.DataFrame,
    per_op: dict[str, dict[str, float]],
    emissions_per_op: dict[str, dict[str, float]],
) -> pd.DataFrame:
    has_ss = ("static_fire_ss", "sel") in lib.noise
    fields = verification_fields(lib, grid, has_ss)
    rows = []
    tgt = targets.set_index("quantity")
    dnl_ext_km = grid.max_extent_km(fields["dnl"], 65.0)
    cd_ext_km = grid.max_extent_km(fields["cdnl"], 60.0)
    rows.append(
        {
            "quantity": "DNL 65 dB contour extent (miles)",
            "ea_value": f"{_t(tgt, 'dnl65_extent_miles', 'value')}-{_t(tgt, 'dnl65_extent_miles', 'value_hi')}",
            "reproduced": dnl_ext_km / 1.609344,
            "basis": "22 day + 3 night launches, 22 day landings of each vehicle, 11 static fires; SEL footprints",
        }
    )
    dig = _digitized_dnl()
    if dig is not None:
        for lvl, r_mi in dig.items():
            rows.append(
                {
                    "quantity": f"DNL {int(lvl)} dB overland contour radius, digitized Figure 21 (miles)",
                    "ea_value": f"{r_mi:.2f}",
                    "reproduced": grid.max_extent_km(fields["dnl"], float(lvl)) / 1.609344,
                    "basis": "same scenario; median radius over the overland sector of the digitized DNL contour",
                }
            )
    rows.append(
        {
            "quantity": "CDNL 60 dB contour extent (miles)",
            "ea_value": f"{_t(tgt, 'cdnl60_extent_miles', 'value')}",
            "reproduced": cd_ext_km / 1.609344,
            "basis": "44 daytime booms; FAA Desk Reference psf-to-CSEL conversion; fitted boom footprints",
        }
    )
    counts = {"launch": 25, "ss_landing": 25, "sh_landing": 25, "static_fire": 11}
    for q, key, unit in (
        ("nox_tons", "nox_tpy", "short tons"),
        ("co_tons", "co_tpy", "short tons"),
        ("co2e_t", "co2e_t", "metric tons"),
    ):
        tot = sum(emissions_per_op[q].get(et, 0.0) * n for et, n in counts.items())
        rows.append(
            {
                "quantity": f"annual {q.split('_')[0].upper()} ({unit})",
                "ea_value": _t(tgt, key, "value"),
                "reproduced": tot,
                "basis": "per-operation values from Table 3/4 divided by the stated activity counts, re-aggregated",
            }
        )
    rows.append(
        {
            "quantity": "deluge water per operation (gallons)",
            "ea_value": f"{_t(tgt, 'deluge_per_operation_gal', 'value')}",
            "reproduced": per_op["deluge_gal"]["launch"],
            "basis": "transcribed value",
        }
    )
    return pd.DataFrame(rows)


def _t(tgt: pd.DataFrame, key: str, col: str) -> str:
    v = tgt.loc[key, col]
    try:
        return f"{float(str(v)):g}"
    except (TypeError, ValueError):
        return str(v)


def _digitized_dnl() -> dict[float, float] | None:
    from .data.digitize import OUT

    if not OUT.exists():
        return None
    d = pd.read_csv(OUT)
    d = d[(d["metric"] == "dnl") & d["r_overland_miles"].notna()]
    return {float(r["level"]): float(r["r_overland_miles"]) for _, r in d.iterrows() if r["level"] in (65.0, 75.0)}


def emissions_per_operation(
    table3: pd.DataFrame, table4: pd.DataFrame, basis: dict[str, float]
) -> dict[str, dict[str, float]]:
    """Per-operation CO, NOx (short tons) and CO2e (t) from the EA tables.

    Launch values include the methane venting line of Table 4 because the
    EA attributes venting to landing attempts, which accompany launches.
    """

    def t3(activity: str, pol: str) -> float:
        row = table3[(table3["activity"] == activity) & (table3["pollutant"] == pol)]
        return float(row["proposed_tpy"].iloc[0])

    def t4(activity: str) -> float:
        return float(table4[table4["activity"] == activity]["proposed_t"].iloc[0])

    nl, nss, nsh, nsf = basis["launch"], basis["ss_landing"], basis["sh_landing"], basis["static_fire"]
    out = {
        "co_tons": {
            "launch": t3("Starship/Super Heavy Launches", "CO") / nl,
            "ss_landing": t3("Starship Landings", "CO") / nss,
            "sh_landing": t3("Super Heavy Landings", "CO") / nsh,
            "static_fire": (t3("Starship Static Test Fires", "CO") + t3("Super Heavy Static Test Fires", "CO")) / nsf,
        },
        "nox_tons": {
            "launch": t3("Starship/Super Heavy Launches", "NOx") / nl,
            "ss_landing": t3("Starship Landings", "NOx") / nss,
            "sh_landing": t3("Super Heavy Landings", "NOx") / nsh,
            "static_fire": (t3("Starship Static Test Fires", "NOx") + t3("Super Heavy Static Test Fires", "NOx")) / nsf,
        },
        "co2e_t": {
            "launch": (t4("Starship/Super Heavy Launches") + t4("Methane Venting")) / nl,
            "ss_landing": t4("Starship Landings") / nss,
            "sh_landing": t4("Super Heavy Landings") / nsh,
            "static_fire": (t4("Starship Static Test Fires") + t4("Super Heavy Static Test Fires")) / nsf,
        },
    }
    return out


# --------------------------------------------------------------------------
# E1-E5: Monte Carlo sweeps
# --------------------------------------------------------------------------


def params_for(base: ScenarioConfig, **over: object) -> ScheduleParams:
    p = base.operations.to_params(base.year)
    closures = over.pop("closures", None)
    p = replace(p, **over)  # type: ignore[arg-type]
    if isinstance(closures, ClosureModel):
        p = replace(p, closures=closures)
    return p


def sweep(
    cases: Sequence[tuple[Mapping[str, float | str], ScheduleParams, ReceptorLevels]],
    spec: MetricSpec,
    n: int,
    seed: np.random.SeedSequence,
    n_workers: int,
) -> pd.DataFrame:
    """Run a list of cases; each case is (labels, params, levels)."""
    frames = []
    for (labels, params, lv), s in zip(cases, seed.spawn(len(cases)), strict=True):
        df = simulate(params, lv, spec, n, s, n_workers)
        for k, v in labels.items():
            df[k] = v
        frames.append(df)
    return pd.concat(frames, ignore_index=True)


def e5_coverage(
    params: ScheduleParams,
    lv: ReceptorLevels,
    receptor: int,
    recovery_days: list[float],
    window_days: list[float],
    season: tuple[float, float],
    disturbance_lamax_db: float,
    n: int,
    seed: np.random.SeedSequence,
) -> pd.DataFrame:
    """Coverage response surface over recovery time and window length.

    The same realizations are reused for every (tau, w) pair so that the
    surface is smooth in the Monte Carlo noise.
    """
    cal = build_calendar(params.year)
    rows = []
    for i, s in enumerate(seed.spawn(n)):
        rng = np.random.default_rng(s)
        real = generate(params, rng, cal)
        t = disturbance_times(real, lv, receptor, disturbance_lamax_db, psf_threshold=1.0)
        u = rng.random(t.shape[0])
        for tau in recovery_days:
            tau_i = -tau * np.log(1.0 - u)  # exponential with mean tau, common random numbers
            cov = coverage(t, tau_i, season, tuple(window_days))
            row = {
                "realization": i,
                "recovery_days": tau,
                "disturbed": cov.disturbed_fraction,
                "longest_gap": cov.longest_gap_days,
                "n_events": int(np.sum((t >= season[0]) & (t < season[1]))),
            }
            for w in window_days:
                row[f"has_window|{w:g}"] = float(cov.has_window[w])
                row[f"n_windows|{w:g}"] = cov.n_windows[w]
            rows.append(row)
    return pd.DataFrame(rows)


# --------------------------------------------------------------------------
# E6: Sobol
# --------------------------------------------------------------------------

SOBOL_OUTPUTS = ("ldn", "cdnl", "na100", "days90", "boom2", "disturbed", "longest_gap", "closure_h", "weekend_avail")


@dataclass
class SobolContext:
    base: ScheduleParams
    lib: FootprintLibrary
    receptors: ReceptorSet
    community_receptor: str
    habitat_receptor: str
    spec: MetricSpec
    factors: list[Factor]
    cadence: float


def _sobol_eval(args: tuple[SobolContext, FloatArray, list[np.random.SeedSequence], int]) -> list[dict[str, float]]:
    ctx, x, seeds, offset = args
    cal = build_calendar(ctx.base.year)
    names = [f.name for f in ctx.factors]
    k_comm = ctx.receptors.names.index(ctx.community_receptor)
    k_hab = ctx.receptors.names.index(ctx.habitat_receptor)
    rows = []
    for i, (row, s) in enumerate(zip(x, seeds, strict=True)):
        v = dict(zip(names, row, strict=True))
        # The random seed is an explicit factor so that Monte Carlo noise is
        # attributed to its own index instead of inflating the total-order
        # indices of the model inputs (common random numbers across the
        # Saltelli matrices).
        rng = np.random.default_rng(int(v["random_seed"])) if "random_seed" in v else np.random.default_rng(s)
        closures = replace(
            ctx.base.closures,
            launch_mean_h=float(v.get("closure_launch_mean_h", ctx.base.closures.launch_mean_h)),
            launch_cv=float(v.get("closure_launch_cv", ctx.base.closures.launch_cv)),
            anomaly_probability=float(v.get("anomaly_probability", ctx.base.closures.anomaly_probability)),
        )
        params = replace(
            ctx.base,
            n_launches=ctx.cadence,
            night_fraction=float(v.get("night_fraction", ctx.base.night_fraction)),
            scrub_probability=float(v.get("scrub_probability", ctx.base.scrub_probability)),
            turnaround_min_days=float(v.get("turnaround_min_days", ctx.base.turnaround_min_days)),
            n_pads=int(np.floor(v.get("n_pads", ctx.base.n_pads))),
            cluster_rate_per_year=float(v.get("cluster_rate_per_year", ctx.base.cluster_rate_per_year)),
            cluster_spread_days=float(v.get("cluster_spread_days", ctx.base.cluster_spread_days)),
            static_fires_per_launch=float(v.get("static_fires_per_launch", ctx.base.static_fires_per_launch)),
            n_external_landings=float(v.get("external_fraction", 0.0)) * ctx.cadence,
            p_ss_landing_local=float(v.get("p_ss_landing_local", ctx.base.p_ss_landing_local)),
            closures=closures,
        )
        lv = receptor_levels(
            ctx.lib,
            ctx.receptors.names,
            ctx.receptors.r_km,
            ctx.receptors.azimuth_deg,
            level_offset_db=float(v.get("level_offset_db", 0.0)),
            boom_scale=float(v.get("boom_scale", 1.0)),
        )
        real = generate(params, rng, cal)
        out: dict[str, float] = {"sample": offset + i}
        out["ldn"] = float(ldn(real, lv)[k_comm])
        out["cdnl"] = float(cdnl(real, lv)[k_comm])
        out["na100"] = float(number_above(real, lv, 100.0)[k_comm])
        out["days90"] = float(exceedance_days(real, lv, 90.0)[k_comm])
        out["boom2"] = float(boom_counts(real, lv, 2.0)[k_comm])
        t = disturbance_times(real, lv, k_hab, ctx.spec.disturbance_lamax_db, psf_threshold=1.0)
        tau = float(v.get("recovery_days", ctx.spec.recovery_days))
        w = float(v.get("window_days", ctx.spec.window_days))
        cov = coverage(t, np.asarray(rng.exponential(tau, size=t.shape[0]), dtype=float), ctx.spec.season, (w,))
        out["disturbed"] = cov.disturbed_fraction
        out["longest_gap"] = cov.longest_gap_days
        acc = access_metrics(real, cal)
        out["closure_h"] = acc.closure_hours_with_anomaly
        out["weekend_avail"] = acc.weekend_daylight_available
        rows.append(out)
    return rows


def e6_sobol(
    ctx: SobolContext, n_base: int, n_boot: int, seed: np.random.SeedSequence, n_workers: int
) -> tuple[pd.DataFrame, pd.DataFrame]:
    x = saltelli_design(ctx.factors, n_base, int(seed.generate_state(1)[0] % (2**31)))
    n = x.shape[0]
    children = seed.spawn(n)
    if n_workers <= 1:
        rows = _sobol_eval((ctx, x, children, 0))
    else:
        bounds = np.linspace(0, n, n_workers * 4 + 1).astype(int)
        tasks = [(ctx, x[a:b], children[a:b], int(a)) for a, b in zip(bounds[:-1], bounds[1:], strict=True)]
        rows = []
        with ProcessPoolExecutor(max_workers=n_workers) as ex:
            for part in ex.map(_sobol_eval, tasks):
                rows.extend(part)
    y = pd.DataFrame(rows).sort_values("sample").reset_index(drop=True)
    tables = []
    for out in SOBOL_OUTPUTS:
        si = sobol_indices(
            ctx.factors, y[out].to_numpy(float), n_boot=n_boot, seed=int(seed.generate_state(1)[0] % (2**31))
        )
        si["output"] = out
        si["cadence"] = ctx.cadence
        tables.append(si)
    return pd.concat(tables, ignore_index=True), y


# --------------------------------------------------------------------------
# E7: gridded fields and contour areas
# --------------------------------------------------------------------------


def expected_counts(df: pd.DataFrame) -> tuple[FloatArray, FloatArray]:
    """Mean day and night counts per event type from a Monte Carlo table."""
    n_day = np.zeros(len(EVENT_TYPES))
    n_night = np.zeros(len(EVENT_TYPES))
    for et in ("launch", "sh_landing", "ss_landing", "static_fire"):
        col = f"n_{et}"
        if col in df:
            n_day[TYPE_CODE[et]] = df[col].mean()
    nl = df["n_launch"].mean()
    fn = df["n_night_launch"].mean() / nl if nl else 0.0
    for et in ("launch", "sh_landing", "ss_landing"):
        tot = n_day[TYPE_CODE[et]]
        n_night[TYPE_CODE[et]] = tot * fn
        n_day[TYPE_CODE[et]] = tot * (1 - fn)
    return n_day, n_night


def grid_fields(
    lib: FootprintLibrary,
    grid: Grid,
    n_day: FloatArray,
    n_night: FloatArray,
    level_offset_db: float = 0.0,
    boom_scale: float = 1.0,
) -> dict[str, FloatArray]:
    sel = np.stack([lib.level(et, "sel", grid.r_km, level_offset_db) for et in EVENT_TYPES])
    lamax = np.stack([lib.level(et, "lamax", grid.r_km, level_offset_db) for et in EVENT_TYPES])
    fields = {"ldn": ldn_from_counts(sel, n_day, n_night)}
    from .footprints import csel_from_psf

    csel = np.stack(
        [
            csel_from_psf(lib.boom(et, grid.r_km, grid.azimuth_deg, boom_scale))
            if et in lib.booms
            else np.full(grid.r_km.shape, -np.inf)
            for et in EVENT_TYPES
        ]
    )
    fields["cdnl"] = ldn_from_counts(csel, n_day, n_night)
    for thr in (90.0, 100.0, 110.0):
        above = (lamax > thr).astype(float)
        fields[f"na{int(thr)}"] = (above * (n_day + n_night)[:, None, None]).sum(axis=0)
        fields[f"na{int(thr)}night"] = (above * n_night[:, None, None]).sum(axis=0)
    for p in (1.0, 2.0, 4.0, 6.0):
        psf = np.stack(
            [
                lib.boom(et, grid.r_km, grid.azimuth_deg, boom_scale) if et in lib.booms else np.zeros(grid.r_km.shape)
                for et in EVENT_TYPES
            ]
        )
        fields[f"boom{p:g}"] = ((psf > p).astype(float) * (n_day + n_night)[:, None, None]).sum(axis=0)
    return fields


# --------------------------------------------------------------------------
# E8: fleet emissions
# --------------------------------------------------------------------------


def e8_emissions(cfg: PipelineConfig, rng: np.random.Generator) -> pd.DataFrame:
    e = cfg.emissions
    inputs = EmissionInputs(
        propellant_per_launch_t=Range(
            e.propellant_per_launch_t.low, e.propellant_per_launch_t.high, e.propellant_per_launch_t.note
        ),
        f_stratosphere=Range(e.f_stratosphere.low, e.f_stratosphere.high, e.f_stratosphere.note),
        ei_g_per_kg={k: Range(v.low, v.high, v.note) for k, v in e.emission_index_g_per_kg.items()},
    )
    sites = [SiteCadence(s.name, s.launches_per_year, s.source) for s in e.sites]
    return inventory(sites, inputs, rng, e.n_samples)
