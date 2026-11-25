"""Single-command reproduction of every number, table and figure input.

    python -m launchcum.pipeline run --config configs/paper.yaml

writes all experiment outputs under the configured results directory, a
``summary.json`` with headline values and ``paper/generated/values.tex``
with LaTeX macros used by the manuscript and the comment.
"""

from __future__ import annotations

import argparse
import json
import re
import sys
import time
from dataclasses import replace
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import xarray as xr

from . import closed_form as cf
from . import experiments as ex
from .data.digitize import OUT as DIGITIZED_CSV
from .data.digitize import as_footprint_points
from .data.transcribed import per_operation_dict, reconcile_all
from .footprints import FootprintLibrary, build_library
from .io import (
    PipelineConfig,
    load_pipeline,
    load_scenario,
    provenance,
    write_csv,
    write_json,
    write_netcdf,
    write_parquet,
)
from .metrics import ReceptorLevels, receptor_levels
from .schedule import ScheduleParams
from .sensitivity import Factor
from .spatial import Grid, ReceptorSet, pad_utm

ROOT = Path(__file__).resolve().parents[2]
GENERATED = ROOT / "paper" / "generated"
SUMMARY_DIR = ROOT / "results" / "summary"

COMMUNITY_EDGE = "Port Isabel (nearest edge)"
COMMUNITY_SPI = "South Padre Island (nearest edge)"
HABITAT = "Piping plover critical habitat (near pad)"
REFUGE = "Lower Rio Grande Valley NWR (near pad)"
PARK = "Boca Chica State Park"
COVERAGE_RECEPTORS = (COMMUNITY_EDGE, COMMUNITY_SPI, HABITAT, REFUGE, PARK)


TEXT_EXTENT_WEIGHT = 0.25


def footprint_points_with_digitized(transcribed: pd.DataFrame) -> pd.DataFrame:
    """Transcribed text points plus digitized contour radii.

    Contour radii digitized from the Appendix B maps are the primary
    footprint data; extents quoted in the EA text are retained at reduced
    weight because several of them disagree with the maps they summarize.
    """
    pts = transcribed.copy()
    if not DIGITIZED_CSV.exists():
        return pts
    dig = as_footprint_points(pd.read_csv(DIGITIZED_CSV))
    has_dig = set(zip(dig["event_type"], dig["metric"], strict=True))
    txt_mask = [
        (et, m) in has_dig and pd.notna(d)
        for et, m, d in zip(pts["event_type"], pts["metric"], pts["distance_miles"], strict=True)
    ]
    pts.loc[txt_mask, "weight"] = pts.loc[txt_mask, "weight"].astype(float) * TEXT_EXTENT_WEIGHT
    return pd.concat([pts, dig], ignore_index=True)


class Runner:
    def __init__(self, config_path: Path) -> None:
        self.config_path = Path(config_path)
        self.cfg: PipelineConfig = load_pipeline(self.config_path)
        self.out = Path(self.cfg.results_dir)
        self.out.mkdir(parents=True, exist_ok=True)
        self.seed = np.random.SeedSequence(self.cfg.seed)
        self.t0 = time.time()
        self.values: dict[str, Any] = {}
        self.summary: dict[str, Any] = {}
        self.scenarios = {Path(p).stem: load_scenario(Path(p)) for p in self.cfg.scenarios}
        self.tables = reconcile_all()
        rec = pd.read_csv(self.cfg.receptors)
        self.pad_x, self.pad_y = pad_utm(self.cfg.pad_lon, self.cfg.pad_lat)
        self.receptors = ReceptorSet.from_lonlat(rec, self.pad_x, self.pad_y)
        self.points = footprint_points_with_digitized(self.tables["footprint_points"])
        self.lib: FootprintLibrary = build_library(self.points, self.receptors.distance_map())
        self.per_op_raw = per_operation_dict(self.tables["per_operation_values"])
        basis = self.per_op_raw["basis_count"]
        self.emissions_per_op = ex.emissions_per_operation(
            self.tables["emission_rates_table3"], self.tables["ghg_table4"], basis
        )
        self.per_operation = {"deluge_gal": self.per_op_raw["deluge_gal"], **self.emissions_per_op}
        self.grid = Grid(self.pad_x, self.pad_y, self.cfg.grid.half_width_km * 1000.0, self.cfg.grid.cell_m)
        thr = self.cfg.thresholds
        self.spec = ex.MetricSpec(
            lamax_thresholds=tuple(thr.lamax_db),
            psf_thresholds=tuple(thr.psf),
            disturbance_lamax_db=thr.disturbance_lamax_db,
            season=tuple(self.cfg.experiments.season),  # type: ignore[arg-type]
            coverage_receptors=tuple(
                self.receptors.names.index(n) for n in COVERAGE_RECEPTORS if n in self.receptors.names
            ),
            per_operation=self.per_operation,
        )
        self.prov_inputs = [
            Path(self.cfg.receptors),
            Path(self.cfg.footprint_points),
            Path(self.cfg.per_operation),
        ] + sorted((ROOT / "data" / "transcribed").glob("*.csv"))

    # ------------------------------------------------------------------
    def prov(self, **extra: Any) -> dict[str, Any]:
        return provenance(
            [self.config_path] + [Path(p) for p in self.cfg.scenarios], self.prov_inputs, self.cfg.seed, extra
        )

    def log(self, msg: str) -> None:
        print(f"[{time.time() - self.t0:8.1f}s] {msg}", flush=True)

    def levels(self, offset_db: float = 0.0, boom_scale: float = 1.0) -> ReceptorLevels:
        return receptor_levels(
            self.lib, self.receptors.names, self.receptors.r_km, self.receptors.azimuth_deg, offset_db, boom_scale
        )

    def base_params(self, scenario: str = "scenario_phase1") -> ScheduleParams:
        return self.scenarios[scenario].operations.to_params(self.scenarios[scenario].year)

    # ------------------------------------------------------------------
    def run_all(self) -> None:
        self.inputs_and_fits()
        self.verification()
        self.scenario_runs()
        self.e1_scaling()
        self.e2_clustering()
        self.e3_night()
        self.e4_closures()
        self.e5_coverage()
        self.e7_footprint_scaling()
        self.e6_sobol()
        self.e8_emissions()
        self.finish()

    # ------------------------------------------------------------------
    def inputs_and_fits(self) -> None:
        self.log("footprint fits")
        res = self.lib.residual_table()
        write_csv(res, self.out / "footprint_fit_residuals.csv", self.prov())
        fits = pd.DataFrame(
            [
                {
                    "event_type": et,
                    "metric": m,
                    "l0_db_at_1km": f.l0,
                    "alpha_db_per_km": f.alpha_per_km,
                    "n_points": f.r_km.size,
                    "rms_residual_db": f.rms_residual_db,
                }
                for (et, m), f in self.lib.noise.items()
            ]
        )
        write_csv(fits, self.out / "footprint_fits.csv", self.prov())
        write_csv(self.receptors.table, self.out / "receptors_used.csv", self.prov())
        write_csv(self.points, self.out / "footprint_points_used.csv", self.prov())
        for (et, m), f in self.lib.noise.items():
            self.values[f"fit_{et}_{m}_l0"] = f.l0
            self.values[f"fit_{et}_{m}_alpha"] = f.alpha_per_km
            self.values[f"fit_{et}_{m}_rms"] = f.rms_residual_db
            self.values[f"fit_{et}_{m}_npts"] = int(f.r_km.size)
        for et, bf in self.lib.booms.items():
            self.values[f"fit_{et}_boom_rms_db"] = (
                float(20 * np.sqrt(np.mean(bf.residual_log10**2))) if bf.residual_log10.size else 0.0
            )
        per_op_rows = [
            {"quantity": q, "event_type": et, "value": v} for q, d in self.per_operation.items() for et, v in d.items()
        ]
        write_csv(pd.DataFrame(per_op_rows), self.out / "per_operation_inputs.csv", self.prov())
        for q, d in self.emissions_per_op.items():
            for et, v in d.items():
                self.values[f"perop_{q}_{et}"] = v
        lv = self.levels()
        lv_tab = pd.DataFrame({"receptor": lv.names, "r_km": lv.r_km, "azimuth_deg": lv.azimuth_deg})
        for et in ("launch", "sh_landing", "ss_landing", "static_fire"):
            from .schedule import TYPE_CODE

            lv_tab[f"lamax_{et}"] = lv.lamax[TYPE_CODE[et]]
            lv_tab[f"sel_{et}"] = lv.sel[TYPE_CODE[et]]
            lv_tab[f"psf_{et}"] = lv.psf[TYPE_CODE[et]]
        write_csv(lv_tab, self.out / "receptor_levels.csv", self.prov())
        self.lv_tab = lv_tab

    # ------------------------------------------------------------------
    def verification(self) -> None:
        self.log("verification against the 2025 Tiered EA")
        tgt = self.tables["verification_targets"]
        vt = ex.verification_table(self.lib, self.grid, tgt, self.per_operation, self.emissions_per_op)
        write_csv(vt, self.out / "verification.csv", self.prov())
        for r in vt.to_dict("records"):
            key = re.sub(r"[^a-z0-9]+", "_", str(r["quantity"]).lower()).strip("_")
            self.values[f"verif_{key}"] = float(r["reproduced"])
        # grid convergence of the verification fields
        fine = Grid(self.pad_x, self.pad_y, self.cfg.grid.half_width_km * 1000.0, self.cfg.grid.convergence_cell_m)
        has_ss = ("static_fire_ss", "sel") in self.lib.noise
        f_coarse = ex.verification_fields(self.lib, self.grid, has_ss)
        f_fine = ex.verification_fields(self.lib, fine, has_ss)
        conv = []
        for name, thr in (("dnl", 65.0), ("cdnl", 60.0)):
            conv.append(
                {
                    "field": name,
                    "threshold": thr,
                    "cell_m_coarse": self.grid.cell_m,
                    "area_coarse_km2": self.grid.area_above(f_coarse[name], thr),
                    "cell_m_fine": fine.cell_m,
                    "area_fine_km2": fine.area_above(f_fine[name], thr),
                }
            )
        conv_df = pd.DataFrame(conv)
        conv_df["relative_difference"] = (conv_df["area_coarse_km2"] - conv_df["area_fine_km2"]) / conv_df[
            "area_fine_km2"
        ]
        write_csv(conv_df, self.out / "grid_convergence_verification.csv", self.prov())
        rel = conv_df["relative_difference"].to_numpy(float)
        area = conv_df["area_coarse_km2"].to_numpy(float)
        self.values["conv_dnl65_rel"] = float(rel[0])
        self.values["conv_cdnl60_rel"] = float(rel[1])
        self.values["verif_dnl65_area_km2"] = float(area[0])
        self.values["verif_cdnl60_area_km2"] = float(area[1])

    # ------------------------------------------------------------------
    def scenario_runs(self) -> None:
        self.log("scenario runs")
        n = self.cfg.monte_carlo.n_realizations
        frames = []
        for (name, sc), s in zip(self.scenarios.items(), self.seed.spawn(len(self.scenarios)), strict=True):
            params = sc.operations.to_params(sc.year)
            lv = self.levels(sc.footprints.level_offset_db, sc.footprints.boom_scale)
            df = ex.simulate(params, lv, self.spec, n, s, self.cfg.monte_carlo.n_workers)
            df["scenario"] = name
            df["hypothetical"] = sc.hypothetical
            frames.append(df)
            self.log(f"  {name}: {len(df)} realizations")
        raw = pd.concat(frames, ignore_index=True)
        write_parquet(raw, self.out / "scenarios_raw.parquet", self.prov())
        summ = ex.summarize(raw, ["scenario"])
        summ.to_csv(self.out / "scenarios_summary.csv")
        self.scenario_raw = raw
        # Monte Carlo standard-error check on headline quantities
        head = [
            c
            for c in raw.columns
            if c.startswith(("ldn|", "na100|", "closure_h", "disturbed|", "days90|"))
            and (c.endswith(COMMUNITY_EDGE) or c.endswith(HABITAT) or "|" not in c)
        ]
        se_rows = []
        for sc_key, g in raw.groupby("scenario"):
            sc_name = str(sc_key)
            for c in head:
                m = float(g[c].mean())
                se = float(g[c].std(ddof=1) / np.sqrt(len(g)))
                se_rows.append(
                    {"scenario": sc_name, "metric": c, "mean": m, "se": se, "se_fraction": abs(se / m) if m else np.nan}
                )
        se_df = pd.DataFrame(se_rows)
        write_csv(se_df, self.out / "monte_carlo_standard_errors.csv", self.prov())
        self.values["mc_max_se_fraction"] = float(np.nanmax(se_df["se_fraction"]))
        self.values["mc_n_realizations"] = n
        for sc_key, g in raw.groupby("scenario"):
            sc_name = str(sc_key)
            tag = sc_name.replace("scenario_", "")
            for metric in ("ldn", "cdnl", "na100", "na90", "days90", "boom2", "boom1"):
                col = f"{metric}|{COMMUNITY_EDGE}"
                if col in g:
                    self.values[f"{tag}_{metric}_portisabel"] = float(g[col].mean())
            for metric in ("disturbed", "longest_gap", "has_window"):
                col = f"{metric}|{HABITAT}"
                if col in g:
                    self.values[f"{tag}_{metric}_habitat"] = float(g[col].mean())
            self.values[f"{tag}_closure_h"] = float(g["closure_h"].mean())
            self.values[f"{tag}_closure_h_q975"] = float(g["closure_h"].quantile(0.975))
            self.values[f"{tag}_weekend_avail"] = float(g["weekend_avail"].mean())
            self.values[f"{tag}_summer_avail"] = float(g["summer_avail"].mean())
            self.values[f"{tag}_deluge_mgal"] = float(g["deluge_gal"].mean() / 1e6)
            self.values[f"{tag}_nox_tons"] = float(g["nox_tons"].mean())
            self.values[f"{tag}_co2e_kt"] = float(g["co2e_t"].mean() / 1e3)
            self.values[f"{tag}_feasible"] = float(g["feasible"].mean())
            self.values[f"{tag}_n_launch"] = float(g["n_launch"].mean())

    # ------------------------------------------------------------------
    def e1_scaling(self) -> None:
        self.log("E1 scaling with cadence")
        e = self.cfg.experiments
        base = replace(self.base_params("scenario_phase1"), model="poisson", weekday_only=False, night_fraction=0.5)
        lv = self.levels()
        cases: list[tuple[dict[str, float | str], ScheduleParams, ReceptorLevels]] = []
        for n_l in e.cadences:
            for frac in e.external_landing_fractions:
                p = replace(base, n_launches=n_l, n_external_landings=frac * n_l)
                cases.append(({"cadence": n_l, "external_fraction": frac}, p, lv))
        raw = ex.sweep(
            cases, self.spec, self.cfg.monte_carlo.n_realizations, self.seed.spawn(1)[0], self.cfg.monte_carlo.n_workers
        )
        write_parquet(raw, self.out / "e1_raw.parquet", self.prov())
        summ = ex.summarize(raw, ["cadence", "external_fraction"])
        summ.to_csv(self.out / "e1_summary.csv")
        # aggregation errors relative to the reference cadence
        ref = e.reference_cadence
        rows = []
        metrics = {
            "deluge_gal": "additive",
            "nox_tons": "additive",
            f"ldn|{COMMUNITY_EDGE}": "energy",
            f"na100|{COMMUNITY_EDGE}": "count",
            f"days90|{COMMUNITY_EDGE}": "exceedance_days",
            f"days90|{COMMUNITY_SPI}": "exceedance_days",
            f"boom2|{COMMUNITY_EDGE}": "count",
            f"disturbed|{HABITAT}": "coverage",
            f"disturbed|{REFUGE}": "coverage",
            f"has_window|{HABITAT}": "window",
            f"longest_gap|{HABITAT}": "gap",
            "closure_h": "closure",
            "weekend_avail": "availability",
        }
        sub0 = raw[(raw["external_fraction"] == 0.0)]
        for metric, cls in metrics.items():
            if metric not in sub0:
                continue
            g = sub0.groupby("cadence")[metric]
            mean = g.mean()
            se = g.std(ddof=1) / np.sqrt(g.size())
            y0 = float(mean.loc[ref])
            for n_l in e.cadences:
                yn = mean.loc[n_l]
                if cls == "energy":
                    pred = y0 + 10 * np.log10(n_l / ref)
                    err = pred - yn
                    rel = np.nan
                else:
                    pred = y0 * n_l / ref
                    err = pred - yn
                    rel = err / yn if yn else np.nan
                rows.append(
                    {
                        "metric": metric,
                        "class": cls,
                        "cadence": n_l,
                        "simulated": yn,
                        "se": se.loc[n_l],
                        "linear_extrapolation": pred,
                        "error": err,
                        "relative_error": rel,
                    }
                )
        err_df = pd.DataFrame(rows)
        write_csv(err_df, self.out / "e1_aggregation_errors.csv", self.prov())
        # bootstrap intervals on the relative error at the top cadence
        rng = np.random.default_rng(self.seed.spawn(1)[0])
        ci_rows = []
        for metric, cls in metrics.items():
            if metric not in sub0 or cls == "energy":
                continue
            a = sub0[sub0["cadence"] == ref][metric].to_numpy(float)
            for n_l in e.cadences:
                b = sub0[sub0["cadence"] == n_l][metric].to_numpy(float)
                boots = []
                for _ in range(400):
                    ya = rng.choice(a, a.size).mean()
                    yb = rng.choice(b, b.size).mean()
                    boots.append((ya * n_l / ref - yb) / yb if yb else np.nan)
                lo, hi = np.nanpercentile(boots, [2.5, 97.5])
                ci_rows.append({"metric": metric, "cadence": n_l, "rel_error_lo": lo, "rel_error_hi": hi})
        ci = pd.DataFrame(ci_rows)
        write_csv(ci, self.out / "e1_aggregation_error_intervals.csv", self.prov())
        self.e1_raw = raw
        top = err_df[err_df["cadence"] == max(e.cadences)].set_index("metric")
        for metric_key, row in top.to_dict("index").items():
            key = re.sub(r"[^a-z0-9]+", "_", str(metric_key).lower()).strip("_")
            self.values[f"e1_relerr_{key}"] = (
                float(row["relative_error"]) if np.isfinite(row["relative_error"]) else float("nan")
            )
            self.values[f"e1_err_{key}"] = float(row["error"])
            self.values[f"e1_sim_{key}"] = float(row["simulated"])
        for n_l in e.cadences:
            tag = f"n{int(n_l)}"
            gc = sub0[sub0["cadence"] == n_l]
            self.values[f"e1_{tag}_disturbed_habitat"] = float(gc[f"disturbed|{HABITAT}"].mean())
            self.values[f"e1_{tag}_ldn_portisabel"] = float(gc[f"ldn|{COMMUNITY_EDGE}"].mean())
            self.values[f"e1_{tag}_days90_spi"] = float(gc[f"days90|{COMMUNITY_SPI}"].mean())
            self.values[f"e1_{tag}_na100_portisabel"] = float(gc[f"na100|{COMMUNITY_EDGE}"].mean())
            self.values[f"e1_{tag}_boom2_portisabel"] = float(gc[f"boom2|{COMMUNITY_EDGE}"].mean())
            self.values[f"e1_{tag}_closure_h"] = float(gc["closure_h"].mean())
        # closed-form coverage check at the habitat receptor
        rate = sub0.groupby("cadence")[f"n_dist|{HABITAT}"].mean() / (e.season[1] - e.season[0])
        phi_cf = cf.poisson_coverage(rate.to_numpy(float), self.spec.recovery_days)
        phi_sim = sub0.groupby("cadence")[f"disturbed|{HABITAT}"].mean().to_numpy(float)
        chk = pd.DataFrame(
            {
                "cadence": rate.index,
                "rate_per_day": rate.to_numpy(float),
                "coverage_closed_form": phi_cf,
                "coverage_simulated": phi_sim,
                "coverage_linear": cf.linear_coverage(rate.to_numpy(float), self.spec.recovery_days),
            }
        )
        write_csv(chk, self.out / "e1_coverage_closed_form_check.csv", self.prov())
        self.values["e1_coverage_cf_maxabsdiff"] = float(np.max(np.abs(phi_cf - phi_sim)))
        # external landings effect at Phase 2 cadence
        top_n = max(e.cadences)
        g0 = raw[(raw["cadence"] == top_n) & (raw["external_fraction"] == 0.0)]
        g1 = raw[(raw["cadence"] == top_n) & (raw["external_fraction"] == 1.0)]
        self.values["e1_ext_boom1_portisabel_0"] = float(g0[f"boom1|{COMMUNITY_EDGE}"].mean())
        self.values["e1_ext_boom1_portisabel_1"] = float(g1[f"boom1|{COMMUNITY_EDGE}"].mean())
        self.values["e1_ext_days90_spi_0"] = float(g0[f"days90|{COMMUNITY_SPI}"].mean())
        self.values["e1_ext_days90_spi_1"] = float(g1[f"days90|{COMMUNITY_SPI}"].mean())
        self.values["e1_ext_ldn_portisabel_0"] = float(g0[f"ldn|{COMMUNITY_EDGE}"].mean())
        self.values["e1_ext_ldn_portisabel_1"] = float(g1[f"ldn|{COMMUNITY_EDGE}"].mean())

    # ------------------------------------------------------------------
    def e2_clustering(self) -> None:
        self.log("E2 temporal model at fixed count")
        e = self.cfg.experiments
        lv = self.levels()
        cases: list[tuple[dict[str, float | str], ScheduleParams, ReceptorLevels]] = []
        for n_l in (100.0, max(e.cadences)):
            for model in e.temporal_models:
                if model == "poisson":
                    p = replace(
                        self.base_params("scenario_phase1"),
                        model="poisson",
                        weekday_only=False,
                        n_launches=n_l,
                        night_fraction=0.5,
                    )
                elif model == "weekday_renewal":
                    p = replace(self.base_params("scenario_phase2"), n_launches=n_l)
                else:
                    p = replace(self.base_params("scenario_phase2_clustered"), n_launches=n_l)
                cases.append(({"cadence": n_l, "model": model}, p, lv))
        raw = ex.sweep(
            cases, self.spec, self.cfg.monte_carlo.n_realizations, self.seed.spawn(1)[0], self.cfg.monte_carlo.n_workers
        )
        write_parquet(raw, self.out / "e2_raw.parquet", self.prov())
        ex.summarize(raw, ["cadence", "model"]).to_csv(self.out / "e2_summary.csv")
        metrics = [
            f"disturbed|{HABITAT}",
            f"longest_gap|{HABITAT}",
            f"has_window|{HABITAT}",
            f"days90|{COMMUNITY_SPI}",
            f"days90|{COMMUNITY_EDGE}",
            "weekend_avail",
            "summer_avail",
            "closure_h",
            f"na100night|{COMMUNITY_EDGE}",
        ]
        rng = np.random.default_rng(self.seed.spawn(1)[0])
        rows = []
        for n_l in sorted(_num(v) for v in raw["cadence"].unique()):
            ref = raw[(raw["cadence"] == n_l) & (raw["model"] == "poisson")]
            for model in e.temporal_models:
                g = raw[(raw["cadence"] == n_l) & (raw["model"] == model)]
                for m in metrics:
                    a = ref[m].to_numpy(float)
                    b = g[m].to_numpy(float)
                    diff = b.mean() - a.mean()
                    boots = [rng.choice(b, b.size).mean() - rng.choice(a, a.size).mean() for _ in range(400)]
                    lo, hi = np.percentile(boots, [2.5, 97.5])
                    rows.append(
                        {
                            "cadence": n_l,
                            "model": model,
                            "metric": m,
                            "mean": b.mean(),
                            "diff_vs_poisson": diff,
                            "ci_lo": lo,
                            "ci_hi": hi,
                            "feasible_fraction": g["feasible"].mean(),
                        }
                    )
        diff_df = pd.DataFrame(rows)
        write_csv(diff_df, self.out / "e2_differences.csv", self.prov())
        top = max(e.cadences)
        for model in e.temporal_models:
            g = raw[(raw["cadence"] == top) & (raw["model"] == model)]
            tag = model.replace("_renewal", "")
            self.values[f"e2_{tag}_disturbed_habitat"] = float(g[f"disturbed|{HABITAT}"].mean())
            self.values[f"e2_{tag}_longest_gap_habitat_h"] = float(g[f"longest_gap|{HABITAT}"].mean() * 24)
            self.values[f"e2_{tag}_days90_spi"] = float(g[f"days90|{COMMUNITY_SPI}"].mean())
            self.values[f"e2_{tag}_weekend_avail"] = float(g["weekend_avail"].mean())
            self.values[f"e2_{tag}_summer_avail"] = float(g["summer_avail"].mean())
            self.values[f"e2_{tag}_has_window_habitat"] = float(g[f"has_window|{HABITAT}"].mean())
            self.values[f"e2_{tag}_feasible"] = float(g["feasible"].mean())
        # Ranking stability: does the ordering of cadences by metric change with the model?
        rank_rows = []
        for m in metrics:
            for model in e.temporal_models:
                gm = raw[raw["model"] == model].groupby("cadence")[m].mean()
                rank_rows.append({"metric": m, "model": model, **{f"n{int(_num(k))}": v for k, v in gm.items()}})
        write_csv(pd.DataFrame(rank_rows), self.out / "e2_cadence_means_by_model.csv", self.prov())
        self.e2_raw = raw

    # ------------------------------------------------------------------
    def e3_night(self) -> None:
        self.log("E3 night fraction")
        e = self.cfg.experiments
        lv = self.levels()
        cases: list[tuple[dict[str, float | str], ScheduleParams, ReceptorLevels]] = []
        for n_l in (float(min(e.cadences)), float(max(e.cadences))):
            for fn in e.night_fractions:
                p = replace(
                    self.base_params("scenario_phase1"),
                    model="poisson",
                    weekday_only=False,
                    n_launches=n_l,
                    night_fraction=fn,
                )
                cases.append(({"cadence": n_l, "night_fraction": fn}, p, lv))
        raw = ex.sweep(
            cases, self.spec, self.cfg.monte_carlo.n_realizations, self.seed.spawn(1)[0], self.cfg.monte_carlo.n_workers
        )
        write_parquet(raw, self.out / "e3_raw.parquet", self.prov())
        ex.summarize(raw, ["cadence", "night_fraction"]).to_csv(self.out / "e3_summary.csv")
        rows = []
        for n_l in sorted(_num(v) for v in raw["cadence"].unique()):
            gm = raw[raw["cadence"] == n_l].groupby("night_fraction")[f"ldn|{COMMUNITY_EDGE}"].mean()
            base = float(gm.loc[0.0])
            for fn_key, val in gm.items():
                fn = _num(fn_key)
                rows.append(
                    {
                        "cadence": n_l,
                        "night_fraction": fn,
                        "ldn_sim": val,
                        "ldn_increment_sim": val - base,
                        "ldn_increment_closed_form": cf.ldn_night_term(0.0, float(fn)),
                    }
                )
        chk = pd.DataFrame(rows)
        write_csv(chk, self.out / "e3_night_closed_form_check.csv", self.prov())
        self.values["e3_night_cf_maxabsdiff_db"] = float(
            np.max(np.abs(chk["ldn_increment_sim"] - chk["ldn_increment_closed_form"]))
        )
        self.values["cf_count_term_25_456_db"] = cf.ldn_count_term(25, 456)
        self.values["cf_night_term_012_05_db"] = cf.ldn_night_term(0.12, 0.5)
        self.values["cf_total_25_456_db"] = cf.ldn_increment(25, 456, 0.12, 0.5)
        self.values["cf_count_term_25_100_db"] = cf.ldn_count_term(25, 100)
        self.values["cf_total_25_100_db"] = cf.ldn_increment(25, 100, 0.12, 0.5)
        # night number-above at the top cadence
        top = raw[raw["cadence"] == max(e.cadences)]
        for fn in (0.0, 0.12, 0.5):
            if fn in set(top["night_fraction"]):
                self.values[f"e3_na100night_portisabel_fn{int(fn * 100):02d}"] = float(
                    top[top["night_fraction"] == fn][f"na100night|{COMMUNITY_EDGE}"].mean()
                )
        # gridded L_dn fields at 25/100/456 with the scenario night shares
        fields = {}
        areas = []
        for n_l, fn, label in (
            (25.0, 0.12, "existing"),
            (100.0, 0.5, "phase1"),
            (456.0, 0.5, "phase2"),
            (456.0, 0.12, "phase2_fn012"),
        ):
            n_day = np.zeros(5)
            n_night = np.zeros(5)
            from .schedule import TYPE_CODE

            for et in ("launch", "sh_landing", "ss_landing"):
                n_day[TYPE_CODE[et]] = n_l * (1 - fn)
                n_night[TYPE_CODE[et]] = n_l * fn
            n_day[TYPE_CODE["static_fire"]] = n_l * self.base_params("scenario_phase1").static_fires_per_launch
            f = ex.grid_fields(self.lib, self.grid, n_day, n_night)
            fields[label] = f
            for thr in (55.0, 60.0, 65.0, 70.0):
                areas.append(
                    {
                        "case": label,
                        "cadence": n_l,
                        "night_fraction": fn,
                        "field": "ldn",
                        "threshold_db": thr,
                        "area_km2": self.grid.area_above(f["ldn"], thr),
                        "max_extent_km": self.grid.max_extent_km(f["ldn"], thr),
                    }
                )
            for thr in (55.0, 60.0, 65.0):
                areas.append(
                    {
                        "case": label,
                        "cadence": n_l,
                        "night_fraction": fn,
                        "field": "cdnl",
                        "threshold_db": thr,
                        "area_km2": self.grid.area_above(f["cdnl"], thr),
                        "max_extent_km": self.grid.max_extent_km(f["cdnl"], thr),
                    }
                )
        area_df = pd.DataFrame(areas)
        write_csv(area_df, self.out / "e3_contour_areas.csv", self.prov())
        for r in area_df.to_dict("records"):
            self.values[f"area_{r['case']}_{r['field']}{int(r['threshold_db'])}_km2"] = float(r["area_km2"])
            self.values[f"extent_{r['case']}_{r['field']}{int(r['threshold_db'])}_km"] = float(r["max_extent_km"])
        lon, lat = self.grid.lonlat()
        ds = xr.Dataset(
            {f"{lab}_{k}": (("y", "x"), v.astype(np.float32)) for lab, f in fields.items() for k, v in f.items()},
            coords={"x": self.grid.x, "y": self.grid.y, "lon": (("y", "x"), lon), "lat": (("y", "x"), lat)},
            attrs={
                "crs": "EPSG:32614",
                "pad_x": self.pad_x,
                "pad_y": self.pad_y,
                "cell_m": self.grid.cell_m,
                "note": "Expected-count fields; L_dn from fitted SEL footprints; CDNL from the FAA Desk Reference form",
            },
        )
        write_netcdf(ds, self.out / "ldn_fields.nc", self.prov())
        # population within the DNL 65 contour from census blocks
        blocks = pd.read_csv(ROOT / "data" / "processed" / "blocks_population.csv")
        from .spatial import lonlat_to_utm

        bx, by = lonlat_to_utm(blocks["lon"].to_numpy(float), blocks["lat"].to_numpy(float))
        br = np.hypot(bx - self.pad_x, by - self.pad_y) / 1000.0
        pop_rows = []
        for label, f in fields.items():
            ext = self.grid.max_extent_km(f["ldn"], 65.0)
            pop_rows.append(
                {
                    "case": label,
                    "ldn65_max_extent_km": ext,
                    "population_within_extent": int(blocks.loc[br <= ext, "pop_2020"].sum()),
                }
            )
            self.values[f"pop_ldn65_{label}"] = int(blocks.loc[br <= ext, "pop_2020"].sum())
        write_csv(pd.DataFrame(pop_rows), self.out / "e3_population_within_ldn65.csv", self.prov())

    # ------------------------------------------------------------------
    def e4_closures(self) -> None:
        self.log("E4 closure caps")
        e = self.cfg.experiments
        lv = self.levels()
        cases: list[tuple[dict[str, float | str], ScheduleParams, ReceptorLevels]] = []
        for dist_name, cc in e.closure_distributions.items():
            for n_l in e.cadences:
                p = replace(
                    self.base_params("scenario_phase1"),
                    model="poisson",
                    weekday_only=False,
                    n_launches=n_l,
                    night_fraction=0.5,
                    closures=cc.to_model(),
                )
                cases.append(({"distribution": dist_name, "cadence": n_l, "mean_h": cc.launch_mean_h}, p, lv))
        # Phase 2 weekday model with the record distribution and the PEA weekend/holiday measures applied
        rec = e.closure_distributions[next(iter(e.closure_distributions))]
        p2 = replace(self.base_params("scenario_phase2"), closures=rec.to_model())
        cases.append(
            (
                {
                    "distribution": "record_weekday_phase2",
                    "cadence": float(max(e.cadences)),
                    "mean_h": rec.launch_mean_h,
                },
                p2,
                lv,
            )
        )
        raw = ex.sweep(
            cases, self.spec, self.cfg.monte_carlo.n_realizations, self.seed.spawn(1)[0], self.cfg.monte_carlo.n_workers
        )
        spec_pea = replace(self.spec, apply_pea_measures=True)
        raw_pea = ex.sweep(
            [
                (
                    {
                        "distribution": "record_weekday_phase2_pea_measures",
                        "cadence": float(max(e.cadences)),
                        "mean_h": rec.launch_mean_h,
                    },
                    p2,
                    lv,
                )
            ],
            spec_pea,
            self.cfg.monte_carlo.n_realizations,
            self.seed.spawn(1)[0],
            self.cfg.monte_carlo.n_workers,
        )
        raw = pd.concat([raw, raw_pea], ignore_index=True)
        write_parquet(raw, self.out / "e4_raw.parquet", self.prov())
        rows = []
        for (dist_key, n_key), g in raw.groupby(["distribution", "cadence"]):
            dist, n_l = str(dist_key), _num(n_key)
            h = g["closure_anom_h"].to_numpy(float)
            h0 = g["closure_h"].to_numpy(float)
            row = {
                "distribution": dist,
                "cadence": n_l,
                "mean_h": g["mean_h"].iloc[0],
                "closure_h_mean": h0.mean(),
                "closure_h_sd": h0.std(ddof=1),
                "closure_anom_h_mean": h.mean(),
                "q975_h": np.quantile(h, 0.975),
                "weekend_avail": g["weekend_avail"].mean(),
                "summer_avail": g["summer_avail"].mean(),
                "n_windows": g["n_closure_windows"].mean(),
            }
            for cap in e.closure_caps_h:
                k = int((h > cap).sum())
                pr = k / h.size
                row[f"p_exceed_{int(cap)}"] = pr
                row[f"p_exceed_{int(cap)}_se"] = np.sqrt(pr * (1 - pr) / h.size)
                row[f"n_exceed_{int(cap)}"] = k
                row[f"budget_per_launch_{int(cap)}_h"] = cap / n_l
            rows.append(row)
        tab = pd.DataFrame(rows)
        write_csv(tab, self.out / "e4_cap_probabilities.csv", self.prov())
        top = max(e.cadences)
        for r in tab.to_dict("records"):
            tag = f"{r['distribution']}_n{int(r['cadence'])}"
            self.values[f"e4_{tag}_closure_h"] = float(r["closure_anom_h_mean"])
            self.values[f"e4_{tag}_weekend_avail"] = float(r["weekend_avail"])
            self.values[f"e4_{tag}_summer_avail"] = float(r["summer_avail"])
            for cap in e.closure_caps_h:
                self.values[f"e4_{tag}_pexceed{int(cap)}"] = float(r[f"p_exceed_{int(cap)}"])
        for cap in e.closure_caps_h:
            self.values[f"budget_per_launch_{int(cap)}_n{int(top)}_h"] = cap / top
            self.values[f"budget_per_launch_{int(cap)}_n100_h"] = cap / 100
        # mean closure hours per launch in the Flight 1-4 record
        self.values["record_mean_launch_closure_h"] = float(
            self.tables["closure_record"]
            .query("flight.str.startswith('Flight')", engine="python")["launch_hours"]
            .mean()
        )
        self.e4_raw = raw

    # ------------------------------------------------------------------
    def e5_coverage(self) -> None:
        self.log("E5 recovery time and window response surfaces")
        e = self.cfg.experiments
        lv = self.levels()
        n = max(self.cfg.monte_carlo.n_realizations // 5, 50)
        frames = []
        k_hab = self.receptors.names.index(HABITAT)
        k_comm = self.receptors.names.index(COMMUNITY_EDGE)
        for n_l in e.cadences:
            for model in ("poisson", "weekday_renewal"):
                if model == "poisson":
                    p = replace(
                        self.base_params("scenario_phase1"),
                        model="poisson",
                        weekday_only=False,
                        n_launches=n_l,
                        night_fraction=0.5,
                    )
                else:
                    p = replace(self.base_params("scenario_phase2"), n_launches=n_l)
                for k, rname in ((k_hab, HABITAT), (k_comm, COMMUNITY_EDGE)):
                    df = ex.e5_coverage(
                        p,
                        lv,
                        k,
                        e.recovery_days,
                        e.window_days,
                        (e.season[0], e.season[1]),
                        self.cfg.thresholds.disturbance_lamax_db,
                        n,
                        self.seed.spawn(1)[0],
                    )
                    df["cadence"] = n_l
                    df["model"] = model
                    df["receptor"] = rname
                    frames.append(df)
        raw = pd.concat(frames, ignore_index=True)
        write_parquet(raw, self.out / "e5_raw.parquet", self.prov())
        g = raw.groupby(["receptor", "model", "cadence", "recovery_days"])
        summ = g.mean(numeric_only=True).drop(columns=["realization"])
        summ["disturbed_se"] = g["disturbed"].std(ddof=1) / np.sqrt(g.size())
        rate = summ["n_events"] / (e.season[1] - e.season[0])
        summ["coverage_closed_form"] = cf.poisson_coverage(
            rate.to_numpy(float), summ.index.get_level_values("recovery_days").to_numpy(float)
        )
        summ["coverage_linear"] = cf.linear_coverage(
            rate.to_numpy(float), summ.index.get_level_values("recovery_days").to_numpy(float)
        )
        summ["lambda_tau"] = rate.to_numpy(float) * summ.index.get_level_values("recovery_days").to_numpy(float)
        summ.to_csv(self.out / "e5_surfaces.csv")
        self.e5_summary = summ
        # regime boundary: smallest recovery time at which the linear extrapolation overstates coverage by >10 percent
        top = max(e.cadences)
        for model in ("poisson", "weekday_renewal"):
            s = summ.xs((HABITAT, model, top), level=("receptor", "model", "cadence"))
            over = (s["coverage_linear"] - s["disturbed"]) / s["disturbed"]
            crossing = s.index[(over > 0.10).to_numpy()]
            self.values[f"e5_{model}_tau_10pct_h"] = float(crossing.min() * 24) if len(crossing) else float("nan")
            self.values[f"e5_{model}_rate_per_day_n{int(top)}"] = float(
                rate.xs((HABITAT, model, top), level=("receptor", "model", "cadence")).iloc[0]
            )
            for tau in (1 / 24.0, 1.0):
                nearest = s.index[int(np.argmin(np.abs(s.index.to_numpy(float) - tau)))]
                r = s.loc[nearest]
                self.values[f"e5_{model}_disturbed_tau{int(round(tau * 24))}h"] = float(r["disturbed"])
                self.values[f"e5_{model}_linear_tau{int(round(tau * 24))}h"] = float(r["coverage_linear"])
                for w in (0.25, 1.0, 7.0):
                    col = f"has_window|{w:g}"
                    if col in s:
                        self.values[f"e5_{model}_haswindow_w{int(round(w * 24))}h_tau{int(round(tau * 24))}h"] = float(
                            r[col]
                        )

    # ------------------------------------------------------------------
    def e7_footprint_scaling(self) -> None:
        self.log("E7 footprint scaling (hypothetical vehicle variants)")
        e = self.cfg.experiments
        rows = []
        from .schedule import TYPE_CODE

        for n_l, fn in ((25.0, 0.12), (456.0, 0.5)):
            n_day = np.zeros(5)
            n_night = np.zeros(5)
            for et in ("launch", "sh_landing", "ss_landing"):
                n_day[TYPE_CODE[et]] = n_l * (1 - fn)
                n_night[TYPE_CODE[et]] = n_l * fn
            n_day[TYPE_CODE["static_fire"]] = n_l * 0.2
            for off in e.footprint_offsets_db:
                for sc in e.boom_scales:
                    f = ex.grid_fields(self.lib, self.grid, n_day, n_night, off, sc)
                    rows.append(
                        {
                            "cadence": n_l,
                            "night_fraction": fn,
                            "level_offset_db": off,
                            "boom_scale": sc,
                            "hypothetical": off != 0 or sc != 1,
                            "ldn65_area_km2": self.grid.area_above(f["ldn"], 65.0),
                            "ldn65_extent_km": self.grid.max_extent_km(f["ldn"], 65.0),
                            "cdnl60_area_km2": self.grid.area_above(f["cdnl"], 60.0),
                            "cdnl60_extent_km": self.grid.max_extent_km(f["cdnl"], 60.0),
                            "na100_ge10_area_km2": self.grid.area_above(f["na100"], 10.0),
                            "boom2_ge10_area_km2": self.grid.area_above(f["boom2"], 10.0),
                            "boom4_ge10_area_km2": self.grid.area_above(f["boom4"], 10.0),
                        }
                    )
        tab = pd.DataFrame(rows)
        write_csv(tab, self.out / "e7_footprint_scaling.csv", self.prov())
        for r in tab.to_dict("records"):
            tag = f"n{int(r['cadence'])}_off{int(r['level_offset_db'])}_bs{int(round(r['boom_scale'] * 10))}"
            self.values[f"e7_{tag}_ldn65_area"] = float(r["ldn65_area_km2"])
            self.values[f"e7_{tag}_cdnl60_area"] = float(r["cdnl60_area_km2"])
            self.values[f"e7_{tag}_boom2_area"] = float(r["boom2_ge10_area_km2"])
        # grid convergence of the baseline Phase 2 field
        fine = Grid(self.pad_x, self.pad_y, self.cfg.grid.half_width_km * 1000.0, self.cfg.grid.convergence_cell_m)
        n_day = np.zeros(5)
        n_night = np.zeros(5)
        for et in ("launch", "sh_landing", "ss_landing"):
            n_day[TYPE_CODE[et]] = 456 * 0.5
            n_night[TYPE_CODE[et]] = 456 * 0.5
        n_day[TYPE_CODE["static_fire"]] = 456 * 0.2
        fc = ex.grid_fields(self.lib, self.grid, n_day, n_night)
        ff = ex.grid_fields(self.lib, fine, n_day, n_night)
        conv = []
        for name, thr in (("ldn", 65.0), ("ldn", 60.0), ("cdnl", 60.0), ("na100", 10.0)):
            ac = self.grid.area_above(fc[name], thr)
            af = fine.area_above(ff[name], thr)
            conv.append(
                {
                    "field": name,
                    "threshold": thr,
                    "area_coarse_km2": ac,
                    "area_fine_km2": af,
                    "relative_difference": (ac - af) / af if af else np.nan,
                }
            )
        conv_df = pd.DataFrame(conv)
        write_csv(conv_df, self.out / "grid_convergence_phase2.csv", self.prov())
        self.values["conv_phase2_ldn65_rel"] = float(conv_df["relative_difference"].to_numpy(float)[0])
        self.values["conv_phase2_max_rel"] = float(np.nanmax(np.abs(conv_df["relative_difference"])))

    # ------------------------------------------------------------------
    def e6_sobol(self) -> None:
        self.log("E6 Sobol sensitivity")
        s = self.cfg.sobol
        factors = [Factor(k, v.low, v.high, v.source_rank, v.note) for k, v in s.factors.items()]
        write_csv(pd.DataFrame([f.__dict__ for f in factors]), self.out / "e6_factor_ranges.csv", self.prov())
        base = replace(self.base_params("scenario_phase2_clustered"), model="campaign", weekday_only=True)
        tables = []
        for n_l, sd in zip(s.cadences, self.seed.spawn(len(s.cadences)), strict=True):
            ctx = ex.SobolContext(
                base=base,
                lib=self.lib,
                receptors=self.receptors,
                community_receptor=COMMUNITY_EDGE,
                habitat_receptor=HABITAT,
                spec=self.spec,
                factors=factors,
                cadence=float(n_l),
            )
            si, y = ex.e6_sobol(ctx, s.n_base, s.n_boot, sd, self.cfg.monte_carlo.n_workers)
            tables.append(si)
            write_parquet(y, self.out / f"e6_samples_n{int(n_l)}.parquet", self.prov())
            self.log(f"  cadence {n_l}: {len(y)} evaluations")
        si_all = pd.concat(tables, ignore_index=True)
        write_csv(si_all, self.out / "e6_sobol_indices.csv", self.prov())
        dom = si_all[si_all["ST"] > s.report_threshold].sort_values(
            ["cadence", "output", "ST"], ascending=[True, True, False]
        )
        write_csv(dom, self.out / "e6_dominant_factors.csv", self.prov())
        for (n_key, out_key), g in si_all.groupby(["cadence", "output"]):
            n_l, out = _num(n_key), str(out_key)
            top = g.sort_values("ST", ascending=False).iloc[0]
            self.values[f"e6_top_{out}_n{int(n_l)}"] = str(top["factor"])
            self.values[f"e6_topst_{out}_n{int(n_l)}"] = float(top["ST"])
            self.values[f"e6_ndom_{out}_n{int(n_l)}"] = int((g["ST"] > s.report_threshold).sum())
        self.values["e6_n_evaluations"] = int(s.n_base * (len(factors) + 2))
        self.values["e6_n_factors"] = len(factors)

    # ------------------------------------------------------------------
    def e8_emissions(self) -> None:
        self.log("E8 fleet stratospheric emissions")
        rng = np.random.default_rng(self.seed.spawn(1)[0])
        tab = ex.e8_emissions(self.cfg, rng)
        write_csv(tab, self.out / "e8_emissions.csv", self.prov())
        for r in tab.to_dict("records"):
            tag = f"{r['site']}_{r['species']}"
            self.values[f"e8_{tag}_median_t"] = float(r["median_t"])
            self.values[f"e8_{tag}_p05_t"] = float(r["p05_t"])
            self.values[f"e8_{tag}_p95_t"] = float(r["p95_t"])
        bc = tab[(tab["site"] == "fleet_total") & (tab["species"] == "BC")].iloc[0]
        for sc in self.cfg.emissions.comparison_scenarios_bc_gg_per_yr:
            self.values[f"e8_bc_fleet_over_{int(sc)}gg_median"] = float(bc["median_t"] / (sc * 1000.0))
        self.values["e8_fleet_launches"] = float(sum(s.launches_per_year for s in self.cfg.emissions.sites))

    # ------------------------------------------------------------------
    def finish(self) -> None:
        self.log("writing summary and LaTeX values")
        self.values["pipeline_name"] = self.cfg.name
        self.values["seed"] = self.cfg.seed
        self.values["n_receptors"] = len(self.receptors.names)
        self.values["grid_cell_m"] = self.cfg.grid.cell_m
        self.values["grid_half_width_km"] = self.cfg.grid.half_width_km
        self.values["runtime_s"] = time.time() - self.t0
        SUMMARY_DIR.mkdir(parents=True, exist_ok=True)
        write_json(self.values, self.out / "summary.json", self.prov())
        if self.cfg.name == "paper":
            write_json(self.values, SUMMARY_DIR / "summary.json", self.prov())
            for name in (
                "verification.csv",
                "e1_aggregation_errors.csv",
                "e2_differences.csv",
                "e4_cap_probabilities.csv",
                "e6_dominant_factors.csv",
                "e8_emissions.csv",
                "footprint_fits.csv",
                "e3_contour_areas.csv",
                "monte_carlo_standard_errors.csv",
                "grid_convergence_phase2.csv",
                "e7_footprint_scaling.csv",
            ):
                src = self.out / name
                if src.exists():
                    (SUMMARY_DIR / name).write_bytes(src.read_bytes())
            GENERATED.mkdir(parents=True, exist_ok=True)
            (GENERATED / "values.tex").write_text(values_tex(self.values), encoding="utf-8")
        else:
            GENERATED.mkdir(parents=True, exist_ok=True)
            (GENERATED / f"values_{self.cfg.name}.tex").write_text(values_tex(self.values), encoding="utf-8")
        self.log("done")


def _num(x: object) -> float:
    return float(str(x))


def macro_name(key: str) -> str:
    parts = re.split(r"[^A-Za-z0-9]+", key)
    name = "val" + "".join(p[:1].upper() + p[1:] for p in parts if p)
    return re.sub(r"\d", lambda m: "ZOTTFFSSEN"[int(m.group())], name)


def fmt(v: Any) -> str:
    if isinstance(v, str):
        return v.replace("_", r"\_")
    if isinstance(v, bool):
        return "1" if v else "0"
    if isinstance(v, int):
        return f"{v:d}"
    if isinstance(v, float):
        if not np.isfinite(v):
            return "n/a"
        a = abs(v)
        if a == 0:
            return "0"
        if a >= 1000:
            return f"{v:,.0f}"
        if a >= 100:
            return f"{v:.0f}"
        if a >= 10:
            return f"{v:.1f}"
        if a >= 1:
            return f"{v:.2f}"
        if a >= 0.01:
            return f"{v:.3f}"
        return f"{v:.2e}"
    return str(v)


def values_tex(values: dict[str, Any]) -> str:
    lines = ["% Generated by launchcum.pipeline; do not edit by hand."]
    for k in sorted(values):
        lines.append(f"\\newcommand{{\\{macro_name(k)}}}{{{fmt(values[k])}}}")
    return "\n".join(lines) + "\n"


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(prog="launchcum")
    sub = ap.add_subparsers(dest="cmd", required=True)
    run = sub.add_parser("run", help="run the full pipeline")
    run.add_argument("--config", type=Path, required=True)
    run.add_argument("--only", nargs="*", default=None, help="subset of stages, e.g. verification e1_scaling")
    args = ap.parse_args(argv)
    r = Runner(args.config)
    if args.only:
        r.inputs_and_fits()
        for stage in args.only:
            getattr(r, stage)()
        r.finish()
    else:
        r.run_all()
    print(
        json.dumps(
            {k: v for k, v in r.values.items() if k.startswith(("verif_", "cf_", "e1_relerr"))}, indent=1, default=str
        )
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
