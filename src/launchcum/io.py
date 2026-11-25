"""Configuration schemas, provenance sidecars and result writers."""

from __future__ import annotations

import hashlib
import json
import platform
import subprocess
from dataclasses import asdict
from datetime import UTC, datetime
from importlib import metadata
from pathlib import Path
from typing import Any, Literal

import numpy as np
import pandas as pd
import xarray as xr
import yaml
from pydantic import BaseModel, ConfigDict, Field, model_validator

from .schedule import ClosureModel, ScheduleParams

PACKAGES_FOR_PROVENANCE = (
    "numpy",
    "scipy",
    "pandas",
    "xarray",
    "geopandas",
    "shapely",
    "pyproj",
    "SALib",
    "matplotlib",
    "pydantic",
    "holidays",
    "pint",
)


class StrictModel(BaseModel):
    model_config = ConfigDict(extra="forbid")


class ClosureConfig(StrictModel):
    launch_mean_h: float = Field(12.9, gt=0)
    launch_cv: float = Field(0.55, ge=0)
    ss_landing_mean_h: float = Field(3.0, ge=0)
    ss_landing_cv: float = Field(0.5, ge=0)
    static_fire_mean_h: float = Field(4.0, ge=0)
    static_fire_cv: float = Field(0.5, ge=0)
    scrub_fraction: float = Field(0.7, ge=0, le=1)
    lead_fraction: float = Field(0.75, ge=0, le=1)
    anomaly_probability: float = Field(0.02, ge=0, le=1)
    anomaly_hours: float = Field(48.0, ge=0)
    external_landing_mean_h: float = Field(3.0, ge=0)
    hypothetical: list[str] = Field(default_factory=list, description="Names of hypothetical parameters")

    def to_model(self) -> ClosureModel:
        d = self.model_dump()
        d.pop("hypothetical")
        return ClosureModel(**d)


class OperationsConfig(StrictModel):
    n_launches: float = Field(gt=0)
    model: Literal["poisson", "weekday_renewal", "campaign"] = "poisson"
    weekday_only: bool = False
    night_fraction: float = Field(0.5, ge=0, le=1)
    n_pads: int = Field(1, ge=1, le=4)
    turnaround_min_days: float = Field(0.5, ge=0)
    gamma_shape: float = Field(2.0, gt=0)
    renewal_slack: float = Field(0.1, ge=0, lt=1)
    cluster_rate_per_year: float = Field(60.0, gt=0)
    cluster_spread_days: float = Field(2.0, gt=0)
    p_sh_landing_local: float = Field(1.0, ge=0, le=1)
    p_ss_landing_local: float = Field(1.0, ge=0, le=1)
    sh_landing_delay_min: float = Field(10.0, ge=0)
    ss_landing_delay_days_min: float = Field(0.04, ge=0)
    ss_landing_delay_days_max: float = Field(3.0, ge=0)
    n_external_landings: float = Field(0.0, ge=0)
    external_weekday_only: bool = False
    static_fires_per_launch: float = Field(0.2, ge=0)
    static_fire_lead_days: tuple[float, float] = (1.0, 14.0)
    scrub_probability: float = Field(0.15, ge=0, lt=1)
    scrub_recycle_days: tuple[float, float] = (1.0, 3.0)
    closures: ClosureConfig = Field(default_factory=ClosureConfig)
    hypothetical: list[str] = Field(default_factory=list)

    @model_validator(mode="after")
    def _check(self) -> OperationsConfig:
        if self.ss_landing_delay_days_max < self.ss_landing_delay_days_min:
            raise ValueError("ss_landing_delay_days_max must be >= min")
        return self

    def to_params(self, year: int) -> ScheduleParams:
        d = self.model_dump()
        d.pop("hypothetical")
        d.pop("closures")
        return ScheduleParams(year=year, closures=self.closures.to_model(), **d)


class FootprintConfig(StrictModel):
    level_offset_db: float = 0.0
    boom_scale: float = Field(1.0, gt=0)
    boom_elongation: float = Field(0.0, ge=0, lt=1)
    hypothetical: list[str] = Field(default_factory=list)


class ScenarioConfig(StrictModel):
    name: str
    description: str = ""
    hypothetical: bool = False
    source: str = ""
    year: int = 2030
    operations: OperationsConfig
    footprints: FootprintConfig = Field(default_factory=FootprintConfig)

    @model_validator(mode="after")
    def _label(self) -> ScenarioConfig:
        if self.hypothetical and "hypothetical" not in self.name:
            raise ValueError("hypothetical scenarios must carry 'hypothetical' in their name")
        return self


class MonteCarloConfig(StrictModel):
    n_realizations: int = Field(10000, ge=10)
    se_target_fraction: float = Field(0.02, gt=0)
    n_workers: int = Field(1, ge=1)


class GridConfig(StrictModel):
    half_width_km: float = Field(30.0, gt=0)
    cell_m: float = Field(250.0, gt=0)
    convergence_cell_m: float = Field(125.0, gt=0)


class ThresholdConfig(StrictModel):
    lamax_db: list[float] = Field(default_factory=lambda: [90.0, 100.0, 110.0, 115.0])
    psf: list[float] = Field(default_factory=lambda: [1.0, 2.0, 4.0, 6.0])
    ldn_db: float = 65.0
    cdnl_db: float = 60.0
    disturbance_lamax_db: float = 90.0


class FactorRange(StrictModel):
    low: float
    high: float
    source_rank: int = Field(ge=0, le=5, description="0 denotes the random seed, which is not a model input")
    note: str = ""


class SobolConfig(StrictModel):
    n_base: int = Field(8192, ge=8)
    n_boot: int = Field(200, ge=10)
    cadences: list[float] = Field(default_factory=lambda: [25.0, 100.0, 456.0])
    factors: dict[str, FactorRange]
    report_threshold: float = 0.05


class SiteConfig(StrictModel):
    name: str
    launches_per_year: float
    source: str


class EmissionsConfig(StrictModel):
    sites: list[SiteConfig]
    propellant_per_launch_t: FactorRange
    f_stratosphere: FactorRange
    emission_index_g_per_kg: dict[str, FactorRange]
    comparison_scenarios_bc_gg_per_yr: list[float] = Field(default_factory=lambda: [10.0, 30.0, 100.0])
    n_samples: int = 20000


def _default_models() -> list[Literal["poisson", "weekday_renewal", "campaign"]]:
    return ["poisson", "weekday_renewal", "campaign"]


class ExperimentsConfig(StrictModel):
    cadences: list[float] = Field(default_factory=lambda: [25.0, 50.0, 100.0, 200.0, 456.0])
    external_landing_fractions: list[float] = Field(default_factory=lambda: [0.0, 0.5, 1.0])
    temporal_models: list[Literal["poisson", "weekday_renewal", "campaign"]] = Field(default_factory=_default_models)
    night_fractions: list[float] = Field(default_factory=lambda: [0.0, 0.1, 0.2, 0.3, 0.4, 0.5])
    closure_caps_h: list[float] = Field(default_factory=lambda: [500.0, 800.0, 2000.0])
    closure_distributions: dict[str, ClosureConfig]
    recovery_days: list[float]
    window_days: list[float]
    season: tuple[float, float] = (59.0, 243.0)
    footprint_offsets_db: list[float] = Field(default_factory=lambda: [0.0, 3.0, 6.0])
    boom_scales: list[float] = Field(default_factory=lambda: [1.0, 1.5, 2.0])
    reference_cadence: float = 25.0


class PipelineConfig(StrictModel):
    name: str
    seed: int
    scenarios: list[Path]
    receptors: Path
    footprint_points: Path
    per_operation: Path
    results_dir: Path
    monte_carlo: MonteCarloConfig = Field(default_factory=MonteCarloConfig)
    grid: GridConfig = Field(default_factory=GridConfig)
    thresholds: ThresholdConfig = Field(default_factory=ThresholdConfig)
    experiments: ExperimentsConfig
    sobol: SobolConfig
    emissions: EmissionsConfig
    pad_lon: float
    pad_lat: float
    pad_source: str


def load_yaml(path: Path) -> dict[str, Any]:
    with open(path, encoding="utf-8") as fh:
        data = yaml.safe_load(fh)
    if not isinstance(data, dict):
        raise ValueError(f"{path} does not contain a mapping")
    return data


def load_scenario(path: Path) -> ScenarioConfig:
    return ScenarioConfig.model_validate(load_yaml(path))


def load_pipeline(path: Path) -> PipelineConfig:
    data = load_yaml(path)
    base = path.parent
    for key in ("scenarios",):
        data[key] = [str((base / p).resolve()) if not Path(p).is_absolute() else p for p in data[key]]
    for key in ("receptors", "footprint_points", "per_operation", "results_dir"):
        p = Path(data[key])
        data[key] = str((base.parent / p).resolve()) if not p.is_absolute() else str(p)
    return PipelineConfig.model_validate(data)


# --------------------------------------------------------------------------
# provenance
# --------------------------------------------------------------------------


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def sha256_text(text: str) -> str:
    return hashlib.sha256(text.encode("utf-8")).hexdigest()


def git_commit() -> str:
    try:
        out = subprocess.run(["git", "rev-parse", "HEAD"], capture_output=True, text=True, check=True)
        return out.stdout.strip()
    except Exception:  # noqa: BLE001 - provenance must never fail the run
        return "unknown"


def package_versions() -> dict[str, str]:
    out = {"python": platform.python_version()}
    for name in PACKAGES_FOR_PROVENANCE:
        try:
            out[name] = metadata.version(name)
        except metadata.PackageNotFoundError:
            out[name] = "missing"
    return out


def provenance(
    config_paths: list[Path], input_paths: list[Path], seed: int, extra: dict[str, Any] | None = None
) -> dict[str, Any]:
    return {
        "created_utc": datetime.now(UTC).isoformat(timespec="seconds"),
        "git_commit": git_commit(),
        "config_hash": sha256_text("".join(Path(p).read_text(encoding="utf-8") for p in config_paths)),
        "config_files": [str(p) for p in config_paths],
        "input_checksums": {str(p): sha256_file(Path(p)) for p in input_paths if Path(p).exists()},
        "package_versions": package_versions(),
        "seed": seed,
        **(extra or {}),
    }


def write_sidecar(path: Path, prov: dict[str, Any]) -> Path:
    side = path.with_suffix(path.suffix + ".json")
    side.write_text(json.dumps(prov, indent=2, sort_keys=True, default=_json_default), encoding="utf-8")
    return side


def _json_default(o: Any) -> Any:
    if isinstance(o, np.generic):
        return o.item()
    if isinstance(o, np.ndarray):
        return o.tolist()
    if isinstance(o, Path):
        return str(o)
    if hasattr(o, "__dataclass_fields__"):
        return asdict(o)
    raise TypeError(f"not serialisable: {type(o)}")


def write_parquet(df: pd.DataFrame, path: Path, prov: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_parquet(path, index=False)
    write_sidecar(path, prov)


def write_csv(df: pd.DataFrame, path: Path, prov: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, index=False)
    write_sidecar(path, prov)


def write_netcdf(ds: xr.Dataset, path: Path, prov: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    ds = ds.copy()
    ds.attrs.update(
        {
            k: json.dumps(v, default=_json_default) if not isinstance(v, str | int | float) else v
            for k, v in prov.items()
        }
    )
    ds.to_netcdf(path)
    write_sidecar(path, prov)


def write_json(obj: Any, path: Path, prov: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps({"provenance": prov, "data": obj}, indent=2, default=_json_default), encoding="utf-8")
