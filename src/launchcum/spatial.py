"""Grid and receptor geometry.

All spatial computations are carried out in UTM zone 14N (EPSG:32614). The
analysis grid is centered on the launch pad and holds distance and azimuth
arrays used to evaluate the isotropic and azimuth-dependent footprints.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd
import shapely
from numpy.typing import NDArray
from pyproj import Transformer
from shapely.geometry.base import BaseGeometry

FloatArray = NDArray[np.float64]
BoolArray = NDArray[np.bool_]

CRS_GEO = "EPSG:4326"
CRS_UTM = "EPSG:32614"

_to_utm = Transformer.from_crs(CRS_GEO, CRS_UTM, always_xy=True)
_to_geo = Transformer.from_crs(CRS_UTM, CRS_GEO, always_xy=True)


def lonlat_to_utm(lon: FloatArray, lat: FloatArray) -> tuple[FloatArray, FloatArray]:
    x, y = _to_utm.transform(np.asarray(lon, float), np.asarray(lat, float))
    return np.asarray(x), np.asarray(y)


def utm_to_lonlat(x: FloatArray, y: FloatArray) -> tuple[FloatArray, FloatArray]:
    lon, lat = _to_geo.transform(np.asarray(x, float), np.asarray(y, float))
    return np.asarray(lon), np.asarray(lat)


@dataclass
class Grid:
    """Regular grid centered on the pad."""

    x0: float
    y0: float
    half_width_m: float
    cell_m: float

    def __post_init__(self) -> None:
        n = int(round(2 * self.half_width_m / self.cell_m))
        edges = -self.half_width_m + self.cell_m * np.arange(n + 1)
        centers = 0.5 * (edges[:-1] + edges[1:])
        self.n = n
        self.x = self.x0 + centers
        self.y = self.y0 + centers
        self.xx, self.yy = np.meshgrid(self.x, self.y)
        dx = self.xx - self.x0
        dy = self.yy - self.y0
        self.r_km = np.hypot(dx, dy) / 1000.0
        self.azimuth_deg = np.degrees(np.arctan2(dx, dy)) % 360.0

    @property
    def cell_area_km2(self) -> float:
        return (self.cell_m / 1000.0) ** 2

    def area_above(self, field: FloatArray, threshold: float) -> float:
        return float(np.sum(field > threshold) * self.cell_area_km2)

    def max_extent_km(self, field: FloatArray, threshold: float) -> float:
        m = field > threshold
        return float(np.max(self.r_km[m])) if m.any() else 0.0

    def mask_from_geometry(self, geom: BaseGeometry) -> BoolArray:
        """Cells whose centers fall inside ``geom`` (UTM coordinates)."""
        return np.asarray(shapely.contains_xy(geom, self.xx, self.yy))

    def lonlat(self) -> tuple[FloatArray, FloatArray]:
        return utm_to_lonlat(self.xx, self.yy)


@dataclass
class ReceptorSet:
    """Named receptor points with distance and azimuth from the pad."""

    table: pd.DataFrame

    @classmethod
    def from_lonlat(cls, df: pd.DataFrame, pad_x: float, pad_y: float) -> ReceptorSet:
        x, y = lonlat_to_utm(df["lon"].to_numpy(float), df["lat"].to_numpy(float))
        out = df.copy()
        out["x"] = x
        out["y"] = y
        out["r_km"] = np.hypot(x - pad_x, y - pad_y) / 1000.0
        out["azimuth_deg"] = np.degrees(np.arctan2(x - pad_x, y - pad_y)) % 360.0
        return cls(out.reset_index(drop=True))

    @property
    def names(self) -> list[str]:
        return [str(n) for n in self.table["name"]]

    @property
    def r_km(self) -> FloatArray:
        return self.table["r_km"].to_numpy(float)

    @property
    def azimuth_deg(self) -> FloatArray:
        return self.table["azimuth_deg"].to_numpy(float)

    def distance_map(self) -> dict[str, float]:
        return dict(zip(self.names, self.r_km.tolist(), strict=True))


def pad_utm(lon: float, lat: float) -> tuple[float, float]:
    x, y = lonlat_to_utm(np.asarray([lon]), np.asarray([lat]))
    return float(x[0]), float(y[0])
