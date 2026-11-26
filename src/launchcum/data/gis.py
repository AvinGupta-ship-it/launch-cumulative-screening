"""Build analysis-ready spatial layers and the receptor table from raw inputs.

Outputs (``data/processed``):

``site_layers.gpkg``
    Layers ``places``, ``blocks``, ``parks``, ``refuges``, ``critical_habitat``
    clipped to the analysis domain, in EPSG:32614.
``receptors.csv``
    Named receptor points with longitude, latitude, kind and derivation.
``blocks_population.csv``
    Populated 2020 census blocks (centroid and population) in the domain.
"""

from __future__ import annotations

import io
import zipfile
from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd
import shapely
from shapely.geometry import Point, box

from ..spatial import CRS_GEO, CRS_UTM, lonlat_to_utm

ROOT = Path(__file__).resolve().parents[3]
RAW = ROOT / "data" / "raw" / "gis"
PROCESSED = ROOT / "data" / "processed"

PLACE_RECEPTORS = {
    "Port Isabel": "Port Isabel",
    "South Padre Island": "South Padre Island",
    "Brownsville": "Brownsville",
    "Laguna Vista": "Laguna Vista",
    "Laguna Heights": "Laguna Heights",
    "Los Fresnos": "Los Fresnos",
}


def _domain(pad_x: float, pad_y: float, half_width_m: float) -> shapely.Polygon:
    return box(pad_x - half_width_m, pad_y - half_width_m, pad_x + half_width_m, pad_y + half_width_m)


def _read_places() -> gpd.GeoDataFrame:
    with zipfile.ZipFile(RAW / "tl_2025_48_place.zip") as zf:
        shp = [n for n in zf.namelist() if n.endswith(".shp")][0]
        parts = {ext: zf.read(shp[:-4] + ext) for ext in (".shp", ".shx", ".dbf", ".prj")}
    tmp = PROCESSED / "_tmp_place"
    tmp.mkdir(parents=True, exist_ok=True)
    for ext, data in parts.items():
        (tmp / f"place{ext}").write_bytes(data)
    g = gpd.read_file(tmp / "place.shp")
    for f in tmp.iterdir():
        f.unlink()
    tmp.rmdir()
    return g.to_crs(CRS_UTM)


def build(pad_lon: float, pad_lat: float, half_width_km: float = 30.0) -> dict[str, Path]:
    PROCESSED.mkdir(parents=True, exist_ok=True)
    px, py = lonlat_to_utm(np.asarray([pad_lon]), np.asarray([pad_lat]))
    pad_x, pad_y = float(px[0]), float(py[0])
    dom = _domain(pad_x, pad_y, half_width_km * 1000.0)
    pad = Point(pad_x, pad_y)
    gpkg = PROCESSED / "site_layers.gpkg"
    if gpkg.exists():
        gpkg.unlink()

    receptors: list[dict[str, object]] = []

    # Places -------------------------------------------------------------
    places = _read_places()
    places["geometry"] = places.geometry.make_valid()
    places = places[places.intersects(dom)].copy()
    places_out = places[["GEOID", "NAME", "NAMELSAD", "INTPTLAT", "INTPTLON", "geometry"]]
    places_out.to_file(gpkg, layer="places", driver="GPKG")
    for name, label in PLACE_RECEPTORS.items():
        row = places[places["NAME"] == name]
        if row.empty:
            continue
        geom = row.geometry.iloc[0]
        receptors.append(
            {
                "name": label,
                "lon": float(row["INTPTLON"].iloc[0]),
                "lat": float(row["INTPTLAT"].iloc[0]),
                "kind": "community",
                "derivation": "TIGER/Line 2025 place internal point",
                "geoid": str(row["GEOID"].iloc[0]),
            }
        )
        near = shapely.ops.nearest_points(pad, geom)[1]
        lon, lat = gpd.GeoSeries([near], crs=CRS_UTM).to_crs(CRS_GEO).iloc[0].coords[0]
        receptors.append(
            {
                "name": f"{label} (nearest edge)",
                "lon": lon,
                "lat": lat,
                "kind": "community_edge",
                "derivation": "point of the TIGER/Line 2025 place polygon nearest the pad",
                "geoid": str(row["GEOID"].iloc[0]),
            }
        )

    # Area hydrography (shorelines of the Gulf, Laguna Madre and bays) -------
    water = gpd.read_file(RAW / "tl_2025_48061_areawater.zip").to_crs(CRS_UTM)
    water["geometry"] = water.geometry.make_valid()
    water = water[water.intersects(dom)].copy()
    water[["HYDROID", "FULLNAME", "MTFCC", "AWATER", "geometry"]].to_file(gpkg, layer="water", driver="GPKG")

    # Census blocks ---------------------------------------------------------
    blocks = gpd.read_file(RAW / "tigerweb_blocks_bbox.geojson").to_crs(CRS_UTM)
    blocks["POP100"] = pd.to_numeric(blocks["POP100"])
    blocks["HU100"] = pd.to_numeric(blocks["HU100"])
    blocks = blocks[blocks.intersects(dom)].copy()
    blocks["cx"] = blocks.geometry.centroid.x
    blocks["cy"] = blocks.geometry.centroid.y
    blocks["r_km"] = np.hypot(blocks["cx"] - pad_x, blocks["cy"] - pad_y) / 1000.0
    blocks[["GEOID", "POP100", "HU100", "AREALAND", "geometry"]].to_file(gpkg, layer="blocks", driver="GPKG")
    pop = blocks[blocks["POP100"] > 0]
    lon, lat = (
        gpd.GeoSeries(gpd.points_from_xy(pop["cx"], pop["cy"]), crs=CRS_UTM).to_crs(CRS_GEO).x,
        gpd.GeoSeries(gpd.points_from_xy(pop["cx"], pop["cy"]), crs=CRS_UTM).to_crs(CRS_GEO).y,
    )
    pd.DataFrame(
        {
            "GEOID": pop["GEOID"].values,
            "pop_2020": pop["POP100"].values,
            "housing_units_2020": pop["HU100"].values,
            "lon": lon.values,
            "lat": lat.values,
            "r_km": pop["r_km"].values,
        }
    ).to_csv(PROCESSED / "blocks_population.csv", index=False)
    village = pop[pop["r_km"] <= 4.0]
    if len(village):
        w = village["POP100"].to_numpy(float)
        vx = float(np.average(village["cx"], weights=w))
        vy = float(np.average(village["cy"], weights=w))
        vlon, vlat = gpd.GeoSeries([Point(vx, vy)], crs=CRS_UTM).to_crs(CRS_GEO).iloc[0].coords[0]
        receptors.append(
            {
                "name": "Boca Chica Village",
                "lon": vlon,
                "lat": vlat,
                "kind": "community",
                "derivation": (
                    f"population-weighted centroid of {len(village)} populated 2020 blocks within 4 km "
                    f"of the pad (population {int(w.sum())})"
                ),
                "geoid": "",
            }
        )

    # Parks -----------------------------------------------------------------
    tp = gpd.read_file(RAW / "tpwd_parks_bocachica.geojson").to_crs(CRS_UTM)
    tp = tp.rename(columns={"LoName": "name"})
    tp["source"] = "TPWD provisional boundary via TxDOT feature service"
    pad_us = gpd.read_file(RAW / "padus_site.geojson").to_crs(CRS_UTM)
    pad_us = pad_us.rename(columns={"Unit_Nm": "name"})
    pad_us["source"] = "USGS PAD-US"
    parks = pd.concat([tp[["name", "source", "geometry"]], pad_us[["name", "source", "geometry"]]], ignore_index=True)
    parks = gpd.GeoDataFrame(parks, geometry="geometry", crs=CRS_UTM)
    parks["geometry"] = parks.geometry.make_valid()
    parks = parks[parks.intersects(dom)]
    parks.to_file(gpkg, layer="parks", driver="GPKG")
    for label, match in (("Boca Chica State Park", "Boca Chica SP"), ("Isla Blanca Park", "Isla Blanca Park")):
        row = parks[parks["name"] == match]
        if row.empty:
            continue
        geom = row.geometry.union_all()
        c = geom.centroid
        lon, lat = gpd.GeoSeries([c], crs=CRS_UTM).to_crs(CRS_GEO).iloc[0].coords[0]
        receptors.append(
            {
                "name": label,
                "lon": lon,
                "lat": lat,
                "kind": "park",
                "derivation": f"centroid of {row['source'].iloc[0]} polygon",
                "geoid": "",
            }
        )

    # Refuges ------------------------------------------------------------------
    ref = pd.concat(
        [
            gpd.read_file(RAW / "nwr_lower_rio_grande_valley.geojson"),
            gpd.read_file(RAW / "nwr_laguna_atascosa.geojson"),
        ],
        ignore_index=True,
    )
    ref = gpd.GeoDataFrame(ref, geometry="geometry", crs=CRS_GEO).to_crs(CRS_UTM)
    ref["geometry"] = ref.geometry.make_valid()
    ref = ref.clip(dom)
    ref[["ORGNAME", "geometry"]].to_file(gpkg, layer="refuges", driver="GPKG")
    for label, org in (
        ("Lower Rio Grande Valley NWR (near pad)", "LOWER RIO GRANDE VALLEY NATIONAL WILDLIFE REFUGE"),
        ("Laguna Atascosa NWR (nearest)", "LAGUNA ATASCOSA NATIONAL WILDLIFE REFUGE"),
    ):
        geom = ref[ref["ORGNAME"] == org].geometry.union_all()
        if geom.is_empty:
            continue
        if "near pad" in label:
            within = geom.intersection(pad.buffer(10_000))
            c = within.centroid if not within.is_empty else shapely.ops.nearest_points(pad, geom)[1]
            der = "centroid of refuge land within 10 km of the pad"
        else:
            c = shapely.ops.nearest_points(pad, geom)[1]
            der = "point of the refuge boundary nearest the pad"
        lon, lat = gpd.GeoSeries([c], crs=CRS_UTM).to_crs(CRS_GEO).iloc[0].coords[0]
        receptors.append({"name": label, "lon": lon, "lat": lat, "kind": "refuge", "derivation": der, "geoid": ""})

    # Critical habitat ---------------------------------------------------------
    ch = pd.concat(
        [gpd.read_file(RAW / "ch_piping_plover_final.geojson"), gpd.read_file(RAW / "ch_red_knot_proposed.geojson")],
        ignore_index=True,
    )
    ch = gpd.GeoDataFrame(ch, geometry="geometry", crs=CRS_GEO).to_crs(CRS_UTM)
    ch["geometry"] = ch.geometry.make_valid()
    ch = ch.clip(dom)
    ch[["comname", "sciname", "status", "fedreg", "pubdate", "geometry"]].to_file(
        gpkg, layer="critical_habitat", driver="GPKG"
    )
    for label, com, how in (
        ("Piping plover critical habitat (near pad)", "Piping Plover", "within15"),
        ("Rufa red knot proposed critical habitat (nearest)", "rufa red knot", "nearest"),
    ):
        geom = ch[ch["comname"] == com].geometry.union_all()
        if geom.is_empty:
            continue
        if how == "within15":
            within = geom.intersection(pad.buffer(15_000))
            c = within.centroid if not within.is_empty else shapely.ops.nearest_points(pad, geom)[1]
            der = "centroid of designated habitat within 15 km of the pad"
        else:
            c = shapely.ops.nearest_points(pad, geom)[1]
            der = "point of the proposed habitat nearest the pad"
        lon, lat = gpd.GeoSeries([c], crs=CRS_UTM).to_crs(CRS_GEO).iloc[0].coords[0]
        receptors.append(
            {"name": label, "lon": lon, "lat": lat, "kind": "critical_habitat", "derivation": der, "geoid": ""}
        )

    rec = pd.DataFrame(receptors)
    rx, ry = lonlat_to_utm(rec["lon"].to_numpy(float), rec["lat"].to_numpy(float))
    rec["r_km"] = np.hypot(rx - pad_x, ry - pad_y) / 1000.0
    rec["azimuth_deg"] = np.degrees(np.arctan2(rx - pad_x, ry - pad_y)) % 360.0
    rec = rec.sort_values("r_km").reset_index(drop=True)
    rec_path = PROCESSED / "receptors.csv"
    rec.to_csv(rec_path, index=False, float_format="%.6f")
    return {"gpkg": gpkg, "receptors": rec_path, "blocks": PROCESSED / "blocks_population.csv"}


def read_layer(name: str) -> gpd.GeoDataFrame:
    return gpd.read_file(PROCESSED / "site_layers.gpkg", layer=name)


def receptor_frame() -> pd.DataFrame:
    return pd.read_csv(PROCESSED / "receptors.csv")


def _unused(_: io.BytesIO) -> None:  # pragma: no cover
    return None
