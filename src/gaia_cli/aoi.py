"""Areas of interest: a vector file (path or URL) holding one Polygon, or an
inline ``bbox``, ``intersects`` Polygon or circle. See dev/stage-file.md#aois.
"""

from __future__ import annotations

import json
from pathlib import Path

import geopandas as gpd
import shapely
import yaml

AOI_FORMS = ({"bbox"}, {"intersects"}, {"latitude", "longitude", "radius_km"})


def parse_aoi(spec, base_dir: Path) -> str | dict:
    """Check an AOI path, URL or inline geometry; resolve relative paths."""
    if isinstance(spec, str):
        if "://" in spec or Path(spec).is_absolute():
            return spec
        path = base_dir / spec
        if not path.exists():
            raise ValueError(f"AOI file {path} not found")
        return str(path)
    if not isinstance(spec, dict) or set(spec) not in AOI_FORMS:
        raise ValueError(
            f"{spec!r} is not an AOI; use a path or URL, {{bbox}}, {{intersects}} "
            "or {latitude, longitude, radius_km}"
        )
    if "bbox" in spec:
        b = spec["bbox"]
        if not (isinstance(b, list) and len(b) == 4 and all(_is_number(v) for v in b)):
            raise ValueError("bbox must be [west, south, east, north]")
        check_lonlat(b[0], b[1])
        check_lonlat(b[2], b[3])
        if b[0] >= b[2] or b[1] >= b[3]:
            raise ValueError("bbox must be [west, south, east, north] with west < east")
    elif "intersects" in spec:
        try:
            geom = shapely.from_geojson(json.dumps(spec["intersects"]))
        except (shapely.errors.ShapelyError, TypeError) as e:
            raise ValueError(f"intersects is not a GeoJSON geometry ({e})") from None
        if geom.geom_type != "Polygon":
            raise ValueError(f"intersects must be a Polygon, not {geom.geom_type}")
    else:
        check_lonlat(spec["longitude"], spec["latitude"])
        r = spec["radius_km"]
        if not _is_number(r) or r <= 0:
            raise ValueError("radius_km must be a positive number")
    return spec


def _is_number(v) -> bool:
    return isinstance(v, (int, float)) and not isinstance(v, bool)


def check_lonlat(lon, lat):
    if not (
        _is_number(lon) and _is_number(lat) and -180 <= lon <= 180 and -90 <= lat <= 90
    ):
        raise ValueError(f"({lon}, {lat}) is not a longitude, latitude pair")


def load_aoi(spec: str | dict) -> gpd.GeoDataFrame:
    """Read a vector file (path or URL) holding one Polygon, or build an
    inline AOI in EPSG:4326."""
    if isinstance(spec, str):
        aoi = gpd.read_file(spec)
        # MultiPolygons and several features aren't supported yet
        if (types := aoi.geom_type.tolist()) != ["Polygon"]:
            raise ValueError(
                f"AOI {spec} must be a single Polygon, not {', '.join(types) or 'empty'}"
            )
        return aoi
    if "bbox" in spec:
        geom = shapely.box(*spec["bbox"])
    elif "intersects" in spec:
        geom = shapely.from_geojson(json.dumps(spec["intersects"]))
    else:
        lon, lat = spec["longitude"], spec["latitude"]
        # Buffer in a local azimuthal equidistant projection so the radius is true
        aeqd = f"+proj=aeqd +lat_0={lat} +lon_0={lon} +units=m"
        circle = gpd.GeoSeries([shapely.Point(0, 0)], crs=aeqd).buffer(
            spec["radius_km"] * 1000, quad_segs=32
        )
        geom = circle.to_crs("EPSG:4326").iloc[0]
    return gpd.GeoDataFrame(geometry=[geom], crs="EPSG:4326")


def aoi_key(spec: str | dict) -> str:
    """A hashable key for *spec*, so each distinct AOI is read once."""
    return spec if isinstance(spec, str) else yaml.safe_dump(spec, sort_keys=True)
