"""
Helper functions shared by the GGOSSS OpenDrift teaching notebook and its solutions notebook.

Everything here is plain numpy / xarray / pyproj code operating on OpenDrift output (the NetCDF file
written by ``o.run(outfile=...)`` or the in-memory ``o.result`` Dataset), so that the diagnostics used
in the exercises (stranded fraction, convex-hull area, centroid separation, coastal sectors, EKE, ...)
are computed the same way everywhere.

Written against OpenDrift 1.14 (``o.result`` is an xarray Dataset with dims (trajectory, time)).
"""
from __future__ import annotations

from datetime import datetime, timedelta
from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr

R_EARTH_KM = 6371.0


# =====================================================================================================
# Forcing preparation
# =====================================================================================================
_U_NAMES = ["u10", "U10", "10u", "ugrd10m", "UGRD_10maboveground", "x_wind"]
_V_NAMES = ["v10", "V10", "10v", "vgrd10m", "VGRD_10maboveground", "y_wind"]


def _find_var(ds, candidates):
    for name in candidates:
        if name in ds.data_vars:
            return name
    for name, da in ds.data_vars.items():  # fall back on CF standard_name
        if da.attrs.get("standard_name") in candidates:
            return name
    return None


def prepare_wind_file(files, out_file, overwrite=False, bbox=None, time_range=None):
    """Merge 10 m wind component file(s) into one CF-compliant file for ``reader_netCDF_CF_generic``.

    ``files`` may be a single file containing both components, or several files (e.g. one file for
    u10 and one for v10, as delivered by the ERA5 ARCO/CDS services). Components are *merged* as
    variables (never concatenated along time: ``cdo cat`` on a u10 file and a v10 file silently
    drops v10, which is what caused the "missing y_wind" error in an earlier version of this notebook).

    The output has variables ``x_wind``/``y_wind`` with the CF standard names OpenDrift looks for,
    a ``time`` coordinate (CDS files call it ``valid_time``) and longitudes in [-180, 180].
    """
    out_file = Path(out_file)
    if out_file.exists() and not overwrite:
        return out_file
    files = [files] if isinstance(files, (str, Path)) else list(files)
    ds = xr.merge([xr.open_dataset(f) for f in files], compat="override", combine_attrs="drop")

    ren = {}
    if "valid_time" in ds.dims or "valid_time" in ds.coords:
        ren["valid_time"] = "time"
    for old, new in [("lon", "longitude"), ("lat", "latitude")]:
        if old in ds.dims or old in ds.coords:
            ren[old] = new
    ds = ds.rename(ren)
    ds = ds.drop_vars([v for v in ["number", "expver", "step", "heightAboveGround", "surface"]
                       if v in ds.variables and v not in ds.dims], errors="ignore")

    uname, vname = _find_var(ds, _U_NAMES), _find_var(ds, _V_NAMES)
    if uname is None or vname is None:
        raise ValueError(f"Could not find both wind components in {files}: variables = {list(ds.data_vars)}")
    ds = ds[[uname, vname]].rename({uname: "x_wind", vname: "y_wind"})

    if float(ds.longitude.max()) > 180:
        ds = ds.assign_coords(longitude=((ds.longitude + 180) % 360) - 180)
    ds = ds.sortby("longitude").sortby("latitude")
    if bbox is not None:
        lon0, lon1, lat0, lat1 = bbox
        ds = ds.sel(longitude=slice(lon0, lon1), latitude=slice(lat0, lat1))
    if time_range is not None:
        ds = ds.sel(time=slice(*time_range))

    ds["x_wind"].attrs = {"standard_name": "x_wind", "long_name": "10 m eastward wind", "units": "m s-1"}
    ds["y_wind"].attrs = {"standard_name": "y_wind", "long_name": "10 m northward wind", "units": "m s-1"}
    ds["longitude"].attrs = {"standard_name": "longitude", "units": "degrees_east", "axis": "X"}
    ds["latitude"].attrs = {"standard_name": "latitude", "units": "degrees_north", "axis": "Y"}
    ds["time"].attrs.update({"standard_name": "time", "axis": "T"})
    ds.attrs = {"Conventions": "CF-1.8", "source_files": ", ".join(str(f) for f in files)}

    enc = {v: {"zlib": True, "complevel": 1, "dtype": "float32"} for v in ["x_wind", "y_wind"]}
    out_file.parent.mkdir(parents=True, exist_ok=True)
    ds.load().to_netcdf(out_file, encoding=enc)
    return out_file


def download_gfs_10m_wind(cycle, lead_hours, bbox, out_file, overwrite=False):
    """Download GFS 0.25 deg 10 m wind (UGRD/VGRD only) for one forecast cycle from the NOAA Open Data
    archive on AWS (``noaa-gfs-bdp-pds``), using HTTP byte ranges so only a few MB are transferred.

    Works both for the latest cycles and for archived cycles (bucket goes back to 2021), which is
    what Exercise 7.1 needs. The legacy NOMADS OPeNDAP ("dods") endpoint is no longer reliable.

    cycle: datetime of the cycle (00/06/12/18 UTC); lead_hours: iterable of forecast hours;
    bbox: (lon_min, lon_max, lat_min, lat_max). Output: CF NetCDF with x_wind/y_wind (valid time).
    """
    import requests
    import tempfile

    out_file = Path(out_file)
    if out_file.exists() and not overwrite:
        return out_file
    base = (f"https://noaa-gfs-bdp-pds.s3.amazonaws.com/gfs.{cycle:%Y%m%d}/{cycle:%H}/atmos/"
            f"gfs.t{cycle:%H}z.pgrb2.0p25.f{{:03d}}")
    lon0, lon1, lat0, lat1 = bbox
    pieces = []
    for lead in lead_hours:
        url = base.format(int(lead))
        idx = requests.get(url + ".idx", timeout=60)
        idx.raise_for_status()
        lines = idx.text.strip().splitlines()
        chunks = []
        for i, line in enumerate(lines):
            if ":10 m above ground:" in line and (":UGRD:" in line or ":VGRD:" in line):
                start = int(line.split(":")[1])
                end = int(lines[i + 1].split(":")[1]) - 1 if i + 1 < len(lines) else ""
                chunks.append((start, end))
        data = b"".join(requests.get(url, headers={"Range": f"bytes={s}-{e}"}, timeout=120).content
                        for s, e in chunks)
        with tempfile.NamedTemporaryFile(suffix=".grib2", delete=False) as tmp:
            tmp.write(data)
        g = xr.open_dataset(tmp.name, engine="cfgrib", backend_kwargs={"indexpath": ""}).load()
        Path(tmp.name).unlink(missing_ok=True)
        g = g.assign_coords(longitude=((g.longitude + 180) % 360) - 180).sortby("longitude").sortby("latitude")
        g = g.sel(longitude=slice(lon0 - 0.5, lon1 + 0.5), latitude=slice(lat0 - 0.5, lat1 + 0.5))
        g = g[["u10", "v10"]].expand_dims(time=[pd.Timestamp(g.valid_time.values)])
        pieces.append(g.drop_vars([c for c in ["step", "valid_time", "heightAboveGround"] if c in g.coords]))
    ds = xr.concat(pieces, dim="time")
    ds.attrs["gfs_cycle"] = f"{cycle:%Y-%m-%d %H} UTC"
    tmp_nc = out_file.with_suffix(".raw.nc")
    ds.to_netcdf(tmp_nc)
    prepare_wind_file(tmp_nc, out_file, overwrite=True)
    tmp_nc.unlink(missing_ok=True)
    with xr.open_dataset(out_file) as check:
        print(f"GFS cycle {cycle:%Y-%m-%d %HZ}: {check.sizes['time']} lead times written to {out_file}")
    return out_file


def croco_gridfile(croco_file, out_file, grid_file=None):
    """Return a small grid file with the masks OpenDrift's ROMS/CROCO reader needs, or None if not needed.

    ``reader_ROMS_native`` requires ``mask_u`` and ``mask_v``; CROCO history files often contain only
    ``mask_rho`` (or no mask at all, in which case give the model grid file ``grid_file``). The u/v masks are
    derived from mask_rho (a u/v point is wet if both adjacent rho points are wet) and written to ``out_file``,
    to be passed as ``reader_ROMS_native.Reader(croco_file, gridfile=out_file)``.
    """
    out_file = Path(out_file)
    with xr.open_dataset(croco_file, decode_times=False) as ds:
        if "mask_u" in ds and "mask_v" in ds:
            return None
    # Reuse an existing mask file if it is newer than the CROCO file (rewriting a file that an open reader is
    # using fails on Windows drives and can corrupt what the reader sees)
    if out_file.exists() and out_file.stat().st_mtime >= Path(croco_file).stat().st_mtime:
        return out_file
    with xr.open_dataset(croco_file, decode_times=False) as ds:
        if "mask_rho" in ds:
            mask_rho = ds.mask_rho.load()
        elif grid_file is not None:
            with xr.open_dataset(grid_file) as g:
                mask_rho = g.mask_rho.load()
        else:
            raise ValueError(f"{croco_file} has no mask_rho: pass the CROCO grid file as grid_file=")
    m = mask_rho.values
    ny, nx = m.shape
    grid = xr.Dataset({
        "mask_u": (("eta_rho", "xi_u"), m[:, :-1] * m[:, 1:]),
        "mask_v": (("eta_v", "xi_rho"), m[:-1, :] * m[1:, :]),
    })
    out_file = Path(out_file)
    out_file.parent.mkdir(parents=True, exist_ok=True)
    grid.to_netcdf(out_file)
    return out_file


def check_coverage(readers, lon, lat, start, end):
    """Print, for each reader, whether it covers a (lon, lat) point and the [start, end] period.

    Useful before a run: OpenDrift silently falls back to default values (e.g. zero current) where no
    reader covers the particles, which gives plausible-looking but meaningless trajectories.
    """
    rows = []
    for r in readers:
        try:
            covered = len(np.atleast_1d(r.covers_positions(np.atleast_1d(lon), np.atleast_1d(lat))[0])) > 0
        except Exception:
            covered = None
        t0, t1 = getattr(r, "start_time", None), getattr(r, "end_time", None)
        time_ok = None if t0 is None else bool(t0 <= start and t1 >= end)
        rows.append(dict(reader=Path(str(r.name)).name, point_covered=covered, period_covered=time_ok,
                         start=t0, end=t1))
    table = pd.DataFrame(rows).set_index("reader")
    print(table.to_string())
    if is_land(lon, lat).any():
        print(f"\nWARNING: ({lon}, {lat}) is on land in the GSHHG coastline used by OpenDrift; elements seeded there are "
              "moved to the nearest ocean position (seed:ocean_only) or deactivated.")
    if not table.point_covered.fillna(True).all() or not table.period_covered.fillna(True).all():
        print("\nWARNING: at least one reader does not cover the release point and/or period. OpenDrift will "
              "use fallback values there - check the data files before interpreting results.")
    return table


# =====================================================================================================
# Seeding helpers
# =====================================================================================================
def seed_oil_mass(o, total_mass_kg, **seed_kwargs):
    """Seed an OpenOil simulation with a prescribed *total* oil mass (kg).

    In OpenDrift >= 1.11 ``OpenOil.seed_elements`` ignores a ``mass_oil=`` argument: the mass per
    element is recomputed from ``m3_per_hour`` x duration x oil density / number of elements.
    Here we seed with the default volume and then rescale the per-element mass so that the sum over
    elements equals ``total_mass_kg`` whatever the oil type (density) and number of elements.
    """
    n_before = o.num_elements_scheduled()
    o.seed_elements(**seed_kwargs)
    sched = o.elements_scheduled
    mass = np.broadcast_to(np.asarray(sched.mass_oil, dtype=float), (o.num_elements_scheduled(),)).copy()
    new = mass[n_before:]
    mass[n_before:] = new * (total_mass_kg / new.sum())
    sched.mass_oil = mass
    return o


def oil_properties(o, temperature_k=288.15):
    """Key properties of the oil record currently selected in an OpenOil object (after seeding)."""
    g = o.oiltype.gnome_oil
    boiling = np.asarray(g.get("boiling_point", []), dtype=float)
    mass_fraction = np.asarray(g.get("mass_fraction", []), dtype=float)
    light = mass_fraction[boiling < 473.15].sum() if boiling.size else np.nan  # boiling point < 200 degC
    return {
        "API gravity [deg]": g.get("api"),
        f"density @ {temperature_k - 273.15:.0f} C [kg/m3]": o.oiltype.density_at_temp(temperature_k),
        f"kinematic viscosity @ {temperature_k - 273.15:.0f} C [cSt]": o.oiltype.kvis_at_temp(temperature_k) * 1e6,
        "pour point [C]": (g.get("pour_point") or np.nan) - 273.15,
        "mass fraction boiling < 200 C [%]": 100 * light,
        "max water content of emulsion [%]": 100 * (g.get("emulsion_water_fraction_max") or np.nan),
    }


# =====================================================================================================
# Reading OpenDrift output
# =====================================================================================================
def as_dataset(sim_or_file):
    """Accept an OpenDrift simulation object, a Dataset or a path, return the output Dataset."""
    if isinstance(sim_or_file, xr.Dataset):
        return sim_or_file
    if isinstance(sim_or_file, (str, Path)):
        return xr.open_dataset(sim_or_file)
    return sim_or_file.result


def status_code(ds, name):
    """Numeric status code of a deactivation reason (e.g. 'stranded'), read from the file metadata.
    Returns None if that status never occurred in the run (e.g. 'stranded' with coastline_action='previous')."""
    meanings = str(ds.status.attrs.get("flag_meanings", "active")).split()
    return meanings.index(name) if name in meanings else None


def filled(ds, var, upto=None):
    """Variable forward-filled along time (OpenDrift writes NaN after an element is deactivated), so that
    element i at time t holds its last valid value - e.g. the stranding position of a stranded element."""
    da = ds[var]
    if upto is not None:
        da = da.isel(time=slice(0, upto + 1))
    # numpy forward fill along the last axis (time): no dependency on bottleneck/numbagg
    da = da.transpose(..., "time")
    v = np.asarray(da.values, dtype=float)
    idx = np.where(np.isfinite(v), np.arange(v.shape[-1]), 0)
    np.maximum.accumulate(idx, axis=-1, out=idx)
    out = np.take_along_axis(v, idx, axis=-1)
    return da.copy(data=out)


def release_time(ds):
    """Time of the first seeded element (first time at which any position is defined)."""
    valid = ds.lon.notnull().any("trajectory").values
    return pd.Timestamp(ds.time.values[np.argmax(valid)])


def time_index(ds, hours_after_release):
    """Index of the output time closest to release + ``hours_after_release`` (clipped to the run)."""
    target = release_time(ds) + pd.Timedelta(hours=hours_after_release)
    idx = int(np.argmin(np.abs(ds.time.values - np.datetime64(target))))
    if abs(pd.Timestamp(ds.time.values[idx]) - target) > pd.Timedelta(hours=1):
        print(f"  note: run output ends at {pd.Timestamp(ds.time.values[-1])}, requested {target}")
    return idx


def state_at(ds, hours_after_release=None):
    """Positions/status of all elements at a given time (default: end of run).

    Deactivated elements keep their last valid position (forward fill). Elements not yet seeded at that
    time are NaN. Returns a DataFrame indexed by trajectory with lon, lat, status (+ z, origin_marker).
    """
    idx = ds.sizes["time"] - 1 if hours_after_release is None else time_index(ds, hours_after_release)
    out = {}
    for var in ["lon", "lat", "status", "z", "origin_marker"]:
        if var in ds:
            out[var] = filled(ds, var, upto=idx).isel(time=-1).values
    return pd.DataFrame(out, index=ds.trajectory.values)


def stranded_mask(ds, hours_after_release=None):
    """Boolean array: element stranded at (or before) the given time."""
    code = status_code(ds, "stranded")
    st = state_at(ds, hours_after_release)
    if code is None:
        return np.zeros(len(st), dtype=bool)
    return (st.status.values == code)


def seeded_mask(ds, hours_after_release=None):
    """Boolean array: element released at (or before) the given time."""
    return np.isfinite(state_at(ds, hours_after_release).lon.values)


def first_impact(ds):
    """Time (h after release) and mean location of the first stranding event, or NaN if none."""
    code = status_code(ds, "stranded")
    if code is None:
        return dict(time_to_first_impact_h=np.nan, first_impact_lon=np.nan, first_impact_lat=np.nan)
    st = filled(ds, "status").values  # (trajectory, time)
    hit = (st == code)
    any_hit = hit.any(axis=0)
    if not any_hit.any():
        return dict(time_to_first_impact_h=np.nan, first_impact_lon=np.nan, first_impact_lat=np.nan)
    it = int(np.argmax(any_hit))
    lon = filled(ds, "lon", upto=it).isel(time=-1).values[hit[:, it]]
    lat = filled(ds, "lat", upto=it).isel(time=-1).values[hit[:, it]]
    t_h = (pd.Timestamp(ds.time.values[it]) - release_time(ds)) / pd.Timedelta(hours=1)
    return dict(time_to_first_impact_h=t_h, first_impact_lon=float(np.mean(lon)),
                first_impact_lat=float(np.mean(lat)))


# =====================================================================================================
# Geometry
# =====================================================================================================
def haversine_km(lon1, lat1, lon2, lat2):
    lon1, lat1, lon2, lat2 = map(np.deg2rad, [lon1, lat1, lon2, lat2])
    a = np.sin((lat2 - lat1) / 2) ** 2 + np.cos(lat1) * np.cos(lat2) * np.sin((lon2 - lon1) / 2) ** 2
    return 2 * R_EARTH_KM * np.arcsin(np.sqrt(a))


def destination_point(lon, lat, bearing_deg, distance_km):
    """Point reached from (lon, lat) travelling distance_km along an initial bearing (deg from North)."""
    lat1, lon1, b = np.deg2rad(lat), np.deg2rad(lon), np.deg2rad(bearing_deg)
    d = distance_km / R_EARTH_KM
    lat2 = np.arcsin(np.sin(lat1) * np.cos(d) + np.cos(lat1) * np.sin(d) * np.cos(b))
    lon2 = lon1 + np.arctan2(np.sin(b) * np.sin(d) * np.cos(lat1), np.cos(d) - np.sin(lat1) * np.sin(lat2))
    return float(np.rad2deg(lon2)), float(np.rad2deg(lat2))


def _local_xy_km(lon, lat, lon0=None, lat0=None):
    import pyproj
    lon0 = np.nanmean(lon) if lon0 is None else lon0
    lat0 = np.nanmean(lat) if lat0 is None else lat0
    proj = pyproj.Proj(proj="aeqd", lon_0=lon0, lat_0=lat0, units="km")
    return proj(np.asarray(lon), np.asarray(lat))


def hull_area_km2(lon, lat):
    """Area (km2) of the convex hull of a set of positions, computed in a local equidistant projection.

    A geometric envelope of one deterministic ensemble at one instant: a convenient proxy for the
    spread of the particles, NOT a probabilistic search area.
    """
    from scipy.spatial import ConvexHull
    lon, lat = np.asarray(lon, float), np.asarray(lat, float)
    ok = np.isfinite(lon) & np.isfinite(lat)
    if ok.sum() < 3:
        return 0.0
    x, y = _local_xy_km(lon[ok], lat[ok])
    try:
        return float(ConvexHull(np.column_stack([x, y])).volume)  # 'volume' of a 2D hull is its area
    except Exception:  # all points (nearly) collinear
        return 0.0


def ensemble_shape(lon, lat):
    """Centroid and principal-axes description of a particle cloud (km): a simple, reproducible way
    of comparing the *shape* of two ensembles (elongation and orientation), beyond their area."""
    lon, lat = np.asarray(lon, float), np.asarray(lat, float)
    ok = np.isfinite(lon) & np.isfinite(lat)
    lon, lat = lon[ok], lat[ok]
    x, y = _local_xy_km(lon, lat)
    cov = np.cov(np.vstack([x, y]))
    evals, evecs = np.linalg.eigh(cov)
    major, minor = np.sqrt(np.maximum(evals[::-1], 0))
    vx, vy = evecs[:, -1]
    orientation = (np.rad2deg(np.arctan2(vx, vy)) + 180) % 180  # bearing of major axis, 0-180 deg from N
    return dict(n=int(ok.sum()), centroid_lon=float(lon.mean()), centroid_lat=float(lat.mean()),
                std_major_km=float(major), std_minor_km=float(minor),
                aspect_ratio=float(major / minor) if minor > 0 else np.inf,
                major_axis_bearing_deg=float(orientation))


# =====================================================================================================
# Coastline helpers (GSHHG full-resolution landmask shipped with OpenDrift)
# =====================================================================================================
def is_land(lon, lat):
    from opendrift.readers.reader_global_landmask import get_mask
    lon, lat = np.atleast_1d(np.asarray(lon, float)), np.atleast_1d(np.asarray(lat, float))
    return np.asarray(get_mask().contains_many(lon, lat), dtype=bool)


def near_coast(lon, lat, distance_km=2.0, n_directions=16):
    """True for positions that have land within ``distance_km`` in at least one of n directions."""
    lon, lat = np.asarray(lon, float), np.asarray(lat, float)
    out = np.zeros(lon.shape, dtype=bool)
    ok = np.isfinite(lon) & np.isfinite(lat)
    for b in np.linspace(0, 360, n_directions, endpoint=False):
        for frac in (0.5, 1.0):
            pts = [destination_point(x, y, b, frac * distance_km) for x, y in zip(lon[ok], lat[ok])]
            if pts:
                px, py = np.array(pts).T
                out[ok] |= is_land(px, py)
    return out


def distance_to_coast_km(lon, lat, max_km=100.0, step_km=0.5, n_directions=36):
    """Approximate distance (km) from a sea point to the nearest land (searching along n bearings)."""
    bearings = np.linspace(0, 360, n_directions, endpoint=False)
    for d in np.arange(step_km, max_km + step_km, step_km):
        pts = np.array([destination_point(lon, lat, b, d) for b in bearings])
        if is_land(pts[:, 0], pts[:, 1]).any():
            return float(d)
    return np.inf


def offshore_point(lon, lat, distance_km=15.0, bearing_deg=None):
    """Move a point ``distance_km`` offshore.

    If ``bearing_deg`` is None, the direction is chosen automatically: opposite to the direction of the
    nearest coast (circular mean of the bearings where land is first found, searching up to 200 km), i.e.
    roughly perpendicular to the local coastline; if that path crosses land, the closest land-free bearing
    is used. Returns (lon, lat, bearing_deg).
    """
    if bearing_deg is not None:
        return (*destination_point(lon, lat, bearing_deg, distance_km), bearing_deg)
    bearings = np.arange(0, 360, 5.0)
    to_coast = None
    for d in np.arange(1.0, 201.0, 1.0):
        pts = np.array([destination_point(lon, lat, b, d) for b in bearings])
        hit = is_land(pts[:, 0], pts[:, 1])
        if hit.any():
            rad = np.deg2rad(bearings[hit])
            to_coast = np.rad2deg(np.arctan2(np.sin(rad).mean(), np.cos(rad).mean())) % 360
            break
    if to_coast is None:
        raise ValueError("No coast within 200 km: 'offshore' is undefined - give bearing_deg explicitly")
    ideal = (to_coast + 180) % 360
    for b in sorted(bearings, key=lambda b: abs((b - ideal + 180) % 360 - 180)):
        path = np.array([destination_point(lon, lat, b, d) for d in np.linspace(0.5, distance_km, 30)])
        if not is_land(path[:, 0], path[:, 1]).any():
            return (*destination_point(lon, lat, b, distance_km), float(b))
    raise ValueError("No land-free direction found - is the release point on land?")


def sector_fractions(ds, sectors, origin=None, hours_after_release=None):
    """Fraction (%) of released elements stranded in each coastal sector.

    sectors: dict name -> (lon_min, lon_max, lat_min, lat_max) boxes enclosing a stretch of coastline.
    origin: optional origin_marker value to select one source. Stranded elements outside every box are
    counted as 'other coast'; elements still afloat as 'at sea (not stranded)'.
    """
    st = state_at(ds, hours_after_release)
    sel = np.isfinite(st.lon.values)
    if origin is not None and "origin_marker" in st:
        sel &= (st.origin_marker.values == origin)
    code = status_code(ds, "stranded")
    stranded = sel & (st.status.values == code) if code is not None else np.zeros(len(st), bool)
    n = sel.sum()
    res, used = {}, np.zeros(len(st), bool)
    for name, (x0, x1, y0, y1) in sectors.items():
        inside = stranded & (st.lon.values >= x0) & (st.lon.values <= x1) & \
                 (st.lat.values >= y0) & (st.lat.values <= y1) & ~used
        used |= inside
        res[name] = 100 * inside.sum() / max(n, 1)
    res["other coast"] = 100 * (stranded & ~used).sum() / max(n, 1)
    res["at sea (not stranded)"] = 100 * (sel & ~stranded).sum() / max(n, 1)
    return pd.Series(res, name=f"n={n}")


# =====================================================================================================
# Wind and ocean diagnostics
# =====================================================================================================
def wind_speed_stats(wind_file, lon, lat, t0, t1, box_deg=0.5):
    """Mean and max 10 m wind speed in a box around (lon, lat) over [t0, t1] (file from prepare_wind_file)."""
    with xr.open_dataset(wind_file) as w:
        w = w.sel(longitude=slice(lon - box_deg, lon + box_deg), latitude=slice(lat - box_deg, lat + box_deg),
                  time=slice(t0, t1))
        spd = np.hypot(w.x_wind, w.y_wind)
        return dict(mean_wind_ms=float(spd.mean()), max_wind_ms=float(spd.max()),
                    mean_wind_cubed_m3s3=float((spd ** 3).mean()), n_times=int(w.sizes["time"]))


def croco_surface_currents(ds, times=None):
    """Surface u, v of a CROCO/ROMS history file interpolated to rho points (time, eta_rho, xi_rho), with
    lon/lat coordinates. Grid-relative components: fine for EKE (u'^2 + v'^2 is rotation invariant).
    ``times``: optional list of datetimes (e.g. ``reader_croco.times``) when the file time axis is not
    CF-decodable (CROCO often stores seconds since a model origin without a proper 'units' attribute)."""
    u = ds.u.isel(s_rho=-1)
    v = ds.v.isel(s_rho=-1)
    u_dims, v_dims = u.dims, v.dims
    xi_u = [d for d in u_dims if d.startswith("xi")][0]
    eta_v = [d for d in v_dims if d.startswith("eta")][0]
    nx, ny = ds.sizes["xi_rho"], ds.sizes["eta_rho"]
    uu = u.values
    vv = v.values
    u_rho = np.full(uu.shape[:-1] + (nx,), np.nan, dtype=float)
    v_rho = np.full(vv.shape[:-2] + (ny, vv.shape[-1]), np.nan, dtype=float)
    u_rho[..., 1:-1] = 0.5 * (uu[..., :-1] + uu[..., 1:])      # u at xi_u -> average onto interior rho points
    v_rho[..., 1:-1, :] = 0.5 * (vv[..., :-1, :] + vv[..., 1:, :])
    if "mask_rho" in ds:
        mask = ds.mask_rho.values == 0
    else:  # land points have identically zero velocity in CROCO output
        mask = np.all(u_rho == 0, axis=0) & np.all(v_rho == 0, axis=0)
    u_rho[..., mask] = np.nan
    v_rho[..., mask] = np.nan
    if times is None:
        times = ds[u.dims[0]].values
        if not np.issubdtype(np.asarray(times).dtype, np.datetime64):
            raise ValueError("CROCO time axis is not decoded to dates: pass times=reader_croco.times")
    times = list(times)
    if times and not isinstance(times[0], (datetime, np.datetime64, pd.Timestamp)):  # e.g. cftime dates
        times = [datetime(t.year, t.month, t.day, t.hour, t.minute, t.second) for t in times]
    times = pd.to_datetime(times)
    coords = dict(time=times, lon=(("eta_rho", "xi_rho"), ds.lon_rho.values),
                  lat=(("eta_rho", "xi_rho"), ds.lat_rho.values))
    dims = ("time", "eta_rho", "xi_rho")
    return xr.DataArray(u_rho, dims=dims, coords=coords), xr.DataArray(v_rho, dims=dims, coords=coords)


def time_mean_eke(u, v):
    """EKE = 1/2 (u'^2 + v'^2) averaged in time, anomalies relative to the time mean (no filtering)."""
    up, vp = u - u.mean("time"), v - v.mean("time")
    return (0.5 * (up ** 2 + vp ** 2)).mean("time")


def circular_diff_deg(a, b):
    """Smallest signed difference a - b between two directions in degrees, in [-180, 180)."""
    return (np.asarray(a) - np.asarray(b) + 180) % 360 - 180


def wind_direction_from(u, v):
    """Meteorological wind direction (deg, direction the wind blows FROM, clockwise from North)."""
    return (np.rad2deg(np.arctan2(-u, -v)) + 360) % 360


# =====================================================================================================
# Plotting
# =====================================================================================================
def map_axes(ax_or_fig=None, extent=None, subplot=(1, 1, 1), title=None):
    """Create a cartopy PlateCarree map axis with the GSHHG coastline used by OpenDrift."""
    import cartopy.crs as ccrs
    import matplotlib.pyplot as plt
    from opendrift.readers.reader_global_landmask import LandmaskFeature
    fig = ax_or_fig if ax_or_fig is not None else plt.figure(figsize=(7, 6))
    ax = fig.add_subplot(*subplot, projection=ccrs.PlateCarree())
    if extent is not None:
        ax.set_extent(extent, crs=ccrs.PlateCarree())
    # scale="f": full-resolution GSHHG polygons bundled with OpenDrift (roaring_landmask) - the same coastline
    # used for stranding, and read locally. Other scales make cartopy download the GSHHG archive (~150 MB).
    ax.add_feature(LandmaskFeature(scale="f"), facecolor="#e8e4d8", edgecolor="0.3", linewidth=0.6, zorder=2)
    gl = ax.gridlines(draw_labels=True, linewidth=0.3, alpha=0.5)
    gl.top_labels = gl.right_labels = False
    if title:
        ax.set_title(title, fontsize=10)
    return ax


def extent_of(*datasets, pad=0.15, points=()):
    """Map extent [lon0, lon1, lat0, lat1] enclosing all positions of the given runs (+ extra points)."""
    lons, lats = [np.asarray([p[0] for p in points])], [np.asarray([p[1] for p in points])]
    for ds in datasets:
        ds = as_dataset(ds)
        lons.append(ds.lon.values.ravel())
        lats.append(ds.lat.values.ravel())
    lon, lat = np.concatenate(lons), np.concatenate(lats)
    lon0, lon1, lat0, lat1 = np.nanmin(lon), np.nanmax(lon), np.nanmin(lat), np.nanmax(lat)
    dx, dy = max(lon1 - lon0, 0.1) * pad, max(lat1 - lat0, 0.1) * pad
    return [lon0 - dx, lon1 + dx, lat0 - dy, lat1 + dy]


def plot_ensemble(ax, ds, hours_after_release=None, color="C0", label=None, tracks=True, max_tracks=300,
                  hull=False, s=4):
    """Plot trajectories (light) and element positions at a given time on a cartopy map axis.
    Stranded elements are drawn as black-edged markers."""
    import cartopy.crs as ccrs
    ds = as_dataset(ds)
    pc = ccrs.PlateCarree()
    if tracks:
        idx = None if hours_after_release is None else time_index(ds, hours_after_release)
        sub = ds.isel(time=slice(0, None if idx is None else idx + 1))
        step = max(1, ds.sizes["trajectory"] // max_tracks)
        ax.plot(sub.lon.values[::step].T, sub.lat.values[::step].T, color=color, lw=0.3, alpha=0.3,
                transform=pc, zorder=3)
    st = state_at(ds, hours_after_release)
    code = status_code(ds, "stranded")
    strd = (st.status.values == code) if code is not None else np.zeros(len(st), bool)
    ax.scatter(st.lon[~strd], st.lat[~strd], s=s, color=color, transform=pc, zorder=4, label=label)
    if strd.any():
        ax.scatter(st.lon[strd], st.lat[strd], s=s + 6, facecolor=color, edgecolor="k", linewidth=0.4,
                   transform=pc, zorder=5)
    if hull:
        from scipy.spatial import ConvexHull
        ok = np.isfinite(st.lon.values)
        pts = np.column_stack([st.lon.values[ok], st.lat.values[ok]])
        if len(pts) >= 3:
            h = ConvexHull(pts)
            v = np.r_[h.vertices, h.vertices[0]]
            ax.plot(pts[v, 0], pts[v, 1], color=color, lw=1.2, ls="--", transform=pc, zorder=5)
    return ax
