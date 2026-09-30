"""
gtools/postprocess.py — CROCO post-processing helpers.

Load and derive fields from CROCO history/average files (croco_his.nc,
croco_avg.nc). Everything is **region-agnostic**: coordinates, mask and
bathymetry are read from the file itself, so the same functions work for any
CROCO configuration (Canary, Gulf of Guinea, Agulhas, ...).

The CROCO vertical-coordinate math (sigma -> depth) is the standard NEW_S_COORD
(Vtransform = 2) formulation; it is adapted from the croco_pytools / somisana
conventions.

Typical use
-----------
    import gtools.postprocess as pp

    ds   = pp.open_history("croco_his.nc")      # xarray Dataset, time decoded
    sst  = pp.surface(ds, "temp")               # top-sigma temperature, (time,eta,xi)
    lon, lat, mask = pp.lonlatmask(ds)          # 2D rho-point coords + land mask
    u, v = pp.surface_uv(ds)                    # surface currents on rho points
    z    = pp.depths(ds)                        # 3D depth of every sigma level (m, negative down)

    # a whole cycled hindcast as one continuous series:
    ds_all = pp.open_run("model-runs/Canary_12", phase="hcast")
"""
from __future__ import annotations
import glob
import os
import numpy as np
import xarray as xr


# ----------------------------------------------------------------------
# Opening files
# ----------------------------------------------------------------------
def open_history(fname, decode_times=True, Yorig=None):
    """Open a single CROCO history/average file as an xarray Dataset.

    CROCO writes 'time' (and 'scrum_time') in seconds since the run's reference
    date (Yorig-01-01). If the file's units are CF-compliant we let xarray decode
    them; otherwise pass Yorig (e.g. 1993 for the hindcast) to build real dates.
    """
    ds = xr.open_dataset(fname, decode_times=False)
    if decode_times and "time" in ds:
        ds = decode_time(ds, Yorig=Yorig)
    return ds


def open_run(run_root, phase="hcast", config=None, Yorig=None):
    """Open all cycles of a cycled run as ONE continuous dataset (concat in time).

    run_root : e.g. ".../hindcast/model-runs/Canary_12"
    phase    : subfolder holding the run, 'hcast' (hindcast) or 'fcst' (forecast)
    Returns a time-sorted concatenation of every cycle's croco_his.nc.
    """
    pattern = os.path.join(run_root, "*", phase, "CROCO_FILES", "croco_his.nc")
    files = sorted(glob.glob(pattern))
    if not files:
        raise FileNotFoundError(f"no croco_his.nc under {pattern}")
    ds_list = [open_history(f, Yorig=Yorig) for f in files]
    # data_vars="all" is xarray's current default; pinning it keeps this
    # behaviour when the default changes in a future version
    ds = xr.concat(ds_list, dim="time", data_vars="all")
    ds = ds.sortby("time")
    # drop duplicate timestamps at cycle joins (last-write-wins)
    _, idx = np.unique(ds["time"].values, return_index=True)
    if len(idx) != ds.sizes["time"]:
        ds = ds.isel(time=np.sort(idx))
    return ds


def decode_time(ds, Yorig=None):
    """Convert CROCO 'time' (seconds since Yorig-01-01) to datetime64.

    Preference order:
      1. CF-compliant units string ('seconds since YYYY-MM-DD...') -> xarray.
      2. Else if Yorig given -> Yorig-01-01 + seconds.
      3. Else leave raw values (and warn). Mirrors somisana handle_time().
    """
    import re
    from datetime import datetime, timedelta

    if "ocean_time" in ds and "time" not in ds:
        ds = ds.rename({"ocean_time": "time"})
    if "time" not in ds:
        return ds

    t = ds["time"]
    vals = np.atleast_1d(t.values)
    if np.issubdtype(np.asarray(vals).dtype, np.datetime64):
        return ds  # already decoded

    units = t.attrs.get("units", "")
    m = re.match(r"seconds since (\d{4}-\d{2}-\d{2}(?: \d{2}:\d{2}:\d{2})?)", units)
    if m:
        try:
            return xr.decode_cf(ds)
        except Exception:
            grp = m.group(1).strip()
            ref = datetime.strptime(grp, "%Y-%m-%d %H:%M:%S" if " " in grp else "%Y-%m-%d")
    elif Yorig is not None:
        ref = datetime(int(Yorig), 1, 1)
    else:
        print("WARNING: could not decode CROCO time (no CF units, no Yorig). "
              "Times left as raw seconds. Pass Yorig=... (e.g. 1993) to fix.")
        return ds

    dt = [ref + timedelta(seconds=float(x)) for x in vals]
    return ds.assign_coords(time=np.array(dt, dtype="datetime64[s]"))


def _decode_time(ds):
    """Backwards-compatible alias (no Yorig)."""
    return decode_time(ds, Yorig=None)


# ----------------------------------------------------------------------
# Coordinates, mask, bathymetry
# ----------------------------------------------------------------------
def lonlatmask(ds, grid="r"):
    """Return (lon2d, lat2d, mask2d) for rho ('r'), u, or v points.

    mask is 1 = ocean, NaN = land (ready to multiply into a field so land
    plots blank).
    """
    if grid == "r":
        lon, lat = ds["lon_rho"].values, ds["lat_rho"].values
        m = ds["mask_rho"].values if "mask_rho" in ds else np.ones_like(lon)
    elif grid == "u":
        lon, lat = ds["lon_u"].values, ds["lat_u"].values
        m = _rho2u(ds["mask_rho"].values) if "mask_rho" in ds else np.ones_like(lon)
    elif grid == "v":
        lon, lat = ds["lon_v"].values, ds["lat_v"].values
        m = _rho2v(ds["mask_rho"].values) if "mask_rho" in ds else np.ones_like(lon)
    else:
        raise ValueError("grid must be 'r', 'u' or 'v'")
    mask = np.where(m > 0, 1.0, np.nan)
    return lon, lat, mask


def bathymetry(ds):
    """Bathymetry h (m, positive down), masked to ocean."""
    _, _, mask = lonlatmask(ds)
    return ds["h"].values * mask


# ----------------------------------------------------------------------
# u,v regridding to rho points
# ----------------------------------------------------------------------
def _u2rho(u):
    """Regrid a u-point field onto the rho grid (extends the xi axis).

    Handles 2D (eta,xi_u), 3D (t|z,eta,xi_u) and 4D (t,z,eta,xi_u).
    Adapted from croco_pytools / somisana u2rho.
    """
    u = np.asarray(u)
    nd = u.ndim
    if nd == 4:
        T, D, Mp, L = u.shape
        ext = np.zeros((T, D, Mp, L + 2))
        ext[:, :, :, 1:-1] = u
        ext[:, :, :, 0] = u[:, :, :, 0]
        ext[:, :, :, -1] = u[:, :, :, -1]
        return 0.5 * (ext[:, :, :, :-1] + ext[:, :, :, 1:])
    elif nd == 3:
        TorD, Mp, L = u.shape
        ext = np.zeros((TorD, Mp, L + 2))
        ext[:, :, 1:-1] = u
        ext[:, :, 0] = u[:, :, 0]
        ext[:, :, -1] = u[:, :, -1]
        return 0.5 * (ext[:, :, :-1] + ext[:, :, 1:])
    else:  # 2D
        Mp, L = u.shape
        ext = np.zeros((Mp, L + 2))
        ext[:, 1:-1] = u
        ext[:, 0] = u[:, 0]
        ext[:, -1] = u[:, -1]
        return 0.5 * (ext[:, :-1] + ext[:, 1:])


def _v2rho(v):
    """Regrid a v-point field onto the rho grid (extends the eta axis).

    Handles 2D (eta_v,xi), 3D (t|z,eta_v,xi) and 4D (t,z,eta_v,xi).
    Adapted from croco_pytools / somisana v2rho.
    """
    v = np.asarray(v)
    nd = v.ndim
    if nd == 4:
        T, D, M, Lp = v.shape
        ext = np.zeros((T, D, M + 2, Lp))
        ext[:, :, 1:-1, :] = v
        ext[:, :, 0, :] = v[:, :, 0, :]
        ext[:, :, -1, :] = v[:, :, -1, :]
        return 0.5 * (ext[:, :, :-1, :] + ext[:, :, 1:, :])
    elif nd == 3:
        TorD, M, Lp = v.shape
        ext = np.zeros((TorD, M + 2, Lp))
        ext[:, 1:-1, :] = v
        ext[:, 0, :] = v[:, 0, :]
        ext[:, -1, :] = v[:, -1, :]
        return 0.5 * (ext[:, :-1, :] + ext[:, 1:, :])
    else:  # 2D
        M, Lp = v.shape
        ext = np.zeros((M + 2, Lp))
        ext[1:-1, :] = v
        ext[0, :] = v[0, :]
        ext[-1, :] = v[-1, :]
        return 0.5 * (ext[:-1, :] + ext[1:, :])


def _rho2u(r):
    return 0.5 * (r[..., :-1] + r[..., 1:])


def _rho2v(r):
    return 0.5 * (r[..., :-1, :] + r[..., 1:, :])


# ----------------------------------------------------------------------
# Surface fields
# ----------------------------------------------------------------------
def surface(ds, var, tindex=None):
    """Surface (top sigma level) of a 3D variable, land-masked.

    var    : 'temp', 'salt', ... (a variable with an s_rho dimension)
    tindex : int time index, or None for all times.
    Returns a DataArray/ndarray shaped (time, eta, xi) or (eta, xi).
    """
    da = ds[var]
    if "s_rho" in da.dims:
        da = da.isel(s_rho=-1)          # top level = surface
    if tindex is not None:
        da = da.isel(time=tindex)
    _, _, mask = lonlatmask(ds)
    return da * mask


def surface_uv(ds, tindex=None):
    """Surface currents (u, v) rotated? No — returned on rho points, grid-aligned.

    Returns (u_rho, v_rho), each (time, eta, xi) or (eta, xi), land-masked.
    NOTE: these are grid-oriented (xi, eta). To get east/north, rotate with the
    grid 'angle' (see rotate_uv).
    """
    u = ds["u"]
    v = ds["v"]
    if "s_rho" in u.dims:
        u = u.isel(s_rho=-1)
        v = v.isel(s_rho=-1)
    if tindex is not None:
        u = u.isel(time=tindex)
        v = v.isel(time=tindex)
    ur = _u2rho(u.values)
    vr = _v2rho(v.values)
    _, _, mask = lonlatmask(ds)
    return ur * mask, vr * mask


def rotate_uv(ds, u_rho, v_rho):
    """Rotate grid-oriented (u,v) on rho points to (east, north) using grid angle."""
    ang = ds["angle"].values
    ue = u_rho * np.cos(ang) - v_rho * np.sin(ang)
    vn = u_rho * np.sin(ang) + v_rho * np.cos(ang)
    return ue, vn


def speed(u, v):
    """Current speed magnitude."""
    return np.sqrt(np.asarray(u) ** 2 + np.asarray(v) ** 2)


# ----------------------------------------------------------------------
# Vertical coordinate: sigma -> depth  (NEW_S_COORD / Vtransform = 2)
# ----------------------------------------------------------------------
def csf(sc, theta_s, theta_b):
    """Sigma stretching function Cs(sc). Adapted from croco_pytools / somisana."""
    one64 = np.float64(1)
    if theta_s > 0.0:
        csrf = (one64 - np.cosh(theta_s * sc)) / (np.cosh(theta_s) - one64)
    else:
        csrf = -(sc ** 2)
    sc1 = csrf + one64
    if theta_b > 0.0:
        Cs = (np.exp(theta_b * sc1) - one64) / (np.exp(theta_b) - one64) - one64
    else:
        Cs = csrf
    return Cs


def z_levels(h, zeta, theta_s, theta_b, hc, N, type="rho", vtransform=2):
    """3D depths (m, negative down) of the sigma levels.

    Vectorised: accepts zeta of shape (M, L) or (T, M, L); returns z of shape
    (T, N, M, L). type = 'rho' or 'w'; vtransform = 1 (OLD) or 2 (NEW).
    Adapted from croco_pytools / somisana z_levels (orig. zlevs.m, P. Penven;
    J. Veitch & G. Fearon).
    """
    zeta = np.asarray(zeta)
    if zeta.ndim == 2:
        zeta = zeta[None, :, :]
    T, M, L = zeta.shape

    if vtransform == 2:
        ds_ = 1.0 / N
        if type == "w":
            sc = np.linspace(-1.0, 0.0, N + 1)
            Cs = csf(sc, theta_s, theta_b)
            N = N + 1
        else:
            sc = ds_ * (np.arange(1, N + 1) - N - 0.5)
            Cs = csf(sc, theta_s, theta_b)
    else:
        if type == "w":
            sc = (np.arange(0, N + 1) - N) / N
            N = N + 1
        else:
            sc = (np.arange(1, N + 1) - N - 0.5) / N
        cff1 = 1.0 / np.sinh(theta_s)
        cff2 = 0.5 / np.tanh(0.5 * theta_s)
        Cs = (1 - theta_b) * cff1 * np.sinh(theta_s * sc) + theta_b * (
            cff2 * np.tanh(theta_s * (sc + 0.5)) - 0.5)

    h = np.where(h == 0, 1e-2, h)
    zeta = np.maximum(zeta, 0.01 - h[None, :, :])
    hb = h[None, :, :]

    if vtransform == 2:
        h2 = hb + hc
        cff = hc * sc[:, None, None] + Cs[:, None, None] * hb
        z = cff * hb / h2 + zeta[:, None, :, :] * (1.0 + cff / h2)
    else:
        hinv = 1.0 / hb
        cff = hc * (sc[:, None, None] - Cs[:, None, None])
        z = cff + Cs[:, None, None] * hb + zeta[:, None, :, :] * (
            1.0 + (cff + Cs[:, None, None] * hb) * hinv)
    return z  # (T, N, M, L)


def depths(ds, tindex=0, type="rho"):
    """3D depth (m, negative down) of every sigma level at one time index.

    Reads the stretching params from the file. Prefers theta_s/theta_b if the
    file stores them; otherwise derives Cs from the file's Cs_rho directly via
    z_levels using theta from global attrs. Returns (N, eta, xi).
    """
    h = ds["h"].values
    zeta = (ds["zeta"].isel(time=tindex).values
            if "time" in ds["zeta"].dims else ds["zeta"].values)
    hc = float(ds["hc"].values)
    N = ds.sizes["s_rho"]
    vt = int(ds["Vtransform"].values) if "Vtransform" in ds else 2
    theta_s = _get_theta(ds, "theta_s")
    theta_b = _get_theta(ds, "theta_b")
    z = z_levels(h, zeta, theta_s, theta_b, hc, N, type=type, vtransform=vt)
    return z[0]  # drop the singleton time axis -> (N, eta, xi)


def _get_theta(ds, name):
    """Fetch theta_s/theta_b from a variable or a global attribute."""
    if name in ds:
        return float(ds[name].values)
    if name in ds.attrs:
        return float(ds.attrs[name])
    # sensible CROCO defaults if truly absent
    return 7.0 if name == "theta_s" else 2.0


# ----------------------------------------------------------------------
# Point time series and vertical sections (basic)
# ----------------------------------------------------------------------
def nearest_index(ds, lon0, lat0):
    """(eta, xi) index of the rho point nearest (lon0, lat0)."""
    lon = ds["lon_rho"].values
    lat = ds["lat_rho"].values
    d = (lon - lon0) ** 2 + (lat - lat0) ** 2
    j, i = np.unravel_index(np.argmin(d), d.shape)
    return int(j), int(i)


def timeseries(ds, var, lon0, lat0, surface_only=True, depth_m=None):
    """Time series of `var` at the point nearest (lon0, lat0).

    Works for raw (temp, salt, u, v, zeta) AND derived (speed, ke, vort, vort_f).
      surface_only=True (default) : top sigma level (or the 2D field for zeta)
      depth_m=<m>                 : interpolate to that true depth at each time
    Returns a labeled DataArray over time (carries CF attrs for auto-labels).
    """
    j, i = nearest_index(ds, lon0, lat0)
    tvals = ds["time"].values
    nt = len(np.atleast_1d(tvals))

    # 2D field (zeta/ssh): straightforward
    if var in ("zeta", "ssh"):
        da = ds["zeta"] if var == "zeta" else ds["zeta"]
        series = da.isel(eta_rho=j, xi_rho=i).values
        out = xr.DataArray(np.atleast_1d(series), dims=("time",),
                           coords={"time": tvals}, name=var)
        _apply_attrs(out, _canon_derived(var))
        if lon0 is not None: out.attrs["lon0"] = lon0
        if lat0 is not None: out.attrs["lat0"] = lat0
        out.attrs["point"] = f"lon={lon0:g}°, lat={lat0:g}°"
        return out

    # 3D field: build per time, take surface or interpolate to depth
    vals = np.full(nt, np.nan)
    for k in range(nt):
        d3 = _var3d(ds, var, tindex=k)            # (s_rho, eta, xi)
        col = d3.isel(eta_rho=j, xi_rho=i).values  # (s_rho,)
        if depth_m is None and surface_only:
            vals[k] = col[-1]                      # top sigma = surface
        elif depth_m is not None:
            z = depths(ds, tindex=k)[:, j, i]      # (s_rho,)
            target = -abs(depth_m)
            if target < np.nanmin(z) or target > np.nanmax(z):
                vals[k] = np.nan                   # shallower than depth
            else:
                order = np.argsort(z)
                vals[k] = np.interp(target, z[order], col[order])
        else:
            vals[k] = col[-1]

    out = xr.DataArray(vals, dims=("time",), coords={"time": tvals}, name=var)
    _apply_attrs(out, _canon_derived(var))
    if lon0 is not None: out.attrs["lon0"] = lon0
    if lat0 is not None: out.attrs["lat0"] = lat0
    out.attrs["point"] = f"lon={lon0:g}°, lat={lat0:g}°"
    if depth_m is not None:
        out.attrs["depth_m"] = depth_m
    return out


# ----------------------------------------------------------------------
# Small helpers
# ----------------------------------------------------------------------
def times(ds):
    """The time coordinate as returned by the file (datetime64 if decoded)."""
    return ds["time"].values


def extent(ds, pad=0.0):
    """[lon_min, lon_max, lat_min, lat_max] of the domain (for map extents)."""
    lon = ds["lon_rho"].values
    lat = ds["lat_rho"].values
    return [float(lon.min()) - pad, float(lon.max()) + pad,
            float(lat.min()) - pad, float(lat.max()) + pad]


# ======================================================================
# Labeled-field extraction (returns DataArrays with CF attrs + coords)
# so the plotting layer can introspect and auto-label.
# ======================================================================
try:
    from gtools.define_attrs import apply_attrs as _apply_attrs
except Exception:
    try:
        from define_attrs import apply_attrs as _apply_attrs
    except Exception:
        def _apply_attrs(da, name, **kw):   # no-op fallback
            return da


def field_map(ds, var, tindex=-1, level=-1):
    """A 2D horizontal field as a labeled DataArray (eta_rho, xi_rho).

    var   : any variable ('temp','salt','zeta',...). 3D vars are sliced at
            sigma `level` (-1 = surface).
    Returns a DataArray carrying lon_rho/lat_rho coords + CF attrs, ready to plot.
    """
    da = ds[var]
    if "time" in da.dims:
        da = da.isel(time=tindex)
    if "s_rho" in da.dims:
        da = da.isel(s_rho=level)
    _, _, mask = lonlatmask(ds)
    da = da * mask
    da = da.assign_coords(lon_rho=ds["lon_rho"], lat_rho=ds["lat_rho"])
    _apply_attrs(da, _canon(var))
    da.name = var
    return da


def profile(ds, var, lon0, lat0, tindex=-1):
    """Vertical profile of `var` at (lon0,lat0): DataArray over depth.

    Works for raw (temp, salt, u, v) AND derived (speed, ke, vort, vort_f)
    fields. Returns a DataArray with a 'depth' coord (m, negative down) + attrs.
    zeta/ssh are 2D and have no profile.
    """
    if var in ("zeta", "ssh"):
        raise ValueError(f"'{var}' is 2D (surface) - no vertical profile")
    j, i = nearest_index(ds, lon0, lat0)
    da3 = _var3d(ds, var, tindex=tindex)       # (s_rho, eta, xi)
    da = da3.isel(eta_rho=j, xi_rho=i)         # (s_rho,)
    z = depths(ds, tindex=tindex)[:, j, i]     # (s_rho,) depths
    da = da.assign_coords(depth=("s_rho", z))
    _apply_attrs(da, _canon_derived(var))
    if lon0 is not None: da.attrs["lon0"] = lon0
    if lat0 is not None: da.attrs["lat0"] = lat0
    da.attrs["point"] = f"lon={lon0:g}°, lat={lat0:g}°"
    da.name = var
    return da


def section(ds, var, lon0, lat0, lon1, lat1, tindex=-1, npts=200):
    """Vertical section of `var` along a straight lon/lat transect.

    Returns a DataArray on (s_rho, points) with coords: distance_km, depth (2D),
    section_lon, section_lat. Simple linear-index interpolation (good for regular
    or gently curvilinear grids); for a great-circle path see the docs.
    """
    # sample points along the transect
    slon = np.linspace(lon0, lon1, npts)
    slat = np.linspace(lat0, lat1, npts)

    if var in ("zeta", "ssh"):
        raise ValueError(f"'{var}' is 2D (surface) - no vertical section")
    da = _var3d(ds, var, tindex=tindex)        # (s_rho, eta, xi), raw or derived
    # nearest grid index for each sample point
    lonr = ds["lon_rho"].values
    latr = ds["lat_rho"].values
    maskr = ds["mask_rho"].values              # 1 water, 0 land
    z3d = depths(ds, tindex=tindex)            # (s_rho, eta, xi)
    N = da.sizes["s_rho"]
    prof = np.full((N, npts), np.nan)
    zsec = np.full((N, npts), np.nan)
    for p in range(npts):
        d = (lonr - slon[p]) ** 2 + (latr - slat[p]) ** 2
        j, i = np.unravel_index(np.argmin(d), d.shape)
        if maskr[j, i] == 0:                    # land -> leave NaN (blank)
            continue
        prof[:, p] = da.isel(eta_rho=j, xi_rho=i).values
        zsec[:, p] = z3d[:, j, i]
    # distance along section (km)
    dist = _haversine_km(lat0, lon0, slat, slon)

    out = xr.DataArray(
        prof, dims=("s_rho", "points"),
        coords={"distance_km": ("points", dist),
                "section_lon": ("points", slon),
                "section_lat": ("points", slat)},
        name=var)
    out = out.assign_coords(depth=(("s_rho", "points"), zsec))
    _apply_attrs(out, _canon_derived(var))
    return out


def hovmoller(ds, var, kind="time_depth", lon0=None, lat0=None,
              lon_line=None, lat_line=None, level=-1):
    """Hovmoller diagram data as a labeled 2D DataArray.

    kind:
      'time_depth' : time vs depth at a point (needs lon0,lat0) -> (time, s_rho)
      'time_lat'   : time vs latitude along a meridian (needs lon0) -> (time, eta_rho)
      'time_lon'   : time vs longitude along a parallel (needs lat0) -> (time, xi_rho)
    For time_lat/time_lon the field is taken at sigma `level` (-1 = surface).
    """
    if kind == "time_depth":
        j, i = nearest_index(ds, lon0, lat0)
        tvals = ds["time"].values
        nt = len(np.atleast_1d(tvals))
        if var in ds and "s_rho" in ds[var].dims and "eta_rho" in ds[var].dims and "xi_rho" in ds[var].dims:
            h = ds[var].isel(eta_rho=j, xi_rho=i)
        else:
            N = ds.sizes["s_rho"]
            mat = np.full((nt, N), np.nan)
            for k in range(nt):
                d3 = _var3d(ds, var, tindex=k)
                mat[k, :] = d3.isel(eta_rho=j, xi_rho=i).values
            h = xr.DataArray(mat, dims=("time", "s_rho"), coords={"time": tvals}, name=var)
        z = depths(ds, tindex=0)[:, j, i]
        h = h.assign_coords(depth=("s_rho", z))
        _apply_attrs(h, _canon_derived(var))
        if lon0 is not None: h.attrs["lon0"] = lon0
        if lat0 is not None: h.attrs["lat0"] = lat0
        if lon0 is not None and lat0 is not None:
            h.attrs["point"] = f"lon={lon0:g}°, lat={lat0:g}°"
        h.name = var
        return h
    da = ds[var] if var in ds else _var3d(ds, var, tindex=-1)
    if kind == "time_lat":
        i = _nearest_col(ds, lon0)                        # fixed xi (meridian)
        if "s_rho" in da.dims:
            da = da.isel(s_rho=level)
        h = da.isel(xi_rho=i)                             # (time, eta_rho)
        h = h.assign_coords(lat=("eta_rho", ds["lat_rho"].values[:, i]))
        _apply_attrs(h, _canon_derived(var))
        if lon0 is not None:
            h.attrs["lon0"] = lon0
            h.attrs["point"] = f"lon={lon0:g}°"
        h.name = var
        return h
    if kind == "time_lon":
        j = _nearest_row(ds, lat0)                        # fixed eta (parallel)
        if "s_rho" in da.dims:
            da = da.isel(s_rho=level)
        h = da.isel(eta_rho=j)                            # (time, xi_rho)
        h = h.assign_coords(lon=("xi_rho", ds["lon_rho"].values[j, :]))
        _apply_attrs(h, _canon_derived(var))
        if lat0 is not None:
            h.attrs["lat0"] = lat0
            h.attrs["point"] = f"lat={lat0:g}°"
        h.name = var
        return h
    raise ValueError("kind must be 'time_depth', 'time_lat' or 'time_lon'")


# ---- small helpers ----
def _canon(var):
    """Map a raw CROCO var name to the attribute-registry key."""
    return {"zeta": "zeta", "temp": "temp", "salt": "salt",
            "u": "u", "v": "v", "w": "w"}.get(var, var)


def _haversine_km(lat0, lon0, lats, lons, R=6371.0):
    lat0r, lon0r = np.radians(lat0), np.radians(lon0)
    latr, lonr = np.radians(lats), np.radians(lons)
    dlat = latr - lat0r
    dlon = lonr - lon0r
    a = np.sin(dlat / 2) ** 2 + np.cos(lat0r) * np.cos(latr) * np.sin(dlon / 2) ** 2
    return R * 2 * np.arcsin(np.sqrt(a))


def _nearest_col(ds, lon0):
    lonmid = ds["lon_rho"].values[ds.sizes["eta_rho"] // 2, :]
    return int(np.argmin(np.abs(lonmid - lon0)))


def _nearest_row(ds, lat0):
    latmid = ds["lat_rho"].values[:, ds.sizes["xi_rho"] // 2]
    return int(np.argmin(np.abs(latmid - lat0)))


def field_at_depth(ds, var, depth_m, tindex=-1):
    """Horizontal field of `var` interpolated to a fixed DEPTH (metres).

    depth_m : positive number = metres below surface (e.g. 30 for 30 m).
    Returns a labeled 2D DataArray (eta_rho, xi_rho); points where the sea floor
    is shallower than depth_m come back NaN (blank) - the physically correct
    'no water at this depth here'.
    """
    target = -abs(float(depth_m))              # depths are negative-down
    da = ds[var]
    if "time" in da.dims:
        da = da.isel(time=tindex)
    vals = da.values                            # (s_rho, eta, xi)
    z = depths(ds, tindex=tindex)               # (s_rho, eta, xi), negative down
    N, eta, xi = vals.shape

    out = np.full((eta, xi), np.nan)
    # linear interpolation in depth at every (eta,xi) column
    for j in range(eta):
        for i in range(xi):
            zc = z[:, j, i]
            vc = vals[:, j, i]
            # skip land / all-nan columns
            if not np.isfinite(vc).any():
                continue
            # only valid if target is within the column's depth range
            zmin, zmax = np.nanmin(zc), np.nanmax(zc)   # zmin=deepest, zmax~surface
            if target < zmin or target > zmax:
                continue                        # deeper than seafloor (or above surface)
            # np.interp needs increasing x -> sort by depth ascending
            order = np.argsort(zc)
            out[j, i] = np.interp(target, zc[order], vc[order])

    _, _, mask = lonlatmask(ds)
    out = out * mask
    da_out = xr.DataArray(out, dims=("eta_rho", "xi_rho"))
    da_out = da_out.assign_coords(lon_rho=ds["lon_rho"], lat_rho=ds["lat_rho"])
    _apply_attrs(da_out, _canon(var))
    da_out.attrs["depth_m"] = depth_m
    da_out.name = var
    return da_out


# ======================================================================
# Velocity at a true depth, and derived vorticity
# ======================================================================
def _var3d_at_depth(ds, values3d, depth_m, tindex=-1):
    """Interpolate a (s_rho, eta, xi) array to a fixed depth (m). Helper."""
    target = -abs(float(depth_m))
    z = depths(ds, tindex=tindex)               # (s_rho, eta, xi) negative down
    N, eta, xi = values3d.shape
    out = np.full((eta, xi), np.nan)
    for j in range(eta):
        for i in range(xi):
            zc = z[:, j, i]; vc = values3d[:, j, i]
            if not np.isfinite(vc).any():
                continue
            if target < np.nanmin(zc) or target > np.nanmax(zc):
                continue
            order = np.argsort(zc)
            out[j, i] = np.interp(target, zc[order], vc[order])
    return out


def uv_at_depth(ds, depth_m, tindex=-1, rotate=True):
    """Currents (u,v) on the rho grid at a true DEPTH (m), optionally rotated to
    east/north. Returns (u, v) 2D arrays; NaN where sea floor shallower than depth.
    """
    u = ds["u"]; v = ds["v"]
    if tindex is not None and "time" in u.dims:
        u = u.isel(time=tindex); v = v.isel(time=tindex)
    ur = _u2rho(u.values)                        # (s_rho, eta, xi)
    vr = _v2rho(v.values)
    u_d = _var3d_at_depth(ds, ur, depth_m, tindex=tindex)
    v_d = _var3d_at_depth(ds, vr, depth_m, tindex=tindex)
    _, _, mask = lonlatmask(ds)
    u_d, v_d = u_d * mask, v_d * mask
    if rotate:
        u_d, v_d = rotate_uv(ds, u_d, v_d)
    return u_d, v_d


def speed_map(ds, tindex=-1, depth_m=None):
    """Current SPEED as a labeled 2D DataArray. Surface if depth_m is None, else
    interpolated to that depth. Ready to plot (carries lon/lat + CF attrs)."""
    if depth_m is None:
        ur, vr = surface_uv(ds, tindex=tindex)
        ue, vn = rotate_uv(ds, ur, vr)
    else:
        ue, vn = uv_at_depth(ds, depth_m, tindex=tindex, rotate=True)
    spd = speed(np.squeeze(ue), np.squeeze(vn))
    da = xr.DataArray(spd, dims=("eta_rho", "xi_rho"))
    da = da.assign_coords(lon_rho=ds["lon_rho"], lat_rho=ds["lat_rho"])
    _apply_attrs(da, "speed")
    if depth_m is not None:
        da.attrs["depth_m"] = depth_m
    da.name = "speed"
    return da


def vorticity(ds, tindex=-1, depth_m=None, normalized=False):
    """Relative vorticity  zeta = dv/dx - du/dy  as a labeled 2D DataArray.

    Surface if depth_m is None, else at that true depth. If normalized=True,
    returns vorticity / f (Rossby-number-like, dimensionless), which is the usual
    way to see eddies (values ~ +/-1). Uses the grid metrics pm, pn.
    Positive = cyclonic in the Northern Hemisphere.
    """
    if depth_m is None:
        ur, vr = surface_uv(ds, tindex=tindex)     # grid-aligned on rho
    else:
        # grid-aligned (not rotated) at depth, for correct d/dx,d/dy on the grid
        ur, vr = uv_at_depth(ds, depth_m, tindex=tindex, rotate=False)
    ur = np.squeeze(ur); vr = np.squeeze(vr)

    # grid spacing (m): pm=1/dx, pn=1/dy on the rho grid
    if "pm" in ds and "pn" in ds:
        pm = ds["pm"].values; pn = ds["pn"].values
    else:
        # fallback: derive from lon/lat (approximate)
        pm = np.gradient(ds["lon_rho"].values, axis=1)
        pn = np.gradient(ds["lat_rho"].values, axis=0)
        pm = 1.0 / (np.abs(pm) * 111000.0 + 1e-9)
        pn = 1.0 / (np.abs(pn) * 111000.0 + 1e-9)

    # central differences on the rho grid, scaled by metrics
    dvdx = np.full_like(vr, np.nan)
    dudy = np.full_like(ur, np.nan)
    dvdx[:, 1:-1] = (vr[:, 2:] - vr[:, :-2]) * 0.5 * pm[:, 1:-1]
    dudy[1:-1, :] = (ur[2:, :] - ur[:-2, :]) * 0.5 * pn[1:-1, :]
    vort = dvdx - dudy

    _, _, mask = lonlatmask(ds)
    name = "vort"
    if normalized and "f" in ds:
        f = ds["f"].values
        vort = vort / np.where(np.abs(f) > 0, f, np.nan)
        name = "vort_f"
    vort = vort * mask

    da = xr.DataArray(vort, dims=("eta_rho", "xi_rho"))
    da = da.assign_coords(lon_rho=ds["lon_rho"], lat_rho=ds["lat_rho"])
    _apply_attrs(da, name)
    if depth_m is not None:
        da.attrs["depth_m"] = depth_m
    da.name = name
    return da


# ======================================================================
# Time-stamp helper + unified field() dispatcher + eddy_view
# ======================================================================
def _time_str(ds, tindex):
    """A human date string for the given time index, or '' if not decodable."""
    try:
        t = np.atleast_1d(ds["time"].values)[tindex]
        if np.issubdtype(np.asarray(t).dtype, np.datetime64):
            return str(np.datetime_as_string(t, unit="h")).replace("T", " ")
    except Exception:
        pass
    return ""


def _stamp(da, ds, tindex, depth_m=None):
    """Attach time_str (+ depth/surface) attrs so plot titles show them."""
    da.attrs["time_str"] = _time_str(ds, tindex)
    if depth_m is not None:
        da.attrs["depth_m"] = depth_m
    else:
        da.attrs["surface"] = 1
    return da


def field(ds, var, depth_m=None, tindex=-1, level=None):
    """Unified scalar-field extractor with depth + time in the attrs.

    depth_m=None -> surface (top sigma) ; depth_m=<m> -> interpolated to depth.
    `level` (sigma index) still available for model-diagnostic slices.
    Returns a labeled 2D DataArray ready for pl.plot().
    """
    if level is not None:
        da = field_map(ds, var, tindex=tindex, level=level)
        da.attrs["time_str"] = _time_str(ds, tindex)
        return da
    if depth_m is None:
        da = field_map(ds, var, tindex=tindex, level=-1)
    else:
        da = field_at_depth(ds, var, depth_m, tindex=tindex)
    return _stamp(da, ds, tindex, depth_m=depth_m)


def eddy_view(ds, base="speed", overlay="vort", depth_m=None, tindex=-1):
    """A combined eddy diagnostic as (base_da, overlay_da).

    base    : 'speed' (default), 'temp', or any scalar var -> the shaded field
    overlay : 'vort' (vorticity contours) or 'uv' (current vectors)
    Returns a tuple you hand to pl.plot_eddy(...) OR use directly:
        base_da carries lon/lat + attrs (shaded);
        overlay is ('vort', vort_da) or ('uv', (u,v)).
    Both stamped with depth + time for the title.
    """
    if base == "speed":
        base_da = speed_map(ds, tindex=tindex, depth_m=depth_m)
    else:
        base_da = field(ds, base, depth_m=depth_m, tindex=tindex)
    _stamp(base_da, ds, tindex, depth_m=depth_m)

    if overlay == "vort":
        ov = vorticity(ds, tindex=tindex, depth_m=depth_m, normalized=True)
        _stamp(ov, ds, tindex, depth_m=depth_m)
        return base_da, ("vort", ov)
    elif overlay == "uv":
        if depth_m is None:
            u, v = surface_uv(ds, tindex=tindex)
            u, v = rotate_uv(ds, u, v)
        else:
            u, v = uv_at_depth(ds, depth_m, tindex=tindex, rotate=True)
        return base_da, ("uv", (np.squeeze(u), np.squeeze(v)))
    raise ValueError("overlay must be 'vort' or 'uv'")


# ======================================================================
# Interior crop — exclude the sponge band along open boundaries
# ======================================================================
def crop_interior(ds, margin_deg=None, margin_cells=None):
    """Return ds with the boundary sponge band removed.

    The sponge (X_SPONGE in croco.in, ~50 km for Canary) damps the mesoscale in
    a band along the open boundaries. For clean plots/statistics, trim it.

      margin_deg   : trim this many degrees off every edge (e.g. 0.5)
      margin_cells : OR trim this many grid cells off every edge (e.g. 6)
    If neither is given, defaults to ~0.5 deg worth of cells.
    """
    eta = ds.sizes["eta_rho"]; xi = ds.sizes["xi_rho"]
    if margin_cells is None:
        if margin_deg is None:
            margin_deg = 0.5
        # estimate cells-per-degree from the grid
        lat = ds["lat_rho"].values
        dlat = np.abs(np.mean(np.diff(lat[:, lat.shape[1] // 2])))
        margin_cells = max(1, int(round(margin_deg / max(dlat, 1e-6))))
    m = int(margin_cells)
    if 2 * m >= min(eta, xi):
        raise ValueError(f"margin {m} too large for grid {eta}x{xi}")
    return ds.isel(eta_rho=slice(m, eta - m),
                   xi_rho=slice(m, xi - m),
                   eta_v=slice(m, eta - m - 1) if "eta_v" in ds.dims else slice(None),
                   xi_u=slice(m, xi - m - 1) if "xi_u" in ds.dims else slice(None),
                   missing_dims="ignore")


def sponge_mask(ds, margin_cells=6):
    """A 2D mask (1 interior, NaN in the sponge band) for shading/annotating."""
    eta = ds.sizes["eta_rho"]; xi = ds.sizes["xi_rho"]
    m = int(margin_cells)
    mask = np.full((eta, xi), np.nan)
    mask[m:eta - m, m:xi - m] = 1.0
    return mask


# ======================================================================
# Derived 3D fields (all sigma levels) so section/profile/timeseries work
# uniformly for temp/salt/u/v/speed/ke/vorticity. zeta stays 2D.
# ======================================================================
def _var3d(ds, field, tindex=-1):
    """Return a 3D DataArray (s_rho, eta_rho, xi_rho) for `field` at one time.

    Raw fields (temp, salt, u, v) are read (u,v regridded to rho). Derived
    fields are computed on every sigma level:
      speed = sqrt(u^2+v^2)   ke = 0.5*(u^2+v^2)   vort = dv/dx - du/dy
    """
    def _pick(v):
        d = ds[v]
        if "time" in d.dims:
            d = d.isel(time=tindex)
        return d

    if field in ("temp", "salt"):
        return _pick(field)

    if field in ("u", "v"):
        raw = _pick(field).values
        rho = _u2rho(raw) if field == "u" else _v2rho(raw)
        return xr.DataArray(rho, dims=("s_rho", "eta_rho", "xi_rho"))

    if field in ("speed", "ke", "vort", "vort_f"):
        u3 = _u2rho(_pick("u").values)          # (s_rho, eta, xi)
        v3 = _v2rho(_pick("v").values)
        # rotate each level to east/north for physical speed/vort
        ang = ds["angle"].values
        cos, sin = np.cos(ang), np.sin(ang)
        N = u3.shape[0]
        out = np.full_like(u3, np.nan)
        if field in ("speed", "ke"):
            for k in range(N):
                ue = u3[k] * cos - v3[k] * sin
                vn = u3[k] * sin + v3[k] * cos
                out[k] = np.sqrt(ue**2 + vn**2) if field == "speed" else 0.5*(ue**2 + vn**2)
        else:  # vorticity per level (grid-aligned d/dx,d/dy with metrics)
            pm = ds["pm"].values if "pm" in ds else None
            pn = ds["pn"].values if "pn" in ds else None
            for k in range(N):
                dvdx = np.full_like(v3[k], np.nan)
                dudy = np.full_like(u3[k], np.nan)
                if pm is not None:
                    dvdx[:, 1:-1] = (v3[k][:, 2:] - v3[k][:, :-2]) * 0.5 * pm[:, 1:-1]
                    dudy[1:-1, :] = (u3[k][2:, :] - u3[k][:-2, :]) * 0.5 * pn[1:-1, :]
                vort = dvdx - dudy
                if field == "vort_f" and "f" in ds:
                    f = ds["f"].values
                    vort = vort / np.where(np.abs(f) > 0, f, np.nan)
                out[k] = vort
        return xr.DataArray(out, dims=("s_rho", "eta_rho", "xi_rho"))

    # fallback: read as-is
    return _pick(field)


def _canon_derived(field):
    """Attribute-registry key for a derived/raw field name."""
    return {"ke": "eke", "vort": "vort", "vort_f": "vort_f",
            "speed": "speed"}.get(field, field)
