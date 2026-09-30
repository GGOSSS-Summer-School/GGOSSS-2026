"""
gtools/plotting.py — CROCO plotting (attribute-driven, generic).

Two layers:
  1. Generic builders — plot any labeled DataArray:
        plot_map(da)          horizontal field on a map
        plot_section(da)      vertical section (distance vs depth)
        plot_profile(da)      single vertical profile (value vs depth)
        plot_hovmoller(da)    2D time-vs-(depth|lat|lon) diagram
  2. Smart wrapper — plot(da) inspects the DataArray's dims and calls the right
     builder automatically.

Labels are read from CF attributes by default (long_name, units) and can be
overridden by keyword (title=, cbar_label=, xlabel=, ylabel=). Pair with the
extractors in gtools.postprocess (field_map, section, profile, hovmoller),
which return labeled DataArrays.

Typical use
-----------
    import gtools.postprocess as pp
    import gtools.plotting    as pl

    ds = pp.open_history("croco_his.nc")

    pl.plot(pp.field_map(ds, "temp"), out="sst.png")                       # map
    pl.plot(pp.section(ds, "temp", -20,21, -17,21), out="sec.png")        # section
    pl.plot(pp.profile(ds, "temp", -19,21), out="prof.png")               # profile
    pl.plot(pp.hovmoller(ds, "temp", "time_lat", lon0=-19), out="hov.png")# Hovmoller
"""
from __future__ import annotations
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

try:
    from gtools.define_attrs import resolve_limits
except Exception:
    try:
        from define_attrs import resolve_limits
    except Exception:
        def resolve_limits(da, cmap=None, vmin=None, vmax=None):
            return (cmap or "viridis", vmin, vmax)

try:
    import cartopy.crs as ccrs
    import cartopy.feature as cfeature
    _HAS_CARTOPY = True
except Exception:
    _HAS_CARTOPY = False


# ----------------------------------------------------------------------
# Shared map dressing + colour-limit convention, used by every map-drawing
# function across gtools (validation.py, validation_satellite.py,
# validation_godae.py, plotting.py itself) so every figure in the toolkit
# looks and behaves the same way.
# ----------------------------------------------------------------------
def add_map_features(ax, coastline_res="10m"):
    """Decorate a cartopy GeoAxes with the standard CROCO map dressing:
    land, lakes, rivers, national/political borders, coastline, and labelled
    lon/lat gridlines (tick labels in degrees).

    No-op (returns False) if cartopy isn't available -- callers should still
    work without it, just without the extra map layers.
    """
    if not _HAS_CARTOPY:
        return False
    ax.add_feature(cfeature.LAND, facecolor="0.85", zorder=3)
    ax.add_feature(cfeature.LAKES, facecolor="#a6cee3", edgecolor="0.4",
                   linewidth=0.4, zorder=3)
    ax.add_feature(cfeature.RIVERS, edgecolor="#4292c6", linewidth=0.5, zorder=4)
    ax.add_feature(cfeature.BORDERS, edgecolor="0.3", linewidth=0.5,
                   linestyle=":", zorder=4)
    ax.coastlines(resolution=coastline_res, linewidth=0.5, zorder=4)
    gl = ax.gridlines(draw_labels=True, linewidth=0.3, color="0.6", alpha=0.4)
    gl.top_labels = False
    gl.right_labels = False
    if hasattr(gl, "geo_labels"):
        # newer cartopy: geo_labels (map-edge labels) triggers the same
        # boundary-polygon code as top/right_labels and can crash with
        # "Points of LinearRing do not form a closed linestring" on some
        # cartopy/shapely combinations -- disable it too, not just top/right.
        gl.geo_labels = False
    return True


def percentile_clim(*fields, pct=(2, 98)):
    """Shared colour-limit convention used across gtools: the (2nd, 98th)
    percentile (by default) of all finite values pooled across one or more
    fields -- robust to the odd outlier pixel, unlike plain nanmin/nanmax.

    Returns (vmin, vmax) as floats, or (None, None) if every field is empty.
    """
    chunks = [np.asarray(f)[np.isfinite(np.asarray(f))].ravel()
             for f in fields if f is not None]
    chunks = [c for c in chunks if c.size]
    if not chunks:
        return None, None
    vals = np.concatenate(chunks)
    lo, hi = np.nanpercentile(vals, pct)
    return float(lo), float(hi)


# ----------------------------------------------------------------------
# attribute helpers — read CF metadata for auto-labels
# ----------------------------------------------------------------------
def _label(da, override=None):
    """A '<long_name> (<units>)' label from attrs, unless overridden."""
    if override is not None:
        return override
    ln = da.attrs.get("long_name", da.name or "")
    un = da.attrs.get("units", "")
    return f"{ln} ({un})" if un else str(ln)


def _short(da, override=None):
    if override is not None:
        return override
    base = da.attrs.get("long_name", da.name or "")
    extra = []
    # point / position stamp
    if "lon0" in da.attrs and "lat0" in da.attrs:
        extra.append(f"lon: {da.attrs['lon0']:g}°, lat: {da.attrs['lat0']:g}°")
    elif "lon0" in da.attrs:
        extra.append(f"lon: {da.attrs['lon0']:g}°")
    elif "lat0" in da.attrs:
        extra.append(f"lat: {da.attrs['lat0']:g}°")
    elif "point" in da.attrs and da.attrs["point"]:
        extra.append(str(da.attrs["point"]))
    # depth stamp
    if "depth_m" in da.attrs:
        extra.append(f"{da.attrs['depth_m']:g} m")
    elif da.attrs.get("surface", 0):
        extra.append("surface")
    # time stamp
    if "time_str" in da.attrs and da.attrs["time_str"]:
        extra.append(da.attrs["time_str"])
    if extra:
        return f"{base}  —  " + "  ".join(extra)
    return base


def _finish(fig, out, dpi=150):
    try:
        fig.tight_layout()
    except Exception:
        pass
    if out:
        fig.savefig(out, dpi=dpi, bbox_inches="tight")
        plt.close(fig)
        return out
    return fig

def get_extend(da, vmin, vmax):
    """
    Determines the appropriate colorbar 'extend' value based on
    whether data values fall outside the [vmin, vmax] range.

    Parameters
    ----------
    da    : array-like  
    vmin  : float       lower bound of the colorbar
    vmax  : float       upper bound of the colorbar

    Returns
    -------
    str : 'neither', 'min', 'max', or 'both'
    """
    data = np.asarray(da.values) if hasattr(da, "values") else np.asarray(da)
    finite = data[np.isfinite(data)]
    below = float(np.nanmin(finite)) < vmin
    above = float(np.nanmax(finite)) > vmax

    if below and above:
        return "both"
    elif below:
        return "min"
    elif above:
        return "max"
    else:
        return "neither"
    
def feature_map(ax):
    import cartopy.feature as cfeature
    """Add features map"""
    ax.add_feature(cfeature.LAND,      facecolor="lightgray", edgecolor="darkgray", linewidth=0.5, zorder=3)
    ax.add_feature(cfeature.COASTLINE, linewidth=0.8,          edgecolor="k",        zorder=4)
    ax.add_feature(cfeature.BORDERS,   edgecolor="gray",       linewidth=0.7,        linestyle="-", zorder=6)
    ax.add_feature(cfeature.RIVERS,    edgecolor="blue",       linewidth=0.6,        zorder=6)
    ax.add_feature(cfeature.LAKES,     facecolor="lightblue",  edgecolor="steelblue",linewidth=0.4, zorder=6)


# ----------------------------------------------------------------------
# 1a. MAP — horizontal field
# ----------------------------------------------------------------------
def plot_map(da, out=None, cmap=None, vmin=None, vmax=None,
             title=None, cbar_label=None, figsize=(8, 8), coastline=True,
             ds=None, isobaths=None, uv=None, uv_kind="current",
             uv_scale=None, uv_skip=None, uv_ref=None):
    """Map a 2D field. Expects lon_rho/lat_rho coords (from pp.field_map).
    cmap/vmin/vmax default to the variable's attributes (override to force).

    Optional overlays (need `ds`, the CROCO dataset, for grid/bathymetry):
      isobaths : list of depths, e.g. [200, 1000, 2000] -> contour lines
      uv       : (u_rho, v_rho) tuple -> current vectors with reference arrow
    """
    cmap, vmin, vmax = resolve_limits(da, cmap, vmin, vmax)
    lon = da["lon_rho"].values
    lat = da["lat_rho"].values
    field = da.values

    if _HAS_CARTOPY and coastline:
        fig = plt.figure(figsize=figsize)
        ax = plt.axes(projection=ccrs.PlateCarree())
        ax.set_extent([lon.min(), lon.max(), lat.min(), lat.max()],
                      crs=ccrs.PlateCarree())
        add_map_features(ax)
        h = ax.pcolormesh(lon, lat, field, cmap=cmap, vmin=vmin, vmax=vmax,
                          shading="auto", transform=ccrs.PlateCarree(), zorder=1)
    else:
        fig, ax = plt.subplots(figsize=figsize)
        h = ax.pcolormesh(lon, lat, field, cmap=cmap, vmin=vmin, vmax=vmax, shading="auto")
        ax.set_xlabel("longitude"); ax.set_ylabel("latitude")

    cb = fig.colorbar(h, ax=ax, shrink=0.8,extend = get_extend(da, vmin, vmax), pad=0.03)
    cb.set_label(_label(da, cbar_label))

    # optional overlays
    if isobaths is not None and ds is not None:
        add_isobaths(ax, ds, isobaths=isobaths)
    if uv is not None and ds is not None:
        u_rho, v_rho = uv
        ext = [lon.min(), lon.max(), lat.min(), lat.max()]
        add_uv(ax, ds, u_rho, v_rho, extents=ext, kind=uv_kind,
               scale=uv_scale, skip=uv_skip, ref_vector=uv_ref)
    # adding feature map to the plot (land,coastline,borders,rivers,lakes)  
    feature_map(ax)
    ax.set_title(title or _short(da))
    return _finish(fig, out)


# ----------------------------------------------------------------------
# 1b. SECTION — vertical slice (distance vs depth)
# ----------------------------------------------------------------------
def plot_section(da, out=None, cmap=None, vmin=None, vmax=None,
                 title=None, cbar_label=None, figsize=(10, 5),contour_colors="k", levels=None):
    """Vertical section. Expects dims (s_rho, points) with 'distance_km' and a
    2D 'depth' coord (from pp.section). cmap/vmin/vmax from attrs by default."""
    cmap, vmin, vmax = resolve_limits(da, cmap, vmin, vmax)
    dist = da["distance_km"].values                       # (points,)
    depth = np.array(da["depth"].values, dtype=float)     # (s_rho, points)
    field = da.values                                     # (s_rho, points)
    # pcolormesh needs finite depth coords; land columns are NaN. Fill the coord
    # (keep the data NaN so land renders blank).
    if not np.isfinite(depth).all():
        col_template = np.nanmean(depth, axis=1, keepdims=True)
        depth = np.where(np.isfinite(depth),
                         depth, np.broadcast_to(col_template, depth.shape))
        if not np.isfinite(depth).all():
            ramp = np.linspace(np.nanmin(depth), 0, depth.shape[0])[:, None]
            depth = np.where(np.isfinite(depth), depth,
                             np.broadcast_to(ramp, depth.shape))
    X = np.tile(dist, (field.shape[0], 1))

    fig, ax = plt.subplots(figsize=figsize)
    h = ax.pcolormesh(X, depth, field, cmap=cmap, vmin=vmin, vmax=vmax, shading="auto")
    cs = ax.contour(
        X, depth, field,
        levels=levels,
        colors=contour_colors,
        linewidths=0.7
    )
    ax.clabel(cs, fontsize=7, inline=True, fmt="%.2f")
    ax.set_xlabel("distance along section (km)")
    ax.set_ylabel("depth (m)")
    cb = fig.colorbar(h, ax=ax, shrink=0.9,extend = get_extend(da, vmin, vmax), pad=0.02)
    cb.set_label(_label(da, cbar_label))
    ax.set_title(title or f"{_short(da)} — section")
    return _finish(fig, out)

# ----------------------------------------------------------------------
# 1c. PROFILE — single vertical profile (value vs depth)
# ----------------------------------------------------------------------
def plot_profile(da, out=None, title=None, xlabel=None, figsize=(4.5, 6),
                 color="C0", marker="o"):
    """Vertical profile. Expects a 'depth' coord (from pp.profile)."""
    depth = da["depth"].values
    val = da.values
    fig, ax = plt.subplots(figsize=figsize)
    ax.plot(val, depth, color=color, marker=marker, ms=3, lw=1.5)
    ax.set_xlabel(_label(da, xlabel))
    ax.set_ylabel("depth (m)")
    ax.set_title(title or _short(da))
    return _finish(fig, out)


# ----------------------------------------------------------------------
# 1d. HOVMOLLER — time vs (depth | lat | lon)
# ----------------------------------------------------------------------
def plot_hovmoller(da, out=None, cmap=None, vmin=None, vmax=None,
                   title=None, cbar_label=None, figsize=(10, 5)):
    """Hovmoller. Expects (time, <second axis>) where the second axis is one of
    s_rho (with 'depth'), eta_rho (with 'lat') or xi_rho (with 'lon')
    (from pp.hovmoller). cmap/vmin/vmax from attrs by default."""
    cmap, vmin, vmax = resolve_limits(da, cmap, vmin, vmax)
    t = da["time"].values
    dims = [d for d in da.dims if d != "time"]
    if not dims:
        raise ValueError("hovmoller needs a time dim and one space/depth dim")
    ydim = dims[0]

    # choose the y coordinate + label from what the extractor attached
    if "depth" in da.coords:
        y = da["depth"].values; ylab = "depth (m)"
    elif "lat" in da.coords:
        y = da["lat"].values; ylab = "latitude"
    elif "lon" in da.coords:
        y = da["lon"].values; ylab = "longitude"
    else:
        y = np.arange(da.sizes[ydim]); ylab = ydim

    field = da.transpose(ydim, "time").values            # (y, time)
    T = np.tile(t, (field.shape[0], 1))
    Y = np.tile(y.reshape(-1, 1), (1, field.shape[1]))

    fig, ax = plt.subplots(figsize=figsize)
    h = ax.pcolormesh(T, Y, field, cmap=cmap, vmin=vmin, vmax=vmax, shading="auto")
    ax.set_xlabel("time"); ax.set_ylabel(ylab)
    cb = fig.colorbar(h, ax=ax, shrink=0.9, extend=get_extend(da, vmin, vmax), pad=0.02)
    cb.set_label(_label(da, cbar_label))
    ax.set_title(title or f"{_short(da)} — Hovmoller")
    fig.autofmt_xdate()
    return _finish(fig, out)


# ----------------------------------------------------------------------
# 2. SMART WRAPPER — auto-detect what to draw
# ----------------------------------------------------------------------
def plot(da, out=None, **kw):
    """Inspect a labeled DataArray and dispatch to the right builder.

    - has lon_rho & lat_rho, no s_rho          -> plot_map
    - has 'distance_km' + 2D 'depth'           -> plot_section
    - has a time dim + one other               -> plot_hovmoller
    - has a 'depth' coord, 1D                   -> plot_profile
    """
    coords = set(da.coords)
    dims = set(da.dims)

    if da.ndim == 1 and "time" in dims and "depth" not in coords:
        return plot_timeseries(da, out=out, **kw)
    if "distance_km" in coords:
        return plot_section(da, out=out, **kw)
    if "time" in dims and len(dims) == 2:
        return plot_hovmoller(da, out=out, **kw)
    if "lon_rho" in coords and "lat_rho" in coords and "s_rho" not in dims:
        return plot_map(da, out=out, **kw)
    if "depth" in coords and da.ndim == 1:
        return plot_profile(da, out=out, **kw)
    # fallback: if 2D with lon/lat treat as map
    if "lon_rho" in coords and "lat_rho" in coords:
        return plot_map(da, out=out, **kw)
    raise ValueError(f"plot(): can't infer plot type from dims={da.dims}, coords={list(da.coords)}")


# ======================================================================
# Composable building blocks (adapted from somisana plotting):
#   add_isobaths, add_uv (with reference vector + auto-scaling), add_colorbar
# Use these to enrich a plot_map axes, or call directly on your own axes.
# ======================================================================
def _uv_autoparams(spd, scale=None, ref_vector=None, skip=None,
                   num_vectors=20, aspect=1.0, kind="current"):
    """Auto-scale vector density + reference size from the speed field.
    kind='current' (m/s ~0.1-1) or 'wind' (m/s ~2-15). Adapted from somisana
    get_uv_params, extended with a wind regime."""
    if skip is None:
        Ny = spd.shape[0]
        skip = max(1, int(Ny / num_vectors * aspect))
    mx = float(np.nanmax(spd)) if np.isfinite(spd).any() else 0.0
    if kind == "wind":
        # wind speeds are ~an order of magnitude larger than currents
        if mx < 5:     s, r = 100.0, 2.0
        elif mx < 10:  s, r = 200.0, 5.0
        elif mx < 20:  s, r = 400.0, 10.0
        else:          s, r = 600.0, 15.0
    else:
        if mx < 0.25:   s, r = 2.5, 0.1
        elif mx < 0.5:  s, r = 5.0, 0.25
        elif mx < 0.75: s, r = 7.5, 0.5
        elif mx < 1.0:  s, r = 10.0, 0.5
        else:           s, r = 15.0, 1.0
    if scale is None:      scale = s
    if ref_vector is None: ref_vector = r
    return scale, skip, ref_vector


def add_isobaths(ax, ds, isobaths=(200, 1000, 2000), color="k",
                 linestyles="dashed", linewidths=1.0, label=True):
    """Overlay bathymetry contours (isobaths) on a map axes. `ds` is a CROCO
    dataset with h/lon_rho/lat_rho. Adapted from somisana plot_isobaths."""
    h = ds["h"].values
    lon = ds["lon_rho"].values
    lat = ds["lat_rho"].values
    kw = {}
    if _HAS_CARTOPY:
        kw["transform"] = ccrs.PlateCarree()
    cs = ax.contour(lon, lat, h, levels=list(isobaths), colors=color,
                    linestyles=linestyles, linewidths=linewidths, zorder=5, **kw)
    if label:
        ax.clabel(cs, fmt="%1.0f", inline=True, fontsize=7)
    return cs


def add_uv(ax, ds, u_rho, v_rho, extents=None, skip=None, scale=None,
           ref_vector=None, color="k", num_vectors=20, kind="current"):
    """Overlay vectors (rho grid) with a labeled reference vector.
    kind='current' or 'wind' picks a sensible magnitude scaling. Adapted from
    somisana plot_uv + get_uv_params."""
    lon = ds["lon_rho"].values
    lat = ds["lat_rho"].values
    u_rho = np.squeeze(np.asarray(u_rho))
    v_rho = np.squeeze(np.asarray(v_rho))
    spd = np.sqrt(u_rho ** 2 + v_rho ** 2)
    if extents is None:
        extents = [float(lon.min()), float(lon.max()), float(lat.min()), float(lat.max())]
    aspect = (extents[3] - extents[2]) / max(1e-9, (extents[1] - extents[0]))
    scale, skip, ref_vector = _uv_autoparams(spd, scale, ref_vector, skip,
                                             num_vectors, aspect, kind=kind)
    kw = dict(scale=scale, color=color, width=0.0035, zorder=6)
    if _HAS_CARTOPY:
        kw["transform"] = ccrs.PlateCarree()
    s = slice(None, None, skip)
    q = ax.quiver(lon[s, s], lat[s, s], u_rho[s, s], v_rho[s, s], **kw)
    # reference vector, top-left
    dx, dy = extents[1] - extents[0], extents[3] - extents[2]
    lx, ly = extents[0] + dx * 0.06, extents[3] - dy * 0.06
    kw2 = dict(scale=scale, color=color, width=0.0035, zorder=7)
    if _HAS_CARTOPY:
        kw2["transform"] = ccrs.PlateCarree()
    ax.quiverkey(
        q,
        X=1.12, 
        Y=1.0,
        U=ref_vector,
        label=f"{ref_vector} m s$^{-1}$",
        labelpos="S",
        coordinates="axes"
    )
    txt_kw = {}
    if _HAS_CARTOPY:
        txt_kw["transform"] = ccrs.PlateCarree()
    return q


# ======================================================================
# Combined eddy view: shaded base field + vorticity contours or current vectors
# ======================================================================
def plot_eddy(base_da, overlay, ds=None, out=None, title=None,
              isobaths=None, uv_kind="current", uv_scale=None, uv_skip=None,
              uv_ref=None, n_contours=8, figsize=(8, 8)):
    """Shaded base field (from eddy_view) + an overlay.

    overlay : ('vort', vort_da)  -> vorticity contours (from vort/f)
              ('uv', (u, v))     -> current vectors
    base_da carries its own cmap/vmin/vmax + depth/time in the title.
    """
    cmap, vmin, vmax = resolve_limits(base_da)
    lon = base_da["lon_rho"].values
    lat = base_da["lat_rho"].values
    field = base_da.values

    if _HAS_CARTOPY:
        fig = plt.figure(figsize=figsize)
        ax = plt.axes(projection=ccrs.PlateCarree())
        ax.set_extent([lon.min(), lon.max(), lat.min(), lat.max()], crs=ccrs.PlateCarree())
        add_map_features(ax)
        h = ax.pcolormesh(lon, lat, field, cmap=cmap, vmin=vmin, vmax=vmax,
                          shading="auto", transform=ccrs.PlateCarree(), zorder=1)
    else:
        fig, ax = plt.subplots(figsize=figsize)
        h = ax.pcolormesh(lon, lat, field, cmap=cmap, vmin=vmin, vmax=vmax, shading="auto")

    cb = fig.colorbar(h, ax=ax, shrink=0.8,extend = get_extend(base_da, vmin, vmax), pad=0.03)
    cb.set_label(_label(base_da))

    kind, data = overlay
    if kind == "vort":
        vd = data
        vlon = vd["lon_rho"].values; vlat = vd["lat_rho"].values
        vmax_c = np.nanmax(np.abs(vd.values))
        levels = np.linspace(-vmax_c, vmax_c, n_contours)
        ckw = {}
        if _HAS_CARTOPY:
            ckw["transform"] = ccrs.PlateCarree()
        cs = ax.contour(vlon, vlat, vd.values, levels=levels,
                        colors="k", linewidths=0.6, alpha=0.7, zorder=5, **ckw)
    elif kind == "uv":
        u, v = data
        ext = [lon.min(), lon.max(), lat.min(), lat.max()]
        add_uv(ax, ds, u, v, extents=ext, kind=uv_kind,
               scale=uv_scale, skip=uv_skip, ref_vector=uv_ref)

    if isobaths is not None and ds is not None:
        add_isobaths(ax, ds, isobaths=isobaths)
    # adding feature map to the plot (land,coastline,borders,rivers,lakes)  
    feature_map(ax)
    ax.set_title(title or _short(base_da))
    return _finish(fig, out)


# ----------------------------------------------------------------------
# TIME SERIES — value vs time at a point
# ----------------------------------------------------------------------
def plot_timeseries(da, out=None, title=None, ylabel=None, figsize=(9, 4),
                    color="C0", marker="o"):
    """Time series. Expects a 1D DataArray over 'time' (from pp.timeseries)."""
    t = da["time"].values
    fig, ax = plt.subplots(figsize=figsize)
    ax.plot(t, da.values, color=color, marker=marker, ms=3, lw=1.4)
    ax.set_xlabel("time")
    ax.set_ylabel(_label(da, ylabel))
    ax.set_title(title or _short(da))
    fig.autofmt_xdate()
    return _finish(fig, out)


# ======================================================================
# Region portrait: grid mesh + bathymetry, side by side (for the gallery)
# ======================================================================
def grid_bathy_map(grd, out=None, title=None, coastline=True, mesh_stride=1,
                   bathy_levels=(0, 20, 50, 100, 200, 300, 500, 1000, 2000, 3000),
                   figsize=(13, 7)):
    """Two-panel 'region portrait' from a CROCO grid file: grid mesh | bathymetry.

    grd: path to croco_grd.nc (or an open dataset).
    Left panel  = the model grid mesh over the coastline.
    Right panel = smoothed bathymetry (h) with a shelf-friendly log-ish scale.
    coastline: overlay a real coastline (needs cartopy; falls back to the land mask).
    mesh_stride: draw every Nth grid line on the left (>1 thins a dense mesh).
    """
    import xarray as xr
    ds = xr.open_dataset(grd) if isinstance(grd, str) else grd
    lon = ds["lon_rho"].values
    lat = ds["lat_rho"].values
    h = ds["h"].values
    mask = ds["mask_rho"].values if "mask_rho" in ds else np.ones_like(h)
    hplot = np.where(mask > 0, h, np.nan)     # blank land

    # cartopy if available
    try:
        import cartopy.crs as ccrs
        # import cartopy.feature as cfeature
        use_cart = coastline
    except Exception:
        use_cart = False
        if coastline:
            print("coastline requested but cartopy not available - grid edge shows the coast")

    extent = [float(np.nanmin(lon)), float(np.nanmax(lon)),
              float(np.nanmin(lat)), float(np.nanmax(lat))]
    subkw = dict(projection=ccrs.PlateCarree()) if use_cart else {}
    fig, (axg, axb) = plt.subplots(1, 2, figsize=figsize, subplot_kw=subkw,
                                   constrained_layout=True)

    # ---- left: grid mesh (OCEAN cells only; land is left blank) ----
    s = max(1, int(mesh_stride))
    tf = dict(transform=ccrs.PlateCarree()) if use_cart else {}
    # break each grid line wherever it crosses land: set land points to NaN so
    # matplotlib lifts the pen there. This keeps the mesh over the computational
    # (ocean) domain and stops it from drawing across the continent.
    lon_o = np.where(mask > 0, lon, np.nan)
    lat_o = np.where(mask > 0, lat, np.nan)
    for j in range(0, lon.shape[0], s):
        axg.plot(lon_o[j, :], lat_o[j, :], color="0.5", lw=0.3, **tf)
    for i in range(0, lon.shape[1], s):
        axg.plot(lon_o[:, i], lat_o[:, i], color="0.5", lw=0.3, **tf)
    # faint fill of the ocean domain so the mesh region reads as "the model"
    axg.pcolormesh(lon, lat, np.where(mask > 0, 1.0, np.nan),
                   cmap="Blues", vmin=0, vmax=3, shading="auto", alpha=0.12, **tf)
    axg.set_title("model grid")

    # ---- right: bathymetry ----
    from matplotlib.colors import BoundaryNorm
    levels = np.array(bathy_levels, dtype=float)
    norm = BoundaryNorm(levels, ncolors=256, extend="max")
    tf = dict(transform=ccrs.PlateCarree()) if use_cart else {}
    hb = axb.pcolormesh(lon, lat, hplot, cmap="YlGnBu", norm=norm, shading="auto", **tf)
    axb.set_title("bathymetry")
    cb = fig.colorbar(hb, ax=axb, shrink=0.85, pad=0.02, ticks=levels)
    cb.set_label("bathymetry (m)")

    for ax in (axg, axb):
        if use_cart:
            # ax.add_feature(cfeature.COASTLINE, linewidth=0.7, zorder=3)
            # adding feature map to the plot (land,coastline,borders,rivers,lakes)  
            feature_map(ax)
            ax.set_extent(extent, crs=ccrs.PlateCarree())
            gl = ax.gridlines(draw_labels=False, linewidth=0.3, color="0.6",
                              linestyle=":", alpha=0.6)
            gl.bottom_labels = True; gl.left_labels = True
            gl.top_labels = False; gl.right_labels = False
        else:
            ax.set_xlim(extent[0], extent[1]); ax.set_ylim(extent[2], extent[3])
            ax.set_aspect("equal", adjustable="box")
            ax.set_xlabel("longitude"); ax.grid(alpha=0.3, ls=":")
        if not use_cart and ax is axg:
            ax.set_ylabel("latitude")

    if title:
        fig.suptitle(title, fontsize=14)
    if isinstance(grd, str):
        ds.close()
    if out:
        fig.savefig(out, dpi=150, bbox_inches="tight")
        plt.close(fig)
        return out
    return fig
