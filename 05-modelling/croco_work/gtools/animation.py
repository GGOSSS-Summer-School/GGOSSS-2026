"""
CROCO animation helper.

Builds an animation of a CROCO field through time, reusing the colour scales,
labels and vector overlays of the static plotting module. Plays inline in a
Jupyter notebook, or writes a GIF/MP4 with ``out=``.

    import gtools.postprocess as pp
    import gtools.animation as anim

    ds = pp.open_history("forecast/scratch/Canary_12/CROCO_FILES/croco_his.nc",
                         Yorig=2000)
    anim.animate(ds, "temp", overlay="wind")                     # inline widget
    anim.animate(ds, "temp", overlay="wind", out="sst.gif")      # write a file
"""
from __future__ import annotations
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import xarray as xr

try:
    from IPython.display import HTML
    _HAS_IPYTHON = True
except Exception:
    _HAS_IPYTHON = False

try:
    import cartopy.crs as ccrs
    import cartopy.feature as cfeature
    _HAS_CARTOPY = True
except Exception:
    _HAS_CARTOPY = False

import gtools.postprocess as pp
import gtools.plotting as pl


# ----------------------------------------------------------------------
# small local helpers (kept here so this module needs nothing new from
# plotting.py)
# ----------------------------------------------------------------------
def _extend_for(da, vmin, vmax):
    """Colorbar 'extend' value: does the data fall outside [vmin, vmax]?"""
    data = np.asarray(da.values) if hasattr(da, "values") else np.asarray(da)
    finite = data[np.isfinite(data)]
    if finite.size == 0:
        return "neither"
    below = float(np.nanmin(finite)) < vmin
    above = float(np.nanmax(finite)) > vmax
    if below and above:
        return "both"
    if below:
        return "min"
    if above:
        return "max"
    return "neither"


def _decorate(ax):
    """Land, coastline and gridlines — matching plot_map's decoration."""
    if not _HAS_CARTOPY:
        return
    ax.add_feature(cfeature.LAND, facecolor="0.85", zorder=3)
    ax.coastlines(resolution="10m", linewidth=0.6, zorder=4)
    gl = ax.gridlines(draw_labels=False, linewidth=0.3, color="0.6", alpha=0.4)
    gl.bottom_labels = True
    gl.left_labels = True
    gl.top_labels = False
    gl.right_labels = False


# ----------------------------------------------------------------------
# main entry point
# ----------------------------------------------------------------------
def animate(ds, var, depth_m=None, overlay=None, uv_depth=None, isobaths=None, tindex_range=None, skip=4,
            scale=None, interval=300, figsize=(8, 7), cmap=None,
            vmin=None, vmax=None, out=None, fps=4, dpi=110):
    """Animate a 2D field through the run.

    Parameters
    ----------
    ds : xarray.Dataset
        CROCO output opened with ``pp.open_history()`` or ``pp.open_run()``.
    var : str
        ``'temp'``, ``'salt'``, ``'zeta'``, ``'speed'``, ``'u'``, ``'v'``, ...
    overlay : str or None
        ``'wind'`` / ``'sustr'`` — surface wind-stress vectors, taken from the
        model's own ``sustr``/``svstr`` (N/m², not 10 m wind speed).
        ``'uv'`` / ``'current'`` — surface current vectors.
        ``None`` — no overlay.
    isobaths : list of float, optional
        Depths (m) to contour, e.g. ``[200, 1000]``.
    tindex_range : (start, end), optional
        Range of time indices. Default: every record.
    skip : int
        Vector subsampling stride (default 4).
    scale : float, optional
        Quiver scale. Smaller values give longer arrows.
    interval : int
        Frame delay in milliseconds, for the inline widget.
    figsize : tuple
        Figure size in inches.
    cmap, vmin, vmax : optional
        Colour overrides. Limits are resolved **once** across the whole series,
        so the scale does not flicker between frames.
    out : str, optional
        Write to this path instead of returning a widget. ``.gif`` uses the
        pillow writer, anything else uses ffmpeg.
    fps : int
        Frames per second when writing a file.
    dpi : int
        Resolution when writing a file.

    Returns
    -------
    str
        The output path, if ``out`` was given.
    IPython.display.HTML
        An inline animation otherwise.
    """
    plt.rcParams["animation.html"] = "jshtml"

    lon2d, lat2d, mask = pp.lonlatmask(ds)
    nt = ds.sizes["time"]
    t_indices = (list(range(nt)) if tindex_range is None
                 else list(range(tindex_range[0], tindex_range[1])))
    if not t_indices:
        raise ValueError("tindex_range selected no frames")

    da_list, u_list, v_list = [], [], []

    for t in t_indices:
        if var == "speed":
            da = pp.speed_map(ds, tindex=t, depth_m=depth_m)
        elif var in ("u", "v"):
            ur, vr = pp.surface_uv(ds, tindex=t)
            ue, vn = pp.rotate_uv(ds, ur, vr)
            vals = ue if var == "u" else vn
            da = xr.DataArray(vals, dims=("eta_rho", "xi_rho"),
                              coords={"lon_rho": ds["lon_rho"],
                                      "lat_rho": ds["lat_rho"]})
            pp._apply_attrs(da, var, rotated=True)
        else:
            da = pp.field(ds, var, depth_m=depth_m, tindex=t)

        pp._stamp(da, ds, tindex=t)
        da_list.append(da)

        if overlay in ("wind", "sustr"):
            if "sustr" not in ds or "svstr" not in ds:
                raise KeyError(
                    "overlay='wind' needs sustr/svstr in the history file. "
                    "Enable them in croco.in's history output block."
                )
            su = pp._u2rho(ds["sustr"].isel(time=t).values)
            sv = pp._v2rho(ds["svstr"].isel(time=t).values)
            se, sn = pp.rotate_uv(ds, su, sv)
            u_list.append(se * mask)
            v_list.append(sn * mask)
        elif overlay in ("uv", "current"):
            if uv_depth is None:
                ur, vr = pp.surface_uv(ds, tindex=t)
                ue, vn = pp.rotate_uv(ds, ur, vr)
            else:
                ue, vn = pp.uv_at_depth(ds, depth_m=uv_depth, tindex=t,
                                        rotate=True)
            u_list.append(ue)
            v_list.append(vn)

    # one colour scale for the whole series, so it does not flicker
    da_concat = xr.concat(da_list, dim="time")
    cmap_r, vmin_r, vmax_r = pl.resolve_limits(da_concat, cmap=cmap,
                                               vmin=vmin, vmax=vmax)
    extend_r = _extend_for(da_concat, vmin_r, vmax_r)

    def _title_for(da_frame):
        time_str = da_frame.attrs.get("time_str", "")
        if var in ("temp", "sst"):
            base = "SST"
            if overlay in ("wind", "sustr"):
                base += " and wind stress"
        elif var in ("zeta", "ssh"):
            base = "Sea Surface Height"
            if overlay in ("uv", "current"):
                base += " and currents"
        elif var == "speed":
            base = "Current speed — surface"
        elif var == "u":
            base = "Zonal current (u) — surface"
        elif var == "v":
            base = "Meridional current (v) — surface"
        else:
            return pl._short(da_frame)
        return f"{base} — {time_str}" if time_str else base

    # ---- figure ------------------------------------------------------
    if _HAS_CARTOPY:
        fig, ax = plt.subplots(figsize=figsize,
                               subplot_kw={"projection": ccrs.PlateCarree()})
        ax.set_extent([lon2d.min(), lon2d.max(), lat2d.min(), lat2d.max()],
                      crs=ccrs.PlateCarree())
        _decorate(ax)
        mesh = ax.pcolormesh(lon2d, lat2d, da_list[0].values, cmap=cmap_r,
                             vmin=vmin_r, vmax=vmax_r, shading="auto",
                             transform=ccrs.PlateCarree(), zorder=1)
    else:
        fig, ax = plt.subplots(figsize=figsize)
        mesh = ax.pcolormesh(lon2d, lat2d, da_list[0].values, cmap=cmap_r,
                             vmin=vmin_r, vmax=vmax_r, shading="auto")
        ax.set_xlabel("longitude")
        ax.set_ylabel("latitude")

    cb = fig.colorbar(mesh, ax=ax, shrink=0.8, pad=0.03, extend=extend_r)
    cb.set_label(pl._label(da_list[0]))

    if isobaths is not None:
        pl.add_isobaths(ax, ds, isobaths=isobaths)

    q_overlay = None
    if overlay is not None and u_list:
        q_color = "white" if var == "speed" else "black"
        q_scale = scale if scale is not None else (
            1.8 if overlay in ("wind", "sustr") else 8.0)
        q_width = 0.004 if var == "speed" else 0.0035
        kw_q = dict(scale=q_scale, color=q_color, width=q_width, zorder=5)
        if _HAS_CARTOPY:
            kw_q["transform"] = ccrs.PlateCarree()
        q_overlay = ax.quiver(lon2d[::skip, ::skip], lat2d[::skip, ::skip],
                              u_list[0][::skip, ::skip],
                              v_list[0][::skip, ::skip], **kw_q)

    title = ax.set_title(_title_for(da_list[0]))

    def update(frame):
        da_f = da_list[frame]
        mesh.set_array(da_f.values.ravel())
        if q_overlay is not None:
            q_overlay.set_UVC(u_list[frame][::skip, ::skip],
                              v_list[frame][::skip, ::skip])
        title.set_text(_title_for(da_f))
        return (mesh, q_overlay, title) if q_overlay is not None else (mesh, title)

    anim_obj = animation.FuncAnimation(fig, update, frames=len(t_indices),
                                       interval=interval, blit=False)

    # ---- output ------------------------------------------------------
    if out:
        writer = "pillow" if str(out).lower().endswith(".gif") else "ffmpeg"
        anim_obj.save(out, writer=writer, fps=fps, dpi=dpi)
        plt.close(fig)
        return out

    plt.close(fig)
    if not _HAS_IPYTHON:
        raise RuntimeError(
            "No IPython available for an inline widget — pass out='file.gif' "
            "to write the animation instead."
        )
    return HTML(anim_obj.to_jshtml())
