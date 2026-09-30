"""gtools/nesting.py — CROCO offline nesting helper.

Turn a CROCO parent output file (e.g. GoG_12, 1/12 deg) into a Mercator/GLORYS-
format ocean file so the standard make_ini / make_bry machinery can build the
ini/bry of a finer child run inside it (e.g. IGOG_36, 1/36 deg).
In the hindcast chain this is automatic: set PARENT=<parent config> in the
child's domain.cfg and step 03 calls parent_to_monthly() (see hindcast/README.md).

The child's crocotools_param.py keeps inputdata='mercator'; this converter makes
the parent CROCO output *look* like a Mercator file:
  temp -> thetao, salt -> so, u -> uo, v -> vo, zeta -> zos
  dims (time, depth, latitude, longitude), depth positive-down (z-levels).

The sigma->z vertical interpolation reuses gtools.postprocess (depths,
field_at_depth, rotate_uv), so the parent's terrain-following levels are mapped
onto standard z-levels the child's make_bry then re-interpolates to child sigma.

Usage
-----
    import gtools.nesting as nest
    nest.croco_to_mercator(
        "hindcast/model-runs/GoG_12/20260601_plain_tides_rivers/hcast/CROCO_FILES/croco_avg.nc",
        "hindcast/scratch/IGOG_36/downloaded_data/PARENT_GoG_12/parent_20260601.nc",
        Yorig=1993)

Then point the child's make_ini / make_bry at the output file.
"""

import numpy as np
import xarray as xr

try:
    import gtools.postprocess as pp
except Exception:
    import postprocess as pp


# Standard GLORYS/Mercator 50-level depth axis (m, positive down).
# Override with depths=... if your product uses a different axis.
MERCATOR_DEPTHS_50 = np.array([
    0.494025, 1.541375, 2.645669, 3.819495, 5.078224, 6.440614, 7.929560,
    9.572997, 11.404999, 13.467140, 15.810070, 18.495559, 21.598820,
    25.211411, 29.444731, 34.434151, 40.344051, 47.373692, 55.764290,
    65.807266, 77.853851, 92.326073, 109.729301, 130.666000, 155.850693,
    186.125595, 222.475204, 266.040314, 318.127411, 380.213013, 453.937714,
    541.088928, 643.566772, 763.333130, 902.339294, 1062.439941, 1245.291016,
    1452.250977, 1684.284058, 1941.892944, 2225.077881, 2533.336426,
    2865.702637, 3220.820801, 3597.031982, 3992.483887, 4405.224121,
    4833.291016, 5274.784180, 5727.917000])


def croco_to_mercator(croco_his, fname_out, Yorig=None,
                      depths=None, tracers=("temp", "salt")):
    """Convert a CROCO history file to a Mercator-format ocean file.

    croco_his : path to the parent CROCO croco_his.nc
    fname_out : output path (a Mercator-like .nc)
    Yorig     : CROCO time origin year (hindcast 1993, forecast 2000)
    depths    : target z-levels (m, positive down); default = Mercator 50-level
    """
    if depths is None:
        depths = MERCATOR_DEPTHS_50
    depths = np.asarray(depths, dtype=float)

    ds = pp.open_history(croco_his, Yorig=Yorig)
    lon = ds["lon_rho"].values          # (eta, xi)
    lat = ds["lat_rho"].values
    time = pp.times(ds)
    nt = len(np.atleast_1d(time))
    nz = depths.size
    eta, xi = lon.shape

    # Mercator uses 1D longitude/latitude axes. CROCO rho grid is (near-)regular
    # for this box, so take the row/column means as the 1D axes.
    lon1d = lon.mean(axis=0)            # (xi,)
    lat1d = lat.mean(axis=1)            # (eta,)

    # allocate outputs (time, depth, lat, lon)
    thetao = np.full((nt, nz, eta, xi), np.nan, dtype="float32")
    so     = np.full((nt, nz, eta, xi), np.nan, dtype="float32")
    uo     = np.full((nt, nz, eta, xi), np.nan, dtype="float32")
    vo     = np.full((nt, nz, eta, xi), np.nan, dtype="float32")
    zos    = np.full((nt, eta, xi), np.nan, dtype="float32")

    maskr = ds["mask_rho"].values.astype(bool)   # True = water

    for k in range(nt):
        # SSH (2D)
        zos[k] = ds["zeta"].isel(time=k).values

        # depth of the TOP sigma cell at this time (negative down), per column
        ztop = pp.depths(ds, tindex=k)[-1]        # (eta, xi), e.g. ~ -0.5 m
        # surface tracers/currents (fill for z-levels above the top cell)
        t_surf = pp.surface(ds, "temp", tindex=k).values
        s_surf = pp.surface(ds, "salt", tindex=k).values
        ur, vr = pp.surface_uv(ds, tindex=k)
        ue_surf, vn_surf = pp.rotate_uv(ds, ur, vr)

        for zi, dm in enumerate(depths):
            th = pp.field_at_depth(ds, "temp", dm, tindex=k).values
            sa = pp.field_at_depth(ds, "salt", dm, tindex=k).values
            u_e, v_n = pp.uv_at_depth(ds, dm, tindex=k, rotate=True)
            u_e, v_n = np.squeeze(u_e), np.squeeze(v_n)
            # a z-level ABOVE the top sigma cell (shallower than |ztop|) and over
            # water -> fill with the surface value (NOT deep sea-floor NaNs)
            above_top = (-abs(dm) > ztop) & maskr
            th[above_top] = t_surf[above_top]
            sa[above_top] = s_surf[above_top]
            u_e[above_top] = ue_surf[above_top]
            v_n[above_top] = vn_surf[above_top]
            thetao[k, zi] = th; so[k, zi] = sa
            uo[k, zi] = u_e; vo[k, zi] = v_n

    out = xr.Dataset(
        {
            "thetao": (("time", "depth", "latitude", "longitude"), thetao),
            "so":     (("time", "depth", "latitude", "longitude"), so),
            "uo":     (("time", "depth", "latitude", "longitude"), uo),
            "vo":     (("time", "depth", "latitude", "longitude"), vo),
            "zos":    (("time", "latitude", "longitude"), zos),
        },
        coords={
            "time": np.atleast_1d(time),
            "depth": depths,
            "latitude": lat1d,
            "longitude": lon1d,
        },
    )
    # CF-ish attributes the reader may look for
    out["depth"].attrs   = {"units": "m", "positive": "down", "long_name": "depth"}
    out["latitude"].attrs  = {"units": "degrees_north"}
    out["longitude"].attrs = {"units": "degrees_east"}
    out["thetao"].attrs = {"units": "degC", "long_name": "sea_water_potential_temperature"}
    out["so"].attrs     = {"units": "psu",  "long_name": "sea_water_salinity"}
    out["uo"].attrs     = {"units": "m s-1","long_name": "eastward_sea_water_velocity"}
    out["vo"].attrs     = {"units": "m s-1","long_name": "northward_sea_water_velocity"}
    out["zos"].attrs    = {"units": "m",    "long_name": "sea_surface_height"}
    out.attrs["title"] = "CROCO parent regridded to Mercator format for offline nesting"
    out.attrs["source"] = croco_his

    enc = {v: {"dtype": "float32", "_FillValue": 9.96921e36} for v in
           ("thetao", "so", "uo", "vo", "zos")}
    out.to_netcdf(fname_out, encoding=enc)
    ds.close()
    print(f"Created Mercator-format parent: {fname_out}")
    print(f"  {nt} time(s), {nz} depths, grid {eta} x {xi}")
    return fname_out


def parent_to_monthly(parent_files, out_dir, Yorig=1993):
    """Offline nesting input for the hindcast chain: convert every parent output
    file (e.g. the spin-up and hindcast croco_avg.nc of the parent's cycles) with
    croco_to_mercator, merge them in time, and write monthly files YYYY_MM.nc --
    the SAME layout as the GLORYS folder. make_ini_hindcast / make_bry_hindcast and
    the cycle driver then build the child's ini/bry from the parent without any
    change (the child's crocotools_param.py keeps inputdata='mercator').

    Daily AVERAGES are the right parent output for a child that runs its own tides:
    averaging removes the parent's tide (hourly snapshots taken once a day would
    alias it), and the child adds the tide itself from TPXO at its boundaries.
    """
    import os
    import pandas as pd
    os.makedirs(out_dir, exist_ok=True)
    parts = []
    for i, f in enumerate(parent_files):
        tmp = os.path.join(out_dir, f"_parent_part{i:02d}.nc")
        croco_to_mercator(f, tmp, Yorig=Yorig)
        parts.append(tmp)
    ds = xr.open_mfdataset(parts, combine="nested", concat_dim="time").sortby("time").load()
    _, keep = np.unique(ds["time"].values, return_index=True)   # overlapping records: keep one
    ds = ds.isel(time=keep)
    t = pd.DatetimeIndex(ds["time"].values)
    written = []
    for ym in sorted(set(zip(t.year, t.month))):
        sub = ds.isel(time=np.where((t.year == ym[0]) & (t.month == ym[1]))[0])
        f = os.path.join(out_dir, f"{ym[0]}_{ym[1]:02d}.nc")
        sub.attrs["title"] = "CROCO parent (daily means) in Mercator/GLORYS format for offline nesting"
        sub.to_netcdf(f + ".tmp")
        os.replace(f + ".tmp", f)
        written.append(f)
        print(f"  {os.path.basename(f)}: {sub.sizes['time']} record(s) "
              f"{str(sub.time.values[0])[:16]} .. {str(sub.time.values[-1])[:16]}")
    ds.close()
    for p in parts:
        os.remove(p)
    return written
