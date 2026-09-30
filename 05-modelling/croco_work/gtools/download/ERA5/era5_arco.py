#!/usr/bin/env python3
"""
era5_arco.py -- ERA5 from the Copernicus Data Store ARCO (Zarr) store, converted to
CROCO online bulk forcing. No CDS request queue: the data are read directly.
(ECMWF reference: ecmwf-training/dss-notebooks, reanalysis-era5-single-levels/arco-access)

Store   https://arco.datastores.ecmwf.int/cadl-arco-geo-002/arco/reanalysis_era5_single_levels/sfc/geoChunked.zarr
        "geo-chunked": each chunk = ~67 000 hours x 4 x 4 points, i.e. made for long
        series over a small area -- the right layout for a regional box over months
        (the time-chunked store holds one global map per chunk: 60x more transfer).
        Measured: t2m, GoG_12 box (+2 deg), May -> Sep 2026 hourly, 340 s, no queue.
Auth    the CDS personal token, read from CDSAPI_KEY or the `key:` line of ~/.cdsapirc
        (sent as an HTTP bearer header; never printed).
Needs   xarray, zarr>=3, fsspec, aiohttp>=3.9

Output  (identical to era5_fast.py / the legacy converter -- what OCEAN/online_get_bulk.F
        reads with ONLINE + ERA_ECMWF): <outputDir>/for_croco/<TAG>_Y<yyyy>M<mm>.nc,
        TAG in T2M Q TP SSR STRD U10M V10M msl (read by CROCO) + SST LSM (informative),
        time in DAYS since Yorig, latitude increasing, missing value 9999.

Fields not in the ARCO store, derived here
  Q    specific humidity from 2 m dewpoint + msl (Bolton 1980), as in era5_fast.py
  SSR  NET surface solar (CROCO uses the SSR file as the net shortwave flux, srflx).
       The store has ssrd (downward, direct + diffuse) and fdir (direct part), so
         SSR = fdir (1 - a_dir(mu)) + (ssrd - fdir) (1 - a_dif)
       with the open-water albedo of Taylor et al. (1996), used by the ECMWF model:
         a_dir = 0.037 / (1.1 mu^1.4 + 0.15),  a_dif = 0.06,
       mu = cosine of the solar zenith angle at the middle of the accumulation hour.
       Over the ocean this matches ERA5's own ssr to a few W m-2; near coasts it
       avoids the land albedo that ERA5's grid-cell ssr mixes in.
  LSM  1 where ERA5 SST is missing (land), 0 elsewhere (daily)
Accumulated fields (tp, strd, ssrd, fdir) are per preceding hour: /3600 -> W m-2,
tp x1000/3600 -> kg m-2 s-1.

Cache   the period extracted per variable is kept in raw/ARCO_<var>_<yyyymmdd>_<yyyymmdd>.nc,
        so a re-run (or another month split) does not download it again.
Env     ERA5_ARCO_URL (store), ERA5_ARCO_PARALLEL (variables at a time, default 3)
"""
from __future__ import annotations

import calendar
import concurrent.futures
import datetime as dt
import os
import time

import numpy as np
import pandas as pd
import xarray as xr

ARCO_URL = os.environ.get(
    "ERA5_ARCO_URL",
    "https://arco.datastores.ecmwf.int/cadl-arco-geo-002/arco/reanalysis_era5_single_levels/sfc/geoChunked.zarr")
N_PARALLEL = int(os.environ.get("ERA5_ARCO_PARALLEL", "3"))
DOMAIN_MARGIN_DEG = 2.0
ARCO_VARS = ["t2m", "d2m", "msl", "u10", "v10", "tp", "strd", "ssrd", "fdir", "sst"]
A_DIF = 0.06


def cds_key():
    key = os.environ.get("CDSAPI_KEY")
    rc = os.path.expanduser("~/.cdsapirc")
    if not key and os.path.exists(rc):
        for line in open(rc):
            if line.startswith("key:"):
                key = line.split(":", 1)[1].strip()
    if not key:
        raise SystemExit("ERA5 ARCO: no CDS token (set CDSAPI_KEY or ~/.cdsapirc 'key:')")
    return key


def open_store():
    import aiohttp
    # aiohttp's default 5-min total timeout per request is too short when the link is busy
    # (e.g. an ocean download in parallel): hundreds of chunk requests share it
    timeout = aiohttp.ClientTimeout(total=None, sock_connect=60, sock_read=600)
    return xr.open_zarr(ARCO_URL, consolidated=True,
                        storage_options={"headers": {"Authorization": f"Bearer {cds_key()}"},
                                         "client_kwargs": {"timeout": timeout}})


def log(msg):
    print(dt.datetime.now().strftime("%H:%M:%S") + "  " + msg, flush=True)


def cos_zenith(times, lat, lon):
    """cos(solar zenith) for times (DatetimeIndex, UTC) x lat x lon (NOAA formulas)."""
    t = pd.DatetimeIndex(times)
    doy = t.dayofyear.values.astype(float)
    hour = t.hour.values + t.minute.values / 60.0
    g = 2 * np.pi / 365.0 * (doy - 1 + (hour - 12) / 24)
    decl = (0.006918 - 0.399912 * np.cos(g) + 0.070257 * np.sin(g) - 0.006758 * np.cos(2 * g)
            + 0.000907 * np.sin(2 * g) - 0.002697 * np.cos(3 * g) + 0.00148 * np.sin(3 * g))
    eqt = 229.18 * (0.000075 + 0.001868 * np.cos(g) - 0.032077 * np.sin(g)
                    - 0.014615 * np.cos(2 * g) - 0.040849 * np.sin(2 * g))       # minutes
    tst = hour[:, None] * 60 + eqt[:, None] + 4 * lon[None, :]                  # true solar time, min
    ha = np.deg2rad(tst / 4 - 180)                                               # (t, lon)
    la = np.deg2rad(lat)[None, :, None]
    d = decl[:, None, None]
    mu = np.sin(la) * np.sin(d) + np.cos(la) * np.cos(d) * np.cos(ha[:, None, :])
    return np.clip(mu, 0.0, 1.0)


def net_solar(ssrd, fdir):
    """SSR (W m-2) from hourly-accumulated ssrd, fdir (J m-2) -- Taylor et al. 1996 albedo."""
    mid = pd.DatetimeIndex(ssrd.time.values) - pd.Timedelta(minutes=30)
    mu = cos_zenith(mid, ssrd.latitude.values, ssrd.longitude.values)
    a_dir = 0.037 / (1.1 * mu ** 1.4 + 0.15)
    sw = ssrd.values / 3600.0
    dr = np.minimum(fdir.values / 3600.0, sw)
    return dr * (1 - a_dir) + (sw - dr) * (1 - A_DIF)


def write_croco(path, tag, units, long_name, lon, lat, days, data, yorig):
    from netCDF4 import Dataset
    if lat[0] > lat[-1]:
        lat, data = lat[::-1], data[:, ::-1, :]
    data = np.where(np.isfinite(data), data, 9999.0)
    with Dataset(path + ".part", "w", format="NETCDF4") as nw:
        nw.createDimension("lon", lon.size)
        nw.createDimension("lat", lat.size)
        nw.createDimension("time", None)
        vlon = nw.createVariable("lon", "f4", ("lon",))
        vlat = nw.createVariable("lat", "f4", ("lat",))
        vt = nw.createVariable("time", "f8", ("time",))
        vd = nw.createVariable(tag, "f4", ("time", "lat", "lon"), zlib=True, complevel=4)
        vlon.long_name, vlon.units = "longitude of RHO-points", "degree_east"
        vlat.long_name, vlat.units = "latitude of RHO-points", "degree_north"
        vt.long_name, vt.units = "Time", f"days since {yorig}-1-1"
        vd.missing_value = 9999.0
        vd.units, vd.long_name = units, long_name
        vlon[:], vlat[:], vt[:], vd[:] = lon, lat, days, data
        nw.history = "era5_arco.py (ERA5 single levels, CDS ARCO store) " + dt.datetime.now().isoformat(timespec="seconds")
    os.replace(path + ".part", path)


def fetch(ds, var, t0, t1, box, raw_dir):
    """One variable over the whole period -> cached NetCDF; returns the DataArray."""
    lon0, lon1, lat0, lat1 = box
    cache = os.path.join(raw_dir, f"ARCO_{var}_{t0:%Y%m%d}_{t1:%Y%m%d}.nc")
    # any earlier extraction of this variable that covers the period (and the box) is reused
    import glob
    for c in sorted(glob.glob(os.path.join(raw_dir, f"ARCO_{var}_*.nc"))):
        try:
            a = xr.open_dataarray(c)
            ok = (pd.Timestamp(a.time.values[0]) <= pd.Timestamp(t0) and pd.Timestamp(a.time.values[-1]) >= pd.Timestamp(t1)
                  and float(a.longitude.min()) <= lon0 + 0.25 and float(a.longitude.max()) >= lon1 - 0.25
                  and float(a.latitude.min()) <= lat0 + 0.25 and float(a.latitude.max()) >= lat1 - 0.25)
            if ok:
                log(f"cached   {os.path.basename(c)}")
                return a.sel(time=slice(t0, t1), latitude=slice(lat0, lat1), longitude=slice(lon0, lon1)).load()
            a.close()
        except Exception:
            pass
    for attempt in range(1, 6):
        try:
            s = time.time()
            log(f"reading  {var} {t0:%Y-%m-%d} -> {t1:%Y-%m-%d %H:%M}" + ("" if attempt == 1 else f" (attempt {attempt})"))
            a = ds[var].sel(time=slice(t0, t1), latitude=slice(lat0, lat1), longitude=slice(lon0, lon1)).load()
            a.encoding = {}
            a.to_netcdf(cache + ".part")
            os.replace(cache + ".part", cache)
            log(f"done     {var}  ({time.time() - s:.0f} s)")
            return a
        except Exception as e:
            # repr(): some network errors (e.g. timeouts) carry an empty message
            log(f"FAILED   {var}: {repr(e)[:200]}")
            if attempt == 5:
                raise
            time.sleep(30 * attempt)


def run(domain, month_start, month_end, output_dir, yorig=1993):
    raw_dir = os.path.join(output_dir, "raw")
    out_dir = os.path.join(output_dir, "for_croco")
    os.makedirs(raw_dir, exist_ok=True)
    os.makedirs(out_dir, exist_ok=True)
    lon0, lon1, lat0, lat1 = (float(x) for x in domain)
    box = (lon0 - DOMAIN_MARGIN_DEG, lon1 + DOMAIN_MARGIN_DEG,
           lat0 - DOMAIN_MARGIN_DEG, lat1 + DOMAIN_MARGIN_DEG)

    ds = open_store()
    store_end = pd.Timestamp(ds.time.values[-1]).to_pydatetime()
    ys, ms = map(int, month_start.split("-"))
    ye, me = map(int, month_end.split("-"))
    t0 = dt.datetime(ys, ms, 1)
    t1 = min(dt.datetime(ye, me, calendar.monthrange(ye, me)[1], 23), store_end)
    # a full-day end only (CROCO needs whole days; a partial last day is dropped)
    if t1.hour != 23:
        t1 = dt.datetime(t1.year, t1.month, t1.day, 23) - dt.timedelta(days=1)
    print(f" ERA5 ARCO store, data up to {store_end:%Y-%m-%d %H:%M}; box {box}; "
          f"period {t0:%Y-%m-%d} -> {t1:%Y-%m-%d %H:%M}")
    if t1 < t0:
        raise SystemExit(" ERA5 ARCO: nothing available in the requested months")

    # months whose for_croco files are all present and reach the period end are skipped
    months = []
    y, m = ys, ms
    while (y, m) <= (ye, me) and dt.datetime(y, m, 1) <= t1:
        months.append((y, m))
        y, m = (y + 1, 1) if m == 12 else (y, m + 1)
    tags = ["T2M", "Q", "TP", "SSR", "STRD", "U10M", "V10M", "msl", "SST", "LSM"]

    def month_done(y, m):
        last = min(dt.datetime(y, m, calendar.monthrange(y, m)[1], 23), t1)
        need = (last - dt.datetime(yorig, 1, 1)).total_seconds() / 86400.0
        for tag in tags:
            p = os.path.join(out_dir, f"{tag}_Y{y}M{m:02d}.nc")
            if not os.path.isfile(p):
                return False
            try:
                with xr.open_dataset(p, decode_times=False) as f:
                    if float(f["time"].values[-1]) < need - (1.0 if tag in ("SST", "LSM") else 0.01):
                        return False
            except Exception:
                return False
        return True

    todo = [(y, m) for (y, m) in months if not month_done(y, m)]
    for (y, m) in months:
        if (y, m) not in todo:
            print(f"  {y}-{m:02d}: all files present -- skipped")
    if not todo:
        return
    p0 = dt.datetime(todo[0][0], todo[0][1], 1)
    p1 = min(dt.datetime(todo[-1][0], todo[-1][1], calendar.monthrange(*todo[-1])[1], 23), t1)

    data = {}
    with concurrent.futures.ThreadPoolExecutor(max_workers=max(1, N_PARALLEL)) as pool:
        futs = {pool.submit(fetch, ds, v, p0, p1, box, raw_dir): v for v in ARCO_VARS}
        for f in concurrent.futures.as_completed(futs):
            data[futs[f]] = f.result()

    lon = data["t2m"].longitude.values.astype("f8")
    lat = data["t2m"].latitude.values.astype("f8")
    ref = dt.datetime(yorig, 1, 1)

    # derived hourly fields over the whole period
    td = data["d2m"].values - 273.15
    es = 611.2 * np.exp(17.67 * td / (td + 243.5))
    q = 0.622 * es / (data["msl"].values - 0.378 * es)
    ssr = net_solar(data["ssrd"], data["fdir"])
    hourly = {
        "T2M":  (data["t2m"].values, "K", "2m_temperature"),
        "Q":    (q, "kg kg-1", "specific_humidity (2 m, from dewpoint and msl)"),
        "TP":   (data["tp"].values * 1000.0 / 3600.0, "kg m-2 s-1", "total_precipitation"),
        "SSR":  (ssr, "W m-2", "surface_net_solar_radiation (ssrd, fdir, Taylor et al. 1996 ocean albedo)"),
        "STRD": (data["strd"].values / 3600.0, "W m-2", "surface_thermal_radiation_downwards"),
        "U10M": (data["u10"].values, "m s-1", "10m_u_component_of_wind"),
        "V10M": (data["v10"].values, "m s-1", "10m_v_component_of_wind"),
        "msl":  (data["msl"].values, "Pa", "mean_sea_level_pressure"),
    }
    times = pd.DatetimeIndex(data["t2m"].time.values)
    sst = data["sst"].values

    for (y, m) in todo:
        sel = (times.year == y) & (times.month == m)
        days = np.array([(t.to_pydatetime() - ref).total_seconds() / 86400.0 for t in times[sel]])
        for tag, (arr, units, ln) in hourly.items():
            write_croco(os.path.join(out_dir, f"{tag}_Y{y}M{m:02d}.nc"), tag, units, ln,
                        lon, lat, days, arr[sel], yorig)
        d00 = sel & (times.hour == 0)                        # SST, LSM: daily at 00:00
        ddays = np.array([(t.to_pydatetime() - ref).total_seconds() / 86400.0 for t in times[d00]])
        write_croco(os.path.join(out_dir, f"SST_Y{y}M{m:02d}.nc"), "SST", "K",
                    "sea_surface_temperature", lon, lat, ddays, sst[d00], yorig)
        write_croco(os.path.join(out_dir, f"LSM_Y{y}M{m:02d}.nc"), "LSM", "(0-1)",
                    "land_sea_mask (1 where ERA5 SST is missing)", lon, lat, ddays,
                    np.isnan(sst[d00]).astype("f4"), yorig)
        print(f"  {y}-{m:02d}: wrote {len(tags)} files, {sel.sum()} hourly records "
              f"(last {times[sel][-1]:%Y-%m-%d %H:%M})")
    n = len([f for f in os.listdir(out_dir) if f.endswith(".nc")])
    print(f"\n ERA5 files conversion done: {n} file(s) in {out_dir}\n")
