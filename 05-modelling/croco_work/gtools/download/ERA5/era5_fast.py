#!/usr/bin/env python3
"""
era5_fast.py -- fast ERA5 (CDS) download + conversion to CROCO online bulk forcing.

Adapted from era5_for_exercise.py (one CDS request per variable per month, all from
the SINGLE-LEVELS product) and made to write exactly what CROCO's online reader
(OCEAN/online_get_bulk.F, keys ONLINE + ERA_ECMWF) expects:

  * file names   <TAG>_Y<yyyy>M<mm>.nc  in <outputDir>/for_croco/, with the tags
                 LSM SST TP STRD SSR T2M Q U10M V10M msl   (msl lower case, as CROCO
                 looks it up; Q -- not Q2M)
  * variables    named like the tag, dims (time, lat, lon), missing_value 9999
  * time         DAYS since Yorig-01-01 (online_get_bulk.F ignores the units
                 attribute and assumes days)
  * latitude     increasing (ERA5 delivers it decreasing)
  * units        tp -> kg m-2 s-1 (x1000/3600), ssr/strd -> W m-2 (/3600), others native

Why it is faster than the legacy ERA5_request.py path
  The legacy path asks for specific humidity on the 1000 hPa PRESSURE-LEVELS product,
  a separate, much slower CDS queue. Here q is derived from 2 m dewpoint + mean sea
  level pressure (Bolton 1980), both single-level fields:
      es(Td) = 611.2 exp(17.67 Td / (Td + 243.5))   [Pa, Td in degC]
      q      = 0.622 es / (p - 0.378 es)            [kg/kg]
  SST and LSM are requested once per day (00:00): CROCO does not need them hourly.

Behaviour
  * current month: only the days already released (today - ERA5_LAG_DAYS, default 6)
  * raw files kept in <outputDir>/raw/ERA5_ecmwf_<VAR>_Y<yyyy>M<mm>.nc (same names as
    the legacy path, so earlier raw downloads are reused); .part while downloading
  * re-running is safe: raw and for_croco files that already cover the available
    days are skipped; a file from an earlier partial month is refreshed
  * ERA5_N_PARALLEL (default 2) requests in flight -- CDS rejects more queued jobs
    per user ("Number queued requests ... temporarily limited"); rejected or failed
    requests are retried with a growing wait

Usage
  python era5_fast.py --domain lonmin,lonmax,latmin,latmax --month_start 2026-06 \
      --month_end 2026-08 --outputDir .../downloaded_data/ERA5 --Yorig 1993
  (normally called by `ggosss26.py download_atmosphere_hindcast`)
"""
from __future__ import annotations

import argparse
import calendar
import concurrent.futures
import datetime as dt
import os
import threading
import time

import numpy as np
from netCDF4 import Dataset, num2date

DATASET = "reanalysis-era5-single-levels"
DOMAIN_MARGIN_DEG = 2.0
N_PARALLEL = int(os.environ.get("ERA5_N_PARALLEL", "2"))
LAG_DAYS = int(os.environ.get("ERA5_LAG_DAYS", "6"))
MAX_TRIES = 8

# raw key -> CDS variable name, hours requested
RAW = {
    "lsm":  ("land_sea_mask", "daily"),
    "sst":  ("sea_surface_temperature", "daily"),
    "tp":   ("total_precipitation", "hourly"),
    "strd": ("surface_thermal_radiation_downwards", "hourly"),
    "ssr":  ("surface_net_solar_radiation", "hourly"),
    "t2m":  ("2m_temperature", "hourly"),
    "d2m":  ("2m_dewpoint_temperature", "hourly"),
    "u10":  ("10m_u_component_of_wind", "hourly"),
    "v10":  ("10m_v_component_of_wind", "hourly"),
    "msl":  ("mean_sea_level_pressure", "hourly"),
}

# CROCO tag -> (raw inputs, units, long name, factor)   ('q' is derived)
OUT = {
    "LSM":  (["lsm"],  "(0-1)",      "land_sea_mask", 1.0),
    "SST":  (["sst"],  "K",          "sea_surface_temperature", 1.0),
    "TP":   (["tp"],   "kg m-2 s-1", "total_precipitation", 1000.0 / 3600.0),
    "STRD": (["strd"], "W m-2",      "surface_thermal_radiation_downwards", 1.0 / 3600.0),
    "SSR":  (["ssr"],  "W m-2",      "surface_net_solar_radiation", 1.0 / 3600.0),
    "T2M":  (["t2m"],  "K",          "2m_temperature", 1.0),
    "Q":    (["d2m", "msl"], "kg kg-1", "specific_humidity (2 m, from dewpoint and msl)", None),
    "U10M": (["u10"],  "m s-1",      "10m_u_component_of_wind", 1.0),
    "V10M": (["v10"],  "m s-1",      "10m_v_component_of_wind", 1.0),
    "msl":  (["msl"],  "Pa",         "mean_sea_level_pressure", 1.0),
}
_COORDS = {"longitude", "latitude", "time", "valid_time", "expver", "number"}

_lock = threading.Lock()


def log(msg):
    with _lock:
        print(dt.datetime.now().strftime("%H:%M:%S") + "  " + msg, flush=True)


def month_range(ms, me):
    y, m = map(int, ms.split("-"))
    ye, mend = map(int, me.split("-"))
    while (y, m) <= (ye, mend):
        yield y, m
        y, m = (y + 1, 1) if m == 12 else (y, m + 1)


def cds_area(domain, margin=DOMAIN_MARGIN_DEG):
    lon0, lon1, lat0, lat1 = (float(x) for x in domain)
    return [min(90.0, lat1 + margin), max(-180.0, lon0 - margin),
            max(-90.0, lat0 - margin), min(180.0, lon1 + margin)]


def available_days(y, m):
    last = dt.date.today() - dt.timedelta(days=LAG_DAYS)
    return [d for d in range(1, calendar.monthrange(y, m)[1] + 1) if dt.date(y, m, d) <= last]


def raw_path(raw_dir, key, y, m):
    return os.path.join(raw_dir, f"ERA5_ecmwf_{key.upper()}_Y{y}M{m:02d}.nc")


def out_path(out_dir, tag, y, m):
    return os.path.join(out_dir, f"{tag}_Y{y}M{m:02d}.nc")


def _last_time(path, yorig=None):
    """Last time record of a file as a datetime (None if unreadable)."""
    try:
        with Dataset(path) as nc:
            tname = "valid_time" if "valid_time" in nc.variables else "time"
            t = nc.variables[tname]
            if len(t) == 0:
                return None
            if yorig is not None:                       # for_croco files: days since Yorig
                return dt.datetime(yorig, 1, 1) + dt.timedelta(days=float(t[-1]))
            d = num2date(t[-1], t.units, getattr(t, "calendar", "standard"))
            return dt.datetime(d.year, d.month, d.day, d.hour)
    except Exception:
        return None


def covers(path, y, m, last_day, yorig=None):
    """True if the file exists and reaches the last available day of the month."""
    if not os.path.isfile(path) or os.path.getsize(path) == 0:
        return False
    t = _last_time(path, yorig)
    return t is not None and t.date() >= dt.date(y, m, last_day)


def request(key, y, m, days, area):
    var, freq = RAW[key]
    hours = ["00:00"] if freq == "daily" else [f"{h:02d}:00" for h in range(24)]
    return {
        "product_type": ["reanalysis"],
        "variable": [var],
        "year": [f"{y}"],
        "month": [f"{m:02d}"],
        "day": [f"{d:02d}" for d in days],
        "time": hours,
        "data_format": "netcdf",
        "download_format": "unarchived",
        "area": area,
    }


def retrieve(task):
    import cdsapi
    key, y, m, days, area, path = task
    name = os.path.basename(path)
    for attempt in range(1, MAX_TRIES + 1):
        try:
            log("request  " + name + ("" if attempt == 1 else f"  (attempt {attempt})"))
            t0 = time.time()
            c = cdsapi.Client(quiet=True, progress=False)
            c.retrieve(DATASET, request(key, y, m, days, area)).download(path + ".part")
            os.replace(path + ".part", path)
            log(f"done     {name}  ({time.time() - t0:.0f} s)")
            return
        except Exception as e:
            log("FAILED   " + name + ": " + (str(e).splitlines() or [repr(e)])[0][:200])
            if attempt == MAX_TRIES:
                raise
            time.sleep(min(60 * attempt, 300))


def read_raw(path):
    """(time datetimes, lon, lat, data[t,lat,lon] float64 with NaN for missing)."""
    with Dataset(path) as nc:
        tname = "valid_time" if "valid_time" in nc.variables else "time"
        t = nc.variables[tname]
        times = num2date(t[:], t.units, getattr(t, "calendar", "standard"))
        lon = np.array(nc.variables["longitude"][:], dtype="f8")
        lat = np.array(nc.variables["latitude"][:], dtype="f8")
        name = [k for k in nc.variables if k not in _COORDS][0]
        v = nc.variables[name]
        arr = v[:]
        if "expver" in v.dimensions:                     # ERA5 / ERA5T merge
            arr = arr.take(0, axis=v.dimensions.index("expver"))
        data = np.ma.filled(np.ma.asarray(arr).astype("f8"), np.nan)
    data = data.reshape(len(times), lat.size, lon.size)
    return times, lon, lat, data


def write_croco(path, tag, units, long_name, lon, lat, days, data, yorig):
    if lat[0] > lat[-1]:                                 # increasing latitude
        lat = lat[::-1]
        data = data[:, ::-1, :]
    data = np.where(np.isfinite(data), data, 9999.0)
    tmp = path + ".part"
    with Dataset(tmp, "w", format="NETCDF4") as nw:
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
        nw.history = "era5_fast.py (ERA5 single levels, CDS) " + dt.datetime.now().isoformat(timespec="seconds")
    os.replace(tmp, path)


def convert(tag, y, m, raw_dir, out_dir, yorig):
    inputs, units, long_name, factor = OUT[tag]
    times, lon, lat, a = read_raw(raw_path(raw_dir, inputs[0], y, m))
    if tag == "Q":
        _, _, _, p = read_raw(raw_path(raw_dir, "msl", y, m))
        td = a - 273.15
        es = 611.2 * np.exp(17.67 * td / (td + 243.5))
        data = 0.622 * es / (p - 0.378 * es)
    else:
        data = a * factor
    t0 = dt.datetime(yorig, 1, 1)
    days = np.array([(dt.datetime(d.year, d.month, d.day, d.hour, d.minute) - t0).total_seconds() / 86400.0
                     for d in times])
    write_croco(out_path(out_dir, tag, y, m), tag, units, long_name, lon, lat, days, data, yorig)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--domain", required=True, help="lonmin,lonmax,latmin,latmax (grid box; +2 deg added)")
    ap.add_argument("--month_start", required=True, help="YYYY-MM")
    ap.add_argument("--month_end", required=True, help="YYYY-MM")
    ap.add_argument("--outputDir", required=True, help="ERA5 base dir (raw/ and for_croco/ inside)")
    ap.add_argument("--Yorig", type=int, default=1993)
    a = ap.parse_args(argv)
    domain = [float(x) for x in a.domain.split(",")] if isinstance(a.domain, str) else a.domain
    run(domain, a.month_start, a.month_end, a.outputDir, a.Yorig)


def run(domain, month_start, month_end, output_dir, yorig=1993):
    raw_dir = os.path.join(output_dir, "raw")
    out_dir = os.path.join(output_dir, "for_croco")
    os.makedirs(raw_dir, exist_ok=True)
    os.makedirs(out_dir, exist_ok=True)
    area = cds_area(domain)
    print(f" ERA5 (single levels) area [N,W,S,E] {area}, months {month_start} .. {month_end}")

    months, tasks = [], []
    for y, m in month_range(month_start, month_end):
        days = available_days(y, m)
        if not days:
            print(f" {y}-{m:02d}: no ERA5 data released yet -- skipped")
            continue
        if len(days) < calendar.monthrange(y, m)[1]:
            print(f" {y}-{m:02d}: ERA5 only up to day {days[-1]} (release delay) -- partial month")
        todo = [t for t in OUT if not covers(out_path(out_dir, t, y, m), y, m, days[-1], yorig)]
        months.append((y, m, days, todo))
        need = sorted({k for t in todo for k in OUT[t][0]})
        for k in need:
            p = raw_path(raw_dir, k, y, m)
            if covers(p, y, m, days[-1]):
                continue
            tasks.append((k, y, m, days, area, p))

    print(f" {len(tasks)} CDS request(s) to do, {N_PARALLEL} in parallel")
    failed = []
    with concurrent.futures.ThreadPoolExecutor(max_workers=max(1, N_PARALLEL)) as pool:
        futs = {pool.submit(retrieve, t): os.path.basename(t[5]) for t in tasks}
        for f in concurrent.futures.as_completed(futs):
            if f.exception() is not None:
                failed.append(futs[f])
    if failed:
        raise SystemExit(" ERA5 requests failed: " + ", ".join(sorted(failed)) +
                         " -- re-run the same command to resume")

    for y, m, days, todo in months:
        for tag in todo:
            convert(tag, y, m, raw_dir, out_dir, yorig)
            print(f"  wrote {os.path.basename(out_path(out_dir, tag, y, m))}")
    n = len([f for f in os.listdir(out_dir) if f.endswith(".nc")])
    print(f"\n ERA5 files conversion done: {n} file(s) in {out_dir}\n")


if __name__ == "__main__":
    main()
