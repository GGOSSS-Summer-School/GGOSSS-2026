"""
Scripts for downloading data from the Copernicus Marine Service (CMEMS)
http://marine.copernicus.eu/services-portfolio/access-to-products/

Authentication note:
    Credentials are NOT hard-coded here. copernicusmarine stores a login in
    ~/.copernicusmarine/ after a one-time `copernicusmarine login`. The helper
    ensure_cmems_login() checks for that login and, if missing, prompts the user
    (or accepts optional usrname/passwd for automated runs). After login, the
    subset commands need no credentials on the command line.
"""
from math import floor, ceil
import numpy as np
from datetime import datetime, timedelta, date
import calendar
import sys, os
from pathlib import Path
import xarray as xr
import subprocess
import time
import threading
import getpass


def is_valid_netcdf_file(file_path):
    # short-circuit on a missing file (avoids HDF5 error spam for files not yet downloaded)
    if not os.path.exists(file_path):
        return False
    try:
        with xr.open_dataset(file_path) as ds:
            return True
    except Exception:
        return False


def netcdf_covers_time_range(file_path, start_date, end_date, tolerance_days=1,
                             max_step_hours=None):
    """True if an existing NetCDF file really covers [start_date, end_date] at the
    expected time step -- not just exists. A previous interrupted run can leave a
    valid file covering one day only; an existence check would skip it forever.
    max_step_hours (e.g. 30 for daily means) also rejects coarser (e.g. monthly) files.
    (from validation.py)"""
    if not is_valid_netcdf_file(file_path):
        return False
    try:
        import pandas as pd
        with xr.open_dataset(file_path) as ds:
            if "time" not in ds:
                return False
            times = pd.to_datetime(np.atleast_1d(ds["time"].values))
            if len(times) == 0:
                return False
            tol = pd.Timedelta(days=tolerance_days)
            if not (times.min() <= pd.Timestamp(start_date) + tol and
                    times.max() >= pd.Timestamp(end_date) - tol):
                return False
            if max_step_hours is not None and len(times) > 1:
                step_h = np.median(np.diff(np.sort(times))).astype("timedelta64[m]") / np.timedelta64(60, "m")
                if step_h > max_step_hours:
                    return False
            return True
    except Exception:
        return False


def ensure_cmems_login(usrname=None, passwd=None):
    """
    Make sure copernicusmarine is authenticated (Python API, no subprocess).
    Credentials are stored by copernicusmarine in ~/.copernicusmarine/, NOT here.
    (from validation.py)
    """
    import copernicusmarine
    try:
        if copernicusmarine.login(check_credentials_valid=True):
            print("CMEMS: already logged in.")
            return
    except Exception:
        pass
    print("CMEMS: no valid login found. Please enter your Copernicus Marine credentials.")
    print("  (Register free at https://data.marine.copernicus.eu/register)")
    if not usrname:
        usrname = input("  CMEMS username: ").strip()
    if not passwd:
        passwd = getpass.getpass("  CMEMS password: ")
    try:
        copernicusmarine.login(username=usrname, password=passwd)
    except Exception as e:
        raise RuntimeError(f"CMEMS login failed. Check your username/password.\n{e}")
    print("CMEMS: login successful (credentials stored for future runs).")


def download_cmems(dataset, varlist, start_date, end_date, domain, depths, outputDir, fname, ver='',
                   expected_max_step_hours=30):
    """
    Download a subset of a CMEMS dataset with the copernicusmarine Python API
    (adapted from validation.py: in-process subset(), coverage-aware skip, retries).
    Assumes ensure_cmems_login() has already run.

    expected_max_step_hours : an existing file covering the dates but at a coarser
        step (e.g. a monthly-mean file where daily is expected) is re-downloaded.
    CMEMS_SERVICE (environment) : force an ARCO store, e.g. arco-geo-series or
        arco-time-series; empty/unset = the toolbox's own choice (default).
    Raises RuntimeError after the last failed attempt (so a pipeline stops).
    """
    import copernicusmarine
    f = os.path.normpath(os.path.join(outputDir, fname))
    if netcdf_covers_time_range(f, start_date, end_date, max_step_hours=expected_max_step_hours):
        print("file already exists and covers the requested range - " + fname)
        return
    elif is_valid_netcdf_file(f):
        print(f"file exists but doesn't cover {start_date.date()} -> {end_date.date()} "
              f"(or is coarser than expected) - re-downloading {fname}")
        os.remove(f)

    subset_kwargs = dict(
        dataset_id=dataset,
        variables=list(varlist),
        minimum_longitude=domain[0], maximum_longitude=domain[1],
        minimum_latitude=domain[2], maximum_latitude=domain[3],
        minimum_depth=depths[0], maximum_depth=depths[1],
        start_datetime=start_date.strftime("%Y-%m-%d 00:00:00"),
        end_datetime=end_date.strftime("%Y-%m-%d 23:59:59"),
        output_directory=os.path.normpath(outputDir),
        output_filename=fname,
        overwrite=True,
        disable_progress_bar=True,
    )
    if ver:
        subset_kwargs["dataset_version"] = ver
    service = os.environ.get("CMEMS_SERVICE", "")
    if service:
        subset_kwargs["service"] = service

    MAX_RETRIES = int(os.environ.get("CMEMS_MAX_RETRIES", "5"))
    RETRY_WAIT = 30
    last = None
    for i in range(MAX_RETRIES):
        t0 = time.time()
        print(f"  {fname}: attempt {i+1} of {MAX_RETRIES} ({dataset}, {start_date.date()} -> {end_date.date()})", flush=True)
        try:
            copernicusmarine.subset(**subset_kwargs)
            if is_valid_netcdf_file(f):
                print(f"  Completed {fname} in {time.time() - t0:.0f} s", flush=True)
                return
            if os.path.exists(f):
                os.unlink(f)
            raise RuntimeError(f"bad NetCDF output: {fname}")
        except Exception as e:
            last = e
            if i < MAX_RETRIES - 1:
                print(f"  Error: {str(e)[:200]}, retrying in {RETRY_WAIT * (i + 1)} s...", flush=True)
                time.sleep(RETRY_WAIT * (i + 1))
    raise RuntimeError(f"CMEMS download failed after {MAX_RETRIES} attempts: {fname}: {last}")


def download_mercator_ops(domain, run_date, hdays, fdays, outputDir, usrname=None, passwd=None):
    """
    Download the operational Mercator ocean output (anfc).

    Credentials are optional and only used to log in the first time:
      - interactive user: omit them -> prompted once by ensure_cmems_login()
      - automated run:     pass usrname/passwd -> no prompt
    """
    # make sure we're authenticated before any downloads (prompts only if needed)
    ensure_cmems_login(usrname, passwd)

    # extend the download range by a day either side so we cover the model run time
    hdays = hdays + 1
    fdays = fdays + 1
    start_date = run_date + timedelta(days=-hdays)
    end_date = run_date + timedelta(days=fdays)

    # hard coding the variable info to extract (variables are in separate files)
    VARIABLES = [
        {
            "name": "so",
            "#": "Salinity in psu",
            "id": "cmems_mod_glo_phy-so_anfc_0.083deg_P1D-m",
            "vars": ["so"],
            "fname": f"mercator_so_{run_date.strftime('%Y%m%d_%H')}.nc"
        },
        {
            "name": "thetao",
            "#": "Temperature in degrees C",
            "id": "cmems_mod_glo_phy-thetao_anfc_0.083deg_P1D-m",
            "vars": ["thetao"],
            "fname": f"mercator_thetao_{run_date.strftime('%Y%m%d_%H')}.nc"
        },
        {
            "name": "zos",
            "#": "SSH in m",
            "id": "cmems_mod_glo_phy_anfc_0.083deg_P1D-m",
            "vars": ["zos"],
            "fname": f"mercator_zos_{run_date.strftime('%Y%m%d_%H')}.nc"
        },
        {
            "name": "uo_vo",
            "#": "uo:Eastward velocity in m/s | vo:Northward velocity in m/s",
            "id": "cmems_mod_glo_phy-cur_anfc_0.083deg_P1D-m",
            "vars": ["uo", "vo"],
            "fname": f"mercator_uo_vo_{run_date.strftime('%Y%m%d_%H')}.nc"
        },
    ]

    # all depths
    depths = [0.493, 5727.918]

    # download the variable files in parallel to save time
    def download_worker(var):
        download_cmems(var["id"], var["vars"], start_date, end_date,
                       domain, depths, outputDir, var["fname"])

    threads = []
    for var in VARIABLES:
        t = threading.Thread(target=download_worker, args=(var,))
        threads.append(t)
        t.start()
    for t in threads:
        t.join()

    # Concatenate the separate NetCDF files
    print("merge NetCDF files")
    output_path = os.path.abspath(os.path.join(outputDir, f"MERCATOR_{run_date.strftime('%Y%m%d_%H')}.nc"))
    datasets = [xr.open_dataset(os.path.join(outputDir, var["fname"])) for var in VARIABLES]
    merged = xr.merge(datasets)
    merged.to_netcdf(output_path, mode="w")
    for ds in datasets:
        ds.close()

    subprocess.call(["chmod", "-R", "775", output_path])


MERCATOR_ANALYSIS = [   # the Mercator global analysis (anfc) is split over 4 datasets
    ("cmems_mod_glo_phy-thetao_anfc_0.083deg_P1D-m", ["thetao"]),
    ("cmems_mod_glo_phy-so_anfc_0.083deg_P1D-m", ["so"]),
    ("cmems_mod_glo_phy-cur_anfc_0.083deg_P1D-m", ["uo", "vo"]),
    ("cmems_mod_glo_phy_anfc_0.083deg_P1D-m", ["zos"]),
]


def download_mercator_monthly(domain, start_date, end_date, depths, outputDir,
                              usrname=None, passwd=None):
    """
    Monthly files YYYY_MM.nc from the Mercator global ANALYSIS (anfc, 1/12 deg,
    daily means) -- the alternative to GLORYS for periods the reanalysis does not
    cover yet (GLORYS lags the present by a few months). Same variable names as
    GLORYS (thetao, so, uo, vo, zos), so the 'mercator' reader of make_ini/make_bry
    uses them unchanged.

    Speed: each of the 4 source datasets is downloaded ONCE for the whole missing
    period (not month by month) from the time-chunked ARCO store, then split into
    monthly files. That store's chunks are 134 days x 2 levels x 32x64 points, so
    a long request wastes little; the map store (1 day x 1 level x 512x2048 points)
    moves 17-60x more data than a regional box needs (measured: ~30 kB/s useful out
    of ~900 kB/s received). CMEMS_SERVICE overrides the store.
    Existing valid monthly files are kept; the whole-period parts (*.period.nc) are
    kept until every month is written, so an interrupted run resumes.
    """
    ensure_cmems_login(usrname, passwd)
    os.makedirs(outputDir, exist_ok=True)

    months = []
    d = datetime(start_date.year, start_date.month, 1)
    while d <= end_date:
        months.append(d)
        d = datetime(d.year + (d.month == 12), d.month % 12 + 1, 1)
    todo = [m for m in months
            if not is_valid_netcdf_file(os.path.join(outputDir, m.strftime('%Y_%m') + '.nc'))]
    for m in months:
        if m not in todo:
            print("file already exists - " + m.strftime('%Y_%m') + ".nc")
    if not todo:
        return
    p0 = todo[0]
    p1 = datetime(todo[-1].year, todo[-1].month,
                  calendar.monthrange(todo[-1].year, todo[-1].month)[1])
    tag = p0.strftime('%Y%m') + '_' + p1.strftime('%Y%m')
    print(f"CMEMS: Mercator analysis {p0.date()} -> {p1.date()} ({len(todo)} month(s)), "
          f"one request per dataset")

    # the whole-period request is always done through the time-chunked store unless
    # the user forces another one
    os.environ.setdefault('CMEMS_SERVICE', 'arco-time-series')
    parts = []
    for ds, vars_ in MERCATOR_ANALYSIS:
        part = f"mercator_{'_'.join(vars_)}_{tag}.period.nc"
        # one request per dataset; the 4 run one after the other (validation.py note:
        # parallel in-process subsets can crash netCDF4/HDF5)
        download_cmems(ds, vars_, p0, p1, domain, depths, outputDir, part)
        parts.append(os.path.join(outputDir, part))

    dss = [xr.open_dataset(p) for p in parts]
    merged = xr.merge(dss, compat="override", join="outer")
    for m0 in todo:
        m1 = datetime(m0.year, m0.month, calendar.monthrange(m0.year, m0.month)[1], 23, 59)
        sub = merged.sel(time=slice(m0, m1))
        if sub.sizes.get('time', 0) == 0:
            print(f"CMEMS: no Mercator records for {m0.strftime('%Y_%m')} -- skipped")
            continue
        sub.attrs["source"] = "Mercator global analysis (CMEMS anfc 1/12 deg), merged by gtools"
        final = os.path.join(outputDir, m0.strftime('%Y_%m') + '.nc')
        sub.to_netcdf(final + '.tmp')
        os.replace(final + '.tmp', final)
        print(f"Completed {m0.strftime('%Y_%m')}.nc (Mercator analysis, {sub.sizes['time']} days)")
    for x in dss:
        x.close()
    for p in parts:
        os.remove(p)


def download_cmems_monthly(dataset, domain, start_date, end_date, varlist, depths,
                           outputDir, usrname=None, passwd=None):
    """
    Monthly files YYYY_MM.nc for any daily CMEMS dataset (GLORYS by default).
    Credentials optional (see download_mercator_ops docstring).

    Speed: the missing months are fetched in ONE request for the whole period from
    the time-chunked ARCO store (chunks span months over a small area, so little is
    wasted), then split into monthly files -- see download_mercator_monthly for the
    measurements. CMEMS_SERVICE overrides the store. A month at the end of the
    dataset (e.g. GLORYS, which lags the present) holds the days that exist.
    """
    ensure_cmems_login(usrname, passwd)
    os.makedirs(outputDir, exist_ok=True)

    months = []
    d = datetime(start_date.year, start_date.month, 1)
    while d <= end_date:
        months.append(d)
        d = datetime(d.year + (d.month == 12), d.month % 12 + 1, 1)
    todo = [m for m in months
            if not is_valid_netcdf_file(os.path.join(outputDir, m.strftime('%Y_%m') + '.nc'))]
    for m in months:
        if m not in todo:
            print("file already exists - " + m.strftime('%Y_%m') + ".nc")
    if not todo:
        return
    p0 = todo[0]
    p1 = datetime(todo[-1].year, todo[-1].month,
                  calendar.monthrange(todo[-1].year, todo[-1].month)[1])
    print(f"CMEMS: {dataset} {p0.date()} -> {p1.date()} ({len(todo)} month(s)), one request")

    os.environ.setdefault('CMEMS_SERVICE', 'arco-time-series')
    ds, pieces = _download_pieces(dataset, varlist, p0, p1, domain, depths, outputDir,
                                  f"{p0:%Y%m}_{p1:%Y%m}")
    last = ds.time.values[-1]
    for m0 in todo:
        m1 = datetime(m0.year, m0.month, calendar.monthrange(m0.year, m0.month)[1], 23, 59)
        sub = ds.sel(time=slice(m0, m1))
        n = sub.sizes.get('time', 0)
        if n == 0:
            print(f"CMEMS: no {dataset} records for {m0:%Y_%m} (dataset ends {str(last)[:10]}) -- skipped")
            continue
        final = os.path.join(outputDir, m0.strftime('%Y_%m') + '.nc')
        sub.to_netcdf(final + '.tmp')
        os.replace(final + '.tmp', final)
        print(f"Completed {m0:%Y_%m}.nc ({n} days"
              + (f", dataset ends {str(last)[:10]})" if n < calendar.monthrange(m0.year, m0.month)[1] else ")"))
    ds.close()
    for p in pieces:
        os.remove(p)


def _piece_worker(args):
    download_cmems(*args)


def _download_with_watchdog(dataset, varlist, t0, t1, domain, depths, outputDir, fname):
    """
    download_cmems in a child process, killed and restarted when it STALLS.

    Why: after a short network drop the toolbox's HTTP reads can hang forever (no read
    timeout), so neither its own retries nor download_cmems' ever start. The watchdog
    watches the bytes the download process reads (network included) and restarts the
    piece after CMEMS_STALL_MIN (default 10) minutes without any. Up to
    CMEMS_MAX_RETRIES (default 5) restarts.
    """
    import glob
    import multiprocessing as mp
    stall = 60 * float(os.environ.get("CMEMS_STALL_MIN", "10"))
    tries = int(os.environ.get("CMEMS_MAX_RETRIES", "5"))
    final = os.path.join(outputDir, fname)
    ctx = mp.get_context("fork")
    for attempt in range(1, tries + 1):
        p = ctx.Process(target=_piece_worker,
                        args=((dataset, varlist, t0, t1, domain, depths, outputDir, fname),))
        p.start()
        last_change, last_state = time.time(), None
        while p.is_alive():
            p.join(20)
            # progress = bytes the download process has READ (network included: /proc/<pid>/io
            # rchar). The output file is no good signal: the toolbox writes it in bursts,
            # sometimes only after tens of minutes of receiving data.
            try:
                with open(f"/proc/{p.pid}/io") as io:
                    state = int([l for l in io if l.startswith("rchar")][0].split()[1])
            except Exception:
                state = tuple((g, os.path.getsize(g), os.path.getmtime(g))
                              for g in glob.glob(final + ".*"))
            if state != last_state:
                last_state, last_change = state, time.time()
            elif time.time() - last_change > stall:
                print(f"  {fname}: no progress for {stall / 60:.0f} min (stalled connection) "
                      f"-- restarting it ({attempt}/{tries})", flush=True)
                p.terminate(); p.join(30)
                if p.is_alive():
                    p.kill(); p.join()
                for g in glob.glob(final + ".*"):
                    os.remove(g)
                break
        if p.exitcode == 0 and is_valid_netcdf_file(final):
            return
        if p.exitcode not in (0, None, -15, -9):
            print(f"  {fname}: download process failed (exit {p.exitcode}), retrying", flush=True)
        time.sleep(30 * attempt)
    raise RuntimeError(f"CMEMS piece {fname} failed after {tries} attempts (stalls or errors)")


# depth bands (m) of the piecewise download: each piece is one variable over one band
DEPTH_BANDS = [(0, 50), (50, 250), (250, 1100), (1100, 3400), (3400, 6000)]


def _download_pieces(dataset, varlist, p0, p1, domain, depths, outputDir, tag):
    """
    Download [p0, p1] as small pieces -- one variable x one depth band each -- and
    return (merged dataset, piece files).

    Why pieces: in the time-chunked ARCO store every chunk holds ONE depth level of ONE
    variable, so splitting by variable and depth band transfers exactly the same data
    as one big request. But a big request (an hour for GLORYS over a month or two) is
    restarted from zero by any network glitch, while a piece is a few minutes: a
    failure repeats only that piece, and a re-run skips the pieces already on disk.
    """
    zmin, zmax = float(depths[0]), float(depths[1])
    jobs = []
    for v in varlist:
        if v in ('zos',):                                   # 2-D: one piece
            jobs.append((v, (zmin, zmax), f"_piece_{tag}_{v}.nc"))
            continue
        for i, (a, b) in enumerate(DEPTH_BANDS):
            lo, hi = max(a, zmin), min(b, zmax)
            if lo < hi:
                jobs.append((v, (lo, hi), f"_piece_{tag}_{v}_z{i}.nc"))
    print(f"CMEMS: {len(jobs)} piece(s) (variable x depth band)")
    pieces = []
    for k, (v, (lo, hi), fname) in enumerate(jobs, 1):
        f = os.path.join(outputDir, fname)
        if is_valid_netcdf_file(f):
            print(f"  [{k}/{len(jobs)}] {fname}: already downloaded")
        else:
            print(f"  [{k}/{len(jobs)}] {v} {lo:g}-{hi:g} m", flush=True)
            # (download_cmems' own "already covers the dates" check is not used here: a
            #  dataset that ends inside the period, like GLORYS, never covers p1)
            _download_with_watchdog(dataset, [v], p0, p1, domain, [lo, hi], outputDir, fname)
        pieces.append(f)
    by_var = {}
    for (v, _, fname), f in zip(jobs, pieces):
        by_var.setdefault(v, []).append(xr.open_dataset(f))
    merged = xr.merge([xr.concat(p, dim="depth").sortby("depth") if "depth" in p[0].dims else p[0]
                       for p in by_var.values()], compat="override", join="outer")
    return merged, pieces
