"""
gtools/validation.py — CROCO model-vs-reference validation. 

Compare CROCO output against the parent product it was downscaled from
(GLORYS reanalysis for hindcasts, Mercator anfc for forecasts), on the CROCO
grid. Compare CROCO output against satellite and some in situ products.

Everything is region-agnostic: coordinates come from the CROCO file, so the same
calls work for Canary, Gulf of Guinea, Agulhas, ...

Design
------
- Load a CROCO surface field (SST, SSH, surface currents) via gtools.postprocess.
- Load the parent field (GLORYS/Mercator: thetao->temp, zos->ssh, uo/vo->u/v).
- Regrid the parent onto the CROCO grid (bilinear, scipy.griddata).
- Compare: difference map + domain statistics (bias, centred RMSE, correlation,
  pattern correlation) \u2014 the statistics adapt somisana's validation.statistics.
- Plot: 3-panel (CROCO | parent | difference) maps, and an SST+wind overlay.

Typical use
-----------
    import gtools.validation as val

    val.compare_sst(
        croco_his="model-runs/Canary_12/20251225/hcast/CROCO_FILES/croco_his.nc",
        parent="downloaded_data/GLORYS/2025_12.nc",
        date="2025-12-27", out="sst_vs_glorys.png")

    val.sst_with_wind(
        croco_his="...croco_his.nc",
        era5_dir="downloaded_data/ERA5/for_croco",
        date="2025-12-27", out="sst_wind.png")
"""
from __future__ import annotations
import sys, os
import subprocess
import glob
from pathlib import Path
import json
import base64
import html as _html
from datetime import datetime, timezone, timedelta, date
import calendar
import time
import threading
import getpass
from math import floor, ceil
import numpy as np
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt

import copernicusmarine

import gtools.postprocess as pp
from gtools.define_attrs import apply_attrs
from gtools.plotting import add_map_features, percentile_clim

try:
    from scipy.interpolate import griddata
    _HAS_SCIPY = True
except Exception:
    _HAS_SCIPY = False

try:
    import cartopy.crs as ccrs
    import cartopy.feature as cfeature
    _HAS_CARTOPY = True
except Exception:
    _HAS_CARTOPY = False


# GLORYS/Mercator (CMEMS) variable names -> our short names
PARENT_VARS = {
    "temp": "thetao",
    "ssh":  "zos",
    "salt": "so",
    "u":    "uo",
    "v":    "vo",
}


# ======================================================================
# downloading Copernicus Marine product
# ======================================================================
def is_valid_netcdf_file(file_path):
    # Short-circuit on a plainly-missing file: this is the common case
    # (checking whether a not-yet-downloaded file already exists), and
    # calling xr.open_dataset() on a missing path makes the underlying
    # HDF5 C library dump a scary-looking (but harmless) diagnostic to
    # stderr before Python ever sees the exception. Skip that entirely.
    if not os.path.exists(file_path):
        return False
    try:
        with xr.open_dataset(file_path) as ds:
            return True
    except:
        return False


def netcdf_covers_time_range(file_path, start_date, end_date, tolerance_days=1,
                             max_step_hours=None):
    """Check whether an existing NetCDF file (a downloaded Mercator variable
    file, or the combined MERCATOR_<cycle>_00.nc reference) actually covers
    [start_date, end_date] AT THE EXPECTED TEMPORAL RESOLUTION, not just
    that it exists and its dates happen to span the right window.

    A bare os.path.exists()/is_valid_netcdf_file() check isn't enough to
    decide whether a re-download is needed: a previous interrupted run, or
    one where a variable download silently failed, can leave a perfectly
    valid, openable NetCDF file that only covers ONE day instead of the
    full cycle window -- and every future run would then see "file exists"
    and skip downloading the rest forever. This checks the file's actual
    time coverage instead.

    A DATE-RANGE-only check still isn't enough on its own: a stale file
    downloaded at a different resolution than currently expected (e.g. a
    leftover hourly PT1H-m file from an earlier version of this pipeline,
    when the resolution used here later changes again) spans the same
    dates but has a different time step -- so it would wrongly pass a
    coverage-only check and never get replaced. Pass max_step_hours (e.g.
    ~30 for the daily-mean P1D-m downloads/reference used here) to also
    reject files whose actual time step is coarser than expected.

    tolerance_days allows for the file's first/last timestep not landing
    exactly on start_date/end_date (e.g. data starting a few hours into
    the requested day).
    """
    if not is_valid_netcdf_file(file_path):
        return False
    try:
        with xr.open_dataset(file_path) as ds:
            if "time" not in ds:
                return False
            times = pd.to_datetime(np.atleast_1d(ds["time"].values))
            if len(times) == 0:
                return False
            tol = pd.Timedelta(days=tolerance_days)
            covers = (times.min() <= pd.Timestamp(start_date) + tol and
                     times.max() >= pd.Timestamp(end_date) - tol)
            if not covers:
                return False
            if max_step_hours is not None and len(times) > 1:
                median_step = np.median(np.diff(np.sort(times))).astype("timedelta64[m]")
                median_step_hours = median_step / np.timedelta64(60, "m")
                if median_step_hours > max_step_hours:
                    return False
            return True
    except Exception:
        return False


def ensure_cmems_login(usrname=None, passwd=None):
    """
    Make sure copernicusmarine is authenticated.
    - If already logged in (valid stored credentials), do nothing.
    - Otherwise prompt for username/password (or use those passed in) and log in once.
    Credentials are stored by copernicusmarine in ~/.copernicusmarine/, NOT in this script.

    This is the ONLY place credentials are used. All download functions rely on the
    stored login, so they take no username/password for the actual data transfer.
    """
    # 1) already logged in? login(check_credentials_valid=True) returns a bool
    #    and does NOT prompt -- the Python-API equivalent of the CLI's
    #    `copernicusmarine login --check-credentials-valid`.
    try:
        if copernicusmarine.login(check_credentials_valid=True):
            print("CMEMS: already logged in.")
            return
    except Exception:
        pass  # not logged in yet, or an invalid/stale credentials file -- fall through to (re)login

    # 2) not logged in -> get credentials (interactive prompt if not supplied)
    print("CMEMS: no valid login found. Please enter your Copernicus Marine credentials.")
    print("  (Register free at https://data.marine.copernicus.eu/register)")
    if not usrname:
        usrname = input("  CMEMS username: ").strip()
    if not passwd:
        passwd = getpass.getpass("  CMEMS password: ")  # hidden input, not echoed

    # 3) store the login for future runs (writes ~/.copernicusmarine credentials)
    try:
        copernicusmarine.login(username=usrname, password=passwd)
    except Exception as e:
        raise RuntimeError(f"CMEMS login failed. Check your username/password.\n{e}")
    print("CMEMS: login successful (credentials stored for future runs).")


def download_cmems(dataset, varlist, start_date, end_date, domain, depths, outputDir, fname, ver='',
                   expected_max_step_hours=None):
    """
    Generic function to download a subset of a CMEMS dataset.
    Assumes ensure_cmems_login() has already run, so no credentials are passed here.

    expected_max_step_hours : if given, an existing file that covers the
        right dates but whose time step is coarser than this (e.g. a stale
        daily-mean file left over from before a product switched to
        hourly) is treated as needing a re-download too, not just an
        outright-missing-dates file. See netcdf_covers_time_range().
    """
    # skip this file if it already exists AND covers the requested date range
    # at the expected resolution (an existence/coverage-only check would
    # wrongly treat a valid-but-incomplete file -- e.g. left over from an
    # interrupted previous run, or a stale file from before a resolution
    # change -- as "already downloaded", permanently skipping the
    # missing/finer-resolution data on every future run)
    f = os.path.normpath(os.path.join(outputDir, fname))
    if netcdf_covers_time_range(f, start_date, end_date, max_step_hours=expected_max_step_hours):
        print("file already exists and covers the requested range - " + fname)
        return
    elif is_valid_netcdf_file(f):
        print(f"file exists but doesn't cover {start_date.date()} -> {end_date.date()} "
             f"(or is at a coarser resolution than expected) - re-downloading {fname}")
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
        overwrite=True,   # the Python API never prompts interactively (unlike the CLI) --
                          # force_download is a deprecated no-op here, so it is omitted
        disable_progress_bar=True,
    )
    if ver:
        subset_kwargs["dataset_version"] = ver

    # allow for a few retries if there was a temporary download error
    MAX_RETRIES = 3
    RETRY_WAIT = 10

    i = 0
    while i < MAX_RETRIES:
        print(f"Attempt {i+1} of {MAX_RETRIES}")
        try:
            copernicusmarine.subset(**subset_kwargs)
            if is_valid_netcdf_file(f):
                print("Completed " + fname)
                return
            else:
                if os.path.exists(f):
                    os.unlink(f)
                raise Exception(f"Mercator download failed (bad NetCDF output): {fname}")
        except Exception as e:
            i += 1
            if i >= MAX_RETRIES:
                print(f"Error: {e}. Giving up after {MAX_RETRIES} attempts - {fname}")
                return
            print(f"Error: {e}, retrying in {RETRY_WAIT} seconds...")
            time.sleep(RETRY_WAIT)
            continue


def download_mercator_ops_sequential(domain, run_date, hdays, fdays, outputDir, usrname=None, passwd=None):
    """Thread-free variant of download_mercator_ops() below: downloads the
    same 4 variables (so, thetao, zos, uo_vo), same file layout, same merge
    step -- but ONE AT A TIME instead of via 4 concurrent threading.Thread
    workers.

    Why this exists: download_mercator_ops() calls copernicusmarine.subset()
    from 4 threads at once. copernicusmarine's download machinery is built
    on asyncio, and the underlying netCDF4/HDF5 C library most conda/pip
    builds ship is NOT thread-safe -- four concurrent subset() calls each
    writing a NetCDF file is a known way to segfault/abort the whole Python
    process, which in Jupyter shows up as "the kernel just died", with no
    Python traceback to catch. This function trades the ~4x parallel
    speed-up for a call path this module already relies on elsewhere
    without issue -- download_cmems() is exactly the same function, just
    called in a plain `for` loop instead of one thread per variable.

    Use this instead of download_mercator_ops() if that one is crashing the
    kernel; everything else (arguments, VARIABLES list, output files,
    merge behaviour, chmod) is identical.
    """
    ensure_cmems_login(usrname, passwd)

    hdays = hdays + 1
    fdays = fdays + 1
    start_date = run_date + timedelta(days=-hdays)
    end_date = run_date + timedelta(days=fdays)

    VARIABLES = [
        {"name": "so", "id": "cmems_mod_glo_phy-so_anfc_0.083deg_P1D-m", "vars": ["so"],
         "fname": f"mercator_so_{run_date.strftime('%Y%m%d_%H')}.nc"},
        {"name": "thetao", "id": "cmems_mod_glo_phy-thetao_anfc_0.083deg_P1D-m", "vars": ["thetao"],
         "fname": f"mercator_thetao_{run_date.strftime('%Y%m%d_%H')}.nc"},
        {"name": "zos", "id": "cmems_mod_glo_phy_anfc_0.083deg_P1D-m", "vars": ["zos"],
         "fname": f"mercator_zos_{run_date.strftime('%Y%m%d_%H')}.nc"},
        {"name": "uo_vo", "id": "cmems_mod_glo_phy-cur_anfc_0.083deg_P1D-m", "vars": ["uo", "vo"],
         "fname": f"mercator_uo_vo_{run_date.strftime('%Y%m%d_%H')}.nc"},
    ]
    depths = [0.493, 5727.918]

    errors = {}
    for var in VARIABLES:
        print(f"-- downloading {var['name']} --")
        try:
            download_cmems(var["id"], var["vars"], start_date, end_date,
                           domain, depths, outputDir, var["fname"],
                           expected_max_step_hours=30)
        except Exception as e:
            errors[var["name"]] = e

    for name, e in errors.items():
        print(f"Mercator download failed for '{name}': {e}")

    available = [var for var in VARIABLES
                if is_valid_netcdf_file(os.path.join(outputDir, var["fname"]))]
    missing = [var["name"] for var in VARIABLES if var not in available]
    if missing:
        print(f"Mercator: {len(missing)}/{len(VARIABLES)} variable file(s) missing "
             f"({missing}) - merging what did download.")
    if not available:
        print("Mercator: no variable files downloaded - nothing to merge.")
        return

    print("merge NetCDF files")
    output_path = os.path.abspath(os.path.join(outputDir, f"MERCATOR_{run_date.strftime('%Y%m%d_%H')}.nc"))
    datasets = [xr.open_dataset(os.path.join(outputDir, var["fname"])) for var in available]
    merged = xr.merge(datasets)
    merged.to_netcdf(output_path, mode="w")
    for ds_ in datasets:
        ds_.close()

    if os.name != "nt":
        try:
            os.chmod(output_path, 0o775)
        except OSError as e:
            print(f"Warning: could not chmod {output_path}: {e}")


def download_mercator_ops(domain, run_date, hdays, fdays, outputDir, usrname=None, passwd=None):
    """
    Download the operational Mercator ocean output (anfc).

    Downloads each variable at the daily-mean resolution (P1D-m) -- the
    hourly product (PT1H-m) is far heavier and isn't needed for this
    validation, so the daily mean is downloaded and used as-is, no
    resampling required.

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

    # download the variable files in parallel to save time; capture any
    # per-thread exception instead of letting it vanish silently (a bare
    # threading.Thread swallows exceptions -- they'd otherwise only show up
    # as a stray missing file when we try to merge below).
    errors = {}

    def download_worker(var):
        try:
            # native P1D-m resolution is 24h; 30h slack catches any
            # drastically coarser/broken stale file without being so tight
            # that minor timestamp jitter causes a false "stale".
            download_cmems(var["id"], var["vars"], start_date, end_date,
                           domain, depths, outputDir, var["fname"],
                           expected_max_step_hours=30)
        except Exception as e:
            errors[var["name"]] = e

    threads = []
    for var in VARIABLES:
        t = threading.Thread(target=download_worker, args=(var,))
        threads.append(t)
        t.start()
    for t in threads:
        t.join()

    for name, e in errors.items():
        print(f"Mercator download failed for '{name}': {e}")

    # Only merge the variable files that actually downloaded -- a missing
    # file here (failed thread, or download_cmems giving up after retries)
    # must not crash xr.open_dataset()/xr.merge() below.
    available = [var for var in VARIABLES
                if is_valid_netcdf_file(os.path.join(outputDir, var["fname"]))]
    missing = [var["name"] for var in VARIABLES if var not in available]
    if missing:
        print(f"Mercator: {len(missing)}/{len(VARIABLES)} variable file(s) missing "
             f"({missing}) - merging what did download.")
    if not available:
        print("Mercator: no variable files downloaded - nothing to merge.")
        return

    print("merge NetCDF files")
    output_path = os.path.abspath(os.path.join(outputDir, f"MERCATOR_{run_date.strftime('%Y%m%d_%H')}.nc"))
    datasets = [xr.open_dataset(os.path.join(outputDir, var["fname"])) for var in available]
    merged = xr.merge(datasets)
    merged.to_netcdf(output_path, mode="w")
    for ds in datasets:
        ds.close()

    if os.name != "nt":   # chmod isn't a thing on Windows -- skip there rather than crash
        try:
            os.chmod(output_path, 0o775)
        except OSError as e:
            print(f"Warning: could not chmod {output_path}: {e}")


def check_cmems_product_available(dataset_id, variables=None, timeout=60):
    """Check whether a CMEMS dataset (product/subdataset) is currently
    reachable on the Copernicus Marine platform, WITHOUT downloading anything.

    Used by the validation notebook/tools so a comparison against a given
    CMEMS product (Mercator forecast, OSTIA, ODYSSEA, in-situ, ...) can be
    skipped cleanly (pass, don't crash) if that product is temporarily
    unavailable, renamed, or the user isn't entitled/logged in.

    Tries the lightweight copernicusmarine.describe(dataset_id=...) call
    (metadata only, no data transfer) via the Python API. Returns True only
    if the dataset_id is found in the catalogue; False for ANY failure
    (not installed, not logged in, network error, dataset renamed/retired,
    ...) -- the caller is expected to treat False as "skip this
    comparison", not raise. describe() raises DatasetNotFound (a subclass
    of Exception) if the id doesn't exist -- caught here like any other
    failure.
    """
    try:
        copernicusmarine.describe(dataset_id=dataset_id)
        return True
    except Exception as e:
        print(f"CMEMS availability check failed for {dataset_id} ({e}) - treating as unavailable")
        return False


# ------------------------------------------------------------------------
# Dataset-ID registry for the products used by the validation notebook
# (gtools/validation.py, validation_satellite.py, validation_godae.py).
# Centralised here so every "is this product available?" / "download this
# product" call agrees on the same CMEMS dataset id.
# ------------------------------------------------------------------------
VALIDATION_DATASETS = {
    # (i) Global Analysis & Forecast physics (the CROCO parent forecast)
    "mercator_forecast": "cmems_mod_glo_phy_anfc_0.083deg_P1D-m",
    # (ii) Satellite SST - two independent products
    "ostia_l4":    "METOFFICE-GLO-SST-L4-NRT-OBS-SST-V2",
    "odyssea_l3s": "IFREMER-GLOB-SST-L3-NRT-OBS_FULL_TIME_SERIE",
    # Satellite SSS - SMOS L4 (gap-filled). NOTE: verify this dataset id
    # against the current CMEMS catalogue (copernicusmarine.describe /
    # marine.copernicus.eu) before relying on it in production -- SSS
    # product ids have changed across catalogue versions more often than
    # the SST ones above.
    "smos_l4_sss": "cmems_obs-mob_glo_phy-sss_nrt_multi_P1D",
    # (iii) CMEMS in-situ TAC (trajectories + profiles), NRT
    "insitu_nrt":  "cmems_obs-ins_glo_phybgcwav_mynrt_na_irr",
}


def dataset_available(name_or_id):
    """Convenience wrapper: look `name_or_id` up in VALIDATION_DATASETS
    first (so callers can pass the short name, e.g. 'ostia_l4'), else treat
    it directly as a CMEMS dataset id, then run check_cmems_product_available.
    """
    dataset_id = VALIDATION_DATASETS.get(name_or_id, name_or_id)
    ok = check_cmems_product_available(dataset_id)
    print(f"  CMEMS product '{name_or_id}' ({dataset_id}): "
          f"{'available' if ok else 'NOT available - skipping'}")
    return ok


SATELLITE_DATASET_KEY = {"OSTIA": "ostia_l4", "ODYSSEA": "odyssea_l3s", "SMOS": "smos_l4_sss"}
def download_satellite_sst(product, domain, start_date, end_date, outputDir,
                           usrname=None, passwd=None):
    """Download a satellite validation product (SST: OSTIA/ODYSSEA, or SSS:
    SMOS), one file per day, into outputDir/<YYYY-MM-DD>.nc -- the layout
    gtools.validation_satellite expects. Despite the name (kept for
    backward compatibility), this works for any product in
    SATELLITE_DATASET_KEY, not just SST.

    product : 'OSTIA', 'ODYSSEA', or 'SMOS' (looked up in VALIDATION_DATASETS
              via SATELLITE_DATASET_KEY).
    domain  : (lon_min, lon_max, lat_min, lat_max).
    Availability-guarded: if the product isn't reachable on the CMEMS
    platform right now, prints a message and returns {} (nothing
    downloaded) rather than raising -- callers should treat that as "skip
    this comparison for this cycle", not a fatal error.

    Returns {date_str ('YYYY-MM-DD'): file_path} for every day that
    downloaded successfully (days that fail --e.g. an ODYSSEA swath gap or
    a SMOS coastal/ice gap-- are skipped individually and simply missing
    from the dict).
    """
    key = SATELLITE_DATASET_KEY[product]
    if not dataset_available(key):
        return {}
    dataset_id = VALIDATION_DATASETS[key]
    ensure_cmems_login(usrname, passwd)
    os.makedirs(outputDir, exist_ok=True)

    files = {}
    day = start_date
    while day <= end_date:
        fname = f"{day.strftime('%Y-%m-%d')}.nc"
        f = os.path.normpath(os.path.join(outputDir, fname))
        if is_valid_netcdf_file(f):
            print(f"  [{product}] {day.date()}: already downloaded - {fname}")
            files[day.strftime("%Y-%m-%d")] = f
            day += timedelta(days=1)
            continue

        try:
            copernicusmarine.subset(
                dataset_id=dataset_id,
                minimum_longitude=domain[0], maximum_longitude=domain[1],
                minimum_latitude=domain[2], maximum_latitude=domain[3],
                start_datetime=f"{day.strftime('%Y-%m-%d')}T00:00:00",
                end_datetime=f"{day.strftime('%Y-%m-%d')}T23:59:59",
                output_directory=os.path.normpath(outputDir),
                output_filename=fname,
                overwrite=True,
                disable_progress_bar=True,
            )
            if is_valid_netcdf_file(f):
                print(f"  [{product}] {day.date()}: downloaded OK - {fname}")
                files[day.strftime("%Y-%m-%d")] = f
            else:
                if os.path.exists(f):
                    os.unlink(f)
                print(f"  [{product}] {day.date()}: no data returned (gap or outage) - skipping")
        except Exception as e:
            print(f"  [{product}] {day.date()}: download failed ({e}) - skipping")
        day += timedelta(days=1)

    return files


# Copernicus Marine In-Situ TAC platform-type codes, as used in the folder
# structure and filename prefix of INSITU_GLO_PHYBGCWAV_DISCRETE_MYNRT_013_030
# (e.g. ".../history/MO/GL_TS_MO_6200069.nc" -> platform type "MO").
# Two categories, matching the file "data type" the platform produces:
#   - "trajectory" -- a moving near-surface track over time (drifters, gliders,
#     ferryboxes, saildrones, vessel tracks) -- what the in-situ validation
#     section treats as "trajectory platforms".
#   - "profile"    -- a depth-resolved cast or fixed-point time series with a
#     vertical dimension (floats, moorings, CTDs, XBTs, thermistor chains) --
#     what the in-situ validation section treats as "profile platforms".
INSITU_PLATFORM_TYPES = {
    "PF": {"label": "Profiling float (Argo)",         "kind": "profile"},
    "MO": {"label": "Fixed mooring",                   "kind": "profile"},
    "XB": {"label": "XBT (expendable bathythermograph)","kind": "profile"},
    "TX": {"label": "Thermistor chain",                 "kind": "profile"},
    "DB": {"label": "Drifting buoy",                    "kind": "trajectory"},
    "DC": {"label": "Drifting buoy (variant/current)",  "kind": "trajectory"},
    "FB": {"label": "Ferrybox (underway, merchant ship)","kind": "trajectory"},
    "TS": {"label": "Thermosalinograph (underway)",     "kind": "trajectory"},
    "VA": {"label": "Vessel of opportunity / other underway", "kind": "trajectory"},
}

def download_insitu(domain, start_date, end_date, outputDir, depth_max=2000,
                    usrname=None, passwd=None, platform_types=None):
    """Download CMEMS in-situ TAC observations (trajectories + profiles),
    one file per platform, into outputDir -- matches gtools.
    validation_godae.load_insitu_many's expected input layout.
 
    domain : (lon_min, lon_max, lat_min, lat_max).
    platform_types : which platform-type codes to download (see
        INSITU_PLATFORM_TYPES) -- defaults to all of them (Argo floats,
        drifting buoys, ferryboxes). Every OTHER platform in the index
        (CTD casts, gliders, moored buoys, voluntary-observing-ship TS,
        wave-parameter files) is excluded BEFORE anything is fetched, so
        nothing out-of-scope is ever downloaded.
    Availability-guarded: returns [] (nothing downloaded, not an error) if
    the in-situ product isn't reachable right now.
 
    IMPORTANT: the GLO in-situ product (cmems_obs-ins_glo_phybgcwav_mynrt_
    na_irr) does NOT support the `subset` service (only some regional
    in-situ products do). This follows Copernicus's own documented
    index-file workflow for this exact dataset instead:
    https://help.marine.copernicus.eu/en/articles/9630028
      1) copernicusmarine.get(..., index_parts=True) downloads 4 small
         index CSVs (history/latest/monthly/platform); we use index_history.
      2) filter that CSV to in-scope platform types (see platform_types)
         overlapping `domain` and [start_date, end_date] (client-side --
         `get` has no geographic subset option).
      3) one batched copernicusmarine.get(..., file_list=...) call to fetch
         exactly the matched files.
 
    Returns the list of downloaded file paths (possibly [] if the index
    can't be parsed or nothing overlaps this cycle -- printed, not raised).
    """
    if platform_types is None:
        platform_types = tuple(INSITU_PLATFORM_TYPES)
    if not dataset_available("insitu_nrt"):
        return []
    dataset_id = VALIDATION_DATASETS["insitu_nrt"]
    ensure_cmems_login(usrname, passwd)
    os.makedirs(outputDir, exist_ok=True)
 
    # ---- step 1: download the 4 index files (history/latest/monthly/platform) ----
    index_dir = os.path.join(outputDir, "_index")
    os.makedirs(index_dir, exist_ok=True)
    try:
        copernicusmarine.get(
            dataset_id=dataset_id,
            index_parts=True,
            output_directory=index_dir,
            overwrite=True,
            disable_progress_bar=True,
        )
    except Exception as e:
        print(f"  [in-situ] could not fetch the platform index ({e}) - skipping")
        return []
 
    history_files = list(Path(index_dir).rglob("index_history.txt"))
    if not history_files:
        print(f"  [in-situ] no index_history.txt found under {index_dir} after "
             f"index_parts=True download - skipping.")
        return []
    index_file = history_files[0]
 
    # ---- step 2: parse (sep=',', skiprows=5 -- per Copernicus's own example
    #      for this exact dataset), keep only in-scope platform types, and
    #      keep platforms overlapping domain+dates ----
    try:
        idx = pd.read_csv(index_file, sep=",", skiprows=5)
        idx.columns = [c.strip() for c in idx.columns]
    except Exception as e:
        print(f"  [in-situ] could not parse {index_file.name} ({e}) - skipping")
        return []
 
    required = ["file_name", "geospatial_lat_min", "geospatial_lat_max",
               "geospatial_lon_min", "geospatial_lon_max",
               "time_coverage_start", "time_coverage_end"]
    if not all(c in idx.columns for c in required):
        print(f"  [in-situ] index file has unexpected columns {list(idx.columns)} - skipping")
        return []
 
    # Exclude out-of-scope platform types FIRST (before any domain/time
    # filtering or download), so nothing outside Argo/drifting-buoy/
    # ferrybox is ever fetched.
    n_index = len(idx)
    idx = idx[idx["file_name"].apply(insitu_platform_type).isin(platform_types)]
    n_excluded = n_index - len(idx)
    if n_excluded:
        print(f"  [in-situ] {n_excluded}/{n_index} platform(s) in the index are "
             f"out-of-scope platform types (only "
             f"{', '.join(f'{k} ({v})' for k, v in INSITU_PLATFORM_TYPES.items() if k in platform_types)} "
             f"are downloaded here) - excluded before download.")
    if idx.empty:
        print(f"  [in-situ] no in-scope platforms in the index - nothing to download.")
        return []
 
    lon_min, lon_max, lat_min, lat_max = domain
    idx["time_coverage_start"] = pd.to_datetime(idx["time_coverage_start"], errors="coerce", utc=True)
    idx["time_coverage_end"] = pd.to_datetime(idx["time_coverage_end"], errors="coerce", utc=True)
    t_start = pd.Timestamp(start_date)
    t_start = t_start.tz_localize("UTC") if t_start.tzinfo is None else t_start
    t_end = pd.Timestamp(end_date)
    t_end = t_end.tz_localize("UTC") if t_end.tzinfo is None else t_end
 
    # overlap (not strict-containment) so platforms/trajectories that only
    # partially cross the domain or window are still included
    overlap = (
        (idx["geospatial_lat_min"] <= lat_max) & (idx["geospatial_lat_max"] >= lat_min) &
        (idx["geospatial_lon_min"] <= lon_max) & (idx["geospatial_lon_max"] >= lon_min) &
        (idx["time_coverage_start"] <= t_end) & (idx["time_coverage_end"] >= t_start)
    )
    matched = idx.loc[overlap]
    if matched.empty:
        print(f"  [in-situ] {len(idx)} platform(s) in the index, none overlap "
             f"this cycle's domain+window - nothing to download.")
        return []
 
    # skip files already downloaded locally (unlike the satellite/Mercator
    # downloaders, this step previously had no such check, so every rerun
    # redownloaded every already-present platform file)
    def _already_local(file_name):
        local_path = os.path.join(outputDir, os.path.basename(file_name))
        return is_valid_netcdf_file(local_path)
 
    already_have = matched["file_name"].apply(_already_local)
    to_fetch = matched.loc[~already_have]
    n_cached = int(already_have.sum())
    if n_cached:
        print(f"  [in-situ] {n_cached}/{len(matched)} platform file(s) already "
             f"downloaded - reusing, not refetching.")
 
    if to_fetch.empty:
        print(f"  [in-situ] all {len(matched)} matched platform file(s) already local.")
    else:
        print(f"  [in-situ] downloading {len(to_fetch)}/{len(matched)} new platform file(s).")
 
        # ---- step 3: one batched download of just the NEW matched files ----
        file_list_path = os.path.join(outputDir, "_insitu_file_list.txt")
        to_fetch["file_name"].to_csv(file_list_path, index=False, header=False)
        try:
            copernicusmarine.get(
                dataset_id=dataset_id,
                file_list=file_list_path,
                output_directory=os.path.normpath(outputDir),
                overwrite=True,
                no_directories=True,
                disable_progress_bar=True,
            )
        except Exception as e:
            print(f"  [in-situ] batch download failed ({e}) - continuing with whatever is cached locally")
 
    files = sorted(str(p) for p in Path(outputDir).glob("*.nc")
                   if is_valid_netcdf_file(str(p)) and insitu_platform_type(p.name) in platform_types)
    print(f"  [in-situ] {len(files)} in-scope platform file(s) available in {outputDir}")
    return files


def download_cmems_monthly(dataset, domain, start_date, end_date, varlist, depths,
                           outputDir, usrname=None, passwd=None):
    """
    Download month by month for any dataset on CMEMS.
    Credentials optional (see download_mercator_ops docstring).
    """
    ensure_cmems_login(usrname, passwd)

    os.makedirs(outputDir, exist_ok=True)

    downloadDate = start_date
    while downloadDate <= end_date:
        print(downloadDate.strftime('%Y-%m'))

        # start and end days of this month
        start_date_download = datetime(downloadDate.year, downloadDate.month, 1)
        day_end = calendar.monthrange(downloadDate.year, downloadDate.month)[1]
        end_date_download = datetime(downloadDate.year, downloadDate.month, day_end)

        # output filename
        fname = str(downloadDate.strftime('%Y_%m')) + '.nc'

        download_cmems(dataset, varlist, start_date_download, end_date_download,
                       domain, depths, outputDir, fname)

        downloadDate = downloadDate + timedelta(days=32)  # ensures we reach next month
        downloadDate = datetime(downloadDate.year, downloadDate.month, 1)


# ======================================================================
# Loading the parent (GLORYS / Mercator)
# ======================================================================
def load_parent(fname, var, date=None, depth_m=None):
    """Load a field from a GLORYS/Mercator file, as (lon2d, lat2d, field2d).

    var     : one of 'temp','ssh','salt','u','v' (mapped to CMEMS names).
    date    : 'YYYY-MM-DD' (nearest time selected); None -> first time.
    depth_m : None -> surface; else nearest parent depth level to depth_m (m).
    """
    ds = xr.open_dataset(fname)
    cmems = PARENT_VARS[var]
    if cmems not in ds:
        raise KeyError(f"{cmems} not in {fname}; has {list(ds.data_vars)}")
    da = ds[cmems]

    # pick time
    if "time" in da.dims:
        da = da.sel(time=np.datetime64(date), method="nearest") if date else da.isel(time=0)
    # depth: surface (index 0) or nearest level to depth_m
    if "depth" in da.dims:
        if depth_m is None:
            da = da.isel(depth=0)
        else:
            da = da.sel(depth=abs(depth_m), method="nearest")

    lon = ds["longitude"].values
    lat = ds["latitude"].values
    lon2d, lat2d = np.meshgrid(lon, lat)
    field = da.values
    ds.close()
    return lon2d, lat2d, field


# ======================================================================
# Regridding parent -> CROCO grid
# ======================================================================
def regrid_to_croco(plon, plat, pfield, croco_ds, method="linear"):
    """Interpolate a parent field (regular grid) onto the CROCO curvilinear grid.

    Returns a 2D array on (eta_rho, xi_rho), land left as NaN via the CROCO mask.
    """
    if not _HAS_SCIPY:
        raise RuntimeError("scipy is required for regridding (pip install scipy)")
    clon, clat, cmask = pp.lonlatmask(croco_ds)
    pts = np.column_stack([plon.ravel(), plat.ravel()])
    vals = pfield.ravel()
    good = np.isfinite(vals)
    out = griddata(pts[good], vals[good], (clon, clat), method=method)
    return out * cmask


# ======================================================================
# Statistics (adapted from somisana validation.statistics)
# ======================================================================
def domain_statistics(model, ref):
    """Compare two 2D fields on the same grid. Returns a dict of stats.

    bias        : mean(model - ref)
    rmse        : root-mean-square difference (total)
    crmse       : centred RMSE (after removing each field's mean)
    corr        : Pearson correlation of the two fields
    model_mean, ref_mean, model_min/max, ref_min/max
    Only points where BOTH are finite are used.
    """
    m = np.asarray(model).ravel()
    r = np.asarray(ref).ravel()
    ok = np.isfinite(m) & np.isfinite(r)
    m, r = m[ok], r[ok]
    if m.size < 2:
        return {k: np.nan for k in
                ("bias", "rmse", "crmse", "corr", "n",
                 "model_mean", "ref_mean", "model_min", "model_max",
                 "ref_min", "ref_max")}
    bias = float(np.mean(m - r))
    rmse = float(np.sqrt(np.mean((m - r) ** 2)))
    ma, ra = m - m.mean(), r - r.mean()
    crmse = float(np.sqrt(np.mean((ma - ra) ** 2)))
    # correlation, guarded against zero variance (constant fields)
    if np.std(m) < 1e-12 or np.std(r) < 1e-12:
        corr = np.nan
    else:
        corr = float(np.corrcoef(m, r)[0, 1])
    return {
        "n": int(m.size),
        "bias": bias, "rmse": rmse, "crmse": crmse, "corr": corr,
        "model_mean": float(m.mean()), "ref_mean": float(r.mean()),
        "model_min": float(m.min()), "model_max": float(m.max()),
        "ref_min": float(r.min()), "ref_max": float(r.max()),
    }


def _print_stats(name, s):
    print(f"  [{name}]  n={s['n']}  bias={s['bias']:+.3f}  "
          f"RMSE={s['rmse']:.3f}  cRMSE={s['crmse']:.3f}  corr={s['corr']:.3f}")


# ======================================================================
# Plot helpers
# ======================================================================
def _panel(ax, lon, lat, field, cmap, vmin, vmax, title, cbar_label):
    if _HAS_CARTOPY:
        h = ax.pcolormesh(lon, lat, field, cmap=cmap, vmin=vmin, vmax=vmax,
                          shading="auto", transform=ccrs.PlateCarree(), zorder=1)
        add_map_features(ax)
    else:
        h = ax.pcolormesh(lon, lat, field, cmap=cmap, vmin=vmin, vmax=vmax, shading="auto")
    ax.set_title(title, fontsize=10)
    cb = plt.colorbar(h, ax=ax, shrink=0.75, pad=0.03)
    cb.set_label(cbar_label, fontsize=8)
    return h


# ======================================================================
# Fixed colour-bar limits for BIAS and RMSE maps, per variable.
#
# Bias maps use the SYMMETRIC range (-limit, +limit) (bias can be either
# sign). RMSE maps use the POSITIVE half of the same range, (0, +limit)
# (RMSE is never negative, so using the same magnitude keeps the two
# panels' colour scales visually comparable at a glance without wasting
# half the RMSE colourbar on values that can't occur).
#
# CROCO/parent panels are NOT covered by this table -- they keep the
# existing percentile_clim (2nd/98th percentile) behaviour, since the
# raw field's range (e.g. 15-28 degC for tropical SST) isn't a fixed
# constant the way an error range plausibly is.
# ======================================================================
BIAS_RMSE_LIMITS = {
    "temp": 1.5,    # degC
    "salt": 1.0,    # PSU
    "speed": 0.3,   # m/s (also used for satellite/GODAE 'speed'/currents)
    "ssh": 0.03,    # m
}


def _save_or_return(fig, out):
    """Shared savefig-with-cartopy-gridliner-fallback, used by every map
    figure in this module. See the comment inside for why the retry
    exists (a known cartopy/shapely edge case, not a bug in this code)."""
    if not fig.get_constrained_layout():
        fig.tight_layout()
    if out:
        try:
            fig.savefig(out, dpi=150, bbox_inches="tight")
        except Exception as e:
            # Defense-in-depth: cartopy's gridliner can still hit degenerate-
            # geometry cases (map-boundary polygon construction failing in
            # shapely) beyond the known geo_labels trigger already disabled
            # in _panel() above. Rather than crash the whole notebook cell,
            # drop every gridliner on this figure and retry once.
            print(f"  warning: first savefig attempt failed ({e}); "
                 f"retrying with gridlines removed...")
            for a in fig.axes:
                for artist in list(getattr(a, "_gridliners", [])):
                    artist.remove()
            fig.savefig(out, dpi=150, bbox_inches="tight")
        plt.close(fig)
        return out
    return fig


def _panel_nocb(ax, lon, lat, field, cmap, vmin, vmax, title):
    """Like _panel(), but WITHOUT its own colorbar -- for figures where
    several panels share one colorbar (added separately by the caller)."""
    if _HAS_CARTOPY:
        h = ax.pcolormesh(lon, lat, field, cmap=cmap, vmin=vmin, vmax=vmax,
                          shading="auto", transform=ccrs.PlateCarree(), zorder=1)
        add_map_features(ax)
    else:
        h = ax.pcolormesh(lon, lat, field, cmap=cmap, vmin=vmin, vmax=vmax, shading="auto")
    ax.set_title(title, fontsize=10)
    return h


def _four_panel(clon, clat, croco, parent, bias_map, rmse_map, title, cbar_label,
                cmap="RdYlBu_r", dcmap="RdBu_r", rcmap="magma_r",
                bias_limit=None, out=None, n_days=None):
    """2x2 comparison figure: [CROCO, parent] top row / [bias, RMSE] bottom
    row. Exactly THREE shared colourbars (not four individual ones): one
    for the top row (CROCO+parent, same scale so the two are directly
    comparable), one for bias, one for RMSE.
 
    bias_limit : if given, bias uses (-bias_limit, +bias_limit) and RMSE
                 uses (0, +bias_limit) -- see BIAS_RMSE_LIMITS. If None,
                 both fall back to a percentile-based range (98th
                 percentile of |bias| / of RMSE).
    n_days     : if given, noted in the RMSE panel's title (e.g. "RMSE
                 (5 days)") so it's clear whether this is a genuine
                 multi-day temporal RMSE or the single-snapshot degenerate
                 case (n_days=1, mathematically RMSE==|bias| for one pair).
    """
    proj = {"projection": ccrs.PlateCarree()} if _HAS_CARTOPY else {}
    fig, axes = plt.subplots(2, 2, figsize=(11, 10), subplot_kw=proj,
                             constrained_layout=True)
 
    vmin, vmax = percentile_clim(croco, parent)
    if bias_limit is not None:
        blim = bias_limit
        rlo, rhi = 0.0, bias_limit
    else:
        finite_bias = bias_map[np.isfinite(bias_map)]
        blim = np.nanpercentile(np.abs(finite_bias), 98) if finite_bias.size else 1.0
        finite_rmse = rmse_map[np.isfinite(rmse_map)]
        rlo, rhi = 0.0, (np.nanpercentile(finite_rmse, 98) if finite_rmse.size else 1.0)
 
    h_topl = _panel_nocb(axes[0, 0], clon, clat, croco, cmap, vmin, vmax, f"CROCO {title}")
    h_topr = _panel_nocb(axes[0, 1], clon, clat, parent, cmap, vmin, vmax, f"Parent {title}")
    h_bias = _panel_nocb(axes[1, 0], clon, clat, bias_map, dcmap, -blim, blim,
                         "Bias (CROCO - parent)")
    rmse_title = "RMSE" + (f" ({n_days} day{'s' if n_days != 1 else ''})" if n_days else "")
    h_rmse = _panel_nocb(axes[1, 1], clon, clat, rmse_map, rcmap, rlo, rhi, rmse_title)
 
    cb_topl = fig.colorbar(h_topl, ax=axes[0, 0], shrink=0.85, pad=0.02, location="right")
    cb_topl.set_label(cbar_label, fontsize=9)
    cb_topr = fig.colorbar(h_topr, ax=axes[0, 1], shrink=0.85, pad=0.02, location="right")
    cb_topr.set_label(cbar_label, fontsize=9)
    cb_bias = fig.colorbar(h_bias, ax=axes[1, 0], shrink=0.85, pad=0.02, location="right")
    cb_bias.set_label(cbar_label, fontsize=9)
    cb_rmse = fig.colorbar(h_rmse, ax=axes[1, 1], shrink=0.85, pad=0.02, location="right")
    cb_rmse.set_label(cbar_label, fontsize=9)
 
    fig.suptitle(title, fontsize=13)
    return _save_or_return(fig, out)
 
 
# ======================================================================
# Public comparisons
# ======================================================================
def _maybe_daily_mean(ds, date, daily_mean):
    """If daily_mean, average CROCO records that fall on `date` into a single
    daily-mean dataset (matching GLORYS/Mercator daily averaging). Returns a ds
    whose time dim is length 1 (the daily mean), so tindex=-1 selects it.
    If date is None, averages ALL records to one daily mean.
    """
    if not daily_mean or "time" not in ds.dims:
        return ds
    if date is not None:
        d = np.datetime64(date, "D")
        tt = pp.times(ds).astype("datetime64[D]")
        sel = (tt == d)
        if sel.any():
            ds = ds.isel(time=np.where(sel)[0])
    # resample to daily means (date=None keeps ALL days; a date keeps just that
    # day -> length 1). resample adds a length-1 time axis to static grid vars,
    # so restore those from the original ds to keep them 2D.
    static = ("lon_rho", "lat_rho", "lon_u", "lat_u", "lon_v", "lat_v",
              "mask_rho", "h", "pm", "pn", "f", "angle", "hc", "Vtransform",
              "theta_s", "theta_b")
    ds_orig = ds
    # keep only numeric time-varying vars for the mean; string/char vars (e.g.
    # 'spherical') and non-time vars break resample().mean()
    import numpy as _np
    keep = [v for v in ds.data_vars
            if "time" in ds[v].dims and _np.issubdtype(ds[v].dtype, _np.number)]
    ds = ds[keep].resample(time="1D").mean()
    ds = ds.transpose("time", ...)
    for v in static:
        if v in ds_orig:
            ds[v] = ds_orig[v]
    ds.attrs = dict(ds_orig.attrs)
    return ds


def _tindex_for_date(ds, date, tindex):
    """If a date is given (and not daily-mean), pick the CROCO record matching it
    so CROCO and the parent are compared on the SAME day. Falls back to the given
    tindex when date is None or no record matches."""
    if date is None:
        return tindex
    try:
        tt = pp.times(ds).astype("datetime64[D]")
        hits = np.where(tt == np.datetime64(date, "D"))[0]
        if hits.size:
            return int(hits[-1])          # last record on that day
    except Exception:
        pass
    return tindex


# ======================================================================
# Multi-day bias/RMSE map accumulation -- shared by compare_sst/sss/ssh/
# currents below. RMSE at a single grid point isn't meaningful from a
# single (CROCO, parent) pair; it's the temporal RMSE over the set of
# `all_dates` given (typically the whole forecast cycle, CYCLE_DAYS in the
# notebook) at each point. If all_dates is a single date (or None), this
# degenerates mathematically to RMSE == |bias| for that one day, which
# is correctly labelled as n_days=1 in the resulting figure rather than
# silently implying a real multi-day statistic.
# ======================================================================
def _field_pair_for_day(ds, parent, var, day, depth_m, tindex, daily_mean, margin_deg):
    """One (croco_2d, parent_on_croco_2d) pair for one day, for var in
    'temp','salt','ssh','speed'. Shared building block for the multi-day
    accumulator and for the single-day panels."""
    ds_use = ds
    if margin_deg is not None:
        ds_use = pp.crop_interior(ds_use, margin_deg=margin_deg)
    ti = tindex
    if daily_mean:
        ds_use = _maybe_daily_mean(ds_use, day, daily_mean)
    else:
        ti = _tindex_for_date(ds_use, day, tindex)

    if var == "ssh":
        _, _, mask = pp.lonlatmask(ds_use)
        croco = ds_use["zeta"].isel(time=ti).values * mask
        plon, plat, pfield = load_parent(parent, "ssh", date=day)
        parent_on_croco = regrid_to_croco(plon, plat, pfield, ds_use)
        croco = croco - np.nanmean(croco)
        parent_on_croco = parent_on_croco - np.nanmean(parent_on_croco)
    elif var == "speed":
        if depth_m is None:
            ur, vr = pp.surface_uv(ds_use, tindex=ti)
            ue, vn = pp.rotate_uv(ds_use, ur, vr)
        else:
            ue, vn = pp.uv_at_depth(ds_use, depth_m, tindex=ti, rotate=True)
        croco = pp.speed(np.squeeze(ue), np.squeeze(vn))
        plon, plat, pu = load_parent(parent, "u", date=day, depth_m=depth_m)
        _, _, pv = load_parent(parent, "v", date=day, depth_m=depth_m)
        pu_c = regrid_to_croco(plon, plat, pu, ds_use)
        pv_c = regrid_to_croco(plon, plat, pv, ds_use)
        parent_on_croco = pp.speed(pu_c, pv_c)
    else:
        if depth_m is None:
            croco = pp.surface(ds_use, var, tindex=ti).values
        else:
            croco = pp.field_at_depth(ds_use, var, depth_m, tindex=ti).values
        plon, plat, pfield = load_parent(parent, var, date=day, depth_m=depth_m)
        parent_on_croco = regrid_to_croco(plon, plat, pfield, ds_use)
    return croco, parent_on_croco


def _multiday_bias_rmse(ds, parent, var, all_dates, depth_m=None, tindex=-1,
                        daily_mean=False, margin_deg=None):
    """Accumulate (croco - parent) across every day in all_dates and return
    (bias_map, rmse_map, n_days) -- bias_map = temporal mean, rmse_map =
    sqrt(temporal mean of squared diff), both per grid point.
    """
    if not all_dates:
        return None, None, 0
    diffs = []
    for day in all_dates:
        try:
            croco, parent_on_croco = _field_pair_for_day(
                ds, parent, var, day, depth_m, tindex, daily_mean, margin_deg)
            diffs.append(croco - parent_on_croco)
        except Exception as e:
            print(f"  (skipping {day} in multi-day RMSE accumulation: {e})")
    if not diffs:
        return None, None, 0
    stack = np.stack(diffs, axis=0)   # (n_days, eta, xi)
    bias_map = np.nanmean(stack, axis=0)
    rmse_map = np.sqrt(np.nanmean(stack ** 2, axis=0))
    return bias_map, rmse_map, len(diffs)


def _compare_grid_field(croco_his, parent, var, label, cbar_label, cmap, date=None,
                        tindex=-1, out=None, Yorig=None, margin_deg=None,
                        daily_mean=False, depth_m=None, all_dates=None):
    """Shared implementation behind compare_sst/sss/ssh/currents: 2x2 figure
    (CROCO, parent, bias, RMSE) + stats, for one day (`date`), with the RMSE
    panel computed across `all_dates` if given (else degenerates to this
    single day, n_days=1). See _three_panel / _multiday_bias_rmse.
    """
    ds = pp.open_history(croco_his, Yorig=Yorig)
    croco, parent_on_croco = _field_pair_for_day(
        ds, parent, var, date, depth_m, tindex, daily_mean, margin_deg)
    clon, clat, _ = pp.lonlatmask(ds)

    diff = croco - parent_on_croco
    rmse = 0.5*(np.sqrt(diff))
    bias_map, rmse_map, n_days = diff, rmse, 1

    s = domain_statistics(croco, parent_on_croco)
    print(f"{label}  CROCO vs parent:"); _print_stats(label, s)
    bias_limit = BIAS_RMSE_LIMITS.get(var)

    cmap = "viridis" if var in ("salt", "speed") else ("Spectral_r" if var == "ssh" else "RdYlBu_r")
    res = _four_panel(clon, clat, croco, parent_on_croco, bias_map, rmse_map,
                      label, cbar_label, cmap=cmap, bias_limit=bias_limit,
                      out=out, n_days=n_days)
    ds.close()
    return res, s


def compare_sst(croco_his, parent, date=None, tindex=-1, out=None, Yorig=None,
                margin_deg=None, daily_mean=False, depth_m=None, all_dates=None):
    """CROCO temperature vs parent on the CROCO grid: 2x2 (CROCO, parent,
    bias, RMSE) map + stats.

    `date` selects the parent field AND the matching CROCO record shown in
    the CROCO/parent/bias panels. `all_dates`, if given (e.g. the whole
    cycle's day list), is used for the RMSE panel (temporal RMSE across
    every day in it) -- otherwise RMSE degenerates to |bias| for `date`
    alone. depth_m=None -> surface (SST); else temperature at that true
    depth. margin_deg trims the sponge band from both.
    """
    dlab = "SST" if depth_m is None else f"temp {depth_m:g} m"
    return _compare_grid_field(croco_his, parent, "temp", dlab, "temperature (\u00b0C)",
                               "RdYlBu_r", date=date, tindex=tindex, out=out, Yorig=Yorig,
                               margin_deg=margin_deg, daily_mean=daily_mean,
                               depth_m=depth_m, all_dates=all_dates)


def compare_sss(croco_his, parent, date=None, tindex=-1, out=None, Yorig=None,
                margin_deg=None, daily_mean=False, depth_m=None, all_dates=None):
    """CROCO salinity vs parent on the CROCO grid: 2x2 map + stats. Mirrors
    compare_sst exactly, for salt/'so' instead of temp/'thetao'.
    depth_m=None -> surface (SSS); else salinity at that true depth.
    """
    dlab = "SSS" if depth_m is None else f"salt {depth_m:g} m"
    return _compare_grid_field(croco_his, parent, "salt", dlab, "salinity (PSU)",
                               "RdYlBu_r", date=date, tindex=tindex, out=out, Yorig=Yorig,
                               margin_deg=margin_deg, daily_mean=daily_mean,
                               depth_m=depth_m, all_dates=all_dates)


def compare_ssh(croco_his, parent, date=None, tindex=-1, out=None, Yorig=None,
                margin_deg=None, daily_mean=False, all_dates=None):
    """CROCO SSH (zeta) vs parent SSH (zos) on the CROCO grid: 2x2 map +
    stats. Compared as anomalies (domain mean removed from both) -- CROCO's
    zeta has no absolute geoid reference, so a raw-level comparison would
    be dominated by that arbitrary offset rather than genuine pattern
    disagreement."""
    return _compare_grid_field(croco_his, parent, "ssh", "SSH anomaly", "SSH' (m)",
                               "viridis", date=date, tindex=tindex, out=out, Yorig=Yorig,
                               margin_deg=margin_deg, daily_mean=daily_mean,
                               depth_m=None, all_dates=all_dates)


def compare_currents(croco_his, parent, date=None, tindex=-1, out=None, skip=4, Yorig=None,
                     margin_deg=None, daily_mean=False, depth_m=None, all_dates=None):
    """CROCO vs parent current SPEED (shaded): 2x2 map + stats.
    depth_m=None -> surface; else compares at that true depth (CROCO
    interpolated from sigma, parent at its nearest depth level)."""
    dlab = "surface" if depth_m is None else f"{depth_m:g} m"
    return _compare_grid_field(croco_his, parent, "speed", f"speed ({dlab})",
                               "speed (m s$^{-1}$)", "viridis", date=date, tindex=tindex,
                               out=out, Yorig=Yorig, margin_deg=margin_deg,
                               daily_mean=daily_mean, depth_m=depth_m, all_dates=all_dates)


def compare_forecast_depth_levels(croco_his, parent, var="temp", date=None,
                                  depths=(0, 120, 300, 1000),
                                  Yorig=None, margin_deg=None, daily_mean=False,
                                  tindex=-1, out=None):
    """CROCO vs Copernicus Marine Forecast, 4x4 figure: one ROW per depth
    level (surface, 120 m, 300 m, 1000 m by default), one COLUMN per panel
    type (CROCO, Copernicus, bias, RMSE) -- the depth-resolved counterpart
    to compare_sst/sss/currents' single-level 2x2 figure.

    var  : 'temp', 'salt', or 'speed' (depths=(0,...) with 0 meaning
           surface, matching depth_m=None elsewhere in this module).
    date : which day's CROCO/parent snapshot this whole figure is for.
           Bias and RMSE are BOTH derived from THIS SAME single day's
           CROCO-vs-Copernicus difference (RMSE reduces to |bias| for a
           single day, n=1) -- so calling this once per day in a loop
           produces a genuinely different bias/RMSE map each time, rather
           than repeating one cycle-wide aggregate across every day.

    Colourbars: each depth ROW gets its OWN shared colourbar for its
    CROCO+Copernicus pair, with vmin/vmax taken from the 2nd/98th
    percentile of THAT ROW's own CROCO+Copernicus values only -- so a
    surface range doesn't wash out the contrast at 1000 m (or vice versa).
    Bias and RMSE each still get ONE shared colourbar across every depth,
    same percentile method as before. 6 colourbars total for 4 depths (one
    per row + one for bias + one for RMSE), not 16.
    """
    if var not in ("temp", "salt", "speed"):
        raise ValueError(f"var={var!r} must be 'temp', 'salt', or 'speed'")
    if date is None:
        raise ValueError("date is required -- bias/RMSE are computed for this single day")

    ds = pp.open_history(croco_his, Yorig=Yorig)
    clon, clat, _ = pp.lonlatmask(ds)

    label = {"temp": "temperature", "salt": "salinity", "speed": "speed"}[var]
    units = {"temp": "\u00b0C", "salt": "PSU", "speed": "m s$^{-1}$"}[var]
    cmap = "viridis" if var in ("salt", "speed") else ("Spectral_r" if var == "ssh" else "RdYlBu_r")
    bias_limit = BIAS_RMSE_LIMITS.get(var)

    rows = []   # per depth: (depth_label, croco, parent, bias_map, rmse_map)
    for d in depths:
        dm = None if d == 0 else d
        try:
            croco, parent_on_croco = _field_pair_for_day(
                ds, parent, var, date, dm, tindex, daily_mean, margin_deg)
        except Exception as e:
            print(f"  depth {d} m: could not build CROCO/parent pair ({e}) - skipping row")
            continue
        # this depth's OWN diff for THIS day only -- not pooled across the
        # cycle, so it changes from one day's call to the next.
        diff = croco - parent_on_croco
        bias_map = diff
        rmse_map = np.abs(diff)   # single-day RMSE == |diff| (n=1)
        dlab = "surface" if d == 0 else f"{d:g} m"
        rows.append((dlab, croco, parent_on_croco, bias_map, rmse_map))
    ds.close()

    if not rows:
        print(f"compare_forecast_depth_levels: no depth level could be compared for {var!r}")
        return None

    if bias_limit is not None:
        blim, rlo, rhi = bias_limit, 0.0, bias_limit
    else:
        all_bias = np.concatenate([r[3][np.isfinite(r[3])].ravel() for r in rows])
        all_rmse = np.concatenate([r[4][np.isfinite(r[4])].ravel() for r in rows])
        blim = np.nanpercentile(np.abs(all_bias), 98) if all_bias.size else 1.0
        rlo, rhi = 0.0, (np.nanpercentile(all_rmse, 98) if all_rmse.size else 1.0)

    proj = {"projection": ccrs.PlateCarree()} if _HAS_CARTOPY else {}
    fig, axes = plt.subplots(len(rows), 4, figsize=(16, 4 * len(rows)),
                             subplot_kw=proj, constrained_layout=True)
    if len(rows) == 1:
        axes = axes[None, :]

    h_bias = h_rmse = None
    for r, (dlab, croco, parent_on_croco, bias_map, rmse_map) in enumerate(rows):
        # this row's OWN colour range -- 2nd/98th percentile of its own
        # CROCO+Copernicus values only, same method as before, just no
        # longer pooled across every depth.
        row_vals = np.concatenate([croco[np.isfinite(croco)].ravel(),
                                   parent_on_croco[np.isfinite(parent_on_croco)].ravel()])
        vmin, vmax = np.nanpercentile(row_vals, [2, 98])

        h_top = _panel_nocb(axes[r, 0], clon, clat, croco, cmap, vmin, vmax,
                            f"CROCO {label}\n{dlab}")
        _panel_nocb(axes[r, 1], clon, clat, parent_on_croco, cmap, vmin, vmax,
                   f"Copernicus {label}\n{dlab}")
        h_bias = _panel_nocb(axes[r, 2], clon, clat, bias_map, "RdBu_r", -blim, blim,
                             f"Bias\n{dlab}")
        h_rmse = _panel_nocb(axes[r, 3], clon, clat, rmse_map, "magma_r", rlo, rhi,
                             f"RMSE\n{dlab}")

        cb_row = fig.colorbar(h_top, ax=[axes[r, 0], axes[r, 1]],
                              shrink=0.85, pad=0.02, location="right")
        cb_row.set_label(f"{label} ({units})", fontsize=9)

    cb_bias = fig.colorbar(h_bias, ax=list(axes[:, 2]), shrink=0.7, pad=0.02, location="right")
    cb_bias.set_label(f"bias ({units})", fontsize=9)
    cb_rmse = fig.colorbar(h_rmse, ax=list(axes[:, 3]), shrink=0.7, pad=0.02, location="right")
    cb_rmse.set_label(f"RMSE ({units})", fontsize=9)

    fig.suptitle(f"CROCO forecast vs Copernicus Marine Forecast -- {label}, {date}", fontsize=14)
    return _save_or_return(fig, out)
    

# ======================================================================
# Domain-wide difference distributions (for boxplots) and boxplot builders
# ======================================================================
def domain_diff(ds, parent, var, day, depth_m=None, tindex=-1, max_points=20000):
    """CROCO-minus-parent difference over the whole domain for one day,
    regridded onto the CROCO grid. SSH is compared as an anomaly (domain
    mean removed from both sides, matching compare_ssh) since CROCO zeta
    has no absolute geoid reference. Returns a flat array of finite
    differences (subsampled if very large) -- the building block for
    bias_boxplot()/rmse_boxplot() below.

    ds : an already-open CROCO history Dataset (pp.open_history(...)).
    """
    croco, parent_on_croco = _field_pair_for_day(
        ds, parent, var, day, depth_m, tindex, daily_mean=False, margin_deg=None)
    diff = (croco - parent_on_croco).ravel()
    diff = diff[np.isfinite(diff)]
    if diff.size > max_points:
        diff = np.random.default_rng(0).choice(diff, max_points, replace=False)
    return diff


def domain_diff_satellite(croco_field, sat_on_croco, max_points=20000):
    """Same as domain_diff(), but for an already-computed (croco, satellite)
    pair on the CROCO grid (see gtools.validation.compare_satellite_grid,
    which does the loading/regridding) -- just flattens+filters+subsamples."""
    diff = (croco_field - sat_on_croco).ravel()
    diff = diff[np.isfinite(diff)]
    if diff.size > max_points:
        diff = np.random.default_rng(0).choice(diff, max_points, replace=False)
    return diff


def bias_boxplot(diffs_by_day, days, ylabel, title, out=None, metric="bias"):
    """Boxplot of domain-wide differences, one box per day.

    diffs_by_day : list of 1D arrays (see domain_diff()), aligned with `days`.
    metric       : 'bias' (signed differences, as given) or 'rmse' (plots
                   |diffs| instead -- the per-pixel error-MAGNITUDE
                   distribution for that day; a true scalar RMSE has no
                   distribution to box, so this is the natural per-pixel
                   analogue: same underlying diffs as the bias boxplot,
                   just unsigned, so day-to-day error MAGNITUDE spread is
                   visible even where signed bias averages out near zero).
    """
    if metric == "rmse":
        diffs_by_day = [np.abs(d) for d in diffs_by_day]
    fig, ax = plt.subplots(figsize=(max(6, 1.0 * len(days)), 4.5))
    positions = np.arange(1, len(days) + 1)
    ax.boxplot(diffs_by_day, positions=positions, showfliers=False)
    ax.set_xticks(positions)
    ax.set_xticklabels(days)
    if metric == "bias":
        ax.axhline(0, color="k", ls="--", lw=1)
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.grid(alpha=0.3)
    plt.setp(ax.get_xticklabels(), rotation=45, ha="right")
    return _save_or_return(fig, out)


def bias_boxplot_multi(diffs_by_day_by_group, days, ylabel, title, out=None,
                       colors=None, metric="bias"):
    """Like bias_boxplot(), but with SEVERAL groups side-by-side per day
    (e.g. OSTIA vs ODYSSEA) -- one narrower box-cluster per day, one colour
    per group, days/groups with no data simply omitted (not errored).

    diffs_by_day_by_group : {group_name: [array_or_empty_per_day, ...]}
    """
    groups = list(diffs_by_day_by_group.keys())
    n = len(groups)
    width = 0.8 / max(n, 1)
    colors = colors or [f"C{i+1}" for i in range(n)]

    fig, ax = plt.subplots(figsize=(max(6, 1.0 * len(days)), 4.5))
    from matplotlib.patches import Patch
    handles = []
    for gi, (name, diffs) in enumerate(diffs_by_day_by_group.items()):
        if metric == "rmse":
            diffs = [np.abs(d) for d in diffs]
        offset = (gi - (n - 1) / 2) * width
        idx = [k for k, d in enumerate(diffs) if d is not None and len(d) > 0]
        if not idx:
            continue
        pos = np.array([k + 1 for k in idx], dtype=float) + offset
        data = [diffs[k] for k in idx]
        bp = ax.boxplot(data, positions=pos, widths=width * 0.9, showfliers=False,
                        patch_artist=True)
        for box in bp["boxes"]:
            box.set_facecolor(colors[gi]); box.set_alpha(0.5)
        handles.append(Patch(facecolor=colors[gi], alpha=0.5, label=name))

    if metric == "bias":
        ax.axhline(0, color="k", ls="--", lw=1)
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.grid(alpha=0.3)
    ax.set_xticks(np.arange(1, len(days) + 1))
    ax.set_xticklabels(days)
    plt.setp(ax.get_xticklabels(), rotation=45, ha="right")
    if handles:
        ax.legend(handles=handles, fontsize=8)
    else:
        ax.text(0.5, 0.5, "no data for this cycle", ha="center", va="center",
                transform=ax.transAxes, fontsize=9, color="gray")
    return _save_or_return(fig, out)

def bias_timeseries(diffs_by_day, days, ylabel, title, out=None, metric="bias"):
    """Mean +/- std timeseries of domain-wide differences, one point per
    day/lead, with a shaded std band (fill_between) -- the timeseries
    analogue of bias_boxplot(): same `diffs_by_day` input (list of 1D
    arrays, see domain_diff()/composite_domain_diff()), aligned with
    `days`, just summarised as mean+/-std instead of a full box.

    metric : 'bias' (signed differences, as given) or 'rmse' (plots the
             mean/std of |diffs| instead, same convention as bias_boxplot).
    Days/leads with no data (empty or None array) show as a gap in the
    line (NaN), not an error or a zero.
    """
    if metric == "rmse":
        diffs_by_day = [np.abs(d) if d is not None else d for d in diffs_by_day]
    means = np.array([np.mean(d) if d is not None and len(d) > 0 else np.nan
                       for d in diffs_by_day])
    stds = np.array([np.std(d) if d is not None and len(d) > 0 else np.nan
                      for d in diffs_by_day])
    x = np.arange(1, len(days) + 1)

    fig, ax = plt.subplots(figsize=(max(6, 1.0 * len(days)), 4.5))
    ax.plot(x, means, "o-", color="C0", lw=1.5, ms=5, label="mean")
    ax.fill_between(x, means - stds, means + stds, color="C0", alpha=0.25,
                     label="+/- 1 std")
    if metric == "bias":
        ax.axhline(0, color="k", ls="--", lw=1)
    ax.set_xticks(x)
    ax.set_xticklabels(days)
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.grid(alpha=0.3)
    ax.legend(fontsize=8)
    plt.setp(ax.get_xticklabels(), rotation=45, ha="right")
    return _save_or_return(fig, out)


def bias_timeseries_multi(diffs_by_day_by_group, days, ylabel, title, out=None,
                          colors=None, metric="bias"):
    """Like bias_timeseries(), but with SEVERAL groups on the same axes
    (e.g. OSTIA vs ODYSSEA) -- one mean line + std band per group,
    days/groups with no data simply gapped (not errored), matching
    bias_boxplot_multi()'s conventions.

    diffs_by_day_by_group : {group_name: [array_or_empty_per_day, ...]}
    """
    groups = list(diffs_by_day_by_group.keys())
    colors = colors or [f"C{i+1}" for i in range(len(groups))]
    x = np.arange(1, len(days) + 1)

    fig, ax = plt.subplots(figsize=(max(6, 1.0 * len(days)), 4.5))
    any_data = False
    for gi, (name, diffs) in enumerate(diffs_by_day_by_group.items()):
        if metric == "rmse":
            diffs = [np.abs(d) if d is not None else d for d in diffs]
        means = np.array([np.mean(d) if d is not None and len(d) > 0 else np.nan
                           for d in diffs])
        stds = np.array([np.std(d) if d is not None and len(d) > 0 else np.nan
                          for d in diffs])
        if np.all(np.isnan(means)):
            continue
        any_data = True
        ax.plot(x, means, "o-", color=colors[gi], lw=1.5, ms=5, label=name)
        ax.fill_between(x, means - stds, means + stds, color=colors[gi], alpha=0.2)

    if metric == "bias":
        ax.axhline(0, color="k", ls="--", lw=1)
    ax.set_xticks(x)
    ax.set_xticklabels(days)
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.grid(alpha=0.3)
    plt.setp(ax.get_xticklabels(), rotation=45, ha="right")
    if any_data:
        ax.legend(fontsize=8)
    else:
        ax.text(0.5, 0.5, "no data for this cycle", ha="center", va="center",
                transform=ax.transAxes, fontsize=9, color="gray")
    return _save_or_return(fig, out)

# ======================================================================
# Pointwise scatter (moved here from the notebook, so every plotting
# function used by 02_validation.ipynb lives in this one module)
# ======================================================================
def scatter_vs_reference(model_field, ref_field, label, units, ax=None, max_points=20000, out=None):
    """Pointwise CROCO-vs-reference scatter with a 1:1 line and the
    domain_statistics() bias/RMSE/corr annotated on the plot."""
    m = np.asarray(model_field).ravel()
    r = np.asarray(ref_field).ravel()
    ok = np.isfinite(m) & np.isfinite(r)
    m, r = m[ok], r[ok]
    if len(m) > max_points:
        idx = np.random.default_rng(0).choice(len(m), max_points, replace=False)
        m, r = m[idx], r[idx]
    s = domain_statistics(m, r)
    made_fig = ax is None
    if made_fig:
        fig, ax = plt.subplots(figsize=(5, 5))
    lo, hi = np.nanpercentile(np.concatenate([m, r]), [1, 99])
    ax.plot([lo, hi], [lo, hi], "k--", lw=1, label="1:1")
    ax.scatter(r, m, s=4, alpha=0.25, color="C0")
    ax.set_xlabel(f"reference {label} ({units})")
    ax.set_ylabel(f"CROCO {label} ({units})")
    ax.set_title(f"{label}: bias={s['bias']:+.2f}  RMSE={s['rmse']:.2f}  corr={s['corr']:.2f}")
    ax.set_xlim(lo, hi); ax.set_ylim(lo, hi)
    ax.set_aspect("equal", adjustable="box")
    ax.legend()
    if made_fig:
        return _save_or_return(fig, out), ax
    return ax


# ======================================================================
# SST + wind overlay (wind from ERA5/GFS for_croco files)
# ======================================================================
def _pick_var(ds, *candidates):
    """Find a data variable by any of `candidates` (case-insensitive), else
    fall back to the only non-coordinate variable. Used by _load_wind so a
    file written by an older/newer version of the ERA5 reconstruction script
    (internal variable named 'U10M' vs 'u10' vs 'u10m') still loads."""
    coord = {"lon", "lat", "longitude", "latitude", "time", "valid_time"}
    by_upper = {v.upper(): v for v in ds.data_vars}
    for c in candidates:
        if c.upper() in by_upper:
            return by_upper[c.upper()]
    rest = [v for v in ds.data_vars if v.lower() not in coord]
    if not rest:
        raise KeyError(
            f"no data variable found in {getattr(ds, 'encoding', {}).get('source', ds)}; "
            f"has {list(ds.data_vars)}"
        )
    return rest[0]


def _load_wind(era5_dir, date, year_month=None):
    """Load U10M/V10M from ERA5 for_croco files for a given date.

    Files are named U10M_Y<year>M<month>.nc etc.; the internal data variable
    may be named U10M, U10, u10m, or u10 depending on which version of the
    ERA5 reconstruction produced the file -- _pick_var() accepts all of them.
    Returns (lon, lat, u, v) at the nearest time to `date`.
    """
    import os
    d = np.datetime64(date)
    y = str(d.astype("datetime64[Y]"))
    m = int(str(d)[5:7])
    tag = f"Y{y}M{m:02d}"
    fu = os.path.join(era5_dir, f"U10M_{tag}.nc")
    fv = os.path.join(era5_dir, f"V10M_{tag}.nc")
    if not os.path.exists(fu) or not os.path.exists(fv):
        raise FileNotFoundError(
            f"wind files not found: {fu} / {fv}\n"
            f"  (expected U10M_{tag}.nc and V10M_{tag}.nc in {era5_dir})"
        )
    dsu = xr.open_dataset(fu); dsv = xr.open_dataset(fv)
    try:
        uvar = _pick_var(dsu, "U10M", "U10")
        vvar = _pick_var(dsv, "V10M", "V10")
        u = dsu[uvar]; v = dsv[vvar]
        if "time" in u.dims:
            u = u.sel(time=d, method="nearest"); v = v.sel(time=d, method="nearest")
        lon_name = "lon" if "lon" in dsu else "longitude"
        lat_name = "lat" if "lat" in dsu else "latitude"
        lon = dsu[lon_name].values
        lat = dsu[lat_name].values
        lon2d, lat2d = np.meshgrid(lon, lat)
        uu, vv = u.values, v.values
    finally:
        dsu.close(); dsv.close()
    return lon2d, lat2d, uu, vv


def sst_with_wind(croco_his, era5_dir, date=None, tindex=-1, out=None,
                  skip=3, title="Canary_12"):
    """CROCO SST shaded, with ERA5 10 m wind vectors on top.

    The upwelling diagnostic: cold coastal SST under upwelling-favourable wind.
    """
    ds = pp.open_history(croco_his)
    sst = pp.surface(ds, "temp", tindex=tindex).values
    clon, clat, _ = pp.lonlatmask(ds)

    # wind (regridded onto CROCO grid so vectors sit on the SST field)
    if date is None:
        date = str(np.datetime_as_string(pp.times(ds)[tindex], unit="D"))
    wlon, wlat, wu, wv = _load_wind(era5_dir, date)
    wu_c = regrid_to_croco(wlon, wlat, wu, ds)
    wv_c = regrid_to_croco(wlon, wlat, wv, ds)

    proj = {"projection": ccrs.PlateCarree()} if _HAS_CARTOPY else {}
    fig, ax = plt.subplots(figsize=(8, 8), subplot_kw=proj)
    h = _panel(ax, clon, clat, sst, "RdYlBu_r",
               np.nanmin(sst), np.nanmax(sst), f"{title} SST + 10 m wind   {date}",
               "SST (\u00b0C)")
    s = slice(None, None, skip)
    qkw = dict(color="k", scale=200, width=0.003, zorder=6)
    if _HAS_CARTOPY:
        qkw["transform"] = ccrs.PlateCarree()
    ax.quiver(clon[s, s], clat[s, s], wu_c[s, s], wv_c[s, s], **qkw)
    fig.tight_layout()
    if out:
        fig.savefig(out, dpi=150, bbox_inches="tight")
        plt.close(fig)
        ds.close()
        return out
    ds.close()
    return fig


# ======================================================================
# Along-cycle validation: how CROCO tracks the parent over lead time
# ======================================================================
def _croco_field(ds, field, tindex, depth_m=None):
    """Return a CROCO 2D field for a given time index.
    field: 'sst'|'temp', 'ssh'|'zeta', or 'speed'.
    depth_m=None -> surface; else at that true depth (temp/speed; ssh is 2D)."""
    if field in ("sst", "temp"):
        if depth_m is None:
            return pp.surface(ds, "temp", tindex=tindex).values
        return pp.field_at_depth(ds, "temp", depth_m, tindex=tindex).values
    if field in ("ssh", "zeta"):
        _, _, mask = pp.lonlatmask(ds)
        return ds["zeta"].isel(time=tindex).values * mask   # 2D, depth n/a
    if field == "speed":
        if depth_m is None:
            ur, vr = pp.surface_uv(ds, tindex=tindex)
            ue, vn = pp.rotate_uv(ds, ur, vr)
        else:
            ue, vn = pp.uv_at_depth(ds, depth_m, tindex=tindex, rotate=True)
        return pp.speed(np.squeeze(ue), np.squeeze(vn))
    raise ValueError(f"unknown field '{field}'")


def _parent_field_on_croco(parent, field, date, croco_ds, depth_m=None):
    """Load the parent field for a date and regrid onto the CROCO grid.
    depth_m=None -> surface; else nearest parent depth level (ssh stays 2D)."""
    if field in ("sst", "temp"):
        plon, plat, pf = load_parent(parent, "temp", date=date, depth_m=depth_m)
    elif field in ("ssh", "zeta"):
        plon, plat, pf = load_parent(parent, "ssh", date=date)   # 2D
    elif field == "speed":
        plon, plat, pu = load_parent(parent, "u", date=date, depth_m=depth_m)
        _, _, pv = load_parent(parent, "v", date=date, depth_m=depth_m)
        pf = np.sqrt(pu ** 2 + pv ** 2)
    else:
        raise ValueError(field)
    return regrid_to_croco(plon, plat, pf, croco_ds)


def cycle_error_growth(croco_his, parent, fields=("sst", "ssh", "speed"), Yorig=None, margin_deg=None,
                       anomaly_ssh=True, daily_mean=False, depth_m=None):
    """Track CROCO-vs-parent statistics along one cycle, per CROCO time record.

    Returns a dict:
       {'lead_days': [...], 'times': [...],
        'sst': {'bias':[...], 'rmse':[...], 'corr':[...]}, 'ssh': {...}, ...}
    Lead day is days since the first record of the cycle.

    daily_mean : if True, average the CROCO records to DAILY means before
      comparing. CROCO is usually sub-daily (e.g. 6-hourly) and carries a diurnal
      cycle, while GLORYS/Mercator are daily averages - comparing them directly
      makes the RMSE oscillate once per day. Daily-averaging removes that aliasing
      and gives a clean growth curve (the physically honest, like-with-like view).
    """
    ds = pp.open_history(croco_his, Yorig=Yorig)
    if margin_deg is not None:
        ds = pp.crop_interior(ds, margin_deg=margin_deg)
    if daily_mean and "time" in ds.dims:
        # group the sub-daily records into daily means (matches GLORYS averaging);
        # use the helper (date=None -> all days) which keeps grid vars 2D
        ds = _maybe_daily_mean(ds, None, True)
    times = pp.times(ds)
    t0 = times[0]
    lead = [(np.datetime64(t) - np.datetime64(t0)) / np.timedelta64(1, "D")
            for t in times]

    out = {"lead_days": lead, "times": [str(np.datetime_as_string(t, unit="h")) for t in times]}
    for f in fields:
        out[f] = {"bias": [], "rmse": [], "crmse": [], "corr": []}

    for k, t in enumerate(times):
        date = str(np.datetime_as_string(t, unit="D"))
        for f in fields:
            croco = _croco_field(ds, f, k, depth_m=depth_m)
            par = _parent_field_on_croco(parent, f, date, ds, depth_m=depth_m)
            a, b = croco, par
            if f in ("ssh", "zeta") and anomaly_ssh:
                a = croco - np.nanmean(croco)
                b = par - np.nanmean(par)
            s = domain_statistics(a, b)
            out[f]["bias"].append(s["bias"])
            out[f]["rmse"].append(s["rmse"])
            out[f]["crmse"].append(s["crmse"])
            out[f]["corr"].append(s["corr"])
    ds.close()
    return out


def plot_error_growth(growth, metric="rmse", out=None, title=None, figsize=(9, 5)):
    """Plot one metric vs lead day for each field, from cycle_error_growth output."""
    lead = growth["lead_days"]
    fields = [k for k in growth if k not in ("lead_days", "times")]
    fig, ax = plt.subplots(figsize=figsize)
    for f in fields:
        ax.plot(lead, growth[f][metric], marker="o", lw=1.8, label=f.upper())
    ax.set_xlabel("lead time (days)")
    ax.set_ylabel(f"{metric.upper()}  (CROCO vs parent)")
    ax.grid(alpha=0.3)
    ax.legend()
    ax.set_title(title or f"Error growth — {metric.upper()} vs lead time")
    fig.tight_layout()
    if out:
        fig.savefig(out, dpi=150, bbox_inches="tight")
        plt.close(fig)
        return out
    return fig


def combined_error_growth(run_root, parent_for_cycle, phase="hcast",
                          fields=("sst", "ssh", "speed"), metric="rmse",
                          out=None, figsize=(10, 6), Yorig=None, margin_deg=None,
                          daily_mean=False, depth_m=None):
    """Overlay error-growth curves for ALL cycles of a cycled run.

    run_root         : .../model-runs/<CONFIG>
    parent_for_cycle : callable(cycle_tag) -> parent filename, OR a single
                       parent filename string used for every cycle.
                       (cycle_tag is the dated folder name, e.g. '20251225'.)
    Produces one subplot per field, with one line per cycle + the mean.
    """
    import glob, os
    cycles = sorted(glob.glob(os.path.join(run_root, "*", phase,
                                           "CROCO_FILES", "croco_his.nc")))
    if not cycles:
        raise FileNotFoundError(f"no croco_his.nc under {run_root}/*/{phase}/CROCO_FILES")

    # gather growth per cycle
    per_cycle = []
    tags = []
    for c in cycles:
        tag = c.split(os.sep)[-4]      # the dated folder
        parent = parent_for_cycle(tag) if callable(parent_for_cycle) else parent_for_cycle
        try:
            g = cycle_error_growth(c, parent, fields=fields,
                                   Yorig=Yorig, margin_deg=margin_deg,
                                   daily_mean=daily_mean, depth_m=depth_m)
            per_cycle.append(g); tags.append(tag)
        except Exception as e:
            print(f"  skip cycle {tag}: {e}")

    if not per_cycle:
        raise RuntimeError("no cycles could be processed")

    # one subplot per field
    fig, axes = plt.subplots(1, len(fields), figsize=figsize, squeeze=False)
    axes = axes[0]
    for ax, f in zip(axes, fields):
        # align on common lead-day length
        Lmin = min(len(g["lead_days"]) for g in per_cycle)
        lead = per_cycle[0]["lead_days"][:Lmin]
        stack = []
        for g, tag in zip(per_cycle, tags):
            y = g[f][metric][:Lmin]
            ax.plot(lead, y, marker="o", lw=1.2, alpha=0.7, label=tag)
            stack.append(y)
        mean = np.nanmean(np.array(stack), axis=0)
        ax.plot(lead, mean, color="k", lw=2.5, label="mean")
        ax.set_title(f.upper())
        ax.set_xlabel("lead time (days)")
        ax.grid(alpha=0.3)
    axes[0].set_ylabel(f"{metric.upper()} (CROCO vs parent)")
    axes[-1].legend(fontsize=8)
    fig.suptitle(f"Error growth across cycles — {metric.upper()}")
    fig.tight_layout()
    if out:
        fig.savefig(out, dpi=150, bbox_inches="tight")
        plt.close(fig)
        return out
    return fig


# --- local helpers for the comparison plots ---
def _finish(fig, out, dpi=150):
    fig.tight_layout()
    if out:
        fig.savefig(out, dpi=dpi, bbox_inches="tight")
        plt.close(fig)
        return out
    return fig


def _haversine_km(lat0, lon0, lats, lons, R=6371.0):
    lat0r, lon0r = np.radians(lat0), np.radians(lon0)
    latr, lonr = np.radians(lats), np.radians(lons)
    dlat = latr - lat0r; dlon = lonr - lon0r
    a = np.sin(dlat/2)**2 + np.cos(lat0r)*np.cos(latr)*np.sin(dlon/2)**2
    return R * 2 * np.arcsin(np.sqrt(a))


# ======================================================================
# CROCO-vs-parent profile / section / time series (overlaid comparisons)
# ======================================================================
def _parent_column(fname, var, lon0, lat0, date=None):
    """Parent vertical profile (value vs depth) at (lon0,lat0).
    Returns (depths_m_negative, values). Depth is negative-down to match CROCO.
    """
    ds = xr.open_dataset(fname)
    cmems = PARENT_VARS[var]
    da = ds[cmems]
    if "time" in da.dims:
        da = da.sel(time=np.datetime64(date), method="nearest") if date else da.isel(time=0)
    da = da.sel(longitude=lon0, latitude=lat0, method="nearest")
    depth = -np.abs(ds["depth"].values)          # negative down
    vals = da.values
    ds.close()
    return depth, vals


def compare_profile(croco_his, parent, var, lon0, lat0, date=None, tindex=-1,
                    out=None, Yorig=None, daily_mean=False):
    """Overlay CROCO and parent vertical profiles of `var` at a point.

    var: 'temp'/'salt'/'u'/'v'/'speed' (mapped to CMEMS for the parent).
    Two lines on one axes (value vs depth).
    """
    ds = pp.open_history(croco_his, Yorig=Yorig)
    if daily_mean:
        ds = _maybe_daily_mean(ds, date, daily_mean)
    cpro = pp.profile(ds, var, lon0, lat0, tindex=tindex)   # CROCO
    cz = cpro["depth"].values
    cv = cpro.values

    # parent profile (map var name; speed from u,v)
    if var == "speed":
        pz, pu = _parent_column(parent, "u", lon0, lat0, date=date)
        _, pv = _parent_column(parent, "v", lon0, lat0, date=date)
        pv_vals = np.sqrt(pu**2 + pv**2)
    else:
        pz, pv_vals = _parent_column(parent, var, lon0, lat0, date=date)

    fig, ax = plt.subplots(figsize=(4.5, 6))
    ax.plot(cv, cz, "o-", color="C0", ms=3, lw=1.6, label="CROCO")
    ax.plot(pv_vals, pz, "s--", color="C3", ms=3, lw=1.6, label="parent")
    lab = cpro.attrs.get("long_name", var)
    un = cpro.attrs.get("units", "")
    ax.set_xlabel(f"{lab} ({un})" if un else lab)
    ax.set_ylabel("depth (m)")
    ax.grid(alpha=0.3); ax.legend()
    ax.set_title(f"{lab} profile  ({lon0:.2f}, {lat0:.2f})")
    ds.close()
    return _finish(fig, out)

def _interp_1d(depth_grid, z_native, v_native):
    """Linear interpolation of one native (depth, value) profile onto a
    common depth_grid, NaN outside the native range (no extrapolation)."""
    z_native, v_native = np.asarray(z_native), np.asarray(v_native)
    good = np.isfinite(z_native) & np.isfinite(v_native)
    z_native, v_native = z_native[good], v_native[good]
    order = np.argsort(z_native)
    z_native, v_native = z_native[order], v_native[order]
    if z_native.size < 2:
        return np.full(len(depth_grid), np.nan)
    return np.interp(depth_grid, z_native, v_native, left=np.nan, right=np.nan)

def _u2rho(u):
    """CROCO u-point field (eta_rho, xi_rho-1) -> rho-point (eta_rho, xi_rho),
    by averaging neighbours (edges just copied, not extrapolated)."""
    out = np.full((u.shape[0], u.shape[1] + 1), np.nan)
    out[:, 1:-1] = 0.5 * (u[:, :-1] + u[:, 1:])
    out[:, 0] = u[:, 0]
    out[:, -1] = u[:, -1]
    return out


def _v2rho(v):
    """CROCO v-point field (eta_rho-1, xi_rho) -> rho-point (eta_rho, xi_rho),
    by averaging neighbours (edges just copied, not extrapolated)."""
    out = np.full((v.shape[0] + 1, v.shape[1]), np.nan)
    out[1:-1, :] = 0.5 * (v[:-1, :] + v[1:, :])
    out[0, :] = v[0, :]
    out[-1, :] = v[-1, :]
    return out
    
def domain_profile(croco_his, parent, var, day, tindex=-1, depth_grid=None,
                   Yorig=None, daily_mean=True, out=None):
    """Full-domain vertical profile of `var`: CROCO domain-mean +/-
    spatial std vs parent domain-mean +/- spatial std, at each depth
    level, for one day -- companion to compare_profile() (single point),
    same two-line overlay but averaged (and std-banded) over the WHOLE
    domain instead of read at one grid cell.

    var must be 'temp', 'salt', or 'speed'. CROCO's native sigma levels
    and the parent's fixed z levels are interpolated onto a common
    `depth_grid` before combining (default 0 to -500 m, 50 levels).
    """
    if var not in ("temp", "salt", "speed"):
        raise ValueError("domain_profile only supports var='temp', 'salt', or 'speed'")
    if depth_grid is None:
        depth_grid = np.linspace(0, -500, 50)

    ds = pp.open_history(croco_his, Yorig=Yorig)
    ti = tindex
    if daily_mean:
        ds = _maybe_daily_mean(ds, day, daily_mean)
    else:
        ti = _tindex_for_date(ds, day, tindex)
    z3 = pp.depths(ds, tindex=int(ti))                # (s_rho, eta, xi)
    _, _, mask = pp.lonlatmask(ds)                    # (eta, xi)
    if var == "speed":
        nlev = z3.shape[0]
        v3 = np.full_like(z3, np.nan)
        for k in range(nlev):
            uk = ds["u"].isel(time=int(ti), s_rho=k).values
            vk = ds["v"].isel(time=int(ti), s_rho=k).values
            ue, vn = pp.rotate_uv(ds, _u2rho(uk), _v2rho(vk))
            v3[k] = pp.speed(ue, vn)
    else:
        v3 = ds[var].isel(time=int(ti)).values         # (s_rho, eta, xi)
    v3m, z3m = v3, z3
    if mask.shape == v3.shape[1:]:
        v3m = np.where(mask, v3, np.nan)
    if mask.shape == z3.shape[1:]:
        z3m = np.where(mask, z3, np.nan)
    cz = np.nanmean(z3m, axis=(1, 2))
    cmean = np.nanmean(v3m, axis=(1, 2))
    cstd = np.nanstd(v3m, axis=(1, 2))
    ds.close()

    ds_p = xr.open_dataset(parent)
    if var == "speed":
        du, dv = ds_p[PARENT_VARS["u"]], ds_p[PARENT_VARS["v"]]
        if "time" in du.dims:
            du = du.sel(time=np.datetime64(day), method="nearest")
            dv = dv.sel(time=np.datetime64(day), method="nearest")
        da = np.sqrt(du**2 + dv**2)
    else:
        da = ds_p[PARENT_VARS[var]]
        if "time" in da.dims:
            da = da.sel(time=np.datetime64(day), method="nearest")
    pz = -np.abs(ds_p["depth"].values)
    dims = [d for d in da.dims if d != "depth"]
    pmean = da.mean(dim=dims, skipna=True).values
    pstd = da.std(dim=dims, skipna=True).values
    ds_p.close()

    cmean_i = _interp_1d(depth_grid, cz, cmean)
    cstd_i  = _interp_1d(depth_grid, cz, cstd)
    pmean_i = _interp_1d(depth_grid, pz, pmean)
    pstd_i  = _interp_1d(depth_grid, pz, pstd)

    fig, ax = plt.subplots(figsize=(4.5, 6))
    ax.plot(cmean_i, depth_grid, "o-", color="C0", ms=3, lw=1.6, label="CROCO")
    ax.fill_betweenx(depth_grid, cmean_i - cstd_i, cmean_i + cstd_i, color="C0", alpha=0.2)
    ax.plot(pmean_i, depth_grid, "s--", color="C3", ms=3, lw=1.6, label="parent")
    ax.fill_betweenx(depth_grid, pmean_i - pstd_i, pmean_i + pstd_i, color="C3", alpha=0.2)
    unit = {"temp": "degC", "salt": "PSU", "speed": "m s$^{-1}$"}[var]
    ax.set_xlabel(f"{var} ({unit})"); ax.set_ylabel("depth (m)")
    ax.grid(alpha=0.3); ax.legend(fontsize=8)
    ax.set_title(f"{var} full-domain profile -- {day}\nmean +/- spatial std")
    return _save_or_return(fig, out)
    
def compare_timeseries(croco_his, parent, var, lon0, lat0, out=None, Yorig=None,
                       surface_only=True, depth_m=None):
    """Overlay CROCO and parent time series of `var` at a point.

    CROCO is sub-daily; the parent (GLORYS/Mercator) is usually a single or few
    daily fields - so the parent line is drawn at its available time(s), often as
    a flat/reference level. Best when the parent covers the CROCO period.
    """
    ds = pp.open_history(croco_his, Yorig=Yorig)
    cts = pp.timeseries(ds, var, lon0, lat0, surface_only=surface_only, depth_m=depth_m)
    ct = cts["time"].values
    cv = cts.values

    # parent series over the same window, at the SAME depth as CROCO
    dsp = xr.open_dataset(parent)
    def _pick_depth(da):
        if "depth" in da.dims:
            if depth_m is None:
                da = da.isel(depth=0)                        # surface
            else:
                da = da.sel(depth=abs(depth_m), method="nearest")  # nearest parent level
        return da
    cmems = PARENT_VARS["u"] if var == "speed" else PARENT_VARS.get(var, var)
    if var == "speed":
        u = _pick_depth(dsp[PARENT_VARS["u"]]); v = _pick_depth(dsp[PARENT_VARS["v"]])
        u = u.sel(longitude=lon0, latitude=lat0, method="nearest")
        v = v.sel(longitude=lon0, latitude=lat0, method="nearest")
        pv_vals = np.sqrt(u.values**2 + v.values**2)
        pt = u["time"].values if "time" in u.dims else None
    else:
        da = _pick_depth(dsp[cmems])
        da = da.sel(longitude=lon0, latitude=lat0, method="nearest")
        pv_vals = np.atleast_1d(da.values)
        pt = da["time"].values if "time" in da.dims else None
    dsp.close()

    fig, ax = plt.subplots(figsize=(9, 4))
    ax.plot(ct, cv, "o-", color="C0", ms=3, lw=1.4, label="CROCO")
    if pt is not None and np.atleast_1d(pt).size == pv_vals.size:
        ax.plot(pt, pv_vals, "s--", color="C3", ms=4, lw=1.4, label="parent")
    else:
        # single parent field -> horizontal reference line
        ax.axhline(float(np.nanmean(pv_vals)), color="C3", ls="--", lw=1.4, label="parent")
    lab = cts.attrs.get("long_name", var); un = cts.attrs.get("units", "")
    ax.set_xlabel("time"); ax.set_ylabel(f"{lab} ({un})" if un else lab)
    ax.grid(alpha=0.3); ax.legend()
    dm = f"  {depth_m:g} m" if depth_m is not None else "  surface"
    ax.set_title(f"{lab}  ({lon0:.2f}, {lat0:.2f}){dm}")
    fig.autofmt_xdate()
    ds.close()
    return _finish(fig, out)


def compare_section(croco_his, parent, var, lon0, lat0, lon1, lat1, date=None,
                    tindex=-1, npts=120, out=None, Yorig=None,
                    depth_max=None, vmin=None, vmax=None, cmap=None,
                    figsize=(13, 5)):
    """Side-by-side vertical sections of `var`: CROCO vs parent, same transect.

    Two panels sharing a colour scale (CROCO | parent). var: temp/salt/speed.
    """
    ds = pp.open_history(croco_his, Yorig=Yorig)
    csec = pp.section(ds, var, lon0, lat0, lon1, lat1, tindex=tindex, npts=npts)
    cdist = csec["distance_km"].values
    cdepth = csec["depth"].values
    cfield = csec.values

    # parent section: sample its 3D field along the same lon/lat line
    dsp = xr.open_dataset(parent)
    slon = np.linspace(lon0, lon1, npts); slat = np.linspace(lat0, lat1, npts)
    pdepth = -np.abs(dsp["depth"].values)
    if var == "speed":
        U = dsp[PARENT_VARS["u"]]; V = dsp[PARENT_VARS["v"]]
        if "time" in U.dims:
            U = U.sel(time=np.datetime64(date), method="nearest") if date else U.isel(time=0)
            V = V.sel(time=np.datetime64(date), method="nearest") if date else V.isel(time=0)
        psec = np.full((len(pdepth), npts), np.nan)
        for p in range(npts):
            u = U.sel(longitude=slon[p], latitude=slat[p], method="nearest").values
            v = V.sel(longitude=slon[p], latitude=slat[p], method="nearest").values
            psec[:, p] = np.sqrt(u**2 + v**2)
    else:
        DA = dsp[PARENT_VARS.get(var, var)]
        if "time" in DA.dims:
            DA = DA.sel(time=np.datetime64(date), method="nearest") if date else DA.isel(time=0)
        psec = np.full((len(pdepth), npts), np.nan)
        for p in range(npts):
            psec[:, p] = DA.sel(longitude=slon[p], latitude=slat[p], method="nearest").values
    dsp.close()

    dist = _haversine_km(lat0, lon0, slat, slon)
    vmin = np.nanmin([np.nanmin(cfield), np.nanmin(psec)])
    vmax = np.nanmax([np.nanmax(cfield), np.nanmax(psec)])
    lab = csec.attrs.get("long_name", var); un = csec.attrs.get("units", "")

    fig, axes = plt.subplots(1, 2, figsize=figsize, sharey=True,
                             constrained_layout=True)
    # pcolormesh needs finite x/y coords; land points have NaN depth. Fill those
    # depths (forward/backward along each column) so coords are finite, but keep
    # the DATA NaN so land still renders blank.
    cdepth_plot = np.array(cdepth, dtype=float)
    # fill NaNs in the depth coord column-wise using nan-free columns as template
    col_template = np.nanmean(cdepth_plot, axis=1, keepdims=True)
    bad = ~np.isfinite(cdepth_plot)
    cdepth_plot = np.where(bad, np.broadcast_to(col_template, cdepth_plot.shape), cdepth_plot)
    # any remaining NaN (fully-empty rows) -> fill with a monotonic ramp
    if not np.isfinite(cdepth_plot).all():
        ramp = np.linspace(np.nanmin(cdepth_plot), 0, cdepth_plot.shape[0])[:, None]
        cdepth_plot = np.where(np.isfinite(cdepth_plot), cdepth_plot,
                               np.broadcast_to(ramp, cdepth_plot.shape))
    X = np.tile(cdist, (cfield.shape[0], 1))
    h0 = axes[0].pcolormesh(X, cdepth_plot, cfield, vmin=vmin, vmax=vmax,
                            cmap=cmap, shading="auto")
    axes[0].set_title(f"CROCO {lab}"); axes[0].set_xlabel("distance (km)"); axes[0].set_ylabel("depth (m)")
    Xp = np.tile(dist, (psec.shape[0], 1))
    Zp = np.tile(pdepth.reshape(-1, 1), (1, npts))
    h1 = axes[1].pcolormesh(Xp, Zp, psec, vmin=vmin, vmax=vmax,
                            cmap=cmap, shading="auto")
    axes[1].set_title(f"Parent {lab}"); axes[1].set_xlabel("distance (km)")
    if depth_max is not None:
        for a in axes:
            a.set_ylim(-abs(depth_max), 0)
    cb = fig.colorbar(h1, ax=axes, shrink=0.85, pad=0.02)
    cb.set_label(f"{lab} ({un})" if un else lab)
    ds.close()
    # not _finish(): it calls tight_layout(), which cannot be combined with
    # the constrained_layout engine used above
    if out:
        fig.savefig(out, dpi=150, bbox_inches="tight")
        plt.close(fig)
        return out
    return fig


# ======================================================================
# Error vs depth: how CROCO-parent agreement changes with depth
# ======================================================================
def error_vs_depth(croco_his, parent, field="speed", depths=(0, 50, 100, 200, 500),
                   date=None, tindex=-1, Yorig=None, margin_deg=None,
                   metric="rmse", out=None, figsize=(6, 5)):
    """Sweep depths and plot the CROCO-vs-parent error for `field` at each.

    field  : 'sst'/'temp' or 'speed' (ssh is 2D, not depth-resolved).
    depths : iterable of depths in m (0 = surface).
    metric : 'rmse', 'bias', 'crmse', or 'corr'.
    Returns (fig_or_path, {depth: stats}). Shows how model-parent agreement
    varies with depth - typically the error is largest at the surface (wind +
    mesoscale) and decreases with depth (slow, large-scale, geostrophic flow).
    """
    ds = pp.open_history(croco_his, Yorig=Yorig)
    if margin_deg is not None:
        ds = pp.crop_interior(ds, margin_deg=margin_deg)
    ti = _tindex_for_date(ds, date, tindex)

    zz, yy, allstats = [], [], {}
    for d in depths:
        dm = None if d == 0 else d
        croco = _croco_field(ds, field, ti, depth_m=dm)
        par = _parent_field_on_croco(parent, field, date, ds, depth_m=dm)
        s = domain_statistics(croco, par)
        allstats[d] = s
        zz.append(d); yy.append(s[metric])
    ds.close()

    fig, ax = plt.subplots(figsize=figsize)
    ax.plot(yy, zz, "o-", lw=1.8)
    ax.invert_yaxis()                       # depth increases downward
    ax.set_ylabel("depth (m)")
    ax.set_xlabel(f"{metric.upper()}  ({field} CROCO vs parent)")
    ax.grid(alpha=0.3)
    ax.set_title(f"{field} error vs depth")
    return _finish(fig, out), allstats


# ======================================================================
# Resolution comparison: two CROCO runs (e.g. 1/25 child vs 1/12 parent)
# side by side, same field, shared colour scale.
# ======================================================================
def _try_cartopy():
    """Return (crs, cfeature) if cartopy is available, else (None, None)."""
    try:
        import cartopy.crs as ccrs
        import cartopy.feature as cfeature
        return ccrs, cfeature
    except Exception:
        return None, None


def _add_coast(ax, ccrs, cfeature, extent):
    """Overlay the standard map dressing (land/lakes/rivers/borders/coastline/
    labelled gridlines) on a cartopy GeoAxes -- see gtools.plotting.add_map_features."""
    ax.set_extent(extent, crs=ccrs.PlateCarree())
    add_map_features(ax)


def compare_resolution(fine_his, coarse_his, var="temp", depth_m=None, tindex=-1,
                       Yorig=None, out=None, cmap=None, vmin=None, vmax=None,
                       labels=("fine (child)", "coarse (parent)"), figsize=(13, 6),
                       normalized=False, coastline=False, diff=False):
    """Two CROCO runs side by side for the same field, shared colour scale.

    Shows whether the finer run resolves more structure than the coarser one.
    var: 'temp','salt','speed','vort','vort_f','zeta' (or any 3D field).
    depth_m: None -> surface; else a true depth.
    normalized: for vorticity, divide by f (vort/f).
    coastline: overlay a real coastline + land + gridlines (needs cartopy; falls
               back to the plain land-mask rendering if cartopy is unavailable).
    diff: add a third panel = fine - (coarse regridded onto the fine grid), i.e.
          what the finer run adds relative to the coarser one (needs scipy).
    Both runs are plotted on their OWN grids; the optional diff panel puts the
    coarse field onto the fine grid first so the subtraction is well defined.
    """
    df = pp.open_history(fine_his,   Yorig=Yorig)
    dc = pp.open_history(coarse_his, Yorig=Yorig)

    def _field(ds):
        if var in ("vort", "vort_f"):
            da = pp.vorticity(ds, depth_m=depth_m, tindex=tindex,
                              normalized=(normalized or var == "vort_f"))
        elif var == "speed":
            da = pp.speed_map(ds, depth_m=depth_m, tindex=tindex)
        elif var == "zeta":
            lon, lat, mask = pp.lonlatmask(ds)
            z = ds["zeta"].isel(time=tindex).values * mask
            return lon, lat, z
        else:
            da = pp.field(ds, var, depth_m=depth_m, tindex=tindex)
        lon, lat, _ = pp.lonlatmask(ds)
        return lon, lat, np.asarray(da.values)

    lonf, latf, ff = _field(df)
    lonc, latc, fc = _field(dc)

    # shared colour scale from both
    allvals = np.concatenate([ff[np.isfinite(ff)].ravel(), fc[np.isfinite(fc)].ravel()])
    if vmin is None: vmin = np.nanpercentile(allvals, 2)
    if vmax is None: vmax = np.nanpercentile(allvals, 98)
    if cmap is None:
        cmap = "Spectral_r" if var in ("vort", "vort_f", "zeta") else \
               ("viridis" if var in ("speed", "salt") else "RdYlBu_r")
        if var in ("vort", "vort_f", "zeta"):   # symmetric about 0
            m = max(abs(vmin), abs(vmax)); vmin, vmax = -m, m

    # shared geographic extent so panels line up (child box may be smaller)
    lon_min = min(np.nanmin(lonf), np.nanmin(lonc))
    lon_max = max(np.nanmax(lonf), np.nanmax(lonc))
    lat_min = min(np.nanmin(latf), np.nanmin(latc))
    lat_max = max(np.nanmax(latf), np.nanmax(latc))
    extent = [lon_min, lon_max, lat_min, lat_max]

    # optional difference: regrid coarse onto the fine grid, subtract
    diff_field = None
    if diff:
        try:
            from scipy.interpolate import griddata
            pts = np.column_stack([lonc[np.isfinite(fc)].ravel(),
                                    latc[np.isfinite(fc)].ravel()])
            vals = fc[np.isfinite(fc)].ravel()
            fc_on_fine = griddata(pts, vals, (lonf, latf), method="linear")
            diff_field = ff - fc_on_fine
            diff_field[~np.isfinite(ff)] = np.nan     # keep land blank
        except Exception as e:
            print(f"diff disabled (need scipy): {e}")
            diff = False

    ccrs, cfeature = (_try_cartopy() if coastline else (None, None))
    use_cartopy = coastline and ccrs is not None
    if coastline and not use_cartopy:
        print("coastline requested but cartopy not available - using plain land mask")

    dlab = "surface" if depth_m is None else f"{depth_m:g} m"
    lab = {"temp": "SST" if depth_m is None else "temp", "vort_f": "vorticity/f",
           "vort": "vorticity", "speed": "speed", "salt": "salinity",
           "zeta": "SSH"}.get(var, var)

    npan = 3 if diff else 2
    subplot_kw = dict(projection=ccrs.PlateCarree()) if use_cartopy else {}
    fig, axes = plt.subplots(1, npan, figsize=(figsize[0]*npan/2, figsize[1]),
                             sharex=not use_cartopy, sharey=not use_cartopy,
                             constrained_layout=True, subplot_kw=subplot_kw)
    panels = [(lonf, latf, ff, labels[0]), (lonc, latc, fc, labels[1])]

    for ax, (lon, lat, f, l) in zip(axes[:2], panels):
        tf = dict(transform=ccrs.PlateCarree()) if use_cartopy else {}
        h = ax.pcolormesh(lon, lat, f, cmap=cmap, vmin=vmin, vmax=vmax, shading="auto", **tf)
        ax.set_title(l, fontsize=12)
        if use_cartopy:
            _add_coast(ax, ccrs, cfeature, extent)
        else:
            ax.set_xlabel("longitude"); ax.set_xlim(lon_min, lon_max)
            ax.set_ylim(lat_min, lat_max); ax.set_aspect("equal", adjustable="box")
    if not use_cartopy:
        axes[0].set_ylabel("latitude")
    # main colorbar (first two panels)
    cb = fig.colorbar(h, ax=axes[:2], location="bottom" if diff else "right",
                      shrink=0.85, pad=0.03, aspect=30)
    cb.set_label(f"{lab} ({dlab})")

    if diff:
        ax = axes[2]
        dm = np.nanpercentile(np.abs(diff_field[np.isfinite(diff_field)]), 98)
        tf = dict(transform=ccrs.PlateCarree()) if use_cartopy else {}
        hd = ax.pcolormesh(lonf, latf, diff_field, cmap="RdBu_r",
                           vmin=-dm, vmax=dm, shading="auto", **tf)
        ax.set_title(f"{labels[0]} − {labels[1]}", fontsize=12)
        if use_cartopy:
            _add_coast(ax, ccrs, cfeature, extent)
        else:
            ax.set_xlabel("longitude"); ax.set_xlim(lon_min, lon_max)
            ax.set_ylim(lat_min, lat_max); ax.set_aspect("equal", adjustable="box")
        cbd = fig.colorbar(hd, ax=ax, location="bottom", shrink=0.85, pad=0.03, aspect=20)
        cbd.set_label(f"Δ{lab}")

    fig.suptitle(f"{lab}  ({dlab})  —  resolution comparison", fontsize=13)
    df.close(); dc.close()
    if out:
        fig.savefig(out, dpi=150, bbox_inches="tight")
        plt.close(fig)
        return out
    return fig


def compare_profile_resolution(fine_his, coarse_his, var, lon0, lat0, tindex=-1,
                               Yorig=None, out=None,
                               labels=("fine (child)", "coarse (parent)")):
    """Overlay vertical profiles from TWO CROCO runs (child vs parent) at a point.

    Unlike compare_profile (CROCO vs Mercator), BOTH inputs are CROCO history
    files, read with pp.profile. Use for nesting: does the finer run resolve a
    different vertical structure at a point?
    var: 'temp','salt','u','v','speed'.
    """
    df = pp.open_history(fine_his,   Yorig=Yorig)
    dc = pp.open_history(coarse_his, Yorig=Yorig)
    pf = pp.profile(df, var, lon0, lat0, tindex=tindex)
    pc = pp.profile(dc, var, lon0, lat0, tindex=tindex)

    fig, ax = plt.subplots(figsize=(4.5, 6))
    ax.plot(pf.values, pf["depth"].values, "o-", color="C0", ms=3, lw=1.6, label=labels[0])
    ax.plot(pc.values, pc["depth"].values, "s--", color="C3", ms=3, lw=1.6, label=labels[1])
    lab = pf.attrs.get("long_name", var); un = pf.attrs.get("units", "")
    ax.set_xlabel(f"{lab} ({un})" if un else lab)
    ax.set_ylabel("depth (m)")
    ax.grid(alpha=0.3); ax.legend()
    ax.set_title(f"{lab} profile  ({lon0:.2f}, {lat0:.2f})")
    df.close(); dc.close()
    return _finish(fig, out)


def compare_section_resolution(fine_his, coarse_his, var, lon0, lat0, lon1, lat1,
                               tindex=-1, Yorig=None, out=None, npts=200,
                               labels=("fine (child)", "coarse (parent)")):
    """Two CROCO vertical sections (child vs parent) along the same transect.

    Both inputs are CROCO history files (read with pp.section), so this works for
    nesting - unlike compare_section, which expects a Mercator parent.
    var: 'temp','salt','u','v','speed','vort'.
    """
    df = pp.open_history(fine_his,   Yorig=Yorig)
    dc = pp.open_history(coarse_his, Yorig=Yorig)
    sf = pp.section(df, var, lon0, lat0, lon1, lat1, tindex=tindex, npts=npts)
    sc = pp.section(dc, var, lon0, lat0, lon1, lat1, tindex=tindex, npts=npts)

    def _prep(sec):
        field = np.asarray(sec.values)                      # (s_rho, npts)
        dist = sec["distance_km"].values
        depth = np.array(sec["depth"].values, dtype=float)  # (s_rho, npts), NaN on land
        # fill NaN depth coords column-wise so pcolormesh has finite coords,
        # but keep the DATA NaN so land renders blank
        col = np.nanmean(depth, axis=1, keepdims=True)
        bad = ~np.isfinite(depth)
        depth = np.where(bad, np.broadcast_to(col, depth.shape), depth)
        if not np.isfinite(depth).all():
            ramp = np.linspace(np.nanmin(depth), 0, depth.shape[0])[:, None]
            depth = np.where(np.isfinite(depth), depth, np.broadcast_to(ramp, depth.shape))
        X = np.tile(dist, (field.shape[0], 1))
        return X, depth, field

    Xf, Zf, Ff = _prep(sf)
    Xc, Zc, Fc = _prep(sc)
    vmin = np.nanmin([np.nanmin(Ff), np.nanmin(Fc)])
    vmax = np.nanmax([np.nanmax(Ff), np.nanmax(Fc)])
    cmap = "viridis" if var in ("speed", "salt") else \
       ("RdBu_r" if var in ("u", "v", "vort", "vort_f") else
        ("Spectral_r" if var == "ssh" else "RdYlBu_r"))
    if var in ("u", "v", "vort", "vort_f"):
        m = max(abs(vmin), abs(vmax)); vmin, vmax = -m, m
    lab = sf.attrs.get("long_name", var); un = sf.attrs.get("units", "")

    fig, axes = plt.subplots(1, 2, figsize=(13, 5), sharey=True, constrained_layout=True)
    for ax, (X, Z, F, l) in zip(axes, [(Xf, Zf, Ff, labels[0]), (Xc, Zc, Fc, labels[1])]):
        h = ax.pcolormesh(X, Z, F, cmap=cmap, vmin=vmin, vmax=vmax, shading="auto")
        ax.set_title(l); ax.set_xlabel("distance (km)")
    axes[0].set_ylabel("depth (m)")
    fig.suptitle(f"{lab} section  ({lon0:.1f},{lat0:.1f})\u2192({lon1:.1f},{lat1:.1f})"
                 f"  \u2014 resolution comparison", fontsize=12)
    cb = fig.colorbar(h, ax=axes, location="right", shrink=0.85, pad=0.02, aspect=30)
    cb.set_label(f"{lab} ({un})" if un else lab)
    df.close(); dc.close()
    if out:
        fig.savefig(out, dpi=150, bbox_inches="tight"); plt.close(fig); return out
    return fig
"""
gtools/validation_godae.py — GODAE OceanView-style validation metrics.

Extends gtools.validation (which does CROCO-vs-parent map comparisons) with
the standard *point/profile* validation protocol used in GODAE OceanView /
CMEMS class-4 intercomparison exercises:

    - CROCO  vs GLORYS12V1 reanalysis          (model-vs-model, "class 1/2")
    - CROCO  vs in-situ observations           (model-vs-obs,  "class 4")
    - GLORYS12V1 vs in-situ observations       (reference skill, for context)

Reference product for in-situ data
-----------------------------------
CMEMS "Global Ocean- In-Situ Observations" (aka INSITU_GLO_PHY),
DOI: https://doi.org/10.48670/moi-00036
(product page: https://data.marine.copernicus.eu/product/INSITU_GLO_PHYBGCWAV_DISCRETE_MYNRT_013_030/description)
Downloaded e.g. via `copernicusmarine` (already a gtools dependency, see
gtools/download/cmems.py) as one NetCDF file per platform/profile, in the
CMEMS in-situ TAC vocabulary (TEMP, PSAL, DEPH, LATITUDE, LONGITUDE, TIME,
plus *_QC quality-flag companions, 1=good).

GODAE-style metrics
--------------------
For a pair of collocated series (model, reference) the standard GODAE
intercomparison scorecard is:

    BIAS   = mean(model - ref)                        (mean error / ME)
    RMSD   = sqrt(mean((model-ref)^2))                 (root-mean-sq diff.)
    URMSD  = sqrt(mean(((model-mean(model))
                        -(ref-mean(ref)))^2))          (unbiased/centred RMSD)
    CORR   = Pearson correlation coefficient
    SI     = URMSD / mean(|ref|)                       (scatter index, %)
    SI_STD = URMSD / std(ref)                          (scatter index, % —
                                                         alternative normalisation,
                                                         more stable for near-
                                                         zero-mean fields like
                                                         SSH anomaly)
    STD_ratio = std(model) / std(ref)                  (for Taylor diagrams)

These are computed globally and, for in-situ profiles, additionally
per standard GODAE depth layer (surface / thermocline / deep), which is how
GODAE OceanView / CMEMS QuID documents typically report skill.

Typical use
-----------
    import gtools.validation_godae as vg

    # 1) CROCO vs GLORYS12V1 (model-vs-model) at the surface, whole domain
    stats = vg.godae_metrics_grid(model_field, ref_field)
    vg.print_scorecard("CROCO vs GLORYS12V1 SST", stats)

    # 2) CROCO (+ Copernicus Marine Forecast for context) vs in-situ profiles
    #    for one forecast cycle
    report = vg.validate_against_insitu(
        croco_his="model-runs/Canary_12/20260724/fcst/CROCO_FILES/croco_his.nc",
        insitu_files=glob.glob("model-runs/Canary_12/20260724/downloaded_data/INSITU/2026-07-2*.nc"),
        glorys_file="model-runs/Canary_12/20260724/downloaded_data/MERCATOR/MERCATOR_20260724_00.nc",
        variables=("temp", "salt"), Yorig=2000)
    vg.print_scorecard_table(report)
    vg.taylor_diagram(report, var="temp", out="taylor_temp.png")
"""


# ======================================================================================
# ======================================================================================
# GODAE OceanView-style validation metrics (formerly gtools/validation_godae.py,
# merged in here so every validation function lives in this one module).
#
# Extends the CROCO-vs-parent MAP comparisons above with the standard
# *point/profile* validation protocol used in GODAE OceanView / CMEMS
# class-4 intercomparison exercises:
#     - CROCO       vs GLORYS12V1/Mercator reanalysis  (model-vs-model, "class 1/2")
#     - CROCO       vs in-situ observations             (model-vs-obs,  "class 4")
#     - GLORYS12V1  vs in-situ observations              (reference skill, for context)
# See GODAE_LAYERS / godae_metrics / validate_against_insitu / taylor_diagram below.
# ======================================================================================
# ======================================================================================


# ======================================================================
# GODAE standard depth layers (m, positive down) — surface/thermocline/deep
# mirrors the layering used in GODAE OceanView / CMEMS class-4 reports.
# ======================================================================
GODAE_LAYERS = {
    "surface (0-10 m)":       (0, 10),
    "mixed layer (10-50 m)":  (10, 50),
    "thermocline (50-300 m)": (50, 300),
    "deep (300-2000 m)":      (300, 2000),
}

# In-situ TAC variable -> our short name, and the model/parent equivalent
INSITU_VARS = {"temp": "TEMP", "salt": "PSAL"}


# ======================================================================
# 1) Core GODAE metrics (works on any two arrays of collocated points)
# ======================================================================
def godae_metrics(model, ref):
    """Compute the standard GODAE OceanView scorecard for collocated 1D arrays.

    model, ref : array-like, same length, collocated in space/time/depth.
    Only points where both are finite are used (NaNs dropped pairwise).

    Returns a dict: n, bias, rmsd, urmsd, corr, si, std_model, std_ref,
    std_ratio, mean_model, mean_ref.
    """
    m = np.asarray(model, dtype=float).ravel()
    r = np.asarray(ref, dtype=float).ravel()
    ok = np.isfinite(m) & np.isfinite(r)
    m, r = m[ok], r[ok]
    if m.size < 2:
        keys = ("n", "bias", "rmsd", "urmsd", "corr", "si", "si_std",
                "std_model", "std_ref", "std_ratio", "mean_model", "mean_ref")
        return {k: (0 if k == "n" else np.nan) for k in keys}

    bias = float(np.mean(m - r))
    rmsd = float(np.sqrt(np.mean((m - r) ** 2)))
    ma, ra = m - m.mean(), r - r.mean()
    urmsd = float(np.sqrt(np.mean((ma - ra) ** 2)))
    std_m, std_r = float(np.std(m)), float(np.std(r))
    corr = float(np.corrcoef(m, r)[0, 1]) if std_m > 1e-12 and std_r > 1e-12 else np.nan
    ref_mean_abs = float(np.mean(np.abs(r)))
    si = float(urmsd / ref_mean_abs * 100.0) if ref_mean_abs > 1e-12 else np.nan
    # si_std: scatter index normalised by std(ref) instead of mean(|ref|).
    # The mean-normalised SI above is the more common GODAE/CMEMS convention,
    # but it is unstable/misleading for near-zero-mean signals (e.g. SSH
    # anomaly, velocity components) where mean(|ref|) can be tiny relative to
    # the field's actual variability — si_std stays well-behaved in that case
    # and is the standard fallback used in the oceanographic literature.
    si_std = float(urmsd / std_r * 100.0) if std_r > 1e-12 else np.nan
    std_ratio = float(std_m / std_r) if std_r > 1e-12 else np.nan

    return {
        "n": int(m.size), "bias": bias, "rmsd": rmsd, "urmsd": urmsd,
        "corr": corr, "si": si, "si_std": si_std, "std_model": std_m, "std_ref": std_r,
        "std_ratio": std_ratio, "mean_model": float(m.mean()), "mean_ref": float(r.mean()),
    }


def godae_metrics_grid(model_field, ref_field):
    """Convenience wrapper: godae_metrics() on two 2D/3D gridded fields."""
    return godae_metrics(model_field, ref_field)


def _time_index_for(ds, date):
    """Nearest CROCO time index for a 'YYYY-MM-DD'-style date, or -1 (last)."""
    if date is None:
        return -1
    times = np.asarray(ds["time"].values)
    target = np.datetime64(date)
    return int(np.argmin(np.abs(times - target)))


def godae_scorecard_croco_vs_glorys(croco_his, parent_file, var, date=None,
                                     depth_m=None, Yorig=None, ssh_anomaly=True):
    """Grid-level GODAE scorecard: CROCO vs its parent GLORYS12V1 (or Mercator
    analysis/forecast) field, on the CROCO grid.

    This is the "class 1/2" (model-vs-model) half of the GODAE OceanView
    protocol: reuses gtools.validation.load_parent + regrid_to_croco for the
    I/O and bilinear regridding (so this stays consistent with gtools'
    existing CROCO-vs-GLORYS maps), then scores the pair with the full GODAE
    metric set (bias/RMSD/uRMSD/corr/SI/std-ratio) instead of just bias/RMSE.

    var      : one of 'temp','ssh','salt','speed' (see validation.PARENT_VARS
               for temp/ssh/salt). NOTE: raw 'u'/'v' are deliberately NOT
               accepted here — CROCO's u,v live on the staggered u/v-point
               grids, not the rho-point grid that regrid_to_croco regrids
               onto, so comparing them directly against the parent's earth-
               relative uo/vo would silently compare mismatched grids /
               conventions. Use 'speed' instead (matches how gtools.
               validation.compare_currents already handles this: CROCO u,v
               are interpolated to rho points and rotated to east/north via
               pp.surface_uv/pp.rotate_uv first, exactly like here).
    date     : 'YYYY-MM-DD' (nearest time on both sides); None -> last CROCO
               step / first parent step (matches load_parent's default)
    depth_m  : None -> surface; else nearest/interpolated depth (m).
               Ignored for var='ssh' (2D field, no vertical dimension).
    ssh_anomaly : only used when var='ssh' (default True). CROCO's zeta has
               no absolute geoid reference, so comparing raw SSH levels
               against the parent's absolute SSH would score an arbitrary
               offset between the two reference levels as spurious 'bias' --
               masking the actual pattern skill this function is meant to
               measure. With ssh_anomaly=True (recommended), the domain
               mean is removed from both fields before scoring, so 'bias'
               reflects a genuine pattern-level offset, not the reference
               mismatch. Set False only if you have already reconciled the
               two reference levels upstream and want the raw comparison.
    Returns a dict from godae_metrics(), plus 'variable', 'date', 'depth_m'.
    """
    if var not in ("temp", "ssh", "salt", "speed"):
        raise ValueError(
            f"var={var!r} not supported here; use 'temp','ssh','salt' or "
            f"'speed' (not raw 'u'/'v' — see this function's docstring)")
    ds = pp.open_history(croco_his, Yorig=Yorig)
    ti = _time_index_for(ds, date)

    if var == "speed":
        if depth_m is None:
            ur, vr = pp.surface_uv(ds, tindex=ti)
            ue, vn = pp.rotate_uv(ds, ur, vr)
        else:
            ue, vn = pp.uv_at_depth(ds, depth_m, tindex=ti, rotate=True)
        model = pp.speed(np.squeeze(ue), np.squeeze(vn))

        plon, plat, pu = load_parent(parent_file, "u", date=date, depth_m=depth_m)
        _, _, pv = load_parent(parent_file, "v", date=date, depth_m=depth_m)
        pu_c = regrid_to_croco(plon, plat, pu, ds, method="linear")
        pv_c = regrid_to_croco(plon, plat, pv, ds, method="linear")
        ref_on_croco = pp.speed(pu_c, pv_c)
    else:
        croco_var = "zeta" if var == "ssh" else var
        d_m = None if var == "ssh" else depth_m
        model_da = pp.field(ds, croco_var, depth_m=d_m, tindex=ti)
        model = np.asarray(model_da.values)

        plon, plat, pfield = load_parent(parent_file, var, date=date, depth_m=depth_m)
        ref_on_croco = regrid_to_croco(plon, plat, pfield, ds, method="linear")

        if var == "ssh" and ssh_anomaly:
            model = model - np.nanmean(model)
            ref_on_croco = ref_on_croco - np.nanmean(ref_on_croco)

    ds.close()
    s = godae_metrics(model, ref_on_croco)
    s.update({"variable": var, "date": date, "depth_m": depth_m})
    return s


def print_scorecard(title, s):
    """Pretty-print one GODAE scorecard dict (see godae_metrics)."""
    print(f"[{title}]")
    print(f"  n={s['n']}  bias={s['bias']:+.3f}  RMSD={s['rmsd']:.3f}  "
          f"uRMSD={s['urmsd']:.3f}  corr={s['corr']:.3f}  "
          f"SI={s['si']:.1f}%  SI_std={s['si_std']:.1f}%  std_ratio={s['std_ratio']:.3f}")


def print_scorecard_table(report):
    """Pretty-print the DataFrame returned by validate_against_insitu()."""
    if report.empty:
        print("No collocated in-situ data found.")
        return
    cols = ["variable", "vs", "layer", "n", "bias", "rmsd", "urmsd", "corr", "si", "si_std"]
    with pd.option_context("display.float_format", "{:.3f}".format,
                            "display.width", 120):
        print(report[cols].to_string(index=False))


# ======================================================================
# 2) CMEMS in-situ TAC reader (doi.org/10.48670/moi-00036)
# ======================================================================
def _qc_to_int(q):
    """Robustly parse one CMEMS QC flag value to an int, or -1 if missing.

    CMEMS in-situ TAC files store QC flags inconsistently across products:
    single-byte strings (b'1'), numpy bytes_ scalars, plain str, or numeric
    dtypes — and missing/fill values show up as '', ' ', 'nan', etc. Decoding
    bytes -> str first (rather than calling int() directly on the raw value)
    avoids a crash on blank/whitespace flags, which int() cannot parse.
    """
    if isinstance(q, (bytes, np.bytes_)):
        q = q.decode("utf-8", errors="ignore")
    s = str(q).strip()
    if s in ("", "nan", "None", "--"):
        return -1
    try:
        return int(float(s))
    except ValueError:
        return -1


def load_insitu(fname, var="temp", qc_flags=(1,)):
    """Read one CMEMS in-situ TAC profile/timeseries file (moi-00036 product).

    Returns a tidy DataFrame with columns: time, lon, lat, depth, value —
    good-quality points only (QC flag in `qc_flags`, CMEMS convention: 1=good).
    Handles both single-profile files (DEPTH dim) and multi-profile/TS files
    (N_PROF x N_LEVELS), which are the two common in-situ TAC layouts.
    """
    cmems_var = INSITU_VARS[var]
    ds = xr.open_dataset(fname)
    if cmems_var not in ds:
        ds.close()
        return pd.DataFrame(columns=["time", "lon", "lat", "depth", "value"])

    val = ds[cmems_var]
    qc_name = f"{cmems_var}_QC"
    depth_name = "DEPH" if "DEPH" in ds else ("DEPTH" if "DEPTH" in ds else None)
    if depth_name is None:
        import warnings
        warnings.warn(
            f"{fname}: no DEPH/DEPTH variable found — assuming depth=0 m for "
            f"every point in this file. Verify this is actually a surface-only "
            f"product (e.g. drifting buoy SST) before trusting the 'surface' "
            f"layer results for it.", stacklevel=2)

    # Broadcast everything to the shape of `val` so we can flatten consistently,
    # whether the file is a single profile (1D depth) or N_PROF x N_LEVELS.
    lat = ds["LATITUDE"].broadcast_like(val)
    lon = ds["LONGITUDE"].broadcast_like(val)
    time = ds["TIME"].broadcast_like(val)
    depth = ds[depth_name].broadcast_like(val) if depth_name else xr.zeros_like(val)

    df = pd.DataFrame({
        "time": np.asarray(time.values).ravel(),
        "lon": np.asarray(lon.values).ravel(),
        "lat": np.asarray(lat.values).ravel(),
        "depth": np.asarray(depth.values).ravel(),
        "value": np.asarray(val.values).ravel(),
    })

    if qc_name in ds:
        qc_raw = np.asarray(ds[qc_name].broadcast_like(val).values).ravel()
        qc = np.array([_qc_to_int(q) for q in qc_raw])
        df = df[np.isin(qc, qc_flags)]

    ds.close()
    df = df.dropna(subset=["value", "lon", "lat", "depth", "time"])
    return df.reset_index(drop=True)


def load_insitu_many(files, var="temp", qc_flags=(1,)):
    """load_insitu() over a list/glob of files, concatenated into one DataFrame."""
    frames = [load_insitu(f, var=var, qc_flags=qc_flags) for f in files]
    frames = [f for f in frames if not f.empty]
    if not frames:
        return pd.DataFrame(columns=["time", "lon", "lat", "depth", "value"])
    return pd.concat(frames, ignore_index=True)


# ======================================================================
# 3) Collocation: sample CROCO / GLORYS at in-situ profile locations
# ======================================================================
def _nearest_time_index(model_times, obs_times):
    """Vectorized nearest-neighbour time index for every obs time.

    model_times, obs_times : arrays of np.datetime64 (any resolution).
    Returns (best_idx, best_dt_hours), both arrays aligned with obs_times.
    Checks both neighbours returned by searchsorted (the value just above
    AND just below each obs time) and picks whichever is actually closer —
    searchsorted alone only gives the insertion point, not the nearer side.
    """
    times_np = np.asarray(model_times)
    obs_np = np.asarray(obs_times)
    n = len(times_np)
    idx = np.clip(np.searchsorted(times_np, obs_np), 0, n - 1)
    idx_left = np.clip(idx - 1, 0, n - 1)
    d_r = np.abs((times_np[idx] - obs_np) / np.timedelta64(1, "h")).astype(float)
    d_l = np.abs((times_np[idx_left] - obs_np) / np.timedelta64(1, "h")).astype(float)
    use_left = d_l < d_r
    best_idx = np.where(use_left, idx_left, idx)
    best_dt = np.where(use_left, d_l, d_r)
    return best_idx, best_dt


def _collocate_model(ds, var, obs, time_tol_hours=6.0):
    """Sample a CROCO history Dataset at each obs (time, lon, lat, depth).

    Nearest neighbour in time (within time_tol_hours, else dropped) and
    horizontal space (nearest rho point; skipped if that point is land);
    linear interpolation in the vertical, using each obs row's own model
    column (correctly rejected as NaN if the obs depth is outside the range
    the model column actually spans, e.g. obs is below the local seafloor).
    Returns an array aligned with `obs` rows (NaN where no match).

    Reads each required CROCO time step's full 3D array only ONCE (not once
    per obs row sharing that time step), which matters for any dataset
    with more than a handful of profiles.
    """
    if obs.empty:
        return np.array([])
    if not np.issubdtype(ds["time"].dtype, np.datetime64):
        raise ValueError(
            f"ds['time'] was not decoded to real dates (dtype={ds['time'].dtype}); "
            f"open the file with pp.open_history(..., Yorig=<reference year>) first")
    model_times = pd.to_datetime(np.atleast_1d(ds["time"].values)).values
    obs_times = pd.to_datetime(obs["time"].values).values
    best_idx, best_dt = _nearest_time_index(model_times, obs_times)

    _, _, mask = pp.lonlatmask(ds)   # 2D, 1=ocean / NaN=land
    out = np.full(len(obs), np.nan)

    for ti in np.unique(best_idx):
        rows = np.where((best_idx == ti) & (best_dt <= time_tol_hours))[0]
        if rows.size == 0:
            continue
        z = pp.depths(ds, tindex=int(ti))              # (s_rho, eta, xi), m, negative down
        v3 = ds[var].isel(time=int(ti)).values          # (s_rho, eta, xi)
        for r in rows:
            row = obs.iloc[r]
            j, i = pp.nearest_index(ds, row["lon"], row["lat"])
            if not np.isfinite(mask[j, i]):
                continue                                 # nearest model point is land -> no match
            zc, vc = z[:, j, i], v3[:, j, i]
            good = np.isfinite(zc) & np.isfinite(vc)
            if good.sum() < 2:
                continue
            depth_pos = -zc[good]                        # positive-down, matches obs['depth']
            order = np.argsort(depth_pos)
            depth_sorted, val_sorted = depth_pos[order], vc[good][order]
            if row["depth"] < depth_sorted[0] or row["depth"] > depth_sorted[-1]:
                continue                                 # obs depth outside this column's range
            out[r] = np.interp(row["depth"], depth_sorted, val_sorted)
    return out


def _collocate_glorys(fname, var, obs):
    """Sample a GLORYS12V1 (or Mercator) file at each obs (time, lon, lat, depth).

    Nearest time step + nearest horizontal grid point + linear interpolation
    in depth using the file's own depth levels. Uses validation.PARENT_VARS
    for the CMEMS variable-name mapping, so both validation modules always
    agree on e.g. 'temp' -> 'thetao'.
    """
    if var not in PARENT_VARS:
        raise KeyError(f"'{var}' has no GLORYS/Mercator mapping; see validation.PARENT_VARS")
    cmems = PARENT_VARS[var]
    ds = xr.open_dataset(fname)
    out = np.full(len(obs), np.nan)
    if cmems not in ds or obs.empty:
        ds.close()
        return out

    da = ds[cmems]
    lon = ds["longitude"].values
    lat = ds["latitude"].values
    depth = ds["depth"].values if "depth" in ds else np.array([0.0])
    has_time = "time" in da.dims
    has_depth = "depth" in da.dims

    if has_time:
        model_times = np.asarray(pd.to_datetime(ds["time"].values).values)
        obs_times = np.asarray(pd.to_datetime(obs["time"].values).values)
        best_idx, _ = _nearest_time_index(model_times, obs_times)
    else:
        best_idx = np.zeros(len(obs), dtype=int)

    # Load the whole (small, regional-subset) field into memory ONCE, with
    # an explicit dimension order -- da.isel(time=, longitude=, latitude=)
    # is dimension-NAME based (order-agnostic) but plain numpy indexing
    # below isn't, so pin the order here rather than assume the file's
    # native layout.
    # Without loading up front, every per-observation extraction re-reads
    # from disk individually -- fine for a handful of points, but on a slow
    # filesystem (network/WSL-mounted drives especially) hundreds of
    # observations means hundreds of separate small disk reads, which is
    # the dominant cost here, not the interpolation itself.
    dim_order = (["time"] if has_time else []) + (["depth"] if has_depth else []) + ["latitude", "longitude"]
    da = da.transpose(*dim_order).load()
    field = da.values

    for pos, (idx, row) in enumerate(obs.iterrows()):
        try:
            xi = int(np.argmin(np.abs(lon - row["lon"])))
            yi = int(np.argmin(np.abs(lat - row["lat"])))
            if has_time and has_depth:
                prof = field[int(best_idx[pos]), :, yi, xi]
            elif has_time:
                prof = field[int(best_idx[pos]), yi, xi]
            elif has_depth:
                prof = field[:, yi, xi]
            else:
                prof = field[yi, xi]
            prof = np.atleast_1d(prof)
            good = np.isfinite(prof)
            if good.sum() < 2:
                continue
            d_good, p_good = depth[good], prof[good]
            order = np.argsort(d_good)
            d_sorted, p_sorted = d_good[order], p_good[order]
            if row["depth"] < d_sorted[0] or row["depth"] > d_sorted[-1]:
                continue
            out[pos] = np.interp(row["depth"], d_sorted, p_sorted)
        except Exception:
            continue
    ds.close()
    return out


def _layer_of(depth_m):
    """GODAE depth-layer label for one obs depth (m, positive down).

    Clamps small negative depths to 0 before bucketing: real in-situ sensors
    (moored buoys, drifters) commonly report slightly negative depth from
    wave motion around a nominally-surface instrument (e.g. -0.3 m) — without
    clamping, `top <= depth_m` fails for every GODAE_LAYERS bucket (all have
    top>=0) and these clearly-surface readings fell through to "below 2000 m",
    which is wrong. Depths below -1 m are left unclamped and simply won't
    match any layer (falls to "below 2000 m" as a catch-all) since they more
    likely indicate a genuine unit/sign error in that observation.
    """
    d = 0.0 if -1.0 <= depth_m < 0.0 else depth_m
    for name, (top, bot) in GODAE_LAYERS.items():
        if top <= d < bot:
            return name
    return "below 2000 m"


# ======================================================================
# 4) Full GODAE-style validation report: CROCO (+ GLORYS) vs in-situ
# ======================================================================
def validate_against_insitu(croco_his, insitu_files, glorys_file=None,
                             variables=("temp", "salt"), Yorig=None,
                             time_tol_hours=6.0, qc_flags=(1,)):
    """Run the full GODAE-style scorecard: CROCO vs in-situ, and (optionally)
    GLORYS12V1 vs in-situ for reference, per variable and per standard depth
    layer (see GODAE_LAYERS).

    croco_his    : path (or glob pattern) to CROCO history file(s)
    insitu_files : list/glob of CMEMS in-situ TAC files (moi-00036 product)
    glorys_file  : optional GLORYS12V1 file covering the same period/area,
                   included in the report as a second "vs" column for context
                   (this is the GODAE OceanView practice of always quoting
                   the driving reanalysis' own skill alongside the downscaled
                   model's skill).
    Returns a tidy DataFrame: variable, vs, layer, n, bias, rmsd, urmsd, corr, si.
    """
    if isinstance(insitu_files, str):
        insitu_files = sorted(glob.glob(insitu_files))
    ds = pp.open_history(croco_his, Yorig=Yorig)
    if not np.issubdtype(ds["time"].dtype, np.datetime64):
        raise ValueError(
            f"{croco_his}: 'time' was not decoded to real dates (dtype="
            f"{ds['time'].dtype}). This usually means the file has no CF "
            f"'units' attribute on time and no Yorig was given -- pass "
            f"Yorig=<reference year> explicitly. Continuing without this "
            f"check would silently misinterpret raw time values as "
            f"nanoseconds-since-1970 and reject every observation as "
            f"outside time_tol_hours, with no indication why.")

    rows = []
    for var in variables:
        obs = load_insitu_many(insitu_files, var=var, qc_flags=qc_flags)
        if obs.empty:
            print(f"  no in-situ '{var}' obs found in the given files")
            continue
        obs["layer"] = obs["depth"].apply(_layer_of)

        model_vals = _collocate_model(ds, var, obs, time_tol_hours=time_tol_hours)
        sources = {"CROCO": model_vals}
        if glorys_file:
            sources["GLORYS12V1"] = _collocate_glorys(glorys_file, var, obs)

        for src_name, vals in sources.items():
            obs_v = obs.assign(model=vals)
            # global (all depths) row
            s = godae_metrics(obs_v["model"], obs_v["value"])
            rows.append({"variable": var, "vs": src_name, "layer": "all", **s})
            # per GODAE depth-layer rows
            for layer, grp in obs_v.groupby("layer"):
                s = godae_metrics(grp["model"], grp["value"])
                rows.append({"variable": var, "vs": src_name, "layer": layer, **s})

    ds.close()
    return pd.DataFrame(rows)


# ======================================================================
# 5) Taylor diagram (standard GODAE OceanView skill-summary plot)
# ======================================================================
def taylor_diagram(report, var, layer="all", out=None):
    """Taylor diagram (std-ratio vs correlation) summarising the scorecard
    for one variable, one depth layer, across every 'vs' source in `report`
    (as produced by validate_against_insitu). Reference = std_ratio=1, corr=1.
    """
    sub = report[(report["variable"] == var) & (report["layer"] == layer)]
    if sub.empty:
        raise ValueError(f"no rows for variable={var!r}, layer={layer!r}")

    # Standard Taylor diagrams only need the [0,90] correlation quadrant for
    # the near-1 correlations typical of ocean temp/salt validation, but a
    # negative correlation would fall at theta>90deg — with thetamax fixed at
    # 90, such a point would silently render just outside/on the sector edge
    # with no indication it's off-scale. Auto-extend to the full [0,180]
    # (correlation -1..1) quadrant whenever a negative correlation is present.
    has_negative_corr = bool((sub["corr"] < 0).any())
    thetamax = 180 if has_negative_corr else 90
    corr_ticks = ([-1.0, -0.5, 0, 0.5, 0.8, 0.9, 0.95, 0.99, 1.0] if has_negative_corr
                  else [0, 0.2, 0.4, 0.6, 0.8, 0.9, 0.95, 0.99, 1.0])

    fig = plt.figure(figsize=(6, 6))
    ax = fig.add_subplot(111, polar=True)
    ax.set_thetamin(0); ax.set_thetamax(thetamax)
    # angle = arccos(correlation); radius = std_ratio
    ax.set_xticks(np.arccos(corr_ticks))
    ax.set_xticklabels([str(c) for c in corr_ticks])
    ax.set_ylabel("normalised std dev (model / ref)")
    finite_ratios = sub["std_ratio"][np.isfinite(sub["std_ratio"])]
    ymax = max(1.6, finite_ratios.max() * 1.2) if len(finite_ratios) else 1.6
    ax.set_ylim(0, ymax)

    # reference point
    ax.plot(0, 1, "k*", ms=16, label="reference")

    for _, row in sub.iterrows():
        if not np.isfinite(row["corr"]) or not np.isfinite(row["std_ratio"]):
            continue
        theta = np.arccos(np.clip(row["corr"], -1, 1))
        ax.plot(theta, row["std_ratio"], "o", ms=10,
                label=f"{row['vs']} (n={row['n']}, RMSD={row['rmsd']:.2f})")

    ax.set_title(f"Taylor diagram — {var}, {layer}", fontsize=12, pad=20)
    ax.legend(loc="upper left", bbox_to_anchor=(1.02, 1.0), fontsize=8)
    fig.tight_layout()
    if out:
        fig.savefig(out, dpi=150, bbox_inches="tight")
        plt.close(fig)
        return out
    return fig


# ======================================================================
# 6) In-situ TRAJECTORY view (drifters / floats / gliders — surface or
#    near-surface platforms whose obs form a track over time)
# ======================================================================
def plot_insitu_trajectory(obs, croco_extent=None, var="temp", out=None, figsize=(8, 7)):
    """Map the in-situ observation locations (as loaded by load_insitu_many),
    coloured by value, so drifting-platform TRACKS (drifters, gliders,
    surface floats) are visible over the CROCO domain — the complement to
    the profile/depth-layer scorecard in validate_against_insitu(), which
    collapses the spatial pattern away.

    obs           : DataFrame from load_insitu()/load_insitu_many() (needs
                    'lon','lat','value'; a 'platform'/'id' column, if
                    present, draws one line per platform instead of just dots).
    croco_extent  : optional (lon_min, lon_max, lat_min, lat_max) to frame the
                    domain (e.g. from `gtools.postprocess.extent(ds)`).
    """
    if obs.empty:
        print("plot_insitu_trajectory: no observations to plot")
        return None
    proj = {"projection": ccrs.PlateCarree()} if _HAS_CARTOPY else {}
    fig, ax = plt.subplots(figsize=figsize, subplot_kw=proj)

    platform_col = next((c for c in ("platform", "platform_code", "id", "wmo")
                         if c in obs.columns), None)
    tf = dict(transform=ccrs.PlateCarree()) if _HAS_CARTOPY else {}
    if platform_col:
        for pid, grp in obs.sort_values("time").groupby(platform_col):
            ax.plot(grp["lon"], grp["lat"], "-", lw=0.8, color="0.5", alpha=0.6, **tf)
    h = ax.scatter(obs["lon"], obs["lat"], c=obs["value"], cmap="RdYlBu_r",
                   s=12, edgecolor="none", **tf)
    cb = fig.colorbar(h, ax=ax, shrink=0.8, pad=0.03)
    cb.set_label(var)

    if _HAS_CARTOPY:
        ax.add_feature(cfeature.LAND, facecolor="0.85", zorder=3)
        ax.coastlines(resolution="10m", linewidth=0.5, zorder=4)
        if croco_extent is not None:
            ax.set_extent(croco_extent, crs=ccrs.PlateCarree())
        gl = ax.gridlines(draw_labels=True, linewidth=0.3, color="0.6", alpha=0.5)
        gl.top_labels = False; gl.right_labels = False
        if hasattr(gl, "geo_labels"):
            gl.geo_labels = False   # see gtools.validation._panel's comment -- same cartopy/shapely crash
    elif croco_extent is not None:
        ax.set_xlim(croco_extent[0], croco_extent[1])
        ax.set_ylim(croco_extent[2], croco_extent[3])

    n_platforms = obs[platform_col].nunique() if platform_col else None
    ttl = f"In-situ {var} trajectories (n={len(obs)} obs"
    ttl += f", {n_platforms} platforms)" if n_platforms else ")"
    ax.set_title(ttl, fontsize=11)
    fig.tight_layout()
    if out:
        try:
            fig.savefig(out, dpi=150, bbox_inches="tight")
        except Exception as e:
            # see the matching comment in gtools.validation._three_panel --
            # same cartopy/shapely gridliner defense-in-depth fallback
            print(f"  warning: first savefig attempt failed ({e}); "
                 f"retrying with gridlines removed...")
            for a in fig.axes:
                for artist in list(getattr(a, "_gridliners", [])):
                    artist.remove()
            fig.savefig(out, dpi=150, bbox_inches="tight")
        plt.close(fig)
        return out
    return fig


# ======================================================================
# 7) In-situ PROFILE overlay: CROCO (+ GLORYS) vs one observed profile
# ======================================================================
def plot_insitu_profile(croco_his, obs, platform_id, var="temp", glorys_file=None,
                        Yorig=None, time_tol_hours=6.0, out=None, figsize=(4.5, 6),
                        platform_col=None):
    """Overlay CROCO's (and optionally GLORYS') vertical profile against ONE
    observed in-situ profile (e.g. one Argo cast, one CTD station) at its
    location/time — the point-level complement to the depth-layer scorecard.

    obs         : DataFrame from load_insitu_many() (needs 'lon','lat','depth',
                  'value','time', plus a platform identifier column).
    platform_id : value in that platform column identifying the single
                  profile to plot (obs is filtered to it).
    """
    platform_col = platform_col or next(
        (c for c in ("platform", "platform_code", "id", "wmo") if c in obs.columns), None)
    prof = obs if platform_col is None else obs[obs[platform_col] == platform_id]
    prof = prof.dropna(subset=["depth", "value"]).sort_values("depth")
    if prof.empty:
        print(f"plot_insitu_profile: no obs found for platform_id={platform_id!r}")
        return None
    lon0, lat0 = float(prof["lon"].iloc[0]), float(prof["lat"].iloc[0])
    date = str(pd.to_datetime(prof["time"].iloc[0]).date())

    fig, ax = plt.subplots(figsize=figsize)
    ax.plot(prof["value"], -prof["depth"], "^-", color="k", ms=4, lw=1.4, label="in-situ")

    ds = pp.open_history(croco_his, Yorig=Yorig)
    try:
        model_times = pd.to_datetime(np.atleast_1d(ds["time"].values)).values
        obs_time = pd.to_datetime([date]).values[0]
        ti = int(np.argmin(np.abs(model_times - obs_time)))
        dt_h = abs((model_times[ti] - obs_time) / np.timedelta64(1, "h"))
        if dt_h > time_tol_hours:
            print(f"  note: nearest CROCO record is {dt_h:.1f} h from the obs "
                  f"(tolerance {time_tol_hours} h) - showing it anyway for context")
        cpro = pp.profile(ds, var, lon0, lat0, tindex=ti)
        ax.plot(cpro.values, cpro["depth"].values, "o-", color="C0", ms=3, lw=1.6, label="CROCO")
        lab = cpro.attrs.get("long_name", var); un = cpro.attrs.get("units", "")
    finally:
        ds.close()

    if glorys_file:
        gz, gv = _parent_column(glorys_file, var, lon0, lat0, date=date)
        ax.plot(gv, gz, "s--", color="C3", ms=3, lw=1.4, label="GLORYS/parent")

    ax.set_xlabel(f"{lab} ({un})" if un else var)
    ax.set_ylabel("depth (m)")
    ax.grid(alpha=0.3); ax.legend()
    ax.set_title(f"{platform_id}  ({lon0:.2f}, {lat0:.2f})  {date}", fontsize=10)
    if out:
        fig.savefig(out, dpi=150, bbox_inches="tight")
        plt.close(fig)
        return out
    return fig



"""
gtools/validation_satellite.py — CROCO forecast vs satellite SST/SSS.

Compares the CROCO forecast against independent Copernicus Marine satellite
products, one file per day:

    - OSTIA  (SST_GLO_SST_L4_NRT_OBSERVATIONS_010_001,
              subdataset METOFFICE-GLO-SST-L4-NRT-OBS-SST-V2)   - SST, L4,
              gap-free, multi-sensor analysis, 0.05deg daily.
    - ODYSSEA (SST_GLO_SST_L3S_NRT_OBSERVATIONS_010_010,
              subdataset IFREMER-GLOB-SST-L3-NRT-OBS_FULL_TIME_SERIE) - SST,
              L3S, single-sensor-type composite, has real data gaps
              (swath/cloud).
    - SMOS   (MULTISCALE_SEA_SURFACE_SALINITY_NRT_015_001 -- verify the
              exact CMEMS dataset id in the catalogue before relying on
              this in production) - SSS, L4, gap-filled, ~25 km / daily.

OSTIA/ODYSSEA are GHRSST-convention files: variable ``analysed_sst``/
``sea_surface_temperature`` in Kelvin. SMOS L4 SSS files carry salinity
directly in PSU (no unit conversion). This module applies each product's
quality/land mask before comparing.

Design mirrors gtools.validation (CROCO-vs-parent), except each day gets
its OWN figure (3 panels: CROCO | satellite | bias), rather than stacking
every day into one tall figure -- easier to read/save/share individually,
especially once a cycle spans more than a couple of days.
    - domain statistics (bias, RMSE, centred RMSE, corr) per day, per product.

Availability
------------
Every public entry point in this module is availability-guarded: use
`gtools.download.cmems.dataset_available('ostia_l4' | 'odyssea_l3s' |
'smos_l4_sss')` first (see the validation notebook), or pass files directly
and this module will simply skip (return None / an empty stats table) any
day whose file is missing or unreadable, rather than raising - a data gap
in an L3S composite (or an SSS retrieval near the coast/ice) is normal, not
an error.

Typical use
-----------
    import gtools.validation_satellite as vs

    files = {"2026-07-24": "model-runs/Canary_12/20260724/downloaded_data/OSTIA/2026-07-24.nc",
             "2026-07-25": "model-runs/Canary_12/20260724/downloaded_data/OSTIA/2026-07-25.nc"}
    figs, stats = vs.compare_satellite_grid(
        "model-runs/Canary_12/20260724/fcst/CROCO_FILES/croco_his.nc",
        files, product="OSTIA", Yorig=2000,
        out_dir="model-runs/Canary_12/validation_20260724")
    # figs == {"2026-07-24": ".../sst_vs_ostia_2026-07-24.png", ...}

    # SMOS SSS works the same way, driven by the same SAT_PRODUCTS registry:
    figs, stats = vs.compare_satellite_grid(croco_his, smos_files, product="SMOS",
                                            out_dir=validation_dir)

# compare_sst_satellite_grid() is kept as a thin backward-compatible alias
# for compare_satellite_grid() (SST products only) -- see below.
"""


# ======================================================================================
# ======================================================================================
# Satellite SST/SSS validation (formerly gtools/validation_satellite.py,
# merged in here so every validation function lives in this one module).
#
# Compares the CROCO forecast against independent Copernicus Marine
# satellite products, one file per day: OSTIA (SST, L4), ODYSSEA (SST,
# L3S), SMOS (SSS, L4). See SAT_PRODUCTS below for the per-product
# variable/mask/QC handling. Figures use the same 2x2 (CROCO, satellite,
# bias, RMSE) layout as compare_sst/sss/ssh/currents above, with the RMSE
# panel computed across every day the product has data for in this cycle.
# ======================================================================================
# ======================================================================================


# ======================================================================
# Product registry (variable names + units, per Copernicus product).
# `crocovar`/`label`/`units`/`cmap` drive which CROCO field this product is
# compared against and how it's labelled -- everything else in this module
# is generic across SST and SSS products.
# ======================================================================
SAT_PRODUCTS = {
    "OSTIA": {
        "dataset_key": "ostia_l4",
        "var": "analysed_sst",
        "qc_var": None,            # L4: gap-filled, effectively no per-pixel QC gate
        "mask_var": "mask",        # 1 = open sea (GHRSST convention: 1=sea,2=land,4=lake,8=ice)
        "units_kelvin": True,
        "crocovar": "temp", "label": "SST", "units": "degC", "cmap": "RdYlBu_r",
    },
    "ODYSSEA": {
        "dataset_key": "odyssea_l3s",
        "var": "sea_surface_temperature",
        "qc_var": "quality_level", # keep >=4 (GHRSST: 0=no data..5=best)
        "qc_good": (4, 5),
        "mask_var": None,
        "units_kelvin": True,
        "crocovar": "temp", "label": "SST", "units": "degC", "cmap": "RdYlBu_r",
    },
    "SMOS": {
        "dataset_key": "smos_l4_sss",
        "var": "sos",              # CMEMS SMOS L4 SSS variable name -- verify
                                    # against the actual product before use
        "qc_var": None,
        "mask_var": None,
        "units_kelvin": False,      # PSU already, no conversion
        "crocovar": "salt", "label": "SSS", "units": "PSU", "cmap": "viridis",
    },
}


# ======================================================================
# Loading
# ======================================================================
def load_satellite_field(fname, product, date=None):
    """Load one satellite file as (lon2d, lat2d, field) in the product's
    natural units (degC for SST products, PSU for SMOS SSS).

    product : one of SAT_PRODUCTS ('OSTIA', 'ODYSSEA', 'SMOS').
    Applies the product's land mask / QC flag before returning; bad/land/
    missing points are NaN. Returns None if the file can't be read or the
    expected variable isn't in it (so callers can skip that day/file).
    """
    spec = SAT_PRODUCTS[product]
    try:
        ds = xr.open_dataset(fname)
    except Exception as e:
        print(f"  [{product}] could not open {fname}: {e}")
        return None
    try:
        if spec["var"] not in ds:
            print(f"  [{product}] '{spec['var']}' not in {fname}; has {list(ds.data_vars)}")
            return None
        da = ds[spec["var"]]
        if "time" in da.dims:
            da = da.sel(time=np.datetime64(date), method="nearest") if date else da.isel(time=0)
        field = np.asarray(da.values, dtype=float)
        if spec["units_kelvin"]:
            field = field - 273.15

        if spec.get("mask_var") and spec["mask_var"] in ds:
            m = np.asarray(ds[spec["mask_var"]].values)
            if m.ndim == field.ndim + 1:      # mask sometimes carries a time dim too
                m = m[0] if date is None else np.squeeze(m)
            field = np.where(m == 1, field, np.nan)
        if spec.get("qc_var") and spec["qc_var"] in ds:
            q = ds[spec["qc_var"]]
            if "time" in q.dims:
                q = q.sel(time=np.datetime64(date), method="nearest") if date else q.isel(time=0)
            q = np.asarray(q.values)
            field = np.where(np.isin(q, spec["qc_good"]), field, np.nan)

        lon_name = "lon" if "lon" in ds else "longitude"
        lat_name = "lat" if "lat" in ds else "latitude"
        lon = ds[lon_name].values
        lat = ds[lat_name].values
        if lon.ndim == 1:
            lon2d, lat2d = np.meshgrid(lon, lat)
        else:
            lon2d, lat2d = lon, lat
        return lon2d, lat2d, field
    finally:
        ds.close()


# Backward-compatible alias (this function was SST-only before SMOS/SSS support).
load_satellite_sst = load_satellite_field


def _croco_field_satellite(croco_his, product, date, Yorig=None, margin_deg=None,
                            daily_mean=True):
    """CROCO surface field (matching `product`'s crocovar) on its native
    grid, for the CROCO record nearest `date` (or the daily mean of that
    date's records if daily_mean=True, matching _field_pair_for_day)."""
    crocovar = SAT_PRODUCTS[product]["crocovar"]
    ds = pp.open_history(croco_his, Yorig=Yorig)
    if margin_deg is not None:
        ds = pp.crop_interior(ds, margin_deg=margin_deg)
    if daily_mean:
        ds = _maybe_daily_mean(ds, date, daily_mean)
        tindex = -1
    else:
        tt = pp.times(ds).astype("datetime64[D]")
        hits = np.where(tt == np.datetime64(date, "D"))[0]
        tindex = int(hits[-1]) if hits.size else -1
    field = pp.surface(ds, crocovar, tindex=tindex).values
    clon, clat, _ = pp.lonlatmask(ds)
    ds.close()
    return clon, clat, field


# ======================================================================
# One figure per day: CROCO | satellite | bias
# ======================================================================
def _day_figure(clon, clat, croco_field, sat_on_croco, bias_map, rmse_map, day_label,
                vmin, vmax, product, n_days=None, figsize=(11, 10)):
    """Build the 2x2 (CROCO | satellite / bias | RMSE) figure for one day.
    bias_map/rmse_map are computed across every day this product has data
    for in the cycle (see compare_satellite_grid) -- same for every day's
    figure, so it reflects the WHOLE cycle's skill, with each day's own
    CROCO/satellite snapshot alongside it for context."""
    spec = SAT_PRODUCTS[product]
    label, units, cmap = spec["label"], spec["units"], spec["cmap"]
    crocovar = spec["crocovar"]
    bias_limit = BIAS_RMSE_LIMITS.get(crocovar)
    if bias_limit is not None:
        blim, rlo, rhi = bias_limit, 0.0, bias_limit
    else:
        finite_bias = bias_map[np.isfinite(bias_map)]
        blim = np.nanpercentile(np.abs(finite_bias), 98) if finite_bias.size else 1.0
        finite_rmse = rmse_map[np.isfinite(rmse_map)]
        rlo, rhi = 0.0, (np.nanpercentile(finite_rmse, 98) if finite_rmse.size else 1.0)

    proj = {"projection": ccrs.PlateCarree()} if _HAS_CARTOPY else {}
    fig, axes = plt.subplots(2, 2, figsize=figsize, subplot_kw=proj,
                             constrained_layout=True)

    h_topl = _panel_nocb(axes[0, 0], clon, clat, croco_field, cmap, vmin, vmax,
                        f"CROCO {label}\n{day_label}")
    h_topr = _panel_nocb(axes[0, 1], clon, clat, sat_on_croco, cmap, vmin, vmax,
               f"{product} {label}\n{day_label}")
    
    h_bias = _panel_nocb(axes[1, 0], clon, clat, bias_map, "RdBu_r", -blim, blim,
                         f"Bias (CROCO - {product})")
    rmse_ttl = f"RMSE ({n_days} day{'s' if n_days != 1 else ''})" if n_days else "RMSE"
    h_rmse = _panel_nocb(axes[1, 1], clon, clat, rmse_map, "magma_r", rlo, rhi, rmse_ttl)

    cb_topl = fig.colorbar(h_topl, ax=axes[0, 0], shrink=0.85, pad=0.02, location="right")
    cb_topl.set_label(units, fontsize=9)
    cb_topr = fig.colorbar(h_topr, ax=axes[0, 1], shrink=0.85, pad=0.02, location="right")
    cb_topr.set_label(units, fontsize=9)
    cb_bias = fig.colorbar(h_bias, ax=axes[1, 0], shrink=0.85, pad=0.02, location="right")
    cb_bias.set_label(units, fontsize=9)
    cb_rmse = fig.colorbar(h_rmse, ax=axes[1, 1], shrink=0.85, pad=0.02, location="right")
    cb_rmse.set_label(units, fontsize=9)

    fig.suptitle(f"CROCO forecast {label} vs {product}", fontsize=13)
    return fig


def compare_satellite_grid(croco_his, sat_files, product, Yorig=None, out_dir=None,
                           margin_deg=None, figsize=(11, 10)):
    """CROCO-vs-satellite comparison, ONE 2x2 FIGURE PER DAY (CROCO |
    satellite / bias | RMSE), rather than stacking every day into a single
    tall figure. Bias/RMSE are each day's OWN comparison (that day's
    croco_field vs. that day's regridded satellite field) -- distinct from
    day to day, not a single cycle-wide map repeated across every figure.

    croco_his : CROCO history file (spans the days requested).
    sat_files : dict {date_str ('YYYY-MM-DD'): satellite_file_path}. Days
                whose file is missing/unreadable are skipped (printed, not
                raised) - this is the expected/normal situation for ODYSSEA
                (L3S has real swath/cloud gaps) and SMOS SSS (coastal/ice
                gaps), and any day CROCO doesn't cover is skipped too.
    product   : 'OSTIA', 'ODYSSEA', or 'SMOS' (see SAT_PRODUCTS).
    out_dir   : if given, each day's figure is saved as
                "<out_dir>/<field>_vs_<product>_<date>.png" (field = 'sst'
                or 'sss', lowercased from the product's label) and closed;
                if None, figures are returned open (not saved).

    Returns (figs, stats_df):
      figs      : dict {date_str: fig_or_saved_path}, one entry per day that
                  could be compared (empty if none could).
      stats_df  : one row per day with bias/rmse/crmse/corr
                  (gtools.validation.domain_statistics), plus a 'MEAN' row
                  across all valid days. Empty if no day could be compared.
    """
    if product not in SAT_PRODUCTS:
        raise ValueError(f"product must be one of {list(SAT_PRODUCTS)}, got {product!r}")
    label = SAT_PRODUCTS[product]["label"]

    rows = []
    panels = []   # (day, clon, clat, croco_field, sat_on_croco)
    for day, fname in sorted(sat_files.items()):
        if not fname or not os.path.exists(fname):
            print(f"  [{product}] {day}: no file - skipping")
            continue
        loaded = load_satellite_field(fname, product, date=day)
        if loaded is None:
            continue
        plon, plat, pfield = loaded
        if not np.isfinite(pfield).any():
            print(f"  [{product}] {day}: all points masked/NaN - skipping")
            continue
        try:
            clon, clat, croco_field = _croco_field_satellite(croco_his, product, day,
                                                    Yorig=Yorig, margin_deg=margin_deg)
        except Exception as e:
            print(f"  [{product}] {day}: CROCO record not available ({e}) - skipping")
            continue

        ds_grid = xr.Dataset({"mask_rho": (("eta_rho", "xi_rho"), np.where(np.isfinite(croco_field), 1.0, np.nan))},
                             coords={"lon_rho": (("eta_rho", "xi_rho"), clon),
                                     "lat_rho": (("eta_rho", "xi_rho"), clat)})
        sat_on_croco = regrid_to_croco(plon, plat, pfield, ds_grid)

        s = domain_statistics(croco_field, sat_on_croco)
        s["date"] = day
        rows.append(s)
        panels.append((day, clon, clat, croco_field, sat_on_croco))
        print(f"  [{product}] {day}:"); _print_stats(label, s)

    stats_df = pd.DataFrame(rows)
    if not stats_df.empty:
        mean_row = stats_df.drop(columns=["date"]).mean(numeric_only=True).to_dict()
        mean_row["date"] = "MEAN"
        stats_df = pd.concat([stats_df, pd.DataFrame([mean_row])], ignore_index=True)

    if not panels:
        print(f"  [{product}] no comparable day found - nothing to plot")
        return {}, stats_df

    # shared colour scale across every day, so different days' figures stay
    # visually comparable (2nd/98th percentile -- see gtools.plotting.percentile_clim)
    vmin, vmax = percentile_clim(*(p[3] for p in panels), *(p[4] for p in panels))

    figs = {}
    for day, clon, clat, croco_field, sat_on_croco in panels:
        # this day's OWN bias/error, not a cycle-wide aggregate -- sqrt of a
        # single squared sample reduces to |diff|, so "RMSE" here is that
        # day's absolute error map (n_days=1, distinct from every other day).
        diff = croco_field - sat_on_croco
        day_bias_map = diff
        day_rmse_map = np.abs(diff)
        fig = _day_figure(clon, clat, croco_field, sat_on_croco, day_bias_map, day_rmse_map,
                          day, vmin, vmax, product, n_days=1, figsize=figsize)
        if out_dir:
            fname = f"{label.lower()}_vs_{product.lower()}_{day}.png"
            path = os.path.join(out_dir, fname)
            fig.savefig(path, dpi=150, bbox_inches="tight")
            plt.close(fig)
            figs[day] = path
        else:
            figs[day] = fig

    return figs, stats_df
    
def compare_sst_satellite_grid(croco_his, sat_files, product, Yorig=None, out=None,
                               margin_deg=None, figsize_per_row=(11, 10)):
    """Backward-compatible SST-only wrapper around compare_satellite_grid().

    Deprecated in favour of compare_satellite_grid() (which also handles
    SMOS SSS and returns one figure per day instead of one stacked figure).
    `out`, if given, is treated as an OUTPUT DIRECTORY (not a single file
    path any more, since there's one figure per day) -- pass out_dir=... to
    compare_satellite_grid() directly to make that explicit in new code.
    """
    figs, stats_df = compare_satellite_grid(croco_his, sat_files, product, Yorig=Yorig,
                                            out_dir=out, margin_deg=margin_deg,
                                            figsize=figsize_per_row)
    return figs, stats_df
"""
gtools/validation_report.py — CROCO HTML validation summary.

Gathers every figure (.png) and statistics table (.csv) written into one
cycle's validation output directory (VALIDATION_DIR, see notebooks/_paths.py
and 02_validation.ipynb) into a single, self-contained, offline-viewable
HTML page -- so a whole cycle's V1 validation can be reviewed (or shared,
e.g. attached to an email or posted somewhere) as ONE file, instead of
opening a folder of a few dozen loose PNGs and CSVs one at a time.

Nothing here re-computes any statistic or re-reads any model/reference
file -- it only reads what 02_validation.ipynb (or validate_all_cycles.sh,
which runs that same notebook per cycle) already wrote to VALIDATION_DIR,
and lays it out. Safe to call multiple times (each call fully rebuilds the
page from whatever is currently in the directory); safe to call on a
directory that's missing some pieces (a product that was unavailable this
cycle, or in-situ having no coverage) -- those sections are simply omitted,
not shown as broken/empty.

Typical use (see 02_validation.ipynb's Section 9)
---------------------------------------------------
    import gtools.validation_report as vr
    html_path = vr.build_html_summary(VALIDATION_DIR, cycle=CYCLE, config=CONFIG)
    print(f"Open {html_path} in a browser to review this cycle.")
"""


# ======================================================================================
# ======================================================================================
# HTML validation summary (formerly gtools/validation_report.py, merged in
# here so every validation function lives in this one module).
#
# Gathers every figure (.png) and statistics table (.csv) written into one
# cycle's VALIDATION_DIR into a single, self-contained, offline-viewable
# HTML page: every image is embedded inline (base64 data URI, not a link
# to a separate file), click-to-enlarge (lightbox), with a "Download"
# link -- see build_html_summary() below.
# ======================================================================================
# ======================================================================================


# ======================================================================
# Figure grouping -- ordered (section title, [filename-prefix, ...]) pairs.
# Every *.png under VALIDATION_DIR is matched against these prefixes, in
# order; the first match wins. Anything matching nothing lands in a final
# "Other figures" catch-all, so a new/renamed figure never silently
# vanishes from the report even if this list falls behind the notebook.
# ======================================================================
FIGURE_GROUPS = [
    ("Bias & RMSE maps -- CROCO vs Copernicus Marine Forecast",
     ["sst_vs_forecast_", "ssh_vs_forecast_", "currents_vs_forecast_", "sss_vs_forecast_"]),
    ("Depth-resolved CROCO vs Copernicus Marine Forecast",
     ["depth_levels_"]),
    ("Vertical profile & error-vs-depth",
     ["profile_", "error_vs_depth_"]),
    ("Scatter plots -- pointwise CROCO vs Copernicus Marine Forecast",
     ["scatter_"]),
    ("Taylor diagrams",
     ["taylor_diagram_"]),
    ("Time series",
     ["timeseries_"]),
    ("Domain-wide bias boxplots (per day)",
     ["boxplot_bias_"]),
    ("Domain-wide RMSE boxplots (per day)",
     ["boxplot_rmse_"]),
    ("Satellite SST -- OSTIA / ODYSSEA",
     ["sst_vs_ostia_", "sst_vs_odyssea_"]),
    ("Satellite SSS -- SMOS",
     ["sss_vs_smos_"]),
    ("In-situ -- trajectory & profile",
     ["insitu_"]),
]

# CSVs are grouped under a readable title the same way; unmatched CSVs get
# their bare filename as the title instead of being dropped.
STATS_TITLES = {
    "sst_vs_ostia_stats.csv":  "Satellite SST vs OSTIA -- daily statistics",
    "sst_vs_odyssea_stats.csv": "Satellite SST vs ODYSSEA -- daily statistics",
    "sss_vs_smos_stats.csv":  "Satellite SSS vs SMOS -- daily statistics",
    "insitu_scorecard.csv":   "In-situ GODAE scorecard (per depth layer)",
}


def _group_for(fname):
    for title, prefixes in FIGURE_GROUPS:
        if any(fname.startswith(p) for p in prefixes):
            return title
    return "Other figures"


def _escape(s):
    return _html.escape(str(s))


def _status_banner(status):
    if status is None:
        return ""
    all_pass = status.get("all_pass")
    if all_pass is None:
        return ""
    css = "pass" if all_pass else "fail"
    label = "PASS" if all_pass else "FAIL"
    rows = ""
    for c in status.get("criteria", []):
        rows += (f"<tr><td>{_escape(c.get('name', ''))}</td>"
                f"<td class=\"{'pass' if c.get('passed') else 'fail'}\">"
                f"{'PASS' if c.get('passed') else 'FAIL'}</td>"
                f"<td>{_escape(c.get('value', ''))}</td></tr>\n")
    table = (f"<table class=\"criteria\"><tr><th>Criterion</th><th>Result</th>"
            f"<th>Value</th></tr>\n{rows}</table>" if rows else "")
    return (f'<div class="banner {css}">Overall V1 status: <b>{label}</b></div>\n{table}')


def _img_data_uri(path):
    """Read a PNG and return it as a data: URI, so the HTML page can embed
    the image bytes directly (<img src="data:image/png;base64,...">)
    instead of a relative link to the file on disk -- the whole point
    being that the resulting .html is a single, self-contained file: it
    still opens and shows every figure even after being copied/emailed/
    moved away from its VALIDATION_DIR, with no PNGs alongside it.
    Returns None (caller falls back to a broken-image-safe placeholder)
    if the file can't be read.
    """
    try:
        with open(path, "rb") as f:
            data = base64.b64encode(f.read()).decode("ascii")
        return f"data:image/png;base64,{data}"
    except Exception:
        return None


def _figure_section(title, files):
    if not files:
        return ""
    cards = ""
    for f in sorted(files):
        uri = _img_data_uri(f)
        fname = _escape(os.path.basename(f))
        if uri:
            # data-full holds the same embedded image; the lightbox (see
            # _LIGHTBOX_JS) reads it on click rather than re-fetching
            # anything -- the whole point of embedding is that the page
            # has no external files to fetch in the first place.
            img_tag = (f'<img src="{uri}" data-full="{uri}" data-name="{fname}" '
                      f'loading="lazy" onclick="sfOpenLightbox(this)">')
        else:
            img_tag = '<div class="warn">(could not read this image)</div>'
        cards += (f'<div class="card">'
                 f'{img_tag}'
                 f'<div class="caption">{fname}</div>'
                 f'</div>\n')
    return f'<h2>{_escape(title)}</h2>\n<div class="gallery">\n{cards}</div>\n'


def _stats_section(csv_path, validation_dir):
    title = STATS_TITLES.get(os.path.basename(csv_path), os.path.basename(csv_path))
    try:
        df = pd.read_csv(csv_path)
    except Exception as e:
        return f'<h3>{_escape(title)}</h3>\n<p class="warn">Could not read this file: {_escape(e)}</p>\n'
    table_html = df.to_html(index=False, classes="stats", float_format=lambda x: f"{x:.4f}")
    return f'<h3>{_escape(title)}</h3>\n{table_html}\n'


_CSS = """
body { font-family: -apple-system, Segoe UI, Helvetica, Arial, sans-serif;
      margin: 0; padding: 0 0 3em 0; background: #f7f7f9; color: #222; }
header { background: #14324d; color: white; padding: 1.2em 2em; }
header h1 { margin: 0 0 0.2em 0; font-size: 1.4em; }
header .meta { font-size: 0.9em; color: #cfe3f7; }
main { max-width: 1200px; margin: 0 auto; padding: 0 2em; }
h2 { border-bottom: 2px solid #14324d; padding-bottom: 0.2em; margin-top: 2em; }
h3 { margin-top: 1.5em; color: #14324d; }
.toc { background: white; border: 1px solid #ddd; border-radius: 6px;
      padding: 1em 1.5em; margin: 1.5em 0; }
.toc a { display: block; padding: 0.15em 0; color: #14324d; text-decoration: none; }
.toc a:hover { text-decoration: underline; }
.banner { padding: 0.8em 1.2em; border-radius: 6px; font-size: 1.1em; margin: 1em 0; }
.banner.pass { background: #e3f6e6; border: 1px solid #2e8b46; color: #1b5e2b; }
.banner.fail { background: #fbe7e7; border: 1px solid #c0392b; color: #7d1f16; }
table.criteria, table.stats { border-collapse: collapse; width: 100%; background: white;
      margin-bottom: 1em; font-size: 0.92em; }
table.criteria th, table.criteria td, table.stats th, table.stats td {
      border: 1px solid #ddd; padding: 0.4em 0.7em; text-align: left; }
table.criteria th, table.stats th { background: #eef2f6; }
td.pass { color: #1b5e2b; font-weight: 600; }
td.fail { color: #7d1f16; font-weight: 600; }
.gallery { display: flex; flex-wrap: wrap; gap: 1em; margin: 1em 0 2em 0; }
.card { background: white; border: 1px solid #ddd; border-radius: 6px; padding: 0.5em;
      width: 340px; box-shadow: 0 1px 3px rgba(0,0,0,0.06); }
.card img { width: 100%; border-radius: 3px; display: block; cursor: zoom-in; }
.card .caption { font-size: 0.78em; color: #555; margin-top: 0.4em; word-break: break-all; }
.warn { color: #b06000; }
footer { text-align: center; color: #888; font-size: 0.85em; margin-top: 3em; }

/* ---- lightbox (click-to-enlarge + download) ---- */
#sf-lightbox { display: none; position: fixed; z-index: 999; inset: 0;
      background: rgba(10,15,25,0.92); text-align: center; padding: 2em; cursor: zoom-out; }
#sf-lightbox.open { display: flex; flex-direction: column; align-items: center; justify-content: center; }
#sf-lightbox img { max-width: 92vw; max-height: 78vh; border-radius: 4px;
      box-shadow: 0 4px 24px rgba(0,0,0,0.5); cursor: default; }
#sf-lightbox-bar { margin-top: 1em; display: flex; gap: 1em; align-items: center; cursor: default; }
#sf-lightbox-bar span { color: #cfe3f7; font-size: 0.9em; }
#sf-lightbox-bar a, #sf-lightbox-close { background: #14324d; color: white; border: none;
      padding: 0.5em 1.1em; border-radius: 5px; text-decoration: none; font-size: 0.9em;
      cursor: pointer; }
#sf-lightbox-bar a:hover, #sf-lightbox-close:hover { background: #1c4468; }
"""

_LIGHTBOX_HTML = (
    '<div id="sf-lightbox" onclick="if(event.target===this) sfCloseLightbox()">'
    '<img id="sf-lightbox-img" src="">'
    '<div id="sf-lightbox-bar">'
    '<span id="sf-lightbox-name"></span>'
    '<a id="sf-lightbox-download" href="" download>Download</a>'
    '<button id="sf-lightbox-close" onclick="sfCloseLightbox()">Close</button>'
    '</div></div>\n')

_LIGHTBOX_JS = """
function sfOpenLightbox(imgEl) {
  var full = imgEl.getAttribute('data-full');
  var name = imgEl.getAttribute('data-name');
  document.getElementById('sf-lightbox-img').src = full;
  document.getElementById('sf-lightbox-name').textContent = name;
  var dl = document.getElementById('sf-lightbox-download');
  dl.href = full;
  dl.setAttribute('download', name);
  document.getElementById('sf-lightbox').classList.add('open');
}
function sfCloseLightbox() {
  document.getElementById('sf-lightbox').classList.remove('open');
}
document.addEventListener('keydown', function(e) {
  if (e.key === 'Escape') sfCloseLightbox();
});
"""


def build_html_summary(validation_dir, cycle=None, config=None, status=None,
                       out_name=None, extra_title=None):
    """Build (or rebuild) a single self-contained HTML summary page for one
    cycle's validation output directory.

    validation_dir : the VALIDATION_DIR for this cycle (see notebooks/
                      _paths.get_validation_dir) -- every *.png/*.csv
                      directly under it is picked up. Every image is
                      embedded inline as a base64 data URI (not a link to
                      the separate PNG file), click-to-enlarge with a
                      Download link (see sfOpenLightbox in the page's own
                      <script>) -- so the resulting .html is fully
                      self-contained: it still shows every figure after
                      being copied/emailed/moved away from this directory,
                      with none of the PNGs alongside it.
    cycle, config   : shown in the header if given; also used to look for
                      VALIDATION_DIR/validation_status.json automatically
                      if `status` isn't passed explicitly (see
                      02_validation.ipynb's Section 5, which writes it).
    status          : optional dict (or None to auto-load
                      validation_status.json if present) with at least
                      'all_pass' and 'criteria' -- rendered as a pass/fail
                      banner at the top of the page. No banner if neither
                      is available (e.g. the Copernicus Marine Forecast was
                      unavailable this cycle, so Section 5 never ran).
    extra_title     : optional extra line appended to the header (e.g. a
                      note about which cycles/products were skipped).

    Returns the path to the written HTML file.
    """
    os.makedirs(validation_dir, exist_ok=True)

    if status is None:
        status_path = os.path.join(validation_dir, "validation_status.json")
        if os.path.exists(status_path):
            try:
                with open(status_path) as f:
                    status = json.load(f)
            except Exception:
                status = None

    pngs = sorted(glob.glob(os.path.join(validation_dir, "*.png")))
    csvs = sorted(glob.glob(os.path.join(validation_dir, "*.csv")))

    grouped = {title: [] for title, _ in FIGURE_GROUPS}
    grouped["Other figures"] = []
    for f in pngs:
        grouped[_group_for(os.path.basename(f))].append(f)

    # ---- header ----
    title_bits = ["CROCO Validation Summary"]
    if config:
        title_bits.append(config)
    if cycle:
        title_bits.append(cycle)
    page_title = " -- ".join(title_bits)
    generated = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M UTC")

    header_html = (
        f'<header><h1>{_escape(page_title)}</h1>'
        f'<div class="meta">Generated {generated}'
        + (f' -- {_escape(extra_title)}' if extra_title else '')
        + f' -- {len(pngs)} figure(s), {len(csvs)} statistics table(s)</div></header>\n')

    # ---- table of contents ----
    toc_entries = [title for title, files in grouped.items() if files]
    if csvs:
        toc_entries.append("Statistics")
    toc_html = ('<div class="toc"><b>Contents</b>\n' +
               "\n".join(f'<a href="#{_escape(t.replace(" ", "_"))}">{_escape(t)}</a>'
                        for t in toc_entries) + "\n</div>\n")

    # ---- body: banner, then each figure group, then stats ----
    body = _status_banner(status)
    body += toc_html
    for title, _ in FIGURE_GROUPS + [("Other figures", [])]:
        files = grouped.get(title, [])
        if not files:
            continue
        section = _figure_section(title, files)
        body += f'<a id="{_escape(title.replace(" ", "_"))}"></a>\n' + section

    if csvs:
        body += f'<a id="Statistics"></a>\n<h2>Statistics</h2>\n'
        for csv_path in csvs:
            body += _stats_section(csv_path, validation_dir)

    if not pngs and not csvs:
        body += '<p class="warn">No figures or statistics found in this directory yet.</p>\n'

    html_doc = (
        "<!DOCTYPE html>\n<html lang=\"en\"><head><meta charset=\"utf-8\">"
        f"<title>{_escape(page_title)}</title><style>{_CSS}</style></head>\n"
        f"<body>\n{header_html}<main>\n{body}"
        f'<footer>CROCO -- DCC process V1 -- '
        f'gtools.validation.build_html_summary()</footer>\n'
        f"</main>\n{_LIGHTBOX_HTML}<script>{_LIGHTBOX_JS}</script>\n</body></html>\n")
    out_name=f"validation_cycle_{cycle}.html"
    out_path = os.path.join(validation_dir, out_name)
    with open(out_path, "w") as f:
        f.write(html_doc)
    return out_path

    


# ======================================================================
# NOTEBOOK-SIMPLIFICATION WRAPPERS (02_validation.ipynb)
# ----------------------------------------------------------------------
# Added so 02_validation.ipynb's cells can each be one function call
# instead of an inline guard + loop + print block. Every function here
# wraps existing building blocks above (compare_sst, godae_scorecard_
# croco_vs_glorys, domain_diff, ...) exactly as the notebook cells used
# to call them by hand -- no comparison logic, thresholds, or file
# naming has changed, only where the loop/guard/print now lives.
# ======================================================================

def check_run(ds):
    """Grid/time sanity check + numerical-stability check on an opened CROCO
    history dataset -- the notebook's Section 1. Prints the same report the
    inline cell used to print, and returns the same `stability_ok` bool
    (used by pass_fail_summary()'s "Numerical stability" criterion).
    """
    print(f"grid          : {ds.sizes['eta_rho']} x {ds.sizes['xi_rho']}  "
          f"({ds.sizes.get('s_rho', '?')} sigma levels)")
    print(f"time steps    : {ds.sizes['time']}")
    print(f"time coverage : {pp.times(ds)[0]}  ->  {pp.times(ds)[-1]}")

    stability_ok = True
    for var in ('zeta', 'temp', 'salt', 'u', 'v'):
        vals = ds[var].values
        n_nan = int(np.isnan(vals).sum())
        n_inf = int(np.isinf(vals).sum())
        ok = (n_nan == 0) and (n_inf == 0)
        stability_ok &= ok
        print(f"  {var:5s}: {'OK  ' if ok else 'FAIL'}  (NaN={n_nan}, Inf={n_inf})")

    print()
    print(f"Numerical stability (FR pass criterion): {'PASS' if stability_ok else 'FAIL'}")
    return stability_ok


def download_references(ds, avail, cycle, reference, main_dir, config, get_paths_mod):
    """Section 1c: resolve the domain/cycle window from `ds`, download the
    Copernicus Marine Forecast reference (if available and not already on
    disk / stale) and the satellite SST (OSTIA, ODYSSEA) and SSS (SMOS)
    products for this cycle, one file per day. Prints exactly the same
    progress messages the inline cell used to.

    `get_paths_mod` is the notebook's already-imported `_paths` module
    (passed in rather than imported here, so this stays decoupled from the
    notebooks/ layout).

    Returns
    -------
    DOMAIN     : (lon_min, lon_max, lat_min, lat_max)
    CYCLE_DAYS : sorted list of 'YYYY-MM-DD' strings, one per calendar day
                 covered by the cycle
    SAT_FILES  : {"OSTIA": {day: path, ...}, "ODYSSEA": {...}, "SMOS": {...}}
    START_DATE, END_DATE : datetime, the cycle's exact time window (also
                 needed by Section 8's optional in-situ download, which
                 wants the real sub-day start/end, not just CYCLE_DAYS'
                 calendar-day strings)
    """
    from datetime import datetime

    clon_full, clat_full, _ = pp.lonlatmask(ds)
    DOMAIN = (float(np.nanmin(clon_full)), float(np.nanmax(clon_full)),
              float(np.nanmin(clat_full)), float(np.nanmax(clat_full)))
    model_times = pd.to_datetime(pp.times(ds))
    START_DATE = model_times[0].to_pydatetime()
    END_DATE = model_times[-1].to_pydatetime()
    CYCLE_DATE = datetime.strptime(cycle, "%Y%m%d")
    print(f"Domain: {DOMAIN}")
    print(f"Cycle window: {START_DATE} -> {END_DATE}")

    CYCLE_DAYS = sorted(set(model_times.strftime('%Y-%m-%d')))

    # ---- (i) Copernicus Marine Forecast (Mercator anfc), combined reference file ----
    print("\n-- Copernicus Marine Forecast --")
    if not avail['mercator_forecast']:
        print("unavailable on the CMEMS platform - Sections 2/2b will be skipped.")
    elif netcdf_covers_time_range(reference, START_DATE, END_DATE, max_step_hours=30):
        print(f"already downloaded and covers the full cycle window at the expected "
              f"resolution: {reference}")
    else:
        if os.path.exists(reference):
            print(f"{reference} exists but either doesn't cover the full cycle window "
                  f"({START_DATE.date()} -> {END_DATE.date()}) or is at a coarser "
                  f"resolution than expected -- re-downloading.")
            os.remove(reference)
        mercator_dir = os.path.dirname(reference)
        fdays = max((END_DATE.date() - CYCLE_DATE.date()).days, 0)
        # download_mercator_ops() downloads its 4 variables via 4 concurrent
        # threads, each calling into copernicusmarine/netCDF4 -- a known
        # kernel-crashing combination (HDF5 isn't thread-safe in most
        # builds). Using the sequential variant here instead; see that
        # function's docstring for the full explanation.
        download_mercator_ops_sequential(DOMAIN, CYCLE_DATE, hdays=0, fdays=fdays, outputDir=mercator_dir)
        if netcdf_covers_time_range(reference, START_DATE, END_DATE, max_step_hours=30):
            print(f"downloaded -> {reference}")
        elif os.path.exists(reference):
            print(f"download ran but {reference} still doesn't cover the full cycle window "
                  f"at the expected resolution -- one or more variable downloads may have "
                  f"failed; check the log above.")
        else:
            print(f"download ran but {reference} wasn't produced - check {mercator_dir} for the actual filename.")

    # ---- (ii) Satellite SST: OSTIA & ODYSSEA, one file per day ----
    print("\n-- Satellite SST --")
    SAT_FILES = {}
    for product in ("OSTIA", "ODYSSEA"):
        sat_dir = get_paths_mod.satellite_dir(main_dir, config, cycle, product)
        SAT_FILES[product] = download_satellite_sst(product, DOMAIN, START_DATE, END_DATE, sat_dir)

    # ---- (ii-b) Satellite SSS: SMOS L4, one file per day ----
    print("\n-- Satellite SSS (SMOS) --")
    if not avail['smos_l4_sss']:
        print("unavailable on the CMEMS platform - Section 7b will be skipped.")
        SAT_FILES['SMOS'] = {}
    else:
        smos_dir = get_paths_mod.satellite_dir(main_dir, config, cycle, "SMOS")
        SAT_FILES['SMOS'] = download_satellite_sst("SMOS", DOMAIN, START_DATE, END_DATE, smos_dir)

    return DOMAIN, CYCLE_DAYS, SAT_FILES, START_DATE, END_DATE


def _validate_maps_by_day(compare_fn, croco_his, reference, days, yorig, validation_dir,
                           prefix, available, skip_msg, **compare_kwargs):
    """Shared body of validate_sst_maps/validate_ssh_maps/validate_current_maps/
    validate_sss_maps: guard, loop `compare_fn` over `days`, print, return
    {day: stats} -- everything Section 2's four near-identical cells used to
    do by hand. Not itself one of the notebook's four public calls, since
    compare_ssh() doesn't take depth_m while the other three do.
    """
    if not (available and os.path.exists(reference)):
        print(skip_msg)
        return {}
    stats_by_day = {}
    for day in days:
        _, stats_by_day[day] = compare_fn(
            croco_his, reference, date=day, Yorig=yorig, all_dates=days,
            out=os.path.join(validation_dir, f'{prefix}_vs_forecast_{day}.png'),
            **compare_kwargs)
    print(f"-> {len(stats_by_day)} figure(s) written, one per day: {days}")
    return stats_by_day


def validate_sst_maps(croco_his, reference, days, yorig, depth_m, validation_dir, available):
    """Section 2 -- SST bias maps, CROCO forecast vs Copernicus Marine Forecast,
    one figure per day. Returns {day: stats_dict}."""
    return _validate_maps_by_day(
        compare_sst, croco_his, reference, days, yorig, validation_dir, 'sst', available,
        "Copernicus Marine Forecast unavailable or REFERENCE missing - skipping this comparison.",
        depth_m=depth_m)


def validate_ssh_maps(croco_his, reference, days, yorig, validation_dir, available):
    """Section 2 -- SSH bias maps, CROCO forecast vs Copernicus Marine Forecast,
    one figure per day. SSH has no depth dimension -- always surface,
    regardless of DEPTH_M. Returns {day: stats_dict}."""
    return _validate_maps_by_day(
        compare_ssh, croco_his, reference, days, yorig, validation_dir, 'ssh', available,
        "Copernicus Marine Forecast unavailable or REFERENCE missing - skipping this comparison.")


def validate_current_maps(croco_his, reference, days, yorig, depth_m, validation_dir, available):
    """Section 2 -- surface-current bias maps, CROCO forecast vs Copernicus
    Marine Forecast, one figure per day. Returns {day: stats_dict}.

    Surface velocities pass criterion is QUALITATIVE (visual consistency
    with expected gyre/coastal-jet circulation) -- inspect the vector maps.
    """
    stats_by_day = _validate_maps_by_day(
        compare_currents, croco_his, reference, days, yorig, validation_dir, 'currents', available,
        "Copernicus Marine Forecast unavailable or REFERENCE missing - skipping this comparison.",
        depth_m=depth_m)
    if stats_by_day:
        print()
        print("Surface velocities pass criterion is QUALITATIVE (visual consistency with")
        print("expected gyre/coastal-jet circulation) -- inspect the vector maps")
    return stats_by_day


def validate_sss_maps(croco_his, reference, days, yorig, depth_m, validation_dir, available):
    """Section 2 -- SSS bias maps, CROCO forecast vs Copernicus Marine Forecast,
    one figure per day. Returns {day: stats_dict}."""
    return _validate_maps_by_day(
        compare_sss, croco_his, reference, days, yorig, validation_dir, 'sss', available,
        "Copernicus Marine Forecast unavailable or REFERENCE missing - skipping this comparison.",
        depth_m=depth_m)


def validate_profiles(ds, croco_his, reference, days, yorig, validation_dir, available,
                       points=None):
    """Section 2b -- point vertical profile + error-vs-depth, CROCO forecast
    vs Copernicus Marine Forecast. One profile figure per (point, day,
    variable in temp/salt/speed), plus one error-vs-depth figure for temp.

    `points`: list of (lon, lat) tuples, or None to auto-pick one coastal-ish
    point (75% across the grid, mid-latitude) -- same default as before.

    Returns (points_used, depth_stats) or (None, None) if skipped.
    """
    if not (available and os.path.exists(reference)):
        print("Copernicus Marine Forecast unavailable or REFERENCE missing - skipping vertical profile/error-vs-depth.")
        return None, None

    clon2, clat2, _ = pp.lonlatmask(ds)
    if points:
        points = list(points)
    else:
        j0p = ds.sizes['eta_rho'] // 2
        i0p = int(0.75 * ds.sizes['xi_rho'])
        points = [(float(clon2[j0p, i0p]), float(clat2[j0p, i0p]))]

    profile_days = sorted(set(pd.to_datetime(pp.times(ds)).strftime('%Y-%m-%d')))
    print(f"Profiling {len(points)} point(s) x {len(profile_days)} day(s): {points}")

    for lon0p, lat0p in points:
        for day in profile_days:
            for var in ('temp', 'salt', 'speed'):
                out_png = os.path.join(
                    validation_dir,
                    f'profile_{var}_vs_forecast_{lon0p:.2f}_{lat0p:.2f}_{day}.png')
                compare_profile(croco_his, reference, var, lon0p, lat0p, date=day,
                                 Yorig=yorig, out=out_png)

    _, depth_stats = error_vs_depth(croco_his, reference, field='temp',
                                     depths=(0, 50, 100, 200, 500), Yorig=yorig)
    return points, depth_stats


def validate_domain_profiles(croco_his, reference, days, yorig, validation_dir, available,
                              variables=('temp', 'salt', 'speed')):
    """Section 2b-bis -- full-domain vertical profile (spatial spread), CROCO
    forecast vs Copernicus Marine Forecast, one figure per (day, variable)."""
    if not (available and os.path.exists(reference)):
        print("Copernicus Marine Forecast unavailable or REFERENCE missing - skipping full-domain vertical profile.")
        return
    for day in days:
        for var in variables:
            out_png = os.path.join(validation_dir, f'profile_domain_{var}_vs_forecast_{day}.png')
            domain_profile(croco_his, reference, var, day, Yorig=yorig, out=out_png)


def validate_depth_levels(croco_his, reference, variable, days, yorig, validation_dir, available,
                           depths=(0, 120, 300, 1000)):
    """Section 2c -- depth-resolved comparison at several levels, CROCO
    forecast vs Copernicus Marine Forecast, one figure per day, for one
    variable ('salt', 'temp' or 'speed'). Called once per variable, same as
    the notebook's three near-identical cells used to (this version prints
    the correct product name for every variable -- the original cells all
    printed "Salinity Forecast" verbatim, even for temp/speed).

    Returns {day: fig}.
    """
    if not (available and os.path.exists(reference)):
        print(f"Copernicus Marine Forecast unavailable or REFERENCE missing - skipping this comparison.")
        return {}
    depth_figs = {}
    for day in days:
        depth_figs[day] = compare_forecast_depth_levels(
            croco_his, reference, var=variable, date=day,
            depths=depths, Yorig=yorig, daily_mean=True,
            out=os.path.join(validation_dir, f'depth_levels_{variable}_vs_forecast_{day}.png'))
    print(f"-> {len(depth_figs)} figure(s) written, one per day: {days}")
    return depth_figs


def scatter_maps_by_day(ds, reference, days, validation_dir):
    """Section 3 -- one 3-panel scatter figure (SSH, SSS, SST vs Copernicus
    Marine Forecast) per day. `ds` must already be open (no availability
    guard here -- callers check AVAIL['mercator_forecast'] and
    os.path.exists(reference) first, same as every other Section-2-family
    cell). Returns {day: out_png}.
    """
    scatter_figs = {}
    for day in days:
        ti = _tindex_for_date(ds, day, -1)
        fig, axes = plt.subplots(1, 3, figsize=(16, 5))

        ssh = ds['zeta'].isel(time=ti).values
        plon, plat, pfield = load_parent(reference, 'ssh', date=day)
        ssh_ref = regrid_to_croco(plon, plat, pfield, ds)
        scatter_vs_reference(ssh, ssh_ref, 'SSH', 'm', ax=axes[0])

        sss = pp.surface(ds, 'salt', tindex=ti).values
        plon, plat, pfield = load_parent(reference, 'salt', date=day)
        sss_ref = regrid_to_croco(plon, plat, pfield, ds)
        scatter_vs_reference(sss, sss_ref, 'SSS', 'PSU', ax=axes[1])

        sst = pp.surface(ds, 'temp', tindex=ti).values
        plon, plat, pfield = load_parent(reference, 'temp', date=day)
        sst_ref = regrid_to_croco(plon, plat, pfield, ds)
        scatter_vs_reference(sst, sst_ref, 'SST', 'degC', ax=axes[2])

        fig.suptitle(day)
        fig.tight_layout()
        out_png = os.path.join(validation_dir, f"scatter_sss_sst_ssh_vs_forecast_{day}.png")
        fig.savefig(out_png, dpi=250, bbox_inches="tight")
        plt.close(fig)
        scatter_figs[day] = out_png
    print(f"-> {len(scatter_figs)} figure(s) written, one per day: {list(scatter_figs)}")
    return scatter_figs


def build_godae_scorecard(croco_his, reference, days, yorig, available,
                           variables=('temp', 'ssh', 'salt', 'speed'), print_table=True):
    """Section 4 -- GODAE scorecard, CROCO forecast vs Copernicus Marine
    Forecast, scored ONE DAY AT A TIME (rather than just the cycle's last
    time step) so the Taylor diagram can plot one figure per day too. SSH's
    arbitrary-geoid-reference issue is handled inside
    godae_scorecard_croco_vs_glorys itself (ssh_anomaly=True by default).

    Returns (report, report_by_day): `report` is every day x every
    variable in one DataFrame; `report_by_day[day]` is just that day's rows.
    """
    if not (available and os.path.exists(reference)):
        print("Copernicus Marine Forecast unavailable or REFERENCE missing - skipping GODAE scorecard (Sections 4/6 below will be skipped too).")
        empty = pd.DataFrame(columns=['variable', 'bias', 'rmsd', 'urmsd', 'corr', 'std_ratio', 'vs', 'layer'])
        return empty, {}

    rows = []
    report_by_day = {}
    for day in days:
        day_rows = []
        for var in variables:
            s = godae_scorecard_croco_vs_glorys(croco_his, reference, var, date=day, Yorig=yorig)
            day_rows.append({**s, 'vs': 'reference', 'layer': 'all'})
        report_by_day[day] = pd.DataFrame(day_rows)
        rows.extend(day_rows)
        if print_table:
            print(f"-- {day} --")
            print_scorecard_table(report_by_day[day])

    report = pd.DataFrame(rows)
    return report, report_by_day


def plot_taylor_diagrams_by_day(report_by_day, validation_dir):
    """Section 4 -- one Taylor diagram per day, all of that day's variables
    (temp/ssh/salt/speed) plotted together on the same polar plot. Distinct
    from taylor_diagram() above, which plots one variable across several
    'vs' sources (used by the in-situ Section 8 scorecard) -- here the
    grouping axis is 'variable', not 'vs'. Returns {day: out_png}.
    """
    if not report_by_day:
        print("No GODAE scorecard available (Copernicus Marine Forecast unavailable) - skipping Taylor diagram.")
        return {}

    taylor_figs = {}
    for day, rep in report_by_day.items():
        fig = plt.figure(figsize=(7, 7))
        ax = fig.add_subplot(111, polar=True)
        has_neg = bool((rep['corr'] < 0).any())
        thetamax = 180 if has_neg else 90
        corr_ticks = ([-1.0, -0.5, 0, 0.5, 0.8, 0.9, 0.95, 0.99, 1.0] if has_neg
                      else [0, 0.2, 0.4, 0.6, 0.8, 0.9, 0.95, 0.99, 1.0])
        ax.set_thetamin(0); ax.set_thetamax(thetamax)
        ax.set_xticks(np.arccos(corr_ticks)); ax.set_xticklabels([str(c) for c in corr_ticks])
        ax.set_rlabel_position(0)
        finite = rep['std_ratio'][np.isfinite(rep['std_ratio'])]
        r_max = max(1.6, finite.max() * 1.2) if len(finite) else 1.6
        ax.set_ylim(0, r_max)
        ax.plot(0, 1, 'k*', ms=16, label='reference')
        for _, row in rep.iterrows():
            theta = np.arccos(np.clip(row['corr'], -1, 1))
            ax.plot(theta, row['std_ratio'], 'o', ms=10,
                    label=f"{row['variable']}  (RMSD={row['rmsd']:.2f})")

        ax.text(0.5, -0.08, 'Normalised standard deviation (CROCO / reference)',
                transform=ax.transAxes, ha='center', va='top', fontsize=10)
        ax.text(np.radians(thetamax / 6), r_max * 1.25, 'Correlation coefficient',
                ha='center', va='center', fontsize=10,
                rotation=90 - thetamax / 2, rotation_mode='anchor')

        ax.set_title(f'Taylor diagram -- CROCO vs reference ({day})', pad=20)
        ax.legend(loc='upper left', bbox_to_anchor=(1.05, 1.0), fontsize=9)
        fig.tight_layout()
        out_png = os.path.join(validation_dir, f"taylor_diagram_{day}.png")
        fig.savefig(out_png, dpi=150, bbox_inches="tight")
        plt.close(fig)
        taylor_figs[day] = out_png
    print(f"-> {len(taylor_figs)} figure(s) written, one per day: {list(taylor_figs)}")
    return taylor_figs


def pass_fail_summary(report, cycle_days, stability_ok, spinup_days=2):
    """Section 5 -- automated pass/fail summary, reading directly from the
    GODAE scorecard `report` (build_godae_scorecard()'s first return value)
    so it can't drift out of sync with the Taylor diagram or with
    run_validation.py's own report.

    The first `spinup_days` day(s) of `cycle_days` are excluded from the
    cycle-mean criteria (every other section still shows every day,
    unfiltered). Prints the same pass/fail table, criteria and worst-day
    lines the inline cell used to.

    Returns (all_pass, eval_days) or (None, cycle_days) if report is empty
    (nothing to score -- reference product was unavailable).
    """
    if report.empty:
        print("Overall V1 status: SKIPPED -- Copernicus Marine Forecast unavailable, "
              "no GODAE scorecard to build the pass/fail summary from.")
        return None, cycle_days

    if len(cycle_days) > spinup_days:
        eval_days = cycle_days[spinup_days:]
    else:
        eval_days = cycle_days[:]
        print(f"Cycle only has {len(cycle_days)} day(s) -- shorter than SPINUP_DAYS="
              f"{spinup_days}, so no day can be excluded; using every day.")
    report_eval = report[report['date'].isin(eval_days)]

    by_var = (report_eval.groupby('variable')[['bias', 'rmsd', 'urmsd', 'corr', 'si', 'si_std', 'std_ratio']]
              .mean().to_dict('index'))

    criteria = [
        ('SST cycle-mean domain-avg RMSD < 0.5 degC',       by_var['temp']['rmsd'] < 0.5,  f"{by_var['temp']['rmsd']:.3f} degC"),
        ('SSH cycle-mean spatial correlation > 0.90',       by_var['ssh']['corr'] > 0.90, f"{by_var['ssh']['corr']:.3f}"),
        ('Salinity cycle-mean domain-avg |bias| < 0.2 PSU', abs(by_var['salt']['bias']) < 0.2, f"{by_var['salt']['bias']:+.3f} PSU"),
        ('Numerical stability (no NaN/Inf)',                stability_ok, 'see section 1'),
    ]

    print(f"Cycle-mean scorecard: {len(eval_days)} of {len(cycle_days)} day(s) "
          f"(excluding {cycle_days[:spinup_days] if len(cycle_days) > spinup_days else []} as spin-up): {eval_days}")
    print()
    print(f"{'Criterion':45s} {'Result':6s}  Value")
    print('-' * 70)
    all_pass = True
    for name, passed, value in criteria:
        all_pass &= passed
        print(f"{name:45s} {'PASS' if passed else 'FAIL':6s}  {value}")
    print('-' * 70)

    if len(eval_days) > 1:
        temp_rows = report_eval[report_eval['variable'] == 'temp']
        ssh_rows  = report_eval[report_eval['variable'] == 'ssh']
        salt_rows = report_eval[report_eval['variable'] == 'salt']
        worst_temp = temp_rows.loc[temp_rows['rmsd'].idxmax()]
        worst_ssh  = ssh_rows.loc[ssh_rows['corr'].idxmin()]
        worst_salt = salt_rows.loc[salt_rows['bias'].abs().idxmax()]
        print()
        print(f"Worst single day -- SST RMSD:   {worst_temp['date']}  ({worst_temp['rmsd']:.3f} degC)")
        print(f"Worst single day -- SSH corr:   {worst_ssh['date']}  ({worst_ssh['corr']:.3f})")
        print(f"Worst single day -- Salt bias:  {worst_salt['date']}  ({worst_salt['bias']:+.3f} PSU)")
    print()
    return all_pass, eval_days


def plot_timeseries_vs_forecast(ds, reference, sat_files, avail, depth_m, validation_dir,
                                 available, profile_points=None):
    """Section 6 -- one figure, 4 stacked subpanels (SSH, temp, salt, current
    speed) at one point: CROCO (daily-mean) vs the Copernicus Marine Forecast
    parent, plus OSTIA/ODYSSEA (temp) and SMOS (salt) where DEPTH_M is
    surface and those products were downloaded (Section 1c).

    `profile_points`: same convention as validate_profiles()'s `points` --
    None picks the same auto coastal-ish point (75% across the grid,
    mid-latitude); otherwise the FIRST (lon, lat) pair is used.

    Returns (lon0, lat0) of the point actually plotted, or (None, None) if
    skipped. NOTE: this also fixes a pre-existing bug where a final
    `for ax in axes: ax.legend(...)` loop silently overwrote every panel's
    carefully-placed per-panel legend with default placement/fontsize=8 --
    that generic re-legend call is removed; each panel's own
    `legend(framealpha=0.5, fontsize=10)` call (already present) is now
    the one that actually takes effect.
    """
    if not (available and os.path.exists(reference)):
        print("Copernicus Marine Forecast unavailable or REFERENCE missing - skipping this comparison.")
        return None, None

    clon, clat, cmask = pp.lonlatmask(ds)

    if not profile_points:
        j0 = ds.sizes['eta_rho'] // 2
        i0 = int(0.75 * ds.sizes['xi_rho'])
        lon0, lat0 = float(clon[j0, i0]), float(clat[j0, i0])
    else:
        lon0, lat0 = profile_points[0]

    dlab = 'surface' if depth_m is None else f'{depth_m:g} m'
    print(f'time series point: ({lon0:.2f}, {lat0:.2f})  --  temp/salt at {dlab}')

    ds_daily = _maybe_daily_mean(ds, date=None, daily_mean=True)
    temp_croco = pp.timeseries(ds_daily, 'temp', lon0, lat0, depth_m=depth_m)
    salt_croco = pp.timeseries(ds_daily, 'salt', lon0, lat0, depth_m=depth_m)
    ssh_croco  = pp.timeseries(ds_daily, 'zeta', lon0, lat0)
    speed_croco = pp.timeseries(ds_daily, 'speed', lon0, lat0, depth_m=depth_m)
    ct = temp_croco['time'].values

    dsp = xr.open_dataset(reference)

    def _parent_point_series(cmems_var, depth_m=None):
        da = dsp[cmems_var]
        if "depth" in da.dims:
            da = da.isel(depth=0) if depth_m is None else da.sel(depth=abs(depth_m), method="nearest")
        da = da.sel(longitude=lon0, latitude=lat0, method="nearest")
        return da['time'].values, np.atleast_1d(da.values)

    pt_temp, pv_temp = _parent_point_series('thetao', depth_m=depth_m)
    pt_salt, pv_salt = _parent_point_series('so', depth_m=depth_m)
    pt_ssh,  pv_ssh  = _parent_point_series('zos')
    pt_u, pv_u = _parent_point_series('uo', depth_m=depth_m)
    pt_v, pv_v = _parent_point_series('vo', depth_m=depth_m)
    pv_speed = np.sqrt(np.asarray(pv_u) ** 2 + np.asarray(pv_v) ** 2)
    dsp.close()

    def _satellite_point_series(product, avail_key):
        if depth_m is not None:
            return [], []
        files = sat_files.get(product, {})
        if not avail.get(avail_key, False) or not files:
            return [], []
        days_, vals_ = [], []
        for day, fname in sorted(files.items()):
            if not fname or not os.path.exists(fname):
                continue
            loaded = load_satellite_field(fname, product, date=day)
            if loaded is None:
                continue
            plon, plat, pfield = loaded
            pfield = np.squeeze(np.asarray(pfield))
            if pfield.shape != plon.shape:
                print(f"  [{product}] {day}: field shape {pfield.shape} doesn't match "
                      f"grid shape {plon.shape} - skipping")
                continue
            d2 = (plon - lon0) ** 2 + (plat - lat0) ** 2
            j, i = np.unravel_index(np.argmin(d2), d2.shape)
            v = pfield[j, i]
            if np.isfinite(v):
                days_.append(np.datetime64(day)); vals_.append(v)
        return days_, vals_

    ostia_t, ostia_v = _satellite_point_series('OSTIA', 'ostia_l4')
    odyssea_t, odyssea_v = _satellite_point_series('ODYSSEA', 'odyssea_l3s')
    smos_t, smos_v = _satellite_point_series('SMOS', 'smos_l4_sss')
    if depth_m is not None:
        print(f"DEPTH_M={depth_m:g} m set -- satellite SST/SSS are surface-only products, "
              "so OSTIA/ODYSSEA/SMOS are omitted below (parent/CROCO panels still shown).")

    fig, axes = plt.subplots(4, 1, figsize=(9, 13), sharex=True)

    axes[0].plot(ct, ssh_croco.values, 'o-', color='C0', ms=3, lw=1.4, label='CROCO)')
    axes[0].plot(pt_ssh, pv_ssh, 's--', color='C3', ms=4, lw=1.4, label='parent')
    axes[0].set_ylabel('SSH (m)')
    axes[0].set_title(f'SSH')
    axes[0].legend(framealpha=0.5, fontsize=10)

    axes[1].plot(ct, temp_croco.values, 'o-', color='C0', ms=3, lw=1.4, label='CROCO')
    axes[1].plot(pt_temp, pv_temp, 's--', color='C3', ms=4, lw=1.4, label='parent')
    if ostia_t:
        axes[1].plot(ostia_t, ostia_v, '^:', color='C1', ms=5, lw=1.2, label='OSTIA')
    if odyssea_t:
        axes[1].plot(odyssea_t, odyssea_v, 'v:', color='C2', ms=5, lw=1.2, label='ODYSSEA')
    axes[1].set_ylabel('temperature (degC)')
    axes[1].set_title(f'temperature  ({dlab})')
    axes[1].legend(framealpha=0.5, fontsize=10)

    axes[2].plot(ct, salt_croco.values, 'o-', color='C0', ms=3, lw=1.4, label='CROCO')
    axes[2].plot(pt_salt, pv_salt, 's--', color='C3', ms=4, lw=1.4, label='parent')
    if smos_t:
        axes[2].plot(smos_t, smos_v, '^:', color='C4', ms=5, lw=1.2, label='SMOS')
    axes[2].set_ylabel('salinity (PSU)')
    axes[2].set_title(f'salinity  ({dlab})')
    axes[2].legend(framealpha=0.5, fontsize=10)

    axes[3].plot(ct, speed_croco.values, 'o-', color='C0', ms=3, lw=1.4, label='CROCO')
    axes[3].plot(pt_u, pv_speed, 's--', color='C3', ms=4, lw=1.4, label='parent')
    axes[3].set_ylabel('speed (m s$^{-1}$)')
    axes[3].set_title(f'current speed  ({dlab})')
    axes[3].set_xlabel('time')
    axes[3].legend(framealpha=0.5, fontsize=10)

    plt.suptitle(f"Location: ({lon0:.2f}, {lat0:.2f}) \n")
    for ax in axes:
        ax.grid(alpha=0.3)
    fig.autofmt_xdate()
    fig.tight_layout()
    fig.savefig(os.path.join(validation_dir, f'timeseries_vs_forecast_{lon0:.2f}E_{lat0:.2f}N.png'),
                dpi=300, bbox_inches='tight')
    return lon0, lat0


def stacked_bias_boxplot(ds, reference, days, depth_m, validation_dir, available):
    """Section 6b -- domain-wide CROCO-minus-parent bias boxplot, one stacked
    figure with 4 subpanels (SSH, temp, salt, speed), one box per day. Also
    writes the RMSE-spread companion figure (|CROCO - parent|). Distinct
    from bias_boxplot()/bias_boxplot_multi() above, which each write ONE
    variable per figure -- this keeps the notebook's original "one file for
    all four variables" output.

    Returns (bias_png, rmse_png) or (None, None) if skipped.
    """
    if not (available and os.path.exists(reference)):
        print("Copernicus Marine Forecast unavailable or REFERENCE missing - skipping this comparison.")
        return None, None

    dlab = 'surface' if depth_m is None else f'{depth_m:g} m'

    ssh_diffs   = [domain_diff(ds, reference, 'ssh',   day) for day in days]
    temp_diffs  = [domain_diff(ds, reference, 'temp',  day, depth_m=depth_m) for day in days]
    salt_diffs  = [domain_diff(ds, reference, 'salt',  day, depth_m=depth_m) for day in days]
    speed_diffs = [domain_diff(ds, reference, 'speed', day, depth_m=depth_m) for day in days]

    def _stacked_boxplot(diffs_list, ylabels, titles, out, metric='bias'):
        data = [[np.abs(d) for d in diffs] for diffs in diffs_list] if metric == 'rmse' else diffs_list
        fig, axes = plt.subplots(len(data), 1, figsize=(max(6, 1.1 * len(days)), 3.3 * len(data)),
                                  sharex=True)
        positions = np.arange(1, len(days) + 1)
        for ax, diffs, ylabel, title in zip(axes, data, ylabels, titles):
            ax.boxplot(diffs, positions=positions, showfliers=False)
            ax.set_xticks(positions); ax.set_xticklabels(days)
            if metric == 'bias':
                ax.axhline(0, color='k', ls='--', lw=1)
            ax.set_ylabel(ylabel); ax.set_title(title); ax.grid(alpha=0.3)
        axes[-1].set_xlabel('day')
        plt.setp(axes[-1].get_xticklabels(), rotation=45, ha='right')
        fig.suptitle("CROCO - parent" if metric == 'bias' else "|CROCO - parent|  (RMSE-spread per day)")
        fig.tight_layout()
        fig.savefig(out, dpi=150, bbox_inches='tight')
        plt.close(fig)
        return out

    diffs_all = [ssh_diffs, temp_diffs, salt_diffs, speed_diffs]
    ylabels = ["SSH' bias (m)", 'temperature bias (degC)', 'salinity bias (PSU)', 'speed bias (m s$^{-1}$)']
    titles = ['SSH anomaly bias', f'temperature bias  ({dlab})', f'salinity bias  ({dlab})',
              f'current speed bias  ({dlab})']

    bias_png = _stacked_boxplot(diffs_all, ylabels, titles,
                                 os.path.join(validation_dir, 'boxplot_bias_vs_forecast.png'), metric='bias')

    rmse_ylabels = [y.replace('bias', '|error|') for y in ylabels]
    rmse_titles = [t.replace('bias', 'RMSE-spread') for t in titles]
    rmse_png = _stacked_boxplot(diffs_all, rmse_ylabels, rmse_titles,
                                 os.path.join(validation_dir, 'boxplot_rmse_vs_forecast.png'), metric='rmse')
    return bias_png, rmse_png


def domain_mean_std_plot(ds, reference, days, depth_m, validation_dir, available):
    """Section 6 (full-domain) -- full-domain mean +/- spatial std, CROCO vs
    parent, one stacked figure with 4 subpanels (SSH, temp, salt, speed),
    one point per day. Returns the output path, or None if skipped.
    """
    if not (available and os.path.exists(reference)):
        print("Copernicus Marine Forecast unavailable or REFERENCE missing - skipping this comparison.")
        return None

    dlab = 'surface' if depth_m is None else f'{depth_m:g} m'

    def _domain_mean_std(var, depth_m=None):
        m_c, s_c, m_p, s_p = [], [], [], []
        for day in days:
            croco2d, parent2d = _field_pair_for_day(
                ds, reference, var, day, depth_m, tindex=-1, daily_mean=True, margin_deg=None)
            m_c.append(np.nanmean(croco2d)); s_c.append(np.nanstd(croco2d))
            m_p.append(np.nanmean(parent2d)); s_p.append(np.nanstd(parent2d))
        return map(np.array, (m_c, s_c, m_p, s_p))

    specs = [('ssh', "SSH' (m)", 'SSH anomaly', None),
             ('temp', 'temperature (degC)', f'temperature ({dlab})', depth_m),
             ('salt', 'salinity (PSU)', f'salinity ({dlab})', depth_m),
             ('speed', 'speed (m s$^{-1}$)', f'current speed ({dlab})', depth_m)]

    fig, axes = plt.subplots(len(specs), 1, figsize=(max(6, 1.1 * len(days)), 3.3 * len(specs)),
                              sharex=True)
    x = np.arange(1, len(days) + 1)
    for ax, (var, ylabel, title, dm) in zip(axes, specs):
        m_c, s_c, m_p, s_p = _domain_mean_std(var, depth_m=dm)
        ax.plot(x, m_c, 'o-', color='C0', lw=1.5, ms=5, label='CROCO')
        ax.fill_between(x, m_c - s_c, m_c + s_c, color='C0', alpha=0.2)
        ax.plot(x, m_p, 's--', color='C3', lw=1.5, ms=5, label='parent')
        ax.fill_between(x, m_p - s_p, m_p + s_p, color='C3', alpha=0.2)
        ax.set_xticks(x); ax.set_xticklabels(days)
        ax.set_ylabel(ylabel); ax.set_title(title); ax.grid(alpha=0.3)
    axes[0].legend(fontsize=8)
    axes[-1].set_xlabel('day')
    plt.setp(axes[-1].get_xticklabels(), rotation=45, ha='right')
    fig.suptitle("CROCO vs parent -- full-domain mean +/- spatial std")
    fig.tight_layout()
    out_png = os.path.join(validation_dir, 'timeseries_mean_domain_vs_forecast.png')
    fig.savefig(out_png, dpi=150, bbox_inches='tight')
    plt.close(fig)
    return out_png


def satellite_bias_boxplots(croco_his, sat_files, avail, days, yorig, validation_dir):
    """Section 6c -- CROCO-minus-satellite domain-wide bias boxplot (per
    day), against OSTIA/ODYSSEA (grouped) and SMOS (separate), plus the
    RMSE-spread companions. Skipped (with a message) if no satellite
    product is available/downloaded for this cycle at all.

    Returns dict of the 4 output paths ({'bias_sst':..., 'rmse_sst':...,
    'bias_sss':..., 'rmse_sss':...}), or {} if skipped.
    """
    sat_avail_key = {"OSTIA": "ostia_l4", "ODYSSEA": "odyssea_l3s", "SMOS": "smos_l4_sss"}
    any_sat_avail = any(avail.get(sat_avail_key[p], False) and sat_files.get(p) for p in sat_avail_key)
    if not any_sat_avail:
        print("No satellite product (OSTIA/ODYSSEA/SMOS) available/downloaded for this cycle - skipping.")
        return {}

    def _sat_diff_for_day(product, day):
        files = sat_files.get(product, {})
        fname = files.get(day)
        if not avail.get(sat_avail_key[product], False) or not fname or not os.path.exists(fname):
            return np.array([])
        loaded = load_satellite_field(fname, product, date=day)
        if loaded is None:
            return np.array([])
        plon, plat, pfield = loaded
        pfield = np.squeeze(np.asarray(pfield))
        if pfield.shape != plon.shape:
            return np.array([])
        try:
            clon, clat, croco_field = _croco_field_satellite(croco_his, product, day, Yorig=yorig)
        except Exception:
            return np.array([])
        grid = xr.Dataset(
            {"mask_rho": (("eta_rho", "xi_rho"), np.where(np.isfinite(croco_field), 1.0, np.nan))},
            coords={"lon_rho": (("eta_rho", "xi_rho"), clon), "lat_rho": (("eta_rho", "xi_rho"), clat)})
        sat_on_croco = regrid_to_croco(plon, plat, pfield, grid)
        return domain_diff_satellite(croco_field, sat_on_croco)

    ostia_diffs   = [_sat_diff_for_day('OSTIA',   day) for day in days]
    odyssea_diffs = [_sat_diff_for_day('ODYSSEA', day) for day in days]
    smos_diffs    = [_sat_diff_for_day('SMOS',    day) for day in days]

    sst_groups = {"OSTIA": ostia_diffs, "ODYSSEA": odyssea_diffs}
    sss_groups = {"SMOS": smos_diffs}

    out = {}
    out['bias_sst'] = bias_boxplot_multi(
        sst_groups, days, 'temperature bias (degC)', 'CROCO - satellite SST bias', colors=['C1', 'C2'],
        out=os.path.join(validation_dir, 'boxplot_bias_sst_vs_satellite.png'))
    out['bias_sss'] = bias_boxplot_multi(
        sss_groups, days, 'salinity bias (PSU)', 'CROCO - satellite SSS bias', colors=['C4'],
        out=os.path.join(validation_dir, 'boxplot_bias_sss_vs_satellite.png'))
    out['rmse_sst'] = bias_boxplot_multi(
        sst_groups, days, 'temperature |error| (degC)', 'CROCO - satellite SST RMSE-spread',
        colors=['C1', 'C2'], metric='rmse',
        out=os.path.join(validation_dir, 'boxplot_rmse_sst_vs_satellite.png'))
    out['rmse_sss'] = bias_boxplot_multi(
        sss_groups, days, 'salinity |error| (PSU)', 'CROCO - satellite SSS RMSE-spread',
        colors=['C4'], metric='rmse',
        out=os.path.join(validation_dir, 'boxplot_rmse_sss_vs_satellite.png'))
    return out


def validate_satellite(croco_his, sat_files, product, avail_flag, yorig, validation_dir,
                        stat_label=None):
    """Sections 7/7b -- CROCO forecast vs one satellite product (OSTIA,
    ODYSSEA, or SMOS): one 2x2 figure per day, plus a CSV of daily
    statistics. Guards on `avail_flag` (AVAIL['ostia_l4'] etc.) and on
    whether any files were actually downloaded for this product+cycle.

    `stat_label` controls the printed "<label> daily ... statistics" header
    and the CSV filename ('sst_vs_<product>_stats.csv' if None); pass e.g.
    'SMOS' for the SSS case so the CSV is named 'sss_vs_smos_stats.csv'
    like the original Section 7b cell did.

    Returns the stats DataFrame (possibly empty), or None if skipped.
    """
    files = sat_files.get(product, {})
    if not avail_flag:
        print(f"{product} unavailable on the CMEMS platform - skipping.")
        return None
    if not files:
        print(f"{product}: nothing downloaded for this cycle (no files/coverage) - skipping.")
        return None

    figs, stats = compare_satellite_grid(croco_his, files, product, Yorig=yorig, out_dir=validation_dir)
    print(f"  -> {len(figs)} figure(s) written, one per day: "
          f"{[os.path.basename(p) for p in figs.values()]}")
    if stats is not None and not stats.empty:
        kind = 'SSS' if product == 'SMOS' else 'SST'
        csv_var = 'sss' if product == 'SMOS' else 'sst'
        print(f"\n{product} daily {kind} statistics:")
        print(stats.to_string(index=False))
        stats.to_csv(os.path.join(validation_dir, f"{csv_var}_vs_{product.lower()}_stats.csv"), index=False)
    return stats

    
