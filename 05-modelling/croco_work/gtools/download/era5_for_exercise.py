#!/usr/bin/env python3
"""
download_atmosphere_hindcast_era5.py -- ERA5 (CDS) download + CROCO bulk-forcing
conversion, standalone reconstruction of ggosss26.py's `download_atmosphere_hindcast`
subcommand.

================================================================================
IMPORTANT -- READ BEFORE USING IN PRODUCTION
================================================================================
This is NOT a copy of CROCO's real `download_atmosphere_hindcast` code --
that implementation was never provided to the author of this script. This file
reproduces the DOCUMENTED, OBSERVABLE behaviour of that command (see
Hindcast_setting_config.md, Step 5b) using the standard `cdsapi` client and
well-established ERA5-to-bulk-forcing formulas. Concretely:

  CONFIRMED from the docs (verbatim), reproduced exactly here:
    - CLI shape: --domain "LON_MIN,LON_MAX,LAT_MIN,LAT_MAX" --month_start
      YYYY-MM --month_end YYYY-MM --outputDir DIR --Yorig YYYY
    - the domain passed in is the GRID BOX, not a pre-padded extent -- this
      tool adds its own 2 deg margin before requesting from CDS
    - 10 variables, logged in this order: lsm, sst, tp, strd, ssr, t2m, q,
      u10, v10, msl
    - one NetCDF file per variable per month, real calendar dates, named
      <VARNAME>_Y%YM%m.nc (e.g. T2M_Y2018M06.nc, LSM_Y2018M06.nc), written to
      outputDir/for_croco/
    - hourly records (recordsperday 24)
    - Yorig propagates into each file's `time` units attribute

  RECONSTRUCTED (standard ERA5/bulk-forcing practice, NOT verified against the
  real CROCO pipeline -- check every one of these against a known-good
  file from your existing for_croco/ output before trusting this in production):
    - the actual CDS request variable names behind "sst"/"tp"/"strd"/"ssr"/
      "t2m"/"q"/"u10"/"v10"/"msl"/"lsm" (see VARIABLES below)
    - "q" (specific humidity) is not a native ERA5 field -- it is DERIVED here
      from 2m_dewpoint_temperature + 2m_temperature + mean_sea_level_pressure
      via the Bolton (1980) saturation-vapour-pressure approximation
    - tp/strd/ssr are ERA5 fields ACCUMULATED OVER THE PRECEDING HOUR; this
      script converts them to instantaneous rates/fluxes (kg m-2 s-1, W m-2)
    - unit conventions for tair (K vs degC) and mslp (Pa vs mb) are
      PARAMETERISED (--tair-units, --pressure-units) rather than hard-coded,
      because CROCO's bulk-flux routine (COARE3/3.5/3.6, or the legacy
      Fairall95 scheme) expects different units across CROCO versions --
      confirm which one your croco.in / get_bulk.F build expects before
      running a real hindcast with this output.

  CROCO BULK-FORCING TAG CONVENTION (for_croco/ filenames and internal
  variable names):
    LSM, SST, TP, STRD, SSR, T2M, Q2M, U10M, V10M, MSL
  The `out` tag in VARIABLES below is used for BOTH the filename
  (<out>_Y<YYYY>M<MM>.nc) and the NetCDF variable name inside the file.
  If your CROCO build expects a different tag for any variable (e.g. Q
  instead of Q2M, or PSFC/SLP instead of MSL), change the `out` field only --
  the `short` field is what maps to the CDS request and the internal
  conversion logic, and must not change.

  NEW-CDS-FORMAT HANDLING (added after the 2024 CDS migration):
    - the new API returns a `valid_time` coordinate instead of `time`, with
      `seconds since 1970-01-01` units. Both are handled here.
    - files may carry a length-2 `expver` dimension (ERA5 vs ERA5T merge);
      the first slice is taken, NaNs in the unused slice are dropped.

  IDEMPOTENT RE-RUNS:
    - `--skip-existing` (default ON) skips any for_croco/<VAR>_Y<YYYY>M<MM>.nc
      that already exists, AND skips re-downloading any _raw/ file already
      present. Safe to re-run after a CDS queue timeout or a partial run.
    - `--no-skip-existing` forces a full re-download and re-write (useful if
      you changed a conversion rule and want to regenerate everything).
================================================================================

Usage (matches the documented CLI):

    python download_atmosphere_hindcast_era5.py \\
        --domain "-10.0,10.0,0.0,10.0" \\
        --month_start 2018-06 --month_end 2018-08 \\
        --outputDir ${HCAST}/downloaded_data/ERA5 \\
        --Yorig 1993

Requires: cdsapi, netCDF4, numpy, python-dateutil
    pip install cdsapi netCDF4 numpy python-dateutil --break-system-packages

Requires a CDS API key in ~/.cdsapirc (see
https://cds.climate.copernicus.eu/how-to-api) -- this script does not manage
credentials for you.
"""
from __future__ import annotations

import argparse
import os
import calendar
from datetime import datetime, timedelta

import numpy as np
from netCDF4 import Dataset, date2num, num2date
from dateutil.relativedelta import relativedelta

try:
    import cdsapi
except ImportError as e:
    raise ImportError(
        "cdsapi is required: pip install cdsapi --break-system-packages\n"
        "and configure ~/.cdsapirc -- see https://cds.climate.copernicus.eu/how-to-api"
    ) from e


# ----------------------------------------------------------------------------
# Variable table.
#   short   : the log/CDS short name. Used to select which CDS variable to
#             request and which conversion branch to run. DO NOT change.
#   out     : the CROCO bulk-forcing tag. Used for BOTH the output filename
#             (<out>_Y<YYYY>M<MM>.nc) and the NetCDF variable name inside
#             the file. Change this if your CROCO build expects a different
#             tag (e.g. Q instead of Q2M, or PSFC instead of MSL).
#   cds_var : the actual CDS/ERA5 request variable name, or None for "q",
#             which is derived rather than requested directly.
#   kind    : 'instant' (value at the timestamp) or 'accum' (accumulated over
#             the preceding hour, needs converting to a rate/flux).
# ----------------------------------------------------------------------------
VARIABLES = [
    {"short": "lsm",  "out": "LSM",  "cds_var": "land_sea_mask",                       "kind": "instant"},
    {"short": "sst",  "out": "SST",  "cds_var": "sea_surface_temperature",             "kind": "instant"},
    {"short": "tp",   "out": "TP",   "cds_var": "total_precipitation",                 "kind": "accum"},
    {"short": "strd", "out": "STRD", "cds_var": "surface_thermal_radiation_downwards", "kind": "accum"},
    {"short": "ssr",  "out": "SSR",  "cds_var": "surface_net_solar_radiation",         "kind": "accum"},
    {"short": "t2m",  "out": "T2M",  "cds_var": "2m_temperature",                      "kind": "instant"},
    {"short": "q",    "out": "Q2M",  "cds_var": None,                                  "kind": "derived"},
    {"short": "u10",  "out": "U10M", "cds_var": "10m_u_component_of_wind",             "kind": "instant"},
    {"short": "v10",  "out": "V10M", "cds_var": "10m_v_component_of_wind",             "kind": "instant"},
    {"short": "msl",  "out": "MSL",  "cds_var": "mean_sea_level_pressure",             "kind": "instant"},
]
# 'q' is derived from these two extra raw downloads (not written to their own
# for_croco/ file -- only the derived Q2M_*.nc is kept, matching the 10-file
# count "10 vars x 3 months" the docs' CHECK step expects).
_HUMIDITY_INPUTS = {
    "d2m": "2m_dewpoint_temperature",
}

CDS_DATASET = "reanalysis-era5-single-levels"
DOMAIN_MARGIN_DEG = 2.0   # "adds its own 2 deg margin" -- Hindcast_setting_config.md 5b

# Coordinate names that are NEVER data variables -- excluded from the
# "pick the only non-coordinate variable" fallback in _read().
_COORD_NAMES = {"longitude", "latitude", "time", "valid_time", "expver", "number"}


def month_range(month_start: str, month_end: str):
    """Yield (year, month) for every calendar month from month_start to
    month_end inclusive, both 'YYYY-MM'."""
    start = datetime.strptime(month_start, "%Y-%m")
    end = datetime.strptime(month_end, "%Y-%m")
    cur = start
    while cur <= end:
        yield cur.year, cur.month
        cur += relativedelta(months=1)


def cds_area(domain: str, margin_deg: float = DOMAIN_MARGIN_DEG):
    """'LON_MIN,LON_MAX,LAT_MIN,LAT_MAX' (the grid box) -> CDS 'area'
    [north, west, south, east], padded by margin_deg on every side.
    ERA5/CDS area values must be within [-180, 180] for longitude and
    [-90, 90] for latitude; this does not handle a domain that straddles the
    antimeridian (+/-180), only the prime meridian (0 deg), consistent with
    ERA5's native -180..180 longitude convention noted in the docs.
    """
    lon_min, lon_max, lat_min, lat_max = (float(x) for x in domain.split(","))
    north = min(90.0, lat_max + margin_deg)
    south = max(-90.0, lat_min - margin_deg)
    west = max(-180.0, lon_min - margin_deg)
    east = min(180.0, lon_max + margin_deg)
    return [north, west, south, east]


def _cds_request_dict(cds_var: str, year: int, month: int, area, time_step_hours: int = 1):
    n_days = calendar.monthrange(year, month)[1]
    return {
        "product_type": "reanalysis",
        "format": "netcdf",
        "variable": cds_var,
        "year": f"{year:04d}",
        "month": f"{month:02d}",
        "day": [f"{d:02d}" for d in range(1, n_days + 1)],
        "time": [f"{h:02d}:00" for h in range(0, 24, time_step_hours)],
        "area": area,   # [north, west, south, east]
    }


def _for_croco_path(output_dir: str, out_tag: str, year: int, month: int) -> str:
    """Path of the final per-variable for_croco NetCDF for a given month."""
    return os.path.join(output_dir, "for_croco", f"{out_tag}_Y{year:04d}M{month:02d}.nc")


def download_raw(client, cds_var: str, year: int, month: int, area, tmp_dir: str,
                 skip_existing: bool = True) -> str:
    """Download one ERA5 variable for one month to a temp file. Returns the
    local path. Skips re-downloading if the file already exists when
    skip_existing is True (safe to re-run after a CDS queue timeout -- see
    the docs' "this is slow, leave it running" note)."""
    os.makedirs(tmp_dir, exist_ok=True)
    out_path = os.path.join(tmp_dir, f"_raw_{cds_var}_{year:04d}{month:02d}.nc")
    if skip_existing and os.path.exists(out_path):
        print(f"  {cds_var} {year:04d}-{month:02d}: raw already present, skipping CDS request")
        return out_path
    print(f"  requesting {cds_var} {year:04d}-{month:02d} from CDS ...")
    client.retrieve(CDS_DATASET, _cds_request_dict(cds_var, year, month, area), out_path)
    return out_path


# ----------------------------------------------------------------------------
# Conversions -- RECONSTRUCTED, see module docstring. Verify against your
# real pipeline's output before production use.
# ----------------------------------------------------------------------------

def accum_to_rate(values_j_or_m: np.ndarray, kind: str, step_seconds: int = 3600) -> np.ndarray:
    """ERA5 hourly-reanalysis accumulated fields (tp, strd, ssr) are
    accumulated over the preceding `step_seconds` (3600 for hourly output).
    Converts to an instantaneous rate/flux:
      - 'precip' : tp is in metres of water -> kg m-2 s-1 (rho_water = 1000 kg/m3)
      - 'radiation' : strd/ssr are in J m-2 -> W m-2
    """
    if kind == "precip":
        return values_j_or_m * 1000.0 / step_seconds
    elif kind == "radiation":
        return values_j_or_m / step_seconds
    raise ValueError(f"unknown accum kind: {kind}")


def specific_humidity_from_dewpoint(d2m_k: np.ndarray, msl_pa: np.ndarray) -> np.ndarray:
    """Specific humidity (kg/kg) at 2 m from dewpoint temperature and mean
    sea level pressure, via Bolton (1980)'s saturation-vapour-pressure
    approximation (accurate to ~0.1% over -35..+35 degC) evaluated AT the
    dewpoint (so it directly gives actual, not saturation, vapour pressure).

        es(Td) = 611.2 * exp(17.67 * (Td_degC) / (Td_degC + 243.5))   [Pa]
        q = 0.622 * es / (p - 0.378 * es)                             [kg/kg]

    d2m_k : 2 m dewpoint temperature, Kelvin
    msl_pa: mean sea level pressure, Pa (used as a surface-pressure stand-in,
            standard practice for ERA5-derived near-surface humidity; if your
            CROCO bulk routine expects humidity referenced to actual surface
            pressure rather than MSLP, adjust here)
    """
    td_c = d2m_k - 273.15
    es_pa = 611.2 * np.exp(17.67 * td_c / (td_c + 243.5))
    q = 0.622 * es_pa / (msl_pa - 0.378 * es_pa)
    return q


def kelvin_to_celsius(t_k: np.ndarray) -> np.ndarray:
    return t_k - 273.15


def pa_to_mb(p_pa: np.ndarray) -> np.ndarray:
    return p_pa / 100.0


# ----------------------------------------------------------------------------
# NetCDF writer -- matches the observed <VAR>_Y%YM%m.nc convention and the
# "Yorig propagates into the time units" behaviour from the docs.
# ----------------------------------------------------------------------------

def write_for_croco(out_path: str, var_name: str, long_name: str, units: str,
                     lon: np.ndarray, lat: np.ndarray, time_hours: np.ndarray,
                     yorig: int, data: np.ndarray):
    """Write one variable's monthly file in a minimal CROCO 'online' bulk-
    forcing layout: dims (time, lat, lon), one data variable named `var_name`
    (which is the CROCO tag, e.g. U10M), a `time` coordinate in hours since
    Yorig-01-01 00:00:00 (so `ncdump -h` shows
    `units = "hours since {yorig}-01-01 00:00:00"`, the same check the docs
    use to confirm Yorig propagated correctly).
    """
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    with Dataset(out_path, "w", format="NETCDF4_CLASSIC") as nc:
        nc.createDimension("lon", lon.size)
        nc.createDimension("lat", lat.size)
        nc.createDimension("time", None)

        v_lon = nc.createVariable("lon", "f8", ("lon",))
        v_lon[:] = lon
        v_lon.units = "degrees_east"

        v_lat = nc.createVariable("lat", "f8", ("lat",))
        v_lat[:] = lat
        v_lat.units = "degrees_north"

        v_time = nc.createVariable("time", "f8", ("time",))
        v_time[:] = time_hours
        v_time.units = f"hours since {yorig:04d}-01-01 00:00:00"
        v_time.calendar = "standard"

        v_data = nc.createVariable(var_name, "f4", ("time", "lat", "lon"),
                                   zlib=True, complevel=4)
        v_data[:] = data
        v_data.units = units
        v_data.long_name = long_name

        nc.history = (f"Generated by download_atmosphere_hindcast_era5.py "
                      f"(reconstruction, see module docstring) on "
                      f"{datetime.utcnow().isoformat()}Z")


# ----------------------------------------------------------------------------
# Main per-month conversion pipeline
# ----------------------------------------------------------------------------

def process_month(client, year: int, month: int, area, output_dir: str,
                   yorig: int, tair_units: str, pressure_units: str,
                   skip_existing: bool = True):
    for_croco_dir = os.path.join(output_dir, "for_croco")
    tmp_dir = os.path.join(output_dir, "_raw")

    # 0) If every final for_croco file for this month already exists, skip the
    #    whole month (no downloads, no conversions).
    if skip_existing:
        missing = [
            v for v in VARIABLES
            if not os.path.exists(_for_croco_path(output_dir, v["out"], year, month))
        ]
        if not missing:
            print(f"  {year:04d}-{month:02d}: all {len(VARIABLES)} for_croco files "
                 f"already present, skipping month entirely")
            return
        if len(missing) < len(VARIABLES):
            print(f"  {year:04d}-{month:02d}: {len(missing)}/{len(VARIABLES)} "
                 f"for_croco files missing -> will fill in: "
                 f"{', '.join(v['out'] for v in missing)}")

    # 1) download every raw variable this month needs (9 direct + d2m for humidity).
    #    Only request raw files that a still-missing output variable depends on.
    if skip_existing:
        needed_shorts = {
            v["short"] for v in VARIABLES
            if not os.path.exists(_for_croco_path(output_dir, v["out"], year, month))
        }
    else:
        needed_shorts = {v["short"] for v in VARIABLES}

    raw_paths = {}
    for v in VARIABLES:
        if v["cds_var"] is None:
            continue
        # 'msl' is also needed to derive 'q', so fetch it if either is missing.
        need_this = (v["short"] in needed_shorts) or (
            v["short"] == "msl" and "q" in needed_shorts
        )
        if not need_this:
            continue
        raw_paths[v["short"]] = download_raw(
            client, v["cds_var"], year, month, area, tmp_dir, skip_existing
        )
    if "q" in needed_shorts:
        raw_paths["d2m"] = download_raw(
            client, _HUMIDITY_INPUTS["d2m"], year, month, area, tmp_dir, skip_existing
        )

    # 2) read coordinates + time from one available file (all requests share the
    #    same area/timestamps, so any file's grid/time applies to all of them).
    #    New-CDS files use `valid_time` with `seconds since 1970-01-01`
    #    instead of `time` with `hours since 1900-01-01` -- decode from
    #    whatever units the file actually declares.
    anchor = raw_paths.get("t2m") or next(iter(raw_paths.values()))
    with Dataset(anchor) as nc0:
        lon = nc0.variables["longitude"][:].astype("f8")
        lat = nc0.variables["latitude"][:].astype("f8")

        tname = "time" if "time" in nc0.variables else "valid_time"
        tvar = nc0.variables[tname]
        t_units = getattr(tvar, "units", "hours since 1900-01-01 00:00:00")
        t_cal = getattr(tvar, "calendar", "standard")
        base_dates = np.array([
            num2date(float(h), t_units, t_cal) for h in tvar[:]
        ])

    time_hours_since_yorig = np.array([
        (d - datetime(yorig, 1, 1)).total_seconds() / 3600.0 for d in base_dates
    ])

    def _read(short):
        """Read one raw file's data array as float64, collapsing the new-CDS
        `expver` merge dimension (ERA5 vs ERA5T) if present, and falling back
        to the only non-coordinate variable if the short name isn't a direct
        match."""
        with Dataset(raw_paths[short]) as nc:
            if short in nc.variables:
                name = short
            else:
                candidates = [k for k in nc.variables if k not in _COORD_NAMES]
                if not candidates:
                    raise KeyError(
                        f"no data variable found in {raw_paths[short]}; "
                        f"variables are {list(nc.variables)}"
                    )
                name = candidates[0]
            v = nc.variables[name]
            arr = v[:].astype("f8")
            if "expver" in v.dimensions:
                arr = arr.take(0, axis=v.dimensions.index("expver"))
        return np.squeeze(arr)

    print(f"  converting + writing for_croco/ files for {year:04d}-{month:02d} ...")

    for v in VARIABLES:
        tag = f"Y{year:04d}M{month:02d}"
        out_path = os.path.join(for_croco_dir, f"{v['out']}_{tag}.nc")

        if skip_existing and os.path.exists(out_path):
            print(f"    {v['out']:5s} -> already exists, skipping")
            continue

        if v["short"] == "q":
            d2m = _read("d2m")
            msl = _read("msl")
            data = specific_humidity_from_dewpoint(d2m, msl)
            units = "kg kg-1"
            long_name = "2 metre specific humidity (derived from dewpoint)"

        elif v["short"] in ("tp", "strd", "ssr"):
            raw = _read(v["short"])
            kind = "precip" if v["short"] == "tp" else "radiation"
            data = accum_to_rate(raw, kind)
            units = "kg m-2 s-1" if kind == "precip" else "W m-2"
            long_name = {
                "tp": "Total precipitation rate",
                "strd": "Surface thermal (longwave) radiation downwards",
                "ssr": "Surface net solar (shortwave) radiation",
            }[v["short"]]

        elif v["short"] == "t2m":
            raw = _read("t2m")
            if tair_units == "celsius":
                data = kelvin_to_celsius(raw)
                units = "degC"
            else:
                data = raw
                units = "K"
            long_name = "2 metre temperature"

        elif v["short"] == "msl":
            raw = _read("msl")
            if pressure_units == "mb":
                data = pa_to_mb(raw)
                units = "mb"
            else:
                data = raw
                units = "Pa"
            long_name = "Mean sea level pressure"

        else:
            data = _read(v["short"])
            units = {
                "lsm": "1 (0-1 land fraction)",
                "sst": "K",
                "u10": "m s-1",
                "v10": "m s-1",
            }[v["short"]]
            long_name = {
                "lsm": "Land-sea mask",
                "sst": "Sea surface temperature",
                "u10": "10 metre U wind component",
                "v10": "10 metre V wind component",
            }[v["short"]]

        write_for_croco(out_path, v["out"], long_name, units, lon, lat,
                        time_hours_since_yorig, yorig, data)
        print(f"    {v['out']:5s} -> {out_path}")

    print(f"  {year:04d}-{month:02d} done "
         f"(LSM, SST, TP, STRD, SSR, T2M, Q2M, U10M, V10M, MSL)")


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--domain", required=True,
                  help="'LON_MIN,LON_MAX,LAT_MIN,LAT_MAX' -- the grid box; "
                       "a 2 deg margin is added automatically before the CDS request")
    p.add_argument("--month_start", required=True, help="YYYY-MM, inclusive")
    p.add_argument("--month_end", required=True, help="YYYY-MM, inclusive")
    p.add_argument("--outputDir", required=True,
                  help="ERA5 files land in <outputDir>/for_croco/")
    p.add_argument("--Yorig", type=int, required=True,
                  help="time-origin year written into every output file's "
                       "time units (e.g. 1993 for a GLORYS-referenced hindcast)")
    p.add_argument("--tair-units", choices=["kelvin", "celsius"], default="kelvin",
                  help="unit for T2M output -- VERIFY against your CROCO bulk-flux "
                       "scheme's expected convention before production use (default: "
                       "kelvin, ERA5's native unit)")
    p.add_argument("--pressure-units", choices=["pa", "mb"], default="pa",
                  help="unit for MSL output -- VERIFY against your CROCO bulk-flux "
                       "scheme's expected convention before production use (default: "
                       "pa, ERA5's native unit)")
    p.add_argument("--skip-existing", dest="skip_existing", action="store_true",
                  default=True,
                  help="skip any for_croco/<VAR>_Y<YYYY>M<MM>.nc that already "
                       "exists, and skip re-downloading raw files already in "
                       "_raw/ (default: ON; makes re-runs idempotent)")
    p.add_argument("--no-skip-existing", dest="skip_existing", action="store_false",
                  help="force a full re-download and re-write, ignoring any "
                       "existing output files")
    args = p.parse_args()

    area = cds_area(args.domain)
    print(f"ERA5 domain (grid box + {DOMAIN_MARGIN_DEG:g} deg margin), "
         f"CDS area [N,W,S,E]: {area}")
    print(f"skip-existing: {args.skip_existing}")

    client = cdsapi.Client()

    months = list(month_range(args.month_start, args.month_end))
    print(f"Downloading {len(months)} month(s): "
         f"{args.month_start} -> {args.month_end}")
    print("This is slow -- CDS queues requests, sometimes for a long time. "
         "Start it and leave it.")

    for year, month in months:
        process_month(client, year, month, area, args.outputDir, args.Yorig,
                      args.tair_units, args.pressure_units,
                      skip_existing=args.skip_existing)

    for_croco_dir = os.path.join(args.outputDir, "for_croco")
    n_files = len([f for f in os.listdir(for_croco_dir) if f.endswith(".nc")]) \
        if os.path.isdir(for_croco_dir) else 0
    print(f"\nERA5 files conversion done.")
    print(f"{n_files} file(s) in {for_croco_dir}  "
         f"(expected {len(VARIABLES)} vars x {len(months)} months = "
         f"{len(VARIABLES) * len(months)})")


if __name__ == "__main__":
    main()

    