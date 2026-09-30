"""
gtools/validation_obs.py — observation references for CROCO validation.

`gtools.validation` compares a CROCO run against the parent product it was
downscaled from. That is a **consistency** check: the parent supplied the
initial and boundary conditions, so agreement is partly guaranteed.

This module adds **independent** references — satellite and observation-based
products from CMEMS — and one entry point that works with any of them:

    import gtools.validation_obs as vo

    ost = vo.download_obs(HIS, "ostia", "~/ggosss26/data/OBS", Yorig=2000)
    vo.compare(HIS, ost, "temp", Yorig=2000)                  # one map
    vo.compare_days(HIS, ost, "temp", days=5, Yorig=2000)     # a day per column
    vo.compare(HIS, ost, "temp", kind="growth", Yorig=2000)   # RMSE vs lead day

The references
--------------
    ostia         OSTIA L4 SST, 0.05deg, gap-free            -> temp (surface)
    duacs         DUACS L4 altimetry, 0.125deg               -> ssh, geostrophic uv
    globcurrent   total surface current, 0.25deg             -> u, v (total)
    armor3d       ARMOR3D L4, 1/8deg, 50 levels              -> temp, salt, ssh, mld
    glorys        the parent product (consistency check)     -> temp, salt, ssh, u, v

Which to use for what
---------------------
    SST map, growth           ostia        finest, gap-free, near-independent
    surface currents          globcurrent  TOTAL current, not geostrophic only
    SSH                       duacs        but see the note on datums below
    profile, section, depth   armor3d      observation-based T/S at depth
    point skill by layer      in-situ      gtools.validation_godae
    "did the downscaling      glorys       consistency, not skill
     stay consistent?"

Every product here is an analysis or a reconstruction, not raw observations.
They carry their own error structure, and none resolves what a 1/12deg model
resolves. Agreement is evidence, not proof; disagreement at scales the
reference cannot see is not necessarily model error.

Two things this module is strict about
--------------------------------------
**Dates.** A CROCO file written without CF time units carries raw seconds, not
dates. Every function here needs real dates — to pick the matching reference
record and to know which days to download — so they raise if the time axis was
not decoded. Pass `Yorig`: 2000 for the forecast track, 1993 for hindcasts.

**Matching.** The reference record is chosen by nearest date, and the gap is
checked. If the nearest record is further away than `max_gap_days` the
comparison raises rather than quietly pairing a Tuesday model field with the
previous Thursday's observation.
"""
from __future__ import annotations

import os
import re
import subprocess

import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap

import gtools.postprocess as pp
import gtools.validation as val


# ======================================================================
# The reference registry
#
# `nrt` and `my` give the near-real-time and multi-year twins where both
# exist — forecasts use nrt, hindcasts use my, the same split as
# Mercator/GLORYS.
#
# `vars` maps our short names to the product's own variable names.
# `offset` and `scale` convert to our units (celsius, metres, m/s).
# `has_depth` says whether depth selection applies.
# ======================================================================
REFERENCES = {
    "ostia": dict(
        nrt="METOFFICE-GLO-SST-L4-NRT-OBS-SST-V2",
        my="METOFFICE-GLO-SST-L4-REP-OBS-SST",
        vars={"temp": "analysed_sst"},
        offset=-273.15,                 # OSTIA ships kelvin
        scale=1.0,
        has_depth=False,
        note="foundation SST — the temperature below the diurnal warm layer, "
             "so a midday model SST can read warmer without being wrong",
    ),
    "odyssea": dict(
        nrt="IFREMER-GLOB-SST-L3-NRT-OBS_FULL_TIME_SERIE",
        my="cmems_obs-sst_glo_phy_my_l3s_P1D-m",
        vars={"temp": "adjusted_sea_surface_temperature"},
        qc={"var": "quality_level", "min": 5},
        offset=-273.15,
        scale=1.0,
        has_depth=False,
        gappy=True,
        note="merged satellite observations, not an analysis: cloud-covered "
             "cells are missing, so coverage varies day to day. More "
             "independent than OSTIA, which assimilates in-situ data \u2014 but "
             "regridding fills those gaps, so check coverage before trusting "
             "a comparison (see the gappy warning in load_reference)",
    ),
    "duacs": dict(
        nrt="cmems_obs-sl_glo_phy-ssh_nrt_allsat-l4-duacs-0.125deg_P1D",
        my=None,
        vars={"ssh": "sla", "u": "ugos", "v": "vgos"},
        offset=0.0, scale=1.0, has_depth=False,
        note="SLA is an anomaly about a mean sea surface, and the geostrophic "
             "velocities exclude Ekman flow — use globcurrent for currents",
    ),
    "globcurrent": dict(
        nrt="cmems_obs-mob_glo_phy-cur_nrt_0.25deg_P1D-m",
        my="cmems_obs-mob_glo_phy-cur_my_0.25deg_P1D-m",
        vars={"u": "uo", "v": "vo"},
        offset=0.0, scale=1.0, has_depth=False,
        note="total surface current: geostrophic + Ekman, which is what CROCO "
             "produces — unlike DUACS, which is geostrophic only",
    ),
    "armor3d": dict(
        nrt="cmems_obs-mob_glo_phy_nrt_0.125deg_P1D-m",
        my="cmems_obs-mob_glo_phy_my_0.125deg_P1D-m",
        vars={"temp": "to", "salt": "so", "ssh": "zo", "mld": "mlotst"},
        offset=0.0, scale=1.0, has_depth=True,
        note="a statistical reconstruction: altimetry projected downward using "
             "Argo-derived covariances. Real skill in the upper ocean where "
             "Argo is dense, less below",
    ),
    "glorys": dict(
        nrt=None, my=None,              # downloaded by ggosss26.py download_ocean
        vars=dict(val.PARENT_VARS),     # temp/ssh/salt/u/v -> thetao/zos/so/uo/vo
        offset=0.0, scale=1.0, has_depth=True,
        note="the product that supplied your boundaries — a consistency check, "
             "not an independent one",
    ),
}

_SOURCE_NAMES = {
    "ostia": "OSTIA", "odyssea": "ODYSSEA", "duacs": "DUACS", "globcurrent": "GlobCurrent",
    "armor3d": "ARMOR3D", "glorys": "GLORYS", "mercator": "Mercator",
}


def source_name(source):
    """A display name for a reference key, for titles."""
    return _SOURCE_NAMES.get(source, source.upper())


# ======================================================================
# Guards
# ======================================================================
def _require_dates(ds, croco_his):
    """Raise unless the CROCO time axis was decoded to real dates.

    Without this the failure surfaces much later and much less clearly: a
    resample raises deep inside pandas, or a download silently requests 1970.
    """
    if not np.issubdtype(ds["time"].dtype, np.datetime64):
        dtype = ds["time"].dtype
        ds.close()
        raise ValueError(
            f"{croco_his}: the time axis was not decoded to real dates "
            f"(dtype={dtype}). The file has no CF 'units' attribute on time, "
            "so a reference year is needed — pass Yorig=2000 for the forecast "
            "track, Yorig=1993 for hindcasts.")


# ======================================================================
# Working out a run's window and domain
# ======================================================================
def run_window(croco_his, Yorig=None, pad_days=1):
    """The date range a CROCO run covers, as ('YYYY-MM-DD', 'YYYY-MM-DD').

    Observation products are daily, so the window is widened to whole days.
    `pad_days` adds a margin on each side — worth having, because a daily mean
    is usually centred at 12:00 while a run starts at 00:00, so the first model
    record needs the previous day's field to bracket it.
    """
    ds = pp.open_history(croco_his, Yorig=Yorig)
    _require_dates(ds, croco_his)
    t = pp.times(ds)
    ds.close()
    t0 = np.datetime64(t[0], "D") - np.timedelta64(pad_days, "D")
    t1 = np.datetime64(t[-1], "D") + np.timedelta64(pad_days, "D")
    return str(t0), str(t1)


def run_domain(croco_his, Yorig=None, margin_deg=0.5):
    """The lon/lat box a CROCO run covers, widened by `margin_deg`.

    The margin matters: regridding onto the CROCO grid interpolates, and
    without reference points outside the model's edge the outermost cells come
    back NaN.
    """
    ds = pp.open_history(croco_his, Yorig=Yorig)
    lon, lat, _ = pp.lonlatmask(ds)
    ds.close()
    return (float(lon.min()) - margin_deg, float(lon.max()) + margin_deg,
            float(lat.min()) - margin_deg, float(lat.max()) + margin_deg)


# ======================================================================
# Download
# ======================================================================
def download_obs(croco_his, source, outdir, Yorig=None, track="nrt",
                 pad_days=1, margin_deg=0.5, force=False, verbose=True):
    """Fetch the observation product covering a run's window and domain.

    croco_his : the run to be validated — its dates and box drive the request
    source    : a key of REFERENCES ('ostia', 'duacs', 'globcurrent', 'armor3d')
    outdir    : where to write; created if absent
    track     : 'nrt' for forecasts, 'my' for hindcasts
    force     : re-download even if the file is already there

    Returns the path to the downloaded file.
    """
    if source not in REFERENCES:
        raise ValueError(f"source must be one of {sorted(REFERENCES)}, got {source!r}")
    spec = REFERENCES[source]
    dataset = spec.get(track)
    if dataset is None:
        other = "my" if track == "nrt" else "nrt"
        raise ValueError(
            f"{source!r} has no {track!r} dataset. "
            + (f"Try track={other!r}." if spec.get(other) else
               "It is not downloaded through this function — see its note."))

    start, end = run_window(croco_his, Yorig=Yorig, pad_days=pad_days)
    lo0, lo1, la0, la1 = run_domain(croco_his, Yorig=Yorig, margin_deg=margin_deg)

    outdir = os.path.expanduser(outdir)
    os.makedirs(outdir, exist_ok=True)
    fname = f"{source}_{start}_{end}.nc"
    path = os.path.join(outdir, fname)

    if os.path.exists(path) and not force:
        if verbose:
            print(f"{path} already exists — pass force=True to re-download")
        return path
    if os.path.exists(path):
        # copernicusmarine writes name_(1).nc rather than overwriting, which
        # leaves you reading the stale file while believing you re-downloaded
        os.remove(path)

    cmd = ["copernicusmarine", "subset",
           "--dataset-id", dataset,
           "--start-datetime", start, "--end-datetime", end,
           "--minimum-longitude", f"{lo0:.4f}", "--maximum-longitude", f"{lo1:.4f}",
           "--minimum-latitude", f"{la0:.4f}", "--maximum-latitude", f"{la1:.4f}",
           "--output-directory", outdir, "--output-filename", fname]
    for v in spec["vars"].values():
        cmd += ["--variable", v]
    if spec.get("qc"):                      # the quality flag, for filtering
        cmd += ["--variable", spec["qc"]["var"]]

    if verbose:
        print(f"{source}: {dataset}")
        print(f"  {start} .. {end}   lon {lo0:.2f}..{lo1:.2f}  lat {la0:.2f}..{la1:.2f}")

    res = subprocess.run(cmd, capture_output=True, text=True)
    if res.returncode != 0:
        # Check if the failure is due to out-of-bounds dates (e.g. ARMOR3D lagging by a few days)
        match = re.search(r"exceed the dataset coordinates \[[^,]+,\s*([\d-]+)", res.stderr + res.stdout)
        if match:
            max_avail = match.group(1)
            if start > max_avail:
                # The entire requested window is in the future relative to dataset
                # Shift the window back so we still download recent valid data
                start_dt = str(np.datetime64(max_avail, "D") - np.timedelta64(7, "D"))
                end_dt = max_avail
                print(f"  WARNING: {source} has no data for this run's window — it "
                      f"ends at {max_avail}. Downloading {start_dt}..{end_dt} "
                      f"instead, which does NOT overlap the run, so comparing "
                      f"against it will raise. ARMOR3D updates weekly.")
            else:
                start_dt = start
                end_dt = max_avail
            if verbose:
                print(f"  note: adjusting {source} date range to available data ({start_dt} .. {end_dt})")
            fname = f"{source}_{start_dt}_{end_dt}.nc"
            path = os.path.join(outdir, fname)
            cmd_adj = list(cmd)
            cmd_adj[cmd_adj.index("--start-datetime") + 1] = start_dt
            cmd_adj[cmd_adj.index("--end-datetime") + 1] = end_dt
            cmd_adj[cmd_adj.index("--output-filename") + 1] = fname
            res2 = subprocess.run(cmd_adj, capture_output=True, text=True)
            if res2.returncode != 0:
                subprocess.run(cmd_adj + ["--force-download"], check=True)
        else:
            # copernicusmarine <2.0 prompted before downloading
            subprocess.run(cmd + ["--force-download"], check=True)

    if verbose and spec.get("note"):
        print(f"  note: {spec['note']}")
    return path


# ======================================================================
# Identifying a reference file
# ======================================================================
def identify(reference):
    """Work out which product a file is: an entry in REFERENCES, or 'croco'.

    Recognised by the variables present, so a file downloaded by hand works
    as well as one from download_obs().
    """
    if isinstance(reference, (list, tuple)):
        return "insitu"
    if isinstance(reference, str) and any(c in reference for c in "*?["):
        return "insitu"

    try:
        d = xr.open_dataset(reference, decode_times=False)
    except Exception:
        return "unknown"
    try:
        names = set(d.variables)
        if "xi_rho" in d.dims:
            return "croco"
        if "analysed_sst" in names:
            return "ostia"
        if "adjusted_sea_surface_temperature" in names:
            return "odyssea"
        if "sla" in names or "adt" in names:
            return "duacs"
        if {"to", "so"} <= names:
            return "armor3d"
        if {"thetao", "so"} <= names:
            return "glorys"
        if {"uo", "vo"} <= names and not ({"thetao", "so"} & names):
            return "globcurrent"
        if {"uo", "vo"} <= names:
            return "glorys"
    finally:
        d.close()
    return "unknown"


# ======================================================================
# Loading a field from any reference, onto the CROCO grid
# ======================================================================
def load_reference(reference, var, croco_ds, date=None, depth_m=None,
                   source=None, max_gap_days=0.5):
    """Read `var` from any supported reference and regrid it onto CROCO's grid.

    The reference record is the one nearest `date`, and the gap is checked:
    if the nearest is further than `max_gap_days` this raises, rather than
    silently pairing the model with an observation from another week.

    Returns a 2D array on (eta_rho, xi_rho), land NaN. Units are converted to
    ours: celsius, metres, m/s.
    """
    source = source or identify(reference)
    if source in ("unknown", "insitu", "croco"):
        raise ValueError(
            f"load_reference cannot handle a {source!r} reference. "
            "CROCO-vs-CROCO goes through gtools.validation.compare_resolution; "
            "in-situ goes through gtools.validation_godae.")

    spec = REFERENCES[source]
    if var == "speed":
        u = load_reference(reference, "u", croco_ds, date, depth_m, source,
                           max_gap_days)
        v = load_reference(reference, "v", croco_ds, date, depth_m, source,
                           max_gap_days)
        return pp.speed(u, v)

    if var not in spec["vars"]:
        raise KeyError(
            f"{source!r} has no {var!r}. It provides "
            f"{sorted(spec['vars'])}." +
            ("  Use armor3d for subsurface temperature and salinity."
             if var in ("temp", "salt") and not spec["has_depth"] else ""))
    if depth_m is not None and not spec["has_depth"]:
        raise ValueError(f"{source!r} is a surface product; depth_m does not apply")

    name = spec["vars"][var]
    d = xr.open_dataset(reference)
    if name not in d:
        d.close()
        raise KeyError(f"{name!r} not in {reference}")
    da = d[name]

    if "time" in da.dims:
        if date is None:
            da = da.isel(time=0)
        else:
            want = np.datetime64(date)
            da = da.sel(time=want, method="nearest")
            got = np.datetime64(da["time"].values, "D")
            gap = abs(float((got - np.datetime64(want, "D"))
                            / np.timedelta64(1, "D")))
            if gap > max_gap_days:
                first = np.datetime_as_string(d["time"].values[0], unit="D")
                last = np.datetime_as_string(d["time"].values[-1], unit="D")
                d.close()
                raise ValueError(
                    f"{os.path.basename(reference)}: the nearest record to "
                    f"{date} is {got}, {gap:.0f} days away. The file covers "
                    f"{first} to {last} — it does not span this run. "
                    "Re-download with download_obs().")
            if gap > 0:
                print(f"  note: {date} matched to the {got} record, "
                      f"{gap:.0f} day away")

    if "depth" in da.dims:
        da = (da.isel(depth=0) if depth_m is None
              else da.sel(depth=abs(depth_m), method="nearest"))

    field = da.values * spec["scale"] + spec["offset"]
    lon = d["longitude"].values
    lat = d["latitude"].values
    d.close()

    if lon.ndim == 1:
        lon, lat = np.meshgrid(lon, lat)

    # A gappy product (L3 observations, not an analysis) has missing cells
    # wherever there was cloud. regrid_to_croco drops non-finite values and
    # then interpolates, so those holes get filled from whatever surrounds
    # them — a field that looks complete but is partly invented. Warn when
    # coverage is poor enough for that to matter.
    if spec.get("gappy"):
        cover = float(np.isfinite(field).sum()) / field.size * 100.0
        if cover < 100.0:
            print(f"  {source}: {cover:.0f}% of cells have data on {date}; "
                  "the rest are interpolated across by the regridding")
        if cover < 50.0:
            import warnings
            warnings.warn(
                f"{source} covers only {cover:.0f}% of the domain on {date}. "
                "More than half the comparison would be against interpolated "
                "values rather than observations.")

    return val.regrid_to_croco(lon, lat, field, croco_ds)


# ======================================================================
# Comparison
# ======================================================================
def compare(croco_his, reference, var="temp", kind="map", method="regrid",
            date=None, tindex=-1, depth_m=None,
            Yorig=None, margin_deg=None, daily_mean=False,
            min_depth=None, max_gap_days=0.5, out=None, verbose=True, **kw):
    """Compare a CROCO run against an observation reference.

    kind : 'map'    three panels — model, reference, difference
           'growth' error metrics against lead time
           'stats'  the statistics alone, no figure

    date : 'YYYY-MM-DD'. If None, the CROCO record given by `tindex` is used
      and the reference is matched to *that record's* date — not to its own
      first record, which would usually be a different day.

    daily_mean : average the CROCO day before comparing. The reference products
      are daily means; CROCO is usually sub-daily and carries a diurnal cycle,
      so without this the difference oscillates once a day.

    Returns (figure, stats) for 'map', the growth structure for 'growth', and
    the stats dict for 'stats'.
    """
    source = identify(reference)
    if source == "croco":
        raise ValueError("that reference is a CROCO file — use "
                         "gtools.validation.compare_resolution")
    if source == "insitu":
        raise ValueError("those are in-situ files — use "
                         "gtools.validation_godae.validate_against_insitu")
    if source == "unknown":
        raise ValueError(f"could not identify the reference {reference!r}")

    # A gappy product (L3 observations) should not be regridded: the holes
    # get filled by interpolation and then counted as agreement. Collocation
    # samples the model at the observed points instead.
    if method == "collocate" or (method == "auto"
                                 and REFERENCES[source].get("gappy")):
        return plot_collocation(croco_his, reference, var=var, date=date,
                                depth_m=depth_m, Yorig=Yorig,
                                daily_mean=daily_mean, min_depth=min_depth,
                                out=out)
    if method not in ("regrid", "auto"):
        raise ValueError("method must be 'regrid', 'collocate' or 'auto'")

    if kind == "growth":
        return growth(croco_his, reference, var=var, depth_m=depth_m,
                      Yorig=Yorig, margin_deg=margin_deg,
                      daily_mean=daily_mean, max_gap_days=max_gap_days,
                      out=out, verbose=verbose, **kw)

    ds = pp.open_history(croco_his, Yorig=Yorig)
    _require_dates(ds, croco_his)
    if margin_deg is not None:
        ds = pp.crop_interior(ds, margin_deg=margin_deg)
    if daily_mean:
        ds = val._maybe_daily_mean(ds, date, daily_mean)
    else:
        tindex = val._tindex_for_date(ds, date, tindex)

    # If no date was given, take it from the CROCO record actually selected,
    # so the reference is matched to the same day rather than to its own
    # first record.
    if date is None:
        date = str(np.datetime_as_string(pp.times(ds)[tindex], unit="D"))

    model = _croco_field(ds, var, tindex, depth_m)
    ref = load_reference(reference, var, ds, date=date, depth_m=depth_m,
                         source=source, max_gap_days=max_gap_days)
    clon, clat, mask = pp.lonlatmask(ds)

    # SSH: the two carry different datums — CROCO's zeta is about its own
    # reference level, DUACS' SLA about a mean sea surface. Only the anomaly
    # is comparable, so remove each field's mean.
    if var == "ssh":
        model = model - np.nanmean(model)
        ref = ref - np.nanmean(ref)

    s = val.domain_statistics(model, ref)

    label = _label(var, depth_m)
    if verbose:
        print(f"{label}  CROCO vs {source_name(source)}   {date}:")
        val._print_stats(label, s)

    if kind == "stats":
        ds.close()
        return s

    fig = three_panel(clon, clat, model, ref, mask, var, source,
                      date=date, depth_m=depth_m, out=out)
    ds.close()
    return fig, s


def growth(croco_his, reference, var="temp", depth_m=None, Yorig=None,
           margin_deg=None, daily_mean=True, max_gap_days=0.5,
           out=None, verbose=True, metric="rmse"):
    """Track the model-vs-reference error along a run, one point per record.

    daily_mean defaults to True here: the observation products are daily means
    while CROCO is usually sub-daily, so comparing raw records makes the error
    oscillate once per day with the diurnal cycle rather than showing growth.

    Returns {'lead_days': [...], 'times': [...], 'bias': [...], 'rmse': [...],
    'crmse': [...], 'corr': [...]}.
    """
    source = identify(reference)
    ds = pp.open_history(croco_his, Yorig=Yorig)
    _require_dates(ds, croco_his)
    if margin_deg is not None:
        ds = pp.crop_interior(ds, margin_deg=margin_deg)
    if daily_mean and "time" in ds.dims:
        ds = val._maybe_daily_mean(ds, None, True)

    times = pp.times(ds)
    lead = [(np.datetime64(t) - np.datetime64(times[0])) / np.timedelta64(1, "D")
            for t in times]
    res = {"lead_days": lead,
           "times": [str(np.datetime_as_string(t, unit="h")) for t in times],
           "bias": [], "rmse": [], "crmse": [], "corr": []}

    for k, t in enumerate(times):
        date = str(np.datetime_as_string(t, unit="D"))
        model = _croco_field(ds, var, k, depth_m)
        ref = load_reference(reference, var, ds, date=date, depth_m=depth_m,
                             source=source, max_gap_days=max_gap_days)
        if var == "ssh":
            model = model - np.nanmean(model)
            ref = ref - np.nanmean(ref)
        s = val.domain_statistics(model, ref)
        for m in ("bias", "rmse", "crmse", "corr"):
            res[m].append(s[m])
    ds.close()

    if verbose:
        print(f"{_label(var, depth_m)}  CROCO vs {source_name(source)}, by lead day:")
        for d, r, b in zip(res["lead_days"], res["rmse"], res["bias"]):
            print(f"   day {d:4.1f}   rmse {r:7.4f}   bias {b:+7.4f}")

    if out:
        fig, ax = plt.subplots(figsize=(9, 5))
        ax.plot(res["lead_days"], res[metric], marker="o", lw=1.8)
        ax.set_xlabel("lead time (days)")
        ax.set_ylabel(f"{metric.upper()}  ({_units(var)})")
        ax.set_title(f"{_label(var, depth_m)} — CROCO vs {source_name(source)}")
        ax.grid(alpha=0.3)
        fig.tight_layout()
        fig.savefig(out, dpi=150, bbox_inches="tight")
        plt.close(fig)
    return res


# ======================================================================
# Composite — many cycles, averaged by lead time
#
# A single forecast tells you how one week went. Averaging several by
# lead time tells you how the system behaves, which is what a
# verification figure is for.
#
# Two things worth being careful about, both of which this does:
#
# RMSE does not average linearly. Taking the mean of three cycles' RMSE
# values is close but wrong, and the error grows when the sample sizes
# differ — which they do here, because a gappy product has a different
# number of valid observations every day. The right quantity is the
# pooled mean square error: sum the squared errors across cycles at each
# lead, divide by the total count, then take the root.
#
# And with a handful of cycles there is no honest confidence interval.
# A bootstrap over three samples estimates nothing. The individual
# cycles are drawn instead, so the spread is visible rather than
# summarised into a band that implies more than is known.
# ======================================================================
def composite(croco_files, references, var="temp", days=None, depth_m=None,
              Yorig=None, daily_mean=True, min_depth=None, method="auto",
              max_gap_days=0.5, out=None, dpi=130, verbose=True):
    """Forecast skill against persistence, averaged across cycles by lead time.

    croco_files : a list of croco_his.nc paths, one per forecast cycle
    references  : one reference path, or a dict {var: path}
    var         : which variable to score

    Returns (composite, per_cycle) where composite is a list of dicts —
    lead_days, rmse_model, rmse_persist, skill, n, n_cycles — and per_cycle
    is what persistence() returned for each run.
    """
    ref = references[var] if isinstance(references, dict) else references

    per_cycle = []
    for f in croco_files:
        try:
            rows = persistence(f, ref, var=var, days=days, depth_m=depth_m,
                               Yorig=Yorig, daily_mean=daily_mean,
                               min_depth=min_depth, method=method,
                               max_gap_days=max_gap_days, verbose=False)
            per_cycle.append(rows)
        except Exception as e:
            print(f"  skipped {os.path.basename(os.path.dirname(f))}: {e}")

    if not per_cycle:
        raise ValueError("no cycle produced usable statistics")

    # pool by lead: sum of squared error and count, then the root of the mean
    leads = sorted({r["lead_days"] for c in per_cycle for r in c})
    comp = []
    for L in leads:
        rows = [r for c in per_cycle for r in c if r["lead_days"] == L]
        n = sum(r["n"] for r in rows)
        if not n:
            continue
        sse_m = sum(r["rmse_model"] ** 2 * r["n"] for r in rows)
        sse_p = sum(r["rmse_persist"] ** 2 * r["n"] for r in rows)
        rm, rp = np.sqrt(sse_m / n), np.sqrt(sse_p / n)
        comp.append({"lead_days": L, "n": n, "n_cycles": len(rows),
                     "rmse_model": float(rm), "rmse_persist": float(rp),
                     "skill": float(1.0 - (rm / rp) ** 2) if rp > 1e-12 else np.nan})

    if verbose:
        source = identify(ref)
        print(f"{_label(var, depth_m)}  CROCO vs {source_name(source)}, "
              f"{len(per_cycle)} cycles pooled:")
        print("   lead  cycles       n    model  persist    skill")
        for r in comp:
            print("   %4.1f  %6d  %7d   %6.3f   %6.3f  %+7.3f"
                  % (r["lead_days"], r["n_cycles"], r["n"],
                     r["rmse_model"], r["rmse_persist"], r["skill"]))

    if out:
        _composite_figure({var: (comp, per_cycle)}, [var], depth_m,
                          identify(ref), out, dpi)
    return comp, per_cycle


def composite_panels(croco_files, references, variables=("speed", "temp", "ssh"),
                     days=None, depth_m=None, Yorig=None, daily_mean=True,
                     min_depth=None, method="auto", max_gap_days=0.5,
                     out=None, dpi=130, verbose=True):
    """One panel per variable: RMSE against lead time, model and persistence.

    references must be a dict {var: path}, since each variable needs its own
    observation product — GlobCurrent for currents, ODYSSEA or OSTIA for SST,
    DUACS for sea level.

    The faint lines behind each pair are the individual cycles.
    """
    if not isinstance(references, dict):
        raise TypeError("composite_panels needs a dict {var: reference path}")

    results, sources = {}, {}
    for v in variables:
        if v not in references:
            print(f"  no reference given for {v!r} — skipped")
            continue
        comp, per = composite(croco_files, references, var=v, days=days,
                              depth_m=depth_m, Yorig=Yorig,
                              daily_mean=daily_mean, min_depth=min_depth,
                              method=method, max_gap_days=max_gap_days,
                              verbose=verbose)
        results[v] = (comp, per)
        sources[v] = identify(references[v])
        if verbose:
            print()

    if not results:
        raise ValueError("nothing to plot")
    if out:
        _composite_figure(results, list(results), depth_m, sources, out, dpi)
    return results


def _composite_figure(results, variables, depth_m, sources, out, dpi):
    """Two rows: RMSE against lead time, and the skill score below."""
    n = len(variables)
    fig, ax = plt.subplots(2, n, figsize=(5.0 * n, 8.5),
                           constrained_layout=True, squeeze=False)

    for k, v in enumerate(variables):
        comp, per = results[v]
        src = sources[v] if isinstance(sources, dict) else sources
        lead = [r["lead_days"] for r in comp]

        # the individual cycles, faint
        for c in per:
            ax[0][k].plot([r["lead_days"] for r in c],
                          [r["rmse_model"] for r in c],
                          color="C0", lw=0.8, alpha=0.35)
            ax[0][k].plot([r["lead_days"] for r in c],
                          [r["rmse_persist"] for r in c],
                          color="C1", lw=0.8, alpha=0.35, ls="--")
            ax[1][k].plot([r["lead_days"] for r in c],
                          [r["skill"] for r in c],
                          color="C2", lw=0.8, alpha=0.35)

        ax[0][k].plot(lead, [r["rmse_model"] for r in comp], "o-",
                      color=col["croco"], lw=2.2, label="CROCO")
        ax[0][k].plot(lead, [r["rmse_persist"] for r in comp], "s--",
                      color="C1", lw=2.2, label="persistence")
        ax[0][k].set_title(f"{_label(v, depth_m)} vs {source_name(src)}")
        ax[0][k].set_ylabel(f"RMSE ({_units(v)})")
        ax[0][k].grid(alpha=0.3)
        if k == 0:
            ax[0][k].legend()

        ax[1][k].axhline(0, color="k", lw=1)
        ax[1][k].plot(lead, [r["skill"] for r in comp], "o-",
                      color="C2", lw=2.2)
        ax[1][k].set_xlabel("lead time (days)")
        ax[1][k].set_ylabel("skill score")
        ax[1][k].grid(alpha=0.3)

    ncyc = len(next(iter(results.values()))[1])
    fig.suptitle(f"Forecast skill against persistence — {ncyc} cycles pooled\n"
                 "faint lines are the individual cycles", fontsize=13)
    fig.savefig(out, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    return out


# ======================================================================
# Three-way skill — GODAE classes 1 and 3 together
#
# The comparison that answers the question the whole project rests on:
# did downscaling improve on the parent?
#
# Three forecasts, scored against the same independent observations:
#
#   CROCO        the downscaled forecast
#   the parent   Mercator's own forecast for the same days — a competitor,
#                not a reference. This is the move that makes the figure
#                mean something: using Mercator as the reference measures
#                whether CROCO tracked its own boundaries, whereas scoring
#                Mercator against observations asks whether the extra
#                resolution bought anything.
#   persistence  CROCO's day-zero field held fixed — the forecast you get
#                for free
#
# All three are sampled at the same points, so the three curves are
# directly comparable rather than each being computed over its own set.
# ======================================================================
def _parent_file_for(croco_his):
    """The Mercator file the driver downloaded for this cycle.

    The driver writes it to <cycle>/downloaded_data/MERCATOR/MERCATOR_<tag>_00.nc
    and croco_his sits at <cycle>/fcst/CROCO_FILES/croco_his.nc, so the cycle
    root is two levels up from the history file's directory.
    """
    import glob as _glob
    cycle = os.path.dirname(os.path.dirname(os.path.dirname(
        os.path.abspath(croco_his))))
    hits = sorted(_glob.glob(os.path.join(cycle, "downloaded_data", "MERCATOR",
                                          "MERCATOR_*.nc")))
    if not hits:
        raise FileNotFoundError(
            f"no Mercator file under {cycle}/downloaded_data/MERCATOR/. "
            "Pass parent_files= explicitly if your layout differs.")
    return hits[-1]


def _sample_at(field, clon, clat, cmask, lon, lat):
    """Interpolate a model field to scattered points, masking dry ground.

    The same masking as collocate(): griddata interpolates across gaps in
    its input rather than returning NaN, so points without wet cells around
    them are dropped instead of being given a value blended from far away.
    """
    from scipy.interpolate import griddata as _griddata

    good = np.isfinite(field)
    tgt = np.column_stack([lon, lat])
    out = _griddata(np.column_stack([clon[good].ravel(), clat[good].ravel()]),
                    field[good].ravel(), tgt, method="linear")
    wet = _griddata(np.column_stack([clon.ravel(), clat.ravel()]),
                    (np.asarray(cmask) > 0).astype(float).ravel(),
                    tgt, method="linear")
    return np.where(np.isfinite(wet) & (wet > 0.99), out, np.nan)


def three_way(croco_his, observations, var="temp", days=None, depth_m=None,
              Yorig=None, daily_mean=True, min_depth=None, parent_file=None,
              max_gap_days=0.5, verbose=True):
    """One cycle: CROCO, the parent and persistence, against the observations.

    Returns a list of dicts, one per lead day, with rmse_croco, rmse_parent
    and rmse_persist plus the counts.
    """
    source = identify(observations)
    if source in ("croco", "insitu", "unknown"):
        raise ValueError(f"three_way needs a gridded observation product, "
                         f"got {source!r}")
    obs_spec = REFERENCES[source]
    parent_file = parent_file or _parent_file_for(croco_his)
    par_spec = REFERENCES["glorys"]

    ds = pp.open_history(croco_his, Yorig=Yorig)
    _require_dates(ds, croco_his)
    if daily_mean and "time" in ds.dims:
        ds = val._maybe_daily_mean(ds, None, True)
    times = pp.times(ds)
    clon, clat, cmask = pp.lonlatmask(ds)
    idx = _select_records(times, days)
    init_field = _croco_field(ds, var, idx[0], depth_m)
    fields = {k: _croco_field(ds, var, k, depth_m) for k in idx}
    ds.close()

    rows = []
    for k in idx:
        date = str(np.datetime_as_string(times[k], unit="D"))

        # the observations, on their own grid, gaps intact
        o, olon, olat = _read_obs(observations, var, date, source,
                                  depth_m, max_gap_days)
        v = np.isfinite(o)
        if not v.any():
            continue
        lo, la = olon[v], olat[v]

        croco = _sample_at(fields[k], clon, clat, cmask, lo, la)
        persist = _sample_at(init_field, clon, clat, cmask, lo, la)
        parent = _read_parent(parent_file, var, date, depth_m, lo, la,
                              max_gap_days)

        ok = (np.isfinite(croco) & np.isfinite(persist) & np.isfinite(parent)
              & np.isfinite(o[v]))
        if min_depth is not None:
            ok &= _depth_at(croco_his, Yorig, lo, la) >= min_depth
        if ok.sum() < 2:
            continue

        oo = o[v][ok]
        if var == "ssh":                    # different datums all round
            oo = oo - np.nanmean(oo)
            cc = croco[ok] - np.nanmean(croco[ok])
            pp_ = persist[ok] - np.nanmean(persist[ok])
            qq = parent[ok] - np.nanmean(parent[ok])
        else:
            cc, pp_, qq = croco[ok], persist[ok], parent[ok]

        sc = godae_metrics(cc, oo)
        sp = godae_metrics(pp_, oo)
        sq = godae_metrics(qq, oo)
        lead = float((np.datetime64(times[k]) - np.datetime64(times[idx[0]]))
                     / np.timedelta64(1, "D"))
        rows.append({"lead_days": lead, "date": date, "n": int(ok.sum()),
                     "rmse_croco": sc["rmsd"], "rmse_parent": sq["rmsd"],
                     "rmse_persist": sp["rmsd"],
                     "bias_croco": sc["bias"], "bias_parent": sq["bias"],
                     "corr_croco": sc["corr"], "corr_parent": sq["corr"]})

    if verbose and rows:
        print(f"{_label(var, depth_m)}  vs {source_name(source)}:")
        print("   lead  date         n     CROCO   parent  persist")
        for r in rows:
            print("   %4.1f  %s  %6d  %6.3f   %6.3f   %6.3f"
                  % (r["lead_days"], r["date"], r["n"], r["rmse_croco"],
                     r["rmse_parent"], r["rmse_persist"]))
    return rows


def _read_obs(fname, var, date, source, depth_m, max_gap_days):
    """The observation field for one day, on its own grid, gaps intact."""
    spec = REFERENCES[source]
    if var == "speed":              # products carry u and v, not the speed
        u, lon, lat = _read_obs(fname, "u", date, source, depth_m, max_gap_days)
        v, _, _ = _read_obs(fname, "v", date, source, depth_m, max_gap_days)
        return np.sqrt(u ** 2 + v ** 2), lon, lat
    d = xr.open_dataset(fname)
    da = d[spec["vars"][var]]
    if "time" in da.dims:
        want = np.datetime64(date)
        da = da.sel(time=want, method="nearest")
        gap = abs(float((np.datetime64(da["time"].values, "D")
                         - np.datetime64(want, "D")) / np.timedelta64(1, "D")))
        if gap > max_gap_days:
            d.close()
            raise ValueError(f"{os.path.basename(fname)}: no record within "
                             f"{max_gap_days} days of {date}")
    if "depth" in da.dims:
        da = (da.isel(depth=0) if depth_m is None
              else da.sel(depth=abs(depth_m), method="nearest"))
    o = da.values * spec["scale"] + spec["offset"]

    qc = spec.get("qc")
    if qc and qc["var"] in d:
        ql = d[qc["var"]]
        if "time" in ql.dims:
            ql = ql.sel(time=da["time"].values, method="nearest")
        o = np.where(ql.values >= qc["min"], o, np.nan)

    lon, lat = d["longitude"].values, d["latitude"].values
    d.close()
    if lon.ndim == 1:
        lon, lat = np.meshgrid(lon, lat)
    return o, lon, lat


def _read_parent(fname, var, date, depth_m, lon, lat, max_gap_days):
    """The parent's forecast for one day, sampled at the observation points."""
    from scipy.interpolate import griddata as _griddata
    spec = REFERENCES["glorys"]

    if var == "speed":
        u = _read_parent(fname, "u", date, depth_m, lon, lat, max_gap_days)
        v = _read_parent(fname, "v", date, depth_m, lon, lat, max_gap_days)
        return np.sqrt(u ** 2 + v ** 2)

    d = xr.open_dataset(fname)
    da = d[spec["vars"][var]]
    if "time" in da.dims:
        want = np.datetime64(date)
        da = da.sel(time=want, method="nearest")
        gap = abs(float((np.datetime64(da["time"].values, "D")
                         - np.datetime64(want, "D")) / np.timedelta64(1, "D")))
        if gap > max_gap_days:
            d.close()
            raise ValueError(f"{os.path.basename(fname)}: the parent has no "
                             f"record within {max_gap_days} days of {date}")
    if "depth" in da.dims:
        da = (da.isel(depth=0) if depth_m is None
              else da.sel(depth=abs(depth_m), method="nearest"))
    f = da.values
    plon, plat = d["longitude"].values, d["latitude"].values
    d.close()
    if plon.ndim == 1:
        plon, plat = np.meshgrid(plon, plat)

    g = np.isfinite(f)
    return _griddata(np.column_stack([plon[g].ravel(), plat[g].ravel()]),
                     f[g].ravel(), np.column_stack([lon, lat]),
                     method="linear")


def _depth_at(croco_his, Yorig, lon, lat):
    """The model's bottom depth at scattered points."""
    from scipy.interpolate import griddata as _griddata
    ds = pp.open_history(croco_his, Yorig=Yorig)
    clon, clat, _ = pp.lonlatmask(ds)
    h = ds["h"].values
    ds.close()
    return _griddata(np.column_stack([clon.ravel(), clat.ravel()]), h.ravel(),
                     np.column_stack([lon, lat]), method="linear")


def skill_panels(croco_files, observations, variables=("temp", "ssh", "speed"),
                 depths=None, days=None, Yorig=None, daily_mean=True,
                 min_depth=None, parent_files=None, max_gap_days=0.5,
                 title=None, colors=None, figsize=None, ylim=None,
                 out=None, dpi=130, verbose=True):
    """One panel per variable: CROCO, the parent and persistence, pooled.

    observations : {var: path} — the independent product for each variable
    depths       : {var: metres} for variables to compare below the surface,
                   e.g. {'speed': 15} for the 15 m current
    parent_files : {croco_his: mercator_path}, or None to find them beside
                   each run

    Squared errors are pooled across cycles at each lead, so the mean is the
    root of the pooled mean square rather than the mean of the RMSEs — those
    differ when the counts do, and a gappy product's counts differ daily.
    """
    depths = depths or {}
    results = {}

    for v in variables:
        if v not in observations:
            print(f"  no observations given for {v!r} — skipped")
            continue
        per_cycle = []
        for f in croco_files:
            pf = (parent_files or {}).get(f)
            try:
                rows = three_way(f, observations[v], var=v, days=days,
                                 depth_m=depths.get(v), Yorig=Yorig,
                                 daily_mean=daily_mean, min_depth=min_depth,
                                 parent_file=pf, max_gap_days=max_gap_days,
                                 verbose=False)
                if rows:
                    per_cycle.append(rows)
            except Exception as e:
                print(f"  {v}: skipped {os.path.basename(os.path.dirname(f))}: {e}")
        if not per_cycle:
            continue

        leads = sorted({r["lead_days"] for c in per_cycle for r in c})
        comp = []
        for L in leads:
            rr = [r for c in per_cycle for r in c if r["lead_days"] == L]
            n = sum(r["n"] for r in rr)
            if not n:
                continue
            row = {"lead_days": L, "n": n, "n_cycles": len(rr)}
            for who in ("croco", "parent", "persist"):
                sse = sum(r[f"rmse_{who}"] ** 2 * r["n"] for r in rr)
                row[f"rmse_{who}"] = float(np.sqrt(sse / n))
            comp.append(row)
        results[v] = (comp, per_cycle, identify(observations[v]))

        if verbose:
            print(f"{_label(v, depths.get(v))}  vs "
                  f"{source_name(identify(observations[v]))}, "
                  f"{len(per_cycle)} cycles pooled:")
            print("   lead       n    CROCO   parent  persist")
            for r in comp:
                print("   %4.1f  %7d   %6.3f   %6.3f   %6.3f"
                      % (r["lead_days"], r["n"], r["rmse_croco"],
                         r["rmse_parent"], r["rmse_persist"]))
            print()

    if not results:
        raise ValueError("nothing to plot")

    if out:
        vs = list(results)
        n = len(vs)
        col = {"croco": "C0", "parent": "C3", "persist": "0.4"}
        col.update(colors or {})
        fig, ax = plt.subplots(1, n, figsize=figsize or (5.0 * n, 4.6),
                               constrained_layout=True, squeeze=False)
        for k, v in enumerate(vs):
            comp, per, src = results[v]
            lead = [r["lead_days"] for r in comp]
            a = ax[0][k]
            for c in per:                       # the individual cycles
                L = [r["lead_days"] for r in c]
                a.plot(L, [r["rmse_croco"] for r in c],
                       color=col["croco"], lw=0.7, alpha=0.3)
                a.plot(L, [r["rmse_parent"] for r in c],
                       color=col["parent"], lw=0.7, alpha=0.3)
                a.plot(L, [r["rmse_persist"] for r in c],
                       color=col["persist"], lw=0.7, alpha=0.3, ls="--")
            a.plot(lead, [r["rmse_croco"] for r in comp], "o-",
                   color="C0", lw=2.2, label="CROCO")
            a.plot(lead, [r["rmse_parent"] for r in comp], "o-",
                   color=col["parent"], lw=2.2, label="Mercator")
            a.plot(lead, [r["rmse_persist"] for r in comp], "s--",
                   color=col["persist"], lw=2.0, label="persistence")
            dm = depths.get(v)
            a.set_title(f"{_label(v, dm)} vs {source_name(src)}")
            a.set_xlabel("lead time (days)")
            a.set_ylabel(f"RMSE ({_units(v)})")
            if ylim and v in ylim:
                a.set_ylim(*ylim[v])
            a.grid(alpha=0.3)
            if k == 0:
                a.legend()
        ncyc = len(next(iter(results.values()))[1])
        fig.suptitle(f"Forecast error against independent observations — "
                     f"{ncyc} cycles pooled\nfaint lines are the individual "
                     "cycles; lead zero is the initial condition", fontsize=12)
        fig.savefig(out, dpi=dpi, bbox_inches="tight")
        plt.close(fig)
    return results


# ======================================================================
# Plotting
#
# These replace val._three_panel, which labels the reference "parent" and
# has no way to know which product it is looking at. Here the source name
# and the date come through into the titles, and the two field panels
# share one colour scale so they are directly comparable.
# ======================================================================
def _coast(ax, lon, lat, mask):
    """Grey land and a coastline, drawn from the model's own mask."""
    ax.pcolormesh(lon, lat, np.where(mask == 0, 1, np.nan),
                  cmap=ListedColormap(["0.85"]), vmin=0, vmax=1, zorder=0)
    ax.contour(lon, lat, mask, levels=[0.5], colors="k", linewidths=0.6)


def three_panel(lon, lat, model, ref, mask, var, source,
                date=None, depth_m=None, out=None, figsize=(17, 6),
                vmin=None, vmax=None, dlim=None, dpi=130):
    """Model, reference and difference, with one shared scale for the fields.

    The two field panels share a colour scale, so a difference visible between
    the panels is a real difference and not a scaling artefact. The difference
    panel gets its own diverging scale, symmetric about zero.
    """
    name = source_name(source)
    label = _label(var, depth_m)
    when = f"  {date}" if date else ""

    both = np.concatenate([model[np.isfinite(model)].ravel(),
                           ref[np.isfinite(ref)].ravel()])
    if vmin is None:
        vmin = float(np.percentile(both, 1))
    if vmax is None:
        vmax = float(np.percentile(both, 99))
    if dlim is None:
        dd = model - ref
        dlim = float(np.nanpercentile(np.abs(dd), 99))

    fig, ax = plt.subplots(1, 3, figsize=figsize, constrained_layout=True)
    for a in ax:
        _coast(a, lon, lat, mask)
        a.set_aspect("equal")
        a.set_xlabel("longitude")
    ax[0].set_ylabel("latitude")

    h0 = ax[0].pcolormesh(lon, lat, model, cmap=_cmap(var), vmin=vmin, vmax=vmax)
    ax[0].set_title(f"CROCO {label}{when}")
    h1 = ax[1].pcolormesh(lon, lat, ref, cmap=_cmap(var), vmin=vmin, vmax=vmax)
    ax[1].set_title(f"{name} {label}{when}")
    fig.colorbar(h1, ax=ax[:2], label=_units(var), shrink=0.85, pad=0.02,
                 extend="both")

    h2 = ax[2].pcolormesh(lon, lat, model - ref, cmap="RdBu_r",
                          vmin=-dlim, vmax=dlim)
    ax[2].set_title(f"CROCO \u2212 {name}")
    fig.colorbar(h2, ax=ax[2], label=_units(var), shrink=0.85, pad=0.02,
                 extend="both")

    if out:
        fig.savefig(out, dpi=dpi, bbox_inches="tight")
        plt.close(fig)
        return out
    return fig


def _select_records(times, days):
    """Work out which record indices to show.

    days may be
        None                every record
        5                   the first five
        (2, 7)              a slice, records 2 to 6
        ('2026-09-01',      a date range, inclusive
         '2026-09-05')
        ['2026-07-11', ...] explicit dates, matched to the nearest record
        [0, 3, 6]           explicit indices

    Returns a list of indices into `times`.
    """
    n = len(times)
    if days is None:
        return list(range(n))
    if isinstance(days, int):
        return list(range(min(days, n)))
    if isinstance(days, tuple) and len(days) == 2 and all(
            isinstance(x, int) for x in days):
        return list(range(*days))[:n]
    if isinstance(days, tuple) and len(days) == 2 and all(
            isinstance(x, str) for x in days):
        # a date range, inclusive: every record between the two
        lo, hi = np.datetime64(days[0], "D"), np.datetime64(days[1], "D")
        idx = [k for k, t in enumerate(times)
               if lo <= np.datetime64(t, "D") <= hi]
        if not idx:
            first = np.datetime_as_string(times[0], unit="D")
            last = np.datetime_as_string(times[-1], unit="D")
            raise ValueError(
                f"no records between {days[0]} and {days[1]}; the run covers "
                f"{first} to {last}.")
        return idx
    if isinstance(days, (list, tuple)):
        if all(isinstance(x, (int, np.integer)) for x in days):
            bad = [k for k in days if k >= n or k < -n]
            if bad:
                raise IndexError(
                    f"record index {bad} out of range; the run has {n} records")
            return [int(k) for k in days]
        # dates: take the nearest record to each, and say how far it was
        want = [np.datetime64(str(d)) for d in days]
        idx = []
        for w in want:
            k = int(np.argmin(np.abs(np.asarray(times) - w)))
            gap = abs(float((np.datetime64(times[k], "D") - np.datetime64(w, "D"))
                            / np.timedelta64(1, "D")))
            if gap > 1.0:
                first = np.datetime_as_string(times[0], unit="D")
                last = np.datetime_as_string(times[-1], unit="D")
                raise ValueError(
                    f"no record near {np.datetime_as_string(w, unit='D')}: the "
                    f"nearest is {np.datetime_as_string(times[k], unit='D')}, "
                    f"{gap:.0f} days away. The run covers {first} to {last}.")
            idx.append(k)
        return idx
    raise TypeError("days must be None, an int, a (start, stop) tuple, "
                    "a list of dates, or a list of indices")


def compare_days(croco_his, reference, var="temp", days=None, depth_m=None,
                 Yorig=None, margin_deg=None, daily_mean=True,
                 rows=("reference", "croco", "difference"),
                 max_gap_days=0.5, vmin=None, vmax=None, dlim=None,
                 cmap=None, dcmap="RdBu_r", title=None, figsize=None,
                 out=None, dpi=130, verbose=True):
    """A grid of comparisons: one column per day, one row per view.

    days : which records to show —
        None                  every record in the run
        5                     the first five
        (2, 7)                records 2 to 6
        ['2026-07-11',        named dates, each matched to the nearest record
         '2026-07-14']
        [0, 3, 6]             record indices

    rows : which views to stack, top to bottom. Any subset and order of
        'reference', 'croco' and 'difference'. The default puts the
        observation first, the model below it, and the difference last — so
        the eye compares the two fields directly and then reads what is left
        over.

    The field rows share one colour scale and the difference row has its own.
    Sharing along and between the field rows is the point: a change between
    days, or between model and reference, then reads as a real change rather
    than as each panel being rescaled to its own range.

    Returns (figure, stats) — or (path, stats) when `out` is given — where
    stats is one dict per column.
    """
    source = identify(reference)
    if source in ("croco", "insitu", "unknown"):
        raise ValueError(f"compare_days needs a gridded reference, got {source!r}")
    name = source_name(source)

    valid = ("reference", "croco", "difference")
    rows = tuple(rows)
    bad = [r for r in rows if r not in valid]
    if bad:
        raise ValueError(f"rows may contain {valid}; got {bad}")
    if not rows:
        raise ValueError("rows must name at least one view")

    ds = pp.open_history(croco_his, Yorig=Yorig)
    _require_dates(ds, croco_his)
    if margin_deg is not None:
        ds = pp.crop_interior(ds, margin_deg=margin_deg)
    if daily_mean and "time" in ds.dims:
        ds = val._maybe_daily_mean(ds, None, True)

    lon, lat, mask = pp.lonlatmask(ds)
    times = pp.times(ds)
    idx = _select_records(times, days)
    if not idx:
        ds.close()
        raise ValueError("no records selected")

    models, refs, dates, stats = [], [], [], []
    for k in idx:
        date = str(np.datetime_as_string(times[k], unit="D"))
        m = _croco_field(ds, var, k, depth_m)
        r = load_reference(reference, var, ds, date=date, depth_m=depth_m,
                           source=source, max_gap_days=max_gap_days)
        if var == "ssh":                       # different datums
            m, r = m - np.nanmean(m), r - np.nanmean(r)
        models.append(m)
        refs.append(r)
        dates.append(date)
        stats.append(val.domain_statistics(m, r))
    ds.close()

    n = len(idx)
    diffs = [m - r for m, r in zip(models, refs)]

    # one scale for every field panel in the figure, one for the differences
    allf = np.concatenate([a[np.isfinite(a)].ravel() for a in models + refs])
    if vmin is None:
        vmin = float(np.percentile(allf, 1))
    if vmax is None:
        vmax = float(np.percentile(allf, 99))
    if dlim is None:
        dlim = float(np.nanpercentile(np.abs(np.concatenate(
            [d[np.isfinite(d)].ravel() for d in diffs])), 99))
    cmap = cmap or _cmap(var)

    label = _label(var, depth_m)
    ylabels = {"reference": f"{name} {label}",
               "croco": f"CROCO {label}",
               "difference": f"CROCO \u2212 {name}"}
    data = {"reference": refs, "croco": models, "difference": diffs}

    nrow = len(rows)
    fig, ax = plt.subplots(nrow, n, figsize=figsize or (3.4 * n, 4.2 * nrow),
                           constrained_layout=True, squeeze=False)

    handles = {}
    for i, row in enumerate(rows):
        for k in range(n):
            a = ax[i][k]
            _coast(a, lon, lat, mask)
            a.set_aspect("equal")
            if row == "difference":
                h = a.pcolormesh(lon, lat, diffs[k], cmap=dcmap,
                                 vmin=-dlim, vmax=dlim)
                a.set_title(f"rmse {stats[k]['rmse']:.2f}   "
                            f"bias {stats[k]['bias']:+.2f}", fontsize=9)
            else:
                h = a.pcolormesh(lon, lat, data[row][k], cmap=cmap,
                                 vmin=vmin, vmax=vmax)
                if i == 0:
                    a.set_title(dates[k], fontsize=10)
            if k:
                a.set_yticklabels([])
            if i < nrow - 1:
                a.set_xticklabels([])
            else:
                a.set_xlabel("longitude")
        ax[i][0].set_ylabel(ylabels[row])
        handles[row] = h

    # one colorbar for the field rows together, one for the difference row
    field_rows = [i for i, r in enumerate(rows) if r != "difference"]
    if field_rows:
        axes = [ax[i][k] for i in field_rows for k in range(n)]
        fig.colorbar(handles[rows[field_rows[0]]], ax=axes, label=_units(var),
                     shrink=0.85, pad=0.02, extend="both")
    if "difference" in rows:
        i = rows.index("difference")
        fig.colorbar(handles["difference"], ax=[ax[i][k] for k in range(n)],
                     label=_units(var), shrink=0.85, pad=0.02, extend="both")

    fig.suptitle(f"CROCO vs {name} — {label}", fontsize=13)

    if verbose:
        print(f"{label}  CROCO vs {name}:")
        for d, s in zip(dates, stats):
            print(f"   {d}   rmse {s['rmse']:7.4f}   bias {s['bias']:+7.4f}"
                  f"   corr {s['corr']:6.3f}")

    if out:
        fig.savefig(out, dpi=dpi, bbox_inches="tight")
        plt.close(fig)
        return out, stats
    return fig, stats


# ======================================================================
# Metrics
#
# The GODAE OceanView scorecard quantities, computed on collocated pairs.
# Self-contained: everything the comparisons here report comes from this
# one function, so the numbers in a map, a scorecard and a scatter are
# the same numbers computed the same way.
# ======================================================================
def godae_metrics(model, ref):
    """The standard scorecard for two collocated 1D arrays.

    model, ref : same length, paired in space and time. Points where either
      is non-finite are dropped pairwise.

    Returns n, bias, rmsd, urmsd, corr, si, si_std, std_model, std_ref,
    std_ratio, mean_model, mean_ref.

        bias   mean(model - ref) — a systematic offset
        rmsd   root-mean-square difference, total error
        urmsd  the same after removing each field's mean, so it measures
               pattern error alone. rmsd**2 = bias**2 + urmsd**2
        corr   Pearson correlation
        si     scatter index: urmsd as a percentage of mean(|ref|). The usual
               GODAE convention, but unstable when the reference has a near-zero
               mean — SSH anomaly, velocity components — so si_std normalises
               by the standard deviation instead and stays well behaved.
    """
    m = np.asarray(model, dtype=float).ravel()
    r = np.asarray(ref, dtype=float).ravel()
    ok = np.isfinite(m) & np.isfinite(r)
    m, r = m[ok], r[ok]

    keys = ("n", "bias", "rmsd", "urmsd", "corr", "si", "si_std",
            "std_model", "std_ref", "std_ratio", "mean_model", "mean_ref")
    if m.size < 2:
        return {k: (0 if k == "n" else np.nan) for k in keys}

    d = m - r
    bias = float(d.mean())
    rmsd = float(np.sqrt((d ** 2).mean()))
    ma, ra = m - m.mean(), r - r.mean()
    urmsd = float(np.sqrt(((ma - ra) ** 2).mean()))
    std_m, std_r = float(m.std()), float(r.std())
    corr = (float(np.corrcoef(m, r)[0, 1])
            if std_m > 1e-12 and std_r > 1e-12 else np.nan)
    ref_abs = float(np.abs(r).mean())
    si = float(urmsd / ref_abs * 100.0) if ref_abs > 1e-12 else np.nan
    si_std = float(urmsd / std_r * 100.0) if std_r > 1e-12 else np.nan

    return {"n": int(m.size), "bias": bias, "rmsd": rmsd, "urmsd": urmsd,
            "corr": corr, "si": si, "si_std": si_std,
            "std_model": std_m, "std_ref": std_r,
            "std_ratio": float(std_m / std_r) if std_r > 1e-12 else np.nan,
            "mean_model": float(m.mean()), "mean_ref": float(r.mean())}


# ======================================================================
# Collocation — for gappy products (L3 observations)
#
# An L3 product carries only the cells that were observed: cloud leaves
# holes, and coverage changes day to day. Regridding such a field onto the
# CROCO grid interpolates across those holes, producing something that
# looks complete but is partly invented, and then counts the invented
# cells as agreement.
#
# Collocation avoids that. Instead of bringing the observations to the
# model grid, it samples the model at each observed cell. The result is a
# pair of 1D arrays — one model value and one observation per valid
# point — which is what the GODAE metrics want anyway.
# ======================================================================
def collocate(croco_his, reference, var="temp", date=None, tindex=-1,
              depth_m=None, Yorig=None, daily_mean=True, min_depth=None,
              max_gap_days=0.5):
    """Sample the model at each valid observation cell.

    Returns a dict with the paired values and where they came from:
        model, obs   1D arrays, same length, one entry per valid observation
        lon, lat     where each pair sits
        date         the day compared
        n_obs        how many observations were valid
        n_grid       how many cells the observation grid has
        coverage     n_obs / n_grid, as a percentage

    min_depth : ignore observations over water shallower than this, in metres.
      Useful where the model and the reference disagree about the coast — see
      the note in compare().
    """
    from scipy.interpolate import griddata as _griddata

    source = identify(reference)
    if source in ("croco", "insitu", "unknown"):
        raise ValueError(f"collocate needs a gridded reference, got {source!r}")
    spec = REFERENCES[source]

    ds = pp.open_history(croco_his, Yorig=Yorig)
    _require_dates(ds, croco_his)
    if daily_mean and "time" in ds.dims:
        ds = val._maybe_daily_mean(ds, date, True)
    else:
        tindex = val._tindex_for_date(ds, date, tindex)
    if date is None:
        date = str(np.datetime_as_string(pp.times(ds)[tindex], unit="D"))

    model_field = _croco_field(ds, var, tindex, depth_m)
    clon, clat, cmask = pp.lonlatmask(ds)
    h = ds["h"].values
    ds.close()

    # the observations, on their own grid, untouched
    name = spec["vars"][var]
    d = xr.open_dataset(reference)
    da = d[name]
    if "time" in da.dims:
        want = np.datetime64(date)
        da = da.sel(time=want, method="nearest")
        got = np.datetime64(da["time"].values, "D")
        gap = abs(float((got - np.datetime64(want, "D")) / np.timedelta64(1, "D")))
        if gap > max_gap_days:
            first = np.datetime_as_string(d["time"].values[0], unit="D")
            last = np.datetime_as_string(d["time"].values[-1], unit="D")
            d.close()
            raise ValueError(
                f"{os.path.basename(reference)}: the nearest record to {date} "
                f"is {got}, {gap:.0f} days away. The file covers {first} to "
                f"{last} — it does not span this run. A silent nearest-match "
                "here would compare one day's model against another day's "
                "observations. Re-download with download_obs().")
        if gap > 0:
            print(f"  note: {date} matched to the {got} record, "
                  f"{gap:.0f} day away")
    if "depth" in da.dims:
        da = (da.isel(depth=0) if depth_m is None
              else da.sel(depth=abs(depth_m), method="nearest"))

    obs = da.values * spec["scale"] + spec["offset"]

    # Quality flag. An L3 product labels each retrieval; anything below the
    # product's own "clear" level is a marginal retrieval — near a cloud edge,
    # at a high satellite zenith angle — and is dropped rather than compared.
    qc = spec.get("qc")
    if qc:
        if qc["var"] in d:
            ql = d[qc["var"]]
            if "time" in ql.dims:
                ql = ql.sel(time=da["time"].values, method="nearest")
            obs = np.where(ql.values >= qc["min"], obs, np.nan)
        else:
            import warnings
            warnings.warn(
                f"{os.path.basename(reference)} has no {qc['var']!r}, so no "
                "quality filtering was applied. Re-download with download_obs() "
                "to include it.")

    olon = d["longitude"].values
    olat = d["latitude"].values
    d.close()
    if olon.ndim == 1:
        olon, olat = np.meshgrid(olon, olat)

    n_grid = obs.size
    valid = np.isfinite(obs)
    n_obs = int(valid.sum())
    if n_obs == 0:
        raise ValueError(f"{source}: no valid observations on {date}")

    # sample the model at the observed points. The model field is complete
    # over its own domain, so this interpolation fills nothing in — unlike
    # regridding the gappy field the other way.
    good = np.isfinite(model_field)
    pts = np.column_stack([clon[good].ravel(), clat[good].ravel()])
    tgt = np.column_stack([olon[valid], olat[valid]])
    model_at_obs = _griddata(pts, model_field[good].ravel(), tgt,
                             method="linear")

    o = obs[valid]
    m = model_at_obs
    lo, la = tgt[:, 0], tgt[:, 1]

    # griddata interpolates across any gap in its input points rather than
    # returning NaN, so an observation that sits where the model has no water
    # still gets a value — blended from wet cells that may be far away and on
    # the wrong side of a coast or a boundary. A row of observations along the
    # model's southern edge came back at 3-7 C in tropical water this way, and
    # those few dozen pairs alone doubled the RMSE.
    #
    # Neither a bounding box nor a nearest-neighbour distance catches this:
    # the points are inside the model's lon/lat span and within a cell of some
    # wet point. What distinguishes them is the mask. Interpolating the
    # land/sea mask to the same targets gives 1.0 only where every surrounding
    # model cell is wet; anywhere near land or past the edge it comes back
    # lower, because dry cells are being blended in.
    allpts = np.column_stack([clon.ravel(), clat.ravel()])
    wet = _griddata(allpts, (np.asarray(cmask) > 0).astype(float).ravel(),
                    tgt, method="linear")
    inside = np.isfinite(m) & np.isfinite(wet) & (wet > 0.99)
    m, o, lo, la = m[inside], o[inside], lo[inside], la[inside]

    if min_depth is not None:
        # the model's bottom depth at each observation point
        hp = _griddata(np.column_stack([clon.ravel(), clat.ravel()]),
                       h.ravel(), (lo, la), method="linear")
        deep = np.isfinite(hp) & (hp >= min_depth)
        m, o, lo, la = m[deep], o[deep], lo[deep], la[deep]

    return {"model": m, "obs": o, "lon": lo, "lat": la, "date": date,
            "n_obs": n_obs, "n_grid": n_grid,
            "coverage": 100.0 * n_obs / n_grid,
            "n_used": int(m.size), "source": source, "var": var,
            "depth_m": depth_m}


def scorecard(croco_his, reference, var="temp", days=None, depth_m=None,
              Yorig=None, daily_mean=True, min_depth=None,
              max_gap_days=0.5, verbose=True):
    """Collocated statistics, one row per day.

    Returns a list of dicts: date, coverage, n, bias, rmsd, urmsd, corr, si.
    The GODAE metrics come from gtools.validation_godae, so these numbers sit
    beside the in-situ scorecard directly.
    """

    ds = pp.open_history(croco_his, Yorig=Yorig)
    _require_dates(ds, croco_his)
    if daily_mean and "time" in ds.dims:
        ds = val._maybe_daily_mean(ds, None, True)
    times = pp.times(ds)
    ds.close()
    idx = _select_records(times, days)

    rows = []
    for k in idx:
        date = str(np.datetime_as_string(times[k], unit="D"))
        c = collocate(croco_his, reference, var=var, date=date,
                      depth_m=depth_m, Yorig=Yorig, daily_mean=daily_mean,
                      min_depth=min_depth, max_gap_days=max_gap_days)
        s = godae_metrics(c["model"], c["obs"])
        rows.append({"date": date, "coverage": c["coverage"],
                     "n_used": c["n_used"], **s})

    if verbose:
        name = source_name(identify(reference))
        print(f"{_label(var, depth_m)}  CROCO vs {name}, collocated:")
        print("   date         cover     n     bias    rmsd   urmsd    corr")
        for r in rows:
            print("   %s  %5.1f%%  %6d  %+6.3f  %6.3f  %6.3f  %6.3f"
                  % (r["date"], r["coverage"], r["n_used"],
                     r["bias"], r["rmsd"], r["urmsd"], r["corr"]))
    return rows


def collocate_days(croco_his, reference, var="temp", days=None, depth_m=None,
                   Yorig=None, daily_mean=True, min_depth=None,
                   rows=("reference", "croco", "difference"),
                   max_gap_days=0.5, vmin=None, vmax=None, dlim=None,
                   cmap=None, dcmap="RdBu_r", title=None,
                   out=None, dpi=130, figsize=None, verbose=True):
    """A grid of collocated comparisons: one column per day, one row per view.

    The counterpart of compare_days() for a **gappy** product. Everything is
    drawn as scattered points at the observed locations, so the gaps stay
    gaps — there is no interpolation of the observations onto the model grid
    and therefore nothing invented where there was cloud.

    days : None (every record), an int (the first N), a (start, stop) tuple,
           a list of dates, or a list of record indices.
    rows : any subset and order of 'reference', 'croco', 'difference'.

    The two field rows share one colour scale and the difference row has its
    own, so a change between days reads as a change rather than as each panel
    being rescaled.

    Returns (figure, stats) — or (path, stats) when `out` is given — with one
    metrics dict per column.
    """
    source = identify(reference)
    if source in ("croco", "insitu", "unknown"):
        raise ValueError(f"collocate_days needs a gridded reference, got {source!r}")
    name = source_name(source)

    valid_rows = ("reference", "croco", "difference")
    rows = tuple(rows)
    bad = [r for r in rows if r not in valid_rows]
    if bad:
        raise ValueError(f"rows may contain {valid_rows}; got {bad}")
    if not rows:
        raise ValueError("rows must name at least one view")

    # which days
    ds = pp.open_history(croco_his, Yorig=Yorig)
    _require_dates(ds, croco_his)
    if daily_mean and "time" in ds.dims:
        ds = val._maybe_daily_mean(ds, None, True)
    times = pp.times(ds)
    ds.close()
    idx = _select_records(times, days)
    if not idx:
        raise ValueError("no records selected")

    cols, stats = [], []
    for k in idx:
        date = str(np.datetime_as_string(times[k], unit="D"))
        c = collocate(croco_his, reference, var=var, date=date,
                      depth_m=depth_m, Yorig=Yorig, daily_mean=daily_mean,
                      min_depth=min_depth, max_gap_days=max_gap_days)
        cols.append(c)
        stats.append(godae_metrics(c["model"], c["obs"]))

    n = len(cols)

    # one scale across every field panel, one across the differences
    allf = np.concatenate([np.concatenate([c["model"], c["obs"]]) for c in cols])
    if vmin is None:
        vmin = float(np.percentile(allf, 1))
    if vmax is None:
        vmax = float(np.percentile(allf, 99))
    alld = np.concatenate([c["model"] - c["obs"] for c in cols])
    if dlim is None:
        dlim = float(np.nanpercentile(np.abs(alld), 99))
    cmap = cmap or _cmap(var)

    label = _label(var, depth_m)
    ylabels = {"reference": f"{name} {label}",
               "croco": f"CROCO {label}",
               "difference": f"CROCO \u2212 {name}"}

    nrow = len(rows)
    if figsize is None:
        figsize = (3.4 * n, 4.2 * nrow)
    fig, ax = plt.subplots(nrow, n, figsize=figsize,
                           constrained_layout=True, squeeze=False)

    # a marker size that keeps the points touching without overlapping
    ms = max(2.0, 40.0 / max(n, 1))

    handles = {}
    for i, row in enumerate(rows):
        for k, c in enumerate(cols):
            a = ax[i][k]
            if row == "difference":
                h = a.scatter(c["lon"], c["lat"], c=c["model"] - c["obs"],
                              s=ms, cmap=dcmap, vmin=-dlim, vmax=dlim)
                a.set_title(f"rmsd {stats[k]['rmsd']:.2f}   "
                            f"bias {stats[k]['bias']:+.2f}", fontsize=9)
            else:
                field = c["obs"] if row == "reference" else c["model"]
                h = a.scatter(c["lon"], c["lat"], c=field, s=ms,
                              cmap=cmap, vmin=vmin, vmax=vmax)
                if i == 0:
                    a.set_title(f"{c['date']}\n{c['coverage']:.0f}% observed",
                                fontsize=10)
            a.set_aspect("equal")
            if k:
                a.set_yticklabels([])
            if i < nrow - 1:
                a.set_xticklabels([])
            else:
                a.set_xlabel("longitude")
        ax[i][0].set_ylabel(ylabels[row])
        handles[row] = h

    field_rows = [i for i, r in enumerate(rows) if r != "difference"]
    if field_rows:
        axes = [ax[i][k] for i in field_rows for k in range(n)]
        fig.colorbar(handles[rows[field_rows[0]]], ax=axes, label=_units(var),
                     shrink=0.85, pad=0.02, extend="both")
    if "difference" in rows:
        i = rows.index("difference")
        fig.colorbar(handles["difference"], ax=[ax[i][k] for k in range(n)],
                     label=_units(var), shrink=0.85, pad=0.02, extend="both")

    fig.suptitle(title or f"CROCO vs {name} \u2014 {label}, collocated",
                 fontsize=13)

    if verbose:
        print(f"{label}  CROCO vs {name}, collocated:")
        print("   date         cover     n     bias    rmsd   urmsd    corr")
        for c, s in zip(cols, stats):
            print("   %s  %5.1f%%  %6d  %+6.3f  %6.3f  %6.3f  %6.3f"
                  % (c["date"], c["coverage"], s["n"], s["bias"],
                     s["rmsd"], s["urmsd"], s["corr"]))

    if out:
        fig.savefig(out, dpi=dpi, bbox_inches="tight")
        plt.close(fig)
        return out, stats
    return fig, stats


def plot_collocation(croco_his, reference, var="temp", date=None,
                     depth_m=None, Yorig=None, daily_mean=True,
                     min_depth=None, vmin=None, vmax=None, cmap=None,
                     title=None, out=None, dpi=130, figsize=(15, 5.5)):
    """Three panels for a gappy product: observations, model, scatter.

    The observation panel shows the gaps as gaps — white where there was
    cloud — rather than hiding them behind an interpolation. The scatter is
    the honest summary: one point per collocated pair, with the metrics
    computed from exactly those pairs.
    """

    c = collocate(croco_his, reference, var=var, date=date, depth_m=depth_m,
                  Yorig=Yorig, daily_mean=daily_mean, min_depth=min_depth)
    s = godae_metrics(c["model"], c["obs"])
    name = source_name(c["source"])
    label = _label(var, depth_m)

    both = np.concatenate([c["model"], c["obs"]])
    if vmin is None:
        vmin = float(np.percentile(both, 1))
    if vmax is None:
        vmax = float(np.percentile(both, 99))
    cmap = cmap or _cmap(var)

    fig, ax = plt.subplots(1, 3, figsize=figsize, constrained_layout=True)

    h0 = ax[0].scatter(c["lon"], c["lat"], c=c["obs"], s=6, cmap=cmap,
                       vmin=vmin, vmax=vmax)
    ax[0].set_title(f"{name} {label}  {c['date']}\n"
                    f"{c['coverage']:.0f}% of the grid observed", fontsize=10)

    ax[1].scatter(c["lon"], c["lat"], c=c["model"], s=6, cmap=cmap,
                  vmin=vmin, vmax=vmax)
    ax[1].set_title(f"CROCO {label}, at those points", fontsize=10)
    fig.colorbar(h0, ax=ax[:2], label=_units(var), shrink=0.85, pad=0.02,
                 extend="both")

    for a in ax[:2]:
        a.set_aspect("equal")
        a.set_xlabel("longitude")
    ax[0].set_ylabel("latitude")

    ax[2].scatter(c["obs"], c["model"], s=4, alpha=0.25, edgecolors="none")
    lims = [vmin, vmax]
    ax[2].plot(lims, lims, "k-", lw=1)
    ax[2].set_xlim(lims); ax[2].set_ylim(lims)
    ax[2].set_aspect("equal")
    ax[2].set_xlabel(f"{name} {_units(var)}")
    ax[2].set_ylabel(f"CROCO {_units(var)}")
    ax[2].set_title(f"n={s['n']}   bias {s['bias']:+.3f}\n"
                    f"rmsd {s['rmsd']:.3f}   corr {s['corr']:.3f}", fontsize=10)
    ax[2].grid(alpha=0.3)

    fig.suptitle(title or f"CROCO vs {name} \u2014 {label}, collocated",
                 fontsize=13)

    if out:
        fig.savefig(out, dpi=dpi, bbox_inches="tight")
        plt.close(fig)
        return out, s
    return fig, s


# ======================================================================
# Persistence — GODAE class 3
#
# The other comparisons answer "how close is the model to the truth?".
# This one answers "is the model worth running?".
#
# Persistence is the forecast you get for free: assume nothing changes,
# and today's state is tomorrow's forecast. A model that cannot beat it
# has not earned its compute. The skill score is
#
#     SS = 1 - MSE_model / MSE_persistence
#
# so SS = 0 means the model matched persistence, SS = 1 means it was
# perfect, and SS < 0 means you would have done better by assuming the
# ocean stood still.
#
# Persistence is a strong baseline at short lead times — the ocean is
# slow, and one day ahead it barely moves — so a low score on day 1 is
# normal and says little. The interesting question is where the curves
# cross: how many days of lead time before the model's advantage
# becomes clear.
# ======================================================================
def persistence(croco_his, reference, var="temp", days=None, depth_m=None,
                Yorig=None, daily_mean=True, min_depth=None, method="auto",
                max_gap_days=0.5, out=None, dpi=130, verbose=True):
    """Forecast skill against persistence, by lead time.

    Persistence here is the **model's own initial state**, held fixed and
    compared against each day's observations — the forecast you would have
    made with no model at all.

    method : 'collocate' compares in observation space (needed for a gappy
             product), 'regrid' brings the reference onto the model grid,
             'auto' picks collocation for anything marked gappy.

    Returns a list of dicts, one per lead day:
        lead_days, date, rmse_model, rmse_persist, skill, n

    where skill = 1 - (rmse_model / rmse_persist)**2. Positive means the
    model beat persistence.
    """
    source = identify(reference)
    if source in ("croco", "insitu", "unknown"):
        raise ValueError(f"persistence needs a gridded reference, got {source!r}")
    name = source_name(source)
    if method == "auto":
        method = "collocate" if REFERENCES[source].get("gappy") else "regrid"

    ds = pp.open_history(croco_his, Yorig=Yorig)
    _require_dates(ds, croco_his)
    if daily_mean and "time" in ds.dims:
        ds = val._maybe_daily_mean(ds, None, True)
    times = pp.times(ds)
    ds.close()
    idx = _select_records(times, days)
    if len(idx) < 2:
        raise ValueError("persistence needs at least two records")

    t0 = idx[0]
    rows = []

    if method == "collocate":
        # The initial field, sampled at each day's observation points. The
        # observation locations move as the cloud moves, so persistence has
        # to be re-sampled per day rather than computed once.
        for k in idx:
            date = str(np.datetime_as_string(times[k], unit="D"))
            now = collocate(croco_his, reference, var=var, date=date,
                            depth_m=depth_m, Yorig=Yorig, daily_mean=daily_mean,
                            min_depth=min_depth, max_gap_days=max_gap_days)
            # the same observations, against the model's day-zero state
            init = _collocate_model_only(croco_his, var, t0, depth_m, Yorig,
                                         daily_mean, now["lon"], now["lat"])
            ok = np.isfinite(init) & np.isfinite(now["model"]) & np.isfinite(now["obs"])
            sm = godae_metrics(now["model"][ok], now["obs"][ok])
            sp = godae_metrics(init[ok], now["obs"][ok])
            rows.append(_skill_row(times[k], times[t0], date, sm, sp))
    else:
        for k in idx:
            date = str(np.datetime_as_string(times[k], unit="D"))
            ds = pp.open_history(croco_his, Yorig=Yorig)
            if daily_mean and "time" in ds.dims:
                ds = val._maybe_daily_mean(ds, None, True)
            now_field = _croco_field(ds, var, k, depth_m)
            init_field = _croco_field(ds, var, t0, depth_m)
            ref = load_reference(reference, var, ds, date=date, depth_m=depth_m,
                                 source=source, max_gap_days=max_gap_days)
            ds.close()
            if var == "ssh":
                now_field = now_field - np.nanmean(now_field)
                init_field = init_field - np.nanmean(init_field)
                ref = ref - np.nanmean(ref)
            sm = godae_metrics(now_field, ref)
            sp = godae_metrics(init_field, ref)
            rows.append(_skill_row(times[k], times[t0], date, sm, sp))

    if verbose:
        print(f"{_label(var, depth_m)}  CROCO vs {name}, against persistence:")
        print("   lead  date         n      model  persist   skill")
        for r in rows:
            print("   %4.1f  %s  %6d  %6.3f   %6.3f  %+7.3f"
                  % (r["lead_days"], r["date"], r["n"], r["rmse_model"],
                     r["rmse_persist"], r["skill"]))
        # a skill of a few thousandths is the two RMSEs agreeing to
        # three decimals, not the model winning
        pos = [r for r in rows[1:] if r["skill"] > 0.01]
        if pos:
            print(f"   the model first beats persistence at lead "
                  f"{pos[0]['lead_days']:.0f} days")
        else:
            print("   the model does not beat persistence in this window")

    if out:
        fig, ax = plt.subplots(1, 2, figsize=(13, 5), constrained_layout=True)
        lead = [r["lead_days"] for r in rows]
        ax[0].plot(lead, [r["rmse_model"] for r in rows], "o-", lw=1.8,
                   label="CROCO")
        ax[0].plot(lead, [r["rmse_persist"] for r in rows], "s--", lw=1.8,
                   label="persistence")
        ax[0].set_xlabel("lead time (days)")
        ax[0].set_ylabel(f"RMSE ({_units(var)})")
        ax[0].set_title(f"{_label(var, depth_m)} vs {name}")
        ax[0].legend(); ax[0].grid(alpha=0.3)

        ax[1].axhline(0, color="k", lw=1)
        ax[1].plot(lead, [r["skill"] for r in rows], "o-", lw=1.8, color="C2")
        ax[1].set_xlabel("lead time (days)")
        ax[1].set_ylabel("skill score")
        ax[1].set_title("1 - (RMSE$_{model}$ / RMSE$_{persistence}$)$^2$\n"
                        "above zero, the model beats persistence", fontsize=10)
        ax[1].grid(alpha=0.3)
        fig.savefig(out, dpi=dpi, bbox_inches="tight")
        plt.close(fig)
        return out, rows
    return rows


def _skill_row(t, t0, date, sm, sp):
    """One lead-day row: the model's error, persistence's, and the score."""
    lead = float((np.datetime64(t) - np.datetime64(t0)) / np.timedelta64(1, "D"))
    rm, rp = sm["rmsd"], sp["rmsd"]
    skill = float(1.0 - (rm / rp) ** 2) if rp > 1e-12 else np.nan
    return {"lead_days": lead, "date": date, "n": sm["n"],
            "rmse_model": rm, "rmse_persist": rp, "skill": skill,
            "bias_model": sm["bias"], "bias_persist": sp["bias"],
            "corr_model": sm["corr"], "corr_persist": sp["corr"]}


def _collocate_model_only(croco_his, var, tindex, depth_m, Yorig, daily_mean,
                          lon, lat):
    """Sample one model record at given points — the persistence field.

    Same masking as collocate(): points without wet model cells around them
    come back NaN rather than being interpolated across.
    """
    from scipy.interpolate import griddata as _griddata

    ds = pp.open_history(croco_his, Yorig=Yorig)
    if daily_mean and "time" in ds.dims:
        ds = val._maybe_daily_mean(ds, None, True)
    field = _croco_field(ds, var, tindex, depth_m)
    clon, clat, cmask = pp.lonlatmask(ds)
    ds.close()

    good = np.isfinite(field)
    tgt = np.column_stack([lon, lat])
    out = _griddata(np.column_stack([clon[good].ravel(), clat[good].ravel()]),
                    field[good].ravel(), tgt, method="linear")
    wet = _griddata(np.column_stack([clon.ravel(), clat.ravel()]),
                    (np.asarray(cmask) > 0).astype(float).ravel(),
                    tgt, method="linear")
    return np.where(np.isfinite(wet) & (wet > 0.99), out, np.nan)


def _croco_field(ds, var, tindex, depth_m):
    """A 2D CROCO field, at the surface or a true depth."""
    _, _, mask = pp.lonlatmask(ds)
    if var == "ssh":
        if depth_m is not None:
            raise ValueError("ssh is two-dimensional; depth_m does not apply")
        return ds["zeta"].isel(time=tindex).values * mask
    if var in ("u", "v"):
        # CROCO stores velocity on the C-grid: u on the eastern cell faces,
        # v on the northern ones, so neither shares the rho grid's shape.
        # Move both to rho points, then rotate from grid-relative to true
        # east/north — the reference products are geographic.
        if depth_m is None:
            ur, vr = pp.surface_uv(ds, tindex=tindex)
            ue, vn = pp.rotate_uv(ds, ur, vr)
        else:
            ue, vn = pp.uv_at_depth(ds, depth_m, tindex=tindex, rotate=True)
        f = np.squeeze(ue) if var == "u" else np.squeeze(vn)
        return f * mask
    if var == "speed":
        if depth_m is None:
            ur, vr = pp.surface_uv(ds, tindex=tindex)
            ue, vn = pp.rotate_uv(ds, ur, vr)
        else:
            ue, vn = pp.uv_at_depth(ds, depth_m, tindex=tindex, rotate=True)
        return pp.speed(np.squeeze(ue), np.squeeze(vn))
    if depth_m is None:
        return pp.surface(ds, var, tindex=tindex).values
    return pp.field_at_depth(ds, var, depth_m, tindex=tindex).values


_LABELS = {"temp": "temperature", "salt": "salinity",
           "u": "eastward velocity", "v": "northward velocity",
           "ssh": "SSH anomaly", "speed": "speed", "mld": "mixed-layer depth"}
_UNITS = {"temp": "temperature (\u00b0C)", "salt": "salinity (PSU)",
          "u": "u (m s$^{-1}$)", "v": "v (m s$^{-1}$)",
          "ssh": "SSH' (m)", "speed": "speed (m s$^{-1}$)", "mld": "MLD (m)"}
_CMAPS = {"temp": "RdYlBu_r", "salt": "viridis",
          "u": "RdBu_r", "v": "RdBu_r",
          "ssh": "viridis", "speed": "viridis", "mld": "viridis"}


def _label(var, depth_m):
    if var == "temp" and depth_m is None:
        return "SST"
    base = _LABELS.get(var, var)
    return base if depth_m is None else f"{base} {depth_m:g} m"


def _units(var):
    return _UNITS.get(var, var)


def _cmap(var):
    return _CMAPS.get(var, "viridis")


# ======================================================================
# What each reference can answer
# ======================================================================
def describe(source=None):
    """Print what a reference provides, or all of them."""
    keys = [source] if source else sorted(REFERENCES)
    for k in keys:
        s = REFERENCES[k]
        depth = "surface and depth" if s["has_depth"] else "surface only"
        print(f"{k}")
        print(f"   variables : {', '.join(sorted(s['vars']))}   ({depth})")
        if s.get("nrt"):
            print(f"   nrt       : {s['nrt']}")
        if s.get("my"):
            print(f"   my        : {s['my']}")
        if s.get("note"):
            print(f"   note      : {s['note']}")
        print()
