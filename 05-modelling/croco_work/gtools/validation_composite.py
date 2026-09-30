"""
gtools/validation_composite.py — CROCO COMPOSITE (multi-cycle)
validation.

02_validation.ipynb / gtools.validation validate ONE forecast cycle,
comparing CROCO against the parent/satellite reference DAY BY CALENDAR
DATE. This module instead merges SEVERAL cycles into one analysis,
indexed by FORECAST LEAD TIME rather than calendar date:

    cycle 20260701's days -> sp1, sp2, fcst1, fcst2, fcst3
    cycle 20260706's days -> sp1, sp2, fcst1, fcst2, fcst3
    cycle 20260711's days -> sp1, sp2, fcst1, fcst2, fcst3, fcst4
    ...
        |
        v  merged/pooled BY LABEL (not by calendar date)
    sp1 : {20260701's day0} + {20260706's day0} + {20260711's day0} + ...
    sp2 : {20260701's day1} + {20260706's day1} + {20260711's day1} + ...
    fcst1 : {20260701's day2} + {20260706's day2} + {20260711's day2} + ...
    fcst2 : ...
    fcst3 : ...
    fcst4 : {20260711's day5} + ... (only cycles that reach this far)

so "how skillful is the forecast on its 3rd day" can be assessed across
every cycle in the composite at once, rather than one cycle's absolute
dates. The first `spinup_days` days of every cycle are labelled 'sp1',
'sp2', ... (CROCO's model-startup transient, same convention as
02_validation.ipynb Section 6) and are kept OUT of the fcst* labels so
they can be excluded from pass/fail scoring the same way.

Every function here is a thin per-lead ACCUMULATOR around the existing
gtools.validation building blocks -- it calls them once per contributing
cycle and pools (concatenates / stacks-and-means) the results, rather than
reimplementing any loading, regridding, statistics or plotting logic. This
guarantees a composite run and the single-cycle 02_validation.ipynb
notebook can never silently disagree on methodology.

Typical use
-----------
    import gtools.postprocess as pp
    import gtools.validation as val
    import gtools.validation_composite as vc

    cycles_info = []
    for cycle in ["20260701", "20260706", "20260711"]:
        croco_his, reference, _ = _paths.get_paths(cycle, config=CONFIG, main_dir=MAIN_DIR)
        days, day_label, label_day = vc.label_cycle_days(croco_his, Yorig=YORIG, spinup_days=2)
        cycles_info.append({"cycle": cycle, "croco_his": croco_his,
                            "reference": reference, "days": days,
                            "day_label": day_label, "label_day": label_day,
                            "sat_files": {"OSTIA": {...}, "ODYSSEA": {...}, "SMOS": {...}}})

    LEADS = vc.all_lead_labels(cycles_info)   # ['sp1','sp2','fcst1','fcst2','fcst3','fcst4']

    vc.compare_composite_field(cycles_info, var="temp", lead="fcst2", Yorig=2000,
                               out="sst_fcst2_composite.png")
"""
from __future__ import annotations
import os
import numpy as np
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt

import gtools.postprocess as pp
import gtools.validation as val
from gtools.plotting import percentile_clim


# ======================================================================
# Lead-time labelling
# ======================================================================
def lead_label(day_index, spinup_days=2):
    """The lead-time label for the day_index'th (0-based) day of a cycle:
    'sp1','sp2',... for the first `spinup_days` days (model spin-up),
    then 'fcst1','fcst2',... for every day after that (genuine forecast
    lead time 1, 2, ...), matching the SPINUP_DAYS convention already used
    in 02_validation.ipynb Section 6."""
    if day_index < spinup_days:
        return f"sp{day_index + 1}"
    return f"fcst{day_index - spinup_days + 1}"


def label_cycle_days(croco_his, Yorig=None, spinup_days=2):
    """Every distinct calendar day this CROCO history file covers (sorted),
    each assigned a lead-time label via lead_label().

    Returns (days, day_to_label, label_to_day):
      days         : sorted ['YYYY-MM-DD', ...] -- same convention as
                     CYCLE_DAYS in 02_validation.ipynb.
      day_to_label : {day: label}
      label_to_day : {label: day} (inverse; unique since days are unique)
    """
    ds = pp.open_history(croco_his, Yorig=Yorig)
    days = sorted({str(d) for d in pp.times(ds).astype("datetime64[D]")})    
    ds.close()
    day_to_label, label_to_day = {}, {}
    for i, day in enumerate(days):
        lbl = lead_label(i, spinup_days=spinup_days)
        day_to_label[day] = lbl
        label_to_day[lbl] = day
    return days, day_to_label, label_to_day


def _lead_sort_key(lbl):
    """Sort key for lead labels: 'sp' before 'fcst', numeric suffix order
    (so 'fcst10' sorts after 'fcst2', unlike a plain string sort)."""
    kind = 0 if lbl.startswith("sp") else 1
    num = int(lbl[2:]) if lbl.startswith("sp") else int(lbl[4:])
    return (kind, num)


def all_lead_labels(cycles_info, fcst_only=False):
    """Every lead-time label present in ANY cycle of cycles_info, sorted
    sp1, sp2, ..., fcst1, fcst2, ... (see _lead_sort_key).

    fcst_only : if True, drop the 'sp*' (spin-up) labels -- use this to
                build the day list a pass/fail summary should be scored
                over, excluding spin-up (mirrors 02_validation.ipynb
                Section 6's SPINUP_DAYS exclusion).
    """
    labels = set()
    for c in cycles_info:
        labels.update(c["label_day"].keys())
    if fcst_only:
        labels = {l for l in labels if l.startswith("fcst")}
    return sorted(labels, key=_lead_sort_key)


def cycles_reaching(cycles_info, lead):
    """The list of cycle names in cycles_info that have a day at `lead`."""
    return [c["cycle"] for c in cycles_info if lead in c["label_day"]]


def composite_domain_series(cycles_info, var, leads, depth_m=None, tindex=-1,
                            Yorig=None, daily_mean=True):
    """Per-lead (croco_values, parent_values), domain-mean of the SAME
    already co-registered val._field_pair_for_day() pair used by
    composite_domain_diff() -- paired per cycle (not pooled bias), so
    CROCO and parent can be plotted as two separate mean+/-std bands."""
    croco_by_lead, parent_by_lead = [], []
    for lead in leads:
        cv, pv = [], []
        for c in cycles_info:
            day = c["label_day"].get(lead)
            if day is None:
                continue
            try:
                ds = pp.open_history(c["croco_his"], Yorig=Yorig)
                croco2d, parent2d = val._field_pair_for_day(
                    ds, c["reference"], var, day, depth_m, tindex, daily_mean, margin_deg=None)
                ds.close()
            except Exception as e:
                print(f"  [{lead}] cycle {c['cycle']} ({day}): domain series failed ({e}) - skipping")
                continue
            cv.append(np.nanmean(croco2d)); pv.append(np.nanmean(parent2d))
        croco_by_lead.append(np.array(cv)); parent_by_lead.append(np.array(pv))
    return croco_by_lead, parent_by_lead


def composite_point_series(cycles_info, var, lon0, lat0, leads, depth_m=None, tindex=-1,
                           Yorig=None, daily_mean=True):
    """Same as composite_domain_series(), but sampled at the nearest CROCO
    grid cell to (lon0, lat0) instead of domain-averaged -- reuses the
    same co-registered field pair, via pp.nearest_index()."""
    croco_by_lead, parent_by_lead = [], []
    for lead in leads:
        cv, pv = [], []
        for c in cycles_info:
            day = c["label_day"].get(lead)
            if day is None:
                continue
            try:
                ds = pp.open_history(c["croco_his"], Yorig=Yorig)
                croco2d, parent2d = val._field_pair_for_day(
                    ds, c["reference"], var, day, depth_m, tindex, daily_mean, margin_deg=None)
                j, i = pp.nearest_index(ds, lon0, lat0)
                ds.close()
            except Exception as e:
                print(f"  [{lead}] cycle {c['cycle']} ({day}): point series failed ({e}) - skipping")
                continue
            cvv, pvv = croco2d[j, i], parent2d[j, i]
            if np.isfinite(cvv) and np.isfinite(pvv):
                cv.append(cvv); pv.append(pvv)
        croco_by_lead.append(np.array(cv)); parent_by_lead.append(np.array(pv))
    return croco_by_lead, parent_by_lead


def composite_satellite_series(cycles_info, product, leads, lon0=None, lat0=None,
                               Yorig=None, margin_deg=None):
    """Per-lead (croco_values, satellite_values) -- domain-mean by default,
    or sampled at the nearest grid point to (lon0, lat0) if given. Reuses
    the SAME loading/regridding as composite_domain_diff_satellite()."""
    croco_by_lead, sat_by_lead = [], []
    for lead in leads:
        cv, sv = [], []
        for c in cycles_info:
            day = c["label_day"].get(lead)
            if day is None:
                continue
            fname = c.get("sat_files", {}).get(product, {}).get(day)
            if not fname or not os.path.exists(fname):
                continue
            loaded = val.load_satellite_field(fname, product, date=day)
            if loaded is None:
                continue
            plon, plat, pfield = loaded
            pfield = np.squeeze(np.asarray(pfield))
            if pfield.shape != plon.shape:
                continue
            try:
                clon, clat, croco_field = val._croco_field_satellite(
                    c["croco_his"], product, day, Yorig=Yorig, margin_deg=margin_deg)
            except Exception as e:
                print(f"  [{product}][{lead}] cycle {c['cycle']} ({day}): "
                     f"CROCO record not available ({e}) - skipping")
                continue
            ds_grid = xr.Dataset(
                {"mask_rho": (("eta_rho", "xi_rho"), np.where(np.isfinite(croco_field), 1.0, np.nan))},
                coords={"lon_rho": (("eta_rho", "xi_rho"), clon),
                        "lat_rho": (("eta_rho", "xi_rho"), clat)})
            sat_on_croco = val.regrid_to_croco(plon, plat, pfield, ds_grid)
            if lon0 is not None and lat0 is not None:
                j, i = np.unravel_index(np.argmin((clon - lon0)**2 + (clat - lat0)**2), clon.shape)
                cvv, svv = croco_field[j, i], sat_on_croco[j, i]
            else:
                cvv, svv = np.nanmean(croco_field), np.nanmean(sat_on_croco)
            if np.isfinite(cvv) and np.isfinite(svv):
                cv.append(cvv); sv.append(svv)
        croco_by_lead.append(np.array(cv)); sat_by_lead.append(np.array(sv))
    return croco_by_lead, sat_by_lead


def composite_two_series_timeseries(series_by_var, leads, n_cycles, out=None,
                                    labels=("CROCO", "parent")):
    """Stacked-panel mean +/- 1 std timeseries, TWO series overlaid per
    panel (CROCO vs parent, or CROCO vs a satellite product), one row per
    variable/product.
    series_by_var : [(ylabel, title, series_a_by_lead, series_b_by_lead), ...]
    """
    fig, axes = plt.subplots(len(series_by_var), 1,
                             figsize=(max(6, 1.1 * len(leads)), 3.3 * len(series_by_var)),
                             sharex=True)
    axes = np.atleast_1d(axes)
    x = np.arange(1, len(leads) + 1)
    colors = ("C0", "C3")
    for ax, (ylabel, title, sa, sb) in zip(axes, series_by_var):
        for s, lab, col in zip((sa, sb), labels, colors):
            m = np.array([np.mean(v) if v is not None and len(v) > 0 else np.nan for v in s])
            sd = np.array([np.std(v) if v is not None and len(v) > 0 else np.nan for v in s])
            ax.plot(x, m, "o-", color=col, lw=1.5, ms=5, label=lab)
            ax.fill_between(x, m - sd, m + sd, color=col, alpha=0.2)
        ax.set_xticks(x); ax.set_xticklabels(leads)
        ax.set_ylabel(ylabel); ax.set_title(title); ax.grid(alpha=0.3)
    axes[0].legend(fontsize=8)
    axes[-1].set_xlabel("lead time")
    plt.setp(axes[-1].get_xticklabels(), rotation=45, ha="right")
    fig.suptitle(f"composite of {n_cycles} cycles")
    fig.tight_layout()
    return val._save_or_return(fig, out)

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
    
def _interp_profile(depth_grid, z_native, v_native):
    z_native, v_native = np.asarray(z_native), np.asarray(v_native)
    good = np.isfinite(z_native) & np.isfinite(v_native)
    z_native, v_native = z_native[good], v_native[good]
    order = np.argsort(z_native)
    z_native, v_native = z_native[order], v_native[order]
    if z_native.size < 2:
        return np.full(len(depth_grid), np.nan)
    return np.interp(depth_grid, z_native, v_native, left=np.nan, right=np.nan)


def _plot_two_profiles(depth_grid, croco_profiles, parent_profiles, used, var, lead,
                       title_suffix, out=None):
    cs, ps = np.array(croco_profiles), np.array(parent_profiles)
    cm, cstd = np.nanmean(cs, axis=0), np.nanstd(cs, axis=0)
    pm, pstd = np.nanmean(ps, axis=0), np.nanstd(ps, axis=0)
    n = len(croco_profiles)

    fig, ax = plt.subplots(figsize=(4.5, 6))
    ax.plot(cm, depth_grid, "o-", color="C0", ms=3, lw=1.6, label="CROCO")
    ax.fill_betweenx(depth_grid, cm - cstd, cm + cstd, color="C0", alpha=0.2)
    ax.plot(pm, depth_grid, "s--", color="C3", ms=3, lw=1.6, label="parent")
    ax.fill_betweenx(depth_grid, pm - pstd, pm + pstd, color="C3", alpha=0.2)
    unit = {"temp": "degC", "salt": "PSU", "speed": "m s$^{-1}$"}[var]
    ax.set_xlabel(f"{var} ({unit})"); ax.set_ylabel("depth (m)")
    ax.grid(alpha=0.3); ax.legend(fontsize=8)
    ax.set_title(f"{var} profile @ {lead}  {title_suffix}\ncomposite of {n} cycle(s): {used}")
    return val._save_or_return(fig, out)


def composite_point_profile(cycles_info, var, lon0, lat0, lead, depth_grid=None,
                            tindex=-1, Yorig=None, daily_mean=False, out=None):
    """CROCO vs parent vertical profile at one point, mean +/- 1 std
    across every cycle reaching `lead` (fill_between for BOTH).
    var must be 'temp', 'salt', or 'speed'."""
    if var not in ("temp", "salt", "speed"):
        raise ValueError("composite_point_profile only supports var='temp', 'salt', or 'speed'")
    if depth_grid is None:
        depth_grid = np.linspace(0, -500, 50)

    croco_profiles, parent_profiles, used = [], [], []
    for c in cycles_info:
        day = c["label_day"].get(lead)
        if day is None:
            continue
        try:
            ds = pp.open_history(c["croco_his"], Yorig=Yorig)
            ti = tindex
            if daily_mean:
                ds = val._maybe_daily_mean(ds, day, daily_mean)
            else:
                ti = val._tindex_for_date(ds, day, tindex)
            cpro = pp.profile(ds, var, lon0, lat0, tindex=ti)   # 'speed' handled natively by pp.profile
            cz, cv = cpro["depth"].values, cpro.values
            ds.close()
        except Exception as e:
            print(f"  [{lead}] cycle {c['cycle']} ({day}): CROCO profile failed ({e}) - skipping")
            continue
        try:
            if var == "speed":
                pz, pu = val._parent_column(c["reference"], "u", lon0, lat0, date=day)
                _, pv = val._parent_column(c["reference"], "v", lon0, lat0, date=day)
                pv_vals = np.sqrt(pu**2 + pv**2)
            else:
                pz, pv_vals = val._parent_column(c["reference"], var, lon0, lat0, date=day)
        except Exception as e:
            print(f"  [{lead}] cycle {c['cycle']} ({day}): parent profile failed ({e}) - skipping")
            continue
        croco_profiles.append(_interp_profile(depth_grid, cz, cv))
        parent_profiles.append(_interp_profile(depth_grid, pz, pv_vals))
        used.append(c["cycle"])
    if not croco_profiles:
        raise ValueError(f"no cycle reaches lead '{lead}' for a profile at ({lon0}, {lat0})")
    return _plot_two_profiles(depth_grid, croco_profiles, parent_profiles, used, var, lead,
                              f"({lon0:.2f}, {lat0:.2f})", out=out)


def composite_domain_profile(cycles_info, var, lead, depth_grid=None, tindex=-1,
                             Yorig=None, daily_mean=False, out=None):
    """CROCO vs parent DOMAIN-MEAN vertical profile, mean +/- 1 std across
    every cycle reaching `lead`. var must be 'temp', 'salt', or 'speed'.
    For 'speed': per sigma/depth level, u/v are combined into a speed
    field FIRST (rotated to east/north via pp.rotate_uv), THEN
    horizontally averaged -- not the speed of the domain-mean u/v."""
    if var not in ("temp", "salt", "speed"):
        raise ValueError("composite_domain_profile only supports var='temp', 'salt', or 'speed'")
    if depth_grid is None:
        depth_grid = np.linspace(0, -500, 50)

    croco_profiles, parent_profiles, used = [], [], []
    for c in cycles_info:
        day = c["label_day"].get(lead)
        if day is None:
            continue
        try:
            ds = pp.open_history(c["croco_his"], Yorig=Yorig)
            ti = tindex
            if daily_mean:
                ds = val._maybe_daily_mean(ds, day, daily_mean)
            else:
                ti = val._tindex_for_date(ds, day, tindex)
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
            cv = np.nanmean(v3m, axis=(1, 2))
            ds.close()
        except Exception as e:
            print(f"  [{lead}] cycle {c['cycle']} ({day}): CROCO domain profile failed ({e}) - skipping")
            continue
        try:
            ds_p = xr.open_dataset(c["reference"])
            if var == "speed":
                du = ds_p[val.PARENT_VARS["u"]]
                dv = ds_p[val.PARENT_VARS["v"]]
                if "time" in du.dims:
                    du = du.sel(time=np.datetime64(day), method="nearest")
                    dv = dv.sel(time=np.datetime64(day), method="nearest")
                da = np.sqrt(du**2 + dv**2)
            else:
                cmems = val.PARENT_VARS[var]
                da = ds_p[cmems]
                if "time" in da.dims:
                    da = da.sel(time=np.datetime64(day), method="nearest")
            pz = -np.abs(ds_p["depth"].values)
            pv = da.mean(dim=[d for d in da.dims if d != "depth"], skipna=True).values
            ds_p.close()
        except Exception as e:
            print(f"  [{lead}] cycle {c['cycle']} ({day}): parent domain profile failed ({e}) - skipping")
            continue
        croco_profiles.append(_interp_profile(depth_grid, cz, cv))
        parent_profiles.append(_interp_profile(depth_grid, pz, pv))
        used.append(c["cycle"])
    if not croco_profiles:
        raise ValueError(f"no cycle reaches lead '{lead}' for a domain profile")
    return _plot_two_profiles(depth_grid, croco_profiles, parent_profiles, used, var, lead,
                              "(full domain)", out=out)
    
# ======================================================================
# Composite CROCO-vs-Copernicus field maps (composite counterpart to
# gtools.validation.compare_sst/sss/ssh/currents), ONE FIGURE PER LEAD,
# composited (mean/RMSE) across every contributing cycle.
# ======================================================================
def compare_composite_field(cycles_info, var, lead, Yorig=None, margin_deg=None,
                            daily_mean=False, depth_m=None, tindex=-1, out=None):
    """CROCO vs Copernicus Marine Forecast at ONE lead-time label (e.g.
    'fcst2'), COMPOSITED across every cycle in cycles_info that reaches
    that lead. 2x2 figure: composite-MEAN CROCO, composite-MEAN
    Copernicus, bias (mean of each contributing cycle's diff), RMSE
    (sqrt of the mean of each cycle's squared diff) -- mathematically the
    SAME accumulation as gtools.validation._multiday_bias_rmse, just one
    sample per CONTRIBUTING CYCLE instead of one sample per day within a
    single cycle.

    var  : 'temp', 'salt', 'ssh', or 'speed' (see
           gtools.validation._field_pair_for_day).
    lead : lead-time label (see label_cycle_days()/all_lead_labels()).
           Cycles that don't reach this lead (shorter forecasts, or a
           cycle still missing its Copernicus reference) are skipped --
           printed, not raised; not every cycle needs to be the same
           length. Cycles whose CROCO grid shape doesn't match the first
           contributing cycle's (e.g. a different CONFIG accidentally
           mixed in) are likewise skipped, not silently mis-stacked.

    Returns (fig_or_path, stats_dict), or None if no cycle could
    contribute at this lead.
    """
    if var not in ("temp", "salt", "ssh", "speed"):
        raise ValueError(f"var={var!r} must be 'temp', 'salt', 'ssh', or 'speed'")
    label = {"temp": "temperature", "salt": "salinity", "ssh": "SSH anomaly",
             "speed": "speed"}[var]
    units = {"temp": "\u00b0C", "salt": "PSU", "ssh": "m", "speed": "m s$^{-1}$"}[var]
    cmap = "viridis" if var in ("salt", "speed") else ("Spectral_r" if var == "ssh" else "RdYlBu_r")

    pairs = []           # (cycle_name, croco, parent_on_croco)
    clon = clat = None
    for c in cycles_info:
        day = c["label_day"].get(lead)
        if day is None:
            continue
        try:
            ds = pp.open_history(c["croco_his"], Yorig=Yorig)
            ds_grid = ds if margin_deg is None else pp.crop_interior(ds, margin_deg=margin_deg)
            this_clon, this_clat, _ = pp.lonlatmask(ds_grid)
            croco, parent_on_croco = val._field_pair_for_day(
                ds, c["reference"], var, day, depth_m, tindex, daily_mean, margin_deg)
            ds.close()
        except Exception as e:
            print(f"  [{lead}] cycle {c['cycle']} ({day}): could not build pair ({e}) - skipping")
            continue
        if clon is None:
            clon, clat = this_clon, this_clat
        elif croco.shape != clon.shape:
            print(f"  [{lead}] cycle {c['cycle']}: grid shape {croco.shape} != "
                 f"{clon.shape} - skipping (different config/grid?)")
            continue
        pairs.append((c["cycle"], croco, parent_on_croco))

    if not pairs:
        print(f"compare_composite_field: no cycle reaches lead={lead!r} for {var!r}")
        return None

    croco_stack = np.stack([p[1] for p in pairs], axis=0)
    parent_stack = np.stack([p[2] for p in pairs], axis=0)
    diff_stack = croco_stack - parent_stack

    croco_mean = np.nanmean(croco_stack, axis=0)
    parent_mean = np.nanmean(parent_stack, axis=0)
    bias_map = np.nanmean(diff_stack, axis=0)
    rmse_map = np.sqrt(np.nanmean(diff_stack ** 2, axis=0))
    n_cycles = len(pairs)
    used = [p[0] for p in pairs]

    s = val.domain_statistics(croco_mean, parent_mean)
    print(f"{label} @ {lead}  (composite of {n_cycles} cycle(s): {used}):")
    val._print_stats(label, s)
    s.update({"variable": var, "lead": lead, "n_cycles": n_cycles, "cycles": used})

    bias_limit = val.BIAS_RMSE_LIMITS.get(var)
    title = f"{label} @ {lead}  ({n_cycles} cycle{'s' if n_cycles != 1 else ''})"
    # n_days intentionally NOT passed through (would mislabel the RMSE
    # panel "N days" when this is actually N CYCLES) -- the cycle count is
    # already in `title`, shown on every panel + the suptitle instead.
    res = val._four_panel(clon, clat, croco_mean, parent_mean, bias_map, rmse_map,
                          title, f"{label} ({units})", cmap=cmap,
                          bias_limit=bias_limit, out=out, n_days=None)
    return res, s


# ======================================================================
# Composite scatter (composite counterpart to Section 3 / scatter_vs_reference)
# ======================================================================
def scatter_composite(cycles_info, lead, Yorig=None, daily_mean=True, out=None):
    """SST + SSH pointwise scatter at ONE lead-time label, with points
    POOLED (concatenated) from every cycle that reaches that lead -- the
    composite counterpart to 02_validation.ipynb Section 3's per-day
    scatter. Reuses gtools.validation._field_pair_for_day for the actual
    CROCO/parent field-building (daily_mean=True averages that day's
    sub-daily CROCO records first, matching compare_composite_field), then
    gtools.validation.scatter_vs_reference for the plotting/stats on the
    pooled arrays.
    """
    sst_m, sst_r, ssh_m, ssh_r, sss_m, sss_r, used = [], [], [], [], [], [], []
    for c in cycles_info:
        day = c["label_day"].get(lead)
        if day is None:
            continue
        try:
            ds = pp.open_history(c["croco_his"], Yorig=Yorig)
            sst, sst_ref = val._field_pair_for_day(
                ds, c["reference"], "temp", day, None, -1, daily_mean, margin_deg=None)
            sss, sss_ref = val._field_pair_for_day(
                ds, c["reference"], "salt", day, None, -1, daily_mean, margin_deg=None)
            ssh, ssh_ref = val._field_pair_for_day(
                ds, c["reference"], "ssh", day, None, -1, daily_mean, margin_deg=None)
            ds.close()
        except Exception as e:
            print(f"  [{lead}] cycle {c['cycle']} ({day}): scatter pair failed ({e}) - skipping")
            continue
        sst_m.append(sst.ravel()); sst_r.append(sst_ref.ravel())
        sss_m.append(sss.ravel()); sss_r.append(sss_ref.ravel())
        ssh_m.append(ssh.ravel()); ssh_r.append(ssh_ref.ravel())
        used.append(c["cycle"])

    if not used:
        print(f"scatter_composite: no cycle reaches lead={lead!r}")
        return None

    fig, axes = plt.subplots(1, 3, figsize=(16, 6))
    val.scatter_vs_reference(np.concatenate(sst_m), np.concatenate(sst_r),
                             "SST", "degC", ax=axes[0])
    val.scatter_vs_reference(np.concatenate(sss_m), np.concatenate(sss_r),
                             "SSS", "PSU", ax=axes[1])
    val.scatter_vs_reference(np.concatenate(ssh_m), np.concatenate(ssh_r),
                             "SSH", "m", ax=axes[2])
    fig.suptitle(f"{lead}  (composite of {len(used)} cycle(s): \n{used})")
    fig.tight_layout()
    return val._save_or_return(fig, out)


# ======================================================================
# Composite GODAE scorecard (composite counterpart to
# gtools.validation.godae_scorecard_croco_vs_glorys)
# ======================================================================
def godae_scorecard_composite(cycles_info, var, lead, depth_m=None, Yorig=None,
                              tindex=-1, daily_mean=False):
    """GODAE scorecard (bias/RMSD/uRMSD/corr/SI/std_ratio -- see
    gtools.validation.godae_metrics) at ONE lead-time label, computed on
    points POOLED from every cycle that reaches that lead.

    Reuses gtools.validation._field_pair_for_day for the exact same
    CROCO/parent field-building logic used by compare_sst/ssh/sss/currents
    and domain_diff (including SSH's domain-mean-removed anomaly
    treatment and speed's rotate-to-earth-relative treatment), then scores
    the pooled points with godae_metrics() -- same metric definitions as
    godae_scorecard_croco_vs_glorys, pooled across cycles at this one lead
    instead of scored for one cycle/day at a time.

    Returns a dict from godae_metrics(), plus 'variable', 'lead',
    'depth_m', 'n_cycles', 'cycles' -- or None if no cycle contributes.
    """
    model_all, ref_all, used = [], [], []
    for c in cycles_info:
        day = c["label_day"].get(lead)
        if day is None:
            continue
        try:
            ds = pp.open_history(c["croco_his"], Yorig=Yorig)
            model, ref_on_croco = val._field_pair_for_day(
                ds, c["reference"], var, day, depth_m, tindex, daily_mean, margin_deg=None)
            ds.close()
        except Exception as e:
            print(f"  [{lead}] cycle {c['cycle']} ({day}): could not score {var!r} ({e}) - skipping")
            continue
        model_all.append(np.asarray(model).ravel())
        ref_all.append(np.asarray(ref_on_croco).ravel())
        used.append(c["cycle"])

    if not used:
        return None
    s = val.godae_metrics(np.concatenate(model_all), np.concatenate(ref_all))
    s.update({"variable": var, "lead": lead, "depth_m": depth_m,
             "n_cycles": len(used), "cycles": used})
    return s


# ======================================================================
# Composite domain-wide differences (for boxplots) -- composite
# counterparts to gtools.validation.domain_diff / domain_diff_satellite.
# Feed the result straight into val.bias_boxplot()/bias_boxplot_multi()
# with LEAD labels instead of calendar days -- no separate composite
# boxplot-drawing function is needed, the existing ones already work on
# any x-axis label.
# ======================================================================
def composite_domain_diff(cycles_info, var, lead, depth_m=None, tindex=-1,
                          Yorig=None, daily_mean=True, max_points=20000):
    """CROCO-minus-parent domain-wide differences at ONE lead-time label,
    POOLED (concatenated) across every cycle that reaches that lead.
    Calls gtools.validation._field_pair_for_day directly (rather than
    domain_diff(), which hardcodes daily_mean=False) once per contributing
    cycle and concatenates the flattened, finite differences, THEN applies
    the max_points subsample -- so the size limit bounds the final pooled
    sample, not each cycle individually.
    """
    parts = []
    for c in cycles_info:
        day = c["label_day"].get(lead)
        if day is None:
            continue
        try:
            ds = pp.open_history(c["croco_his"], Yorig=Yorig)
            croco, parent_on_croco = val._field_pair_for_day(
                ds, c["reference"], var, day, depth_m, tindex, daily_mean, margin_deg=None)
            ds.close()
        except Exception as e:
            print(f"  [{lead}] cycle {c['cycle']} ({day}): domain_diff failed ({e}) - skipping")
            continue
        d = (croco - parent_on_croco).ravel()
        d = d[np.isfinite(d)]
        parts.append(d)
    if not parts:
        return np.array([])
    diff = np.concatenate(parts)
    if diff.size > max_points:
        diff = np.random.default_rng(0).choice(diff, max_points, replace=False)
    return diff


def composite_domain_diff_satellite(cycles_info, product, lead, Yorig=None,
                                    margin_deg=None, max_points=20000):
    """CROCO-minus-satellite domain-wide differences at ONE lead-time
    label, pooled across every cycle that has a downloaded `product` file
    for that lead's day. Composite counterpart to
    gtools.validation.domain_diff_satellite(); does its own
    loading/regridding (mirroring compare_satellite_grid's inner loop)
    since domain_diff_satellite itself takes an already-computed pair.
    """
    parts = [] 
    for c in cycles_info:
        day = c["label_day"].get(lead)
        if day is None:
            continue
        fname = c.get("sat_files", {}).get(product, {}).get(day)
        if not fname or not os.path.exists(fname):
            continue
        loaded = val.load_satellite_field(fname, product, date=day)
        if loaded is None:
            continue
        plon, plat, pfield = loaded
        pfield = np.squeeze(np.asarray(pfield))
        if pfield.shape != plon.shape:
            continue
        try:
            clon, clat, croco_field = val._croco_field_satellite(
                c["croco_his"], product, day, Yorig=Yorig, margin_deg=margin_deg, daily_mean=True)
        except Exception as e:
            print(f"  [{product}][{lead}] cycle {c['cycle']} ({day}): "
                 f"CROCO record not available ({e}) - skipping")
            continue
        ds_grid = xr.Dataset(
            {"mask_rho": (("eta_rho", "xi_rho"), np.where(np.isfinite(croco_field), 1.0, np.nan))},
            coords={"lon_rho": (("eta_rho", "xi_rho"), clon),
                    "lat_rho": (("eta_rho", "xi_rho"), clat)})
        sat_on_croco = val.regrid_to_croco(plon, plat, pfield, ds_grid)
        parts.append(val.domain_diff_satellite(croco_field, sat_on_croco, max_points=max_points))
    if not parts:
        return np.array([])
    diff = np.concatenate(parts)
    if diff.size > max_points:
        diff = np.random.default_rng(0).choice(diff, max_points, replace=False)
    return diff


# ======================================================================
# Composite CROCO-vs-satellite field maps (composite counterpart to
# gtools.validation.compare_satellite_grid), ONE FIGURE PER LEAD.
# ======================================================================
def compare_composite_satellite_grid(cycles_info, product, lead, Yorig=None,
                                     margin_deg=None, out=None, figsize=(11, 10)):
    """CROCO-vs-satellite comparison at ONE lead-time label, composited
    (pixel-stacked mean/RMSE) across every cycle that has BOTH a CROCO
    record and a downloaded `product` file for that lead's day. One 2x2
    figure per lead (CROCO composite mean | satellite composite mean /
    bias | RMSE) -- mirrors gtools.validation.compare_satellite_grid's
    per-day figure, composited across cycles instead of drawn per single
    cycle/day. Reuses val._day_figure for the actual plotting.

    Returns (fig_or_path, stats_dict), or (None, None) if no cycle could
    contribute.
    """
    if product not in val.SAT_PRODUCTS:
        raise ValueError(f"product must be one of {list(val.SAT_PRODUCTS)}, got {product!r}")
    spec = val.SAT_PRODUCTS[product]
    label = spec["label"]

    pairs = []   # (cycle_name, clon, clat, croco_field, sat_on_croco)
    ref_clon = ref_clat = None
    for c in cycles_info:
        day = c["label_day"].get(lead)
        if day is None:
            continue
        fname = c.get("sat_files", {}).get(product, {}).get(day)
        if not fname or not os.path.exists(fname):
            continue
        loaded = val.load_satellite_field(fname, product, date=day)
        if loaded is None:
            continue
        plon, plat, pfield = loaded
        pfield = np.squeeze(np.asarray(pfield))
        if pfield.shape != plon.shape:
            continue
        try:
            clon, clat, croco_field = val._croco_field_satellite(
                c["croco_his"], product, day, Yorig=Yorig, margin_deg=margin_deg, daily_mean=True)
        except Exception as e:
            print(f"  [{product}][{lead}] cycle {c['cycle']} ({day}): "
                 f"CROCO record not available ({e}) - skipping")
            continue
        if ref_clon is None:
            ref_clon, ref_clat = clon, clat
        elif croco_field.shape != ref_clon.shape:
            print(f"  [{product}][{lead}] cycle {c['cycle']}: grid shape mismatch - skipping")
            continue
        ds_grid = xr.Dataset(
            {"mask_rho": (("eta_rho", "xi_rho"), np.where(np.isfinite(croco_field), 1.0, np.nan))},
            coords={"lon_rho": (("eta_rho", "xi_rho"), clon),
                    "lat_rho": (("eta_rho", "xi_rho"), clat)})
        sat_on_croco = val.regrid_to_croco(plon, plat, pfield, ds_grid)
        pairs.append((c["cycle"], croco_field, sat_on_croco))

    if not pairs:
        print(f"compare_composite_satellite_grid: no cycle has [{product}] data at lead={lead!r}")
        return None, None

    croco_stack = np.stack([p[1] for p in pairs], axis=0)
    sat_stack = np.stack([p[2] for p in pairs], axis=0)
    diff_stack = croco_stack - sat_stack

    croco_mean = np.nanmean(croco_stack, axis=0)
    sat_mean = np.nanmean(sat_stack, axis=0)
    bias_map = np.nanmean(diff_stack, axis=0)
    rmse_map = np.sqrt(np.nanmean(diff_stack ** 2, axis=0))
    n_cycles = len(pairs)
    used = [p[0] for p in pairs]

    s = val.domain_statistics(croco_mean, sat_mean)
    print(f"[{product}] {label} @ {lead}  (composite of {n_cycles} cycle(s): {used}):")
    val._print_stats(label, s)
    s.update({"variable": product, "lead": lead, "n_cycles": n_cycles, "cycles": used})

    vmin, vmax = percentile_clim(croco_mean, sat_mean)
    day_label = f"{lead}  ({n_cycles} cycle{'s' if n_cycles != 1 else ''})"
    fig = val._day_figure(ref_clon, ref_clat, croco_mean, sat_mean, bias_map, rmse_map,
                          day_label, vmin, vmax, product, n_days=None, figsize=figsize)
    if out:
        fig.savefig(out, dpi=150, bbox_inches="tight")
        plt.close(fig)
        return out, s
    return fig, s

def build_html_summary_composite(composite_dir, composite_id=None, config=None,
                                 cycles=None, status=None, extra_title=None):
    """Composite-aware wrapper around gtools.validation.build_html_summary():
    gathers every *.png/*.csv already written into COMPOSITE_DIR by this
    notebook (bias/RMSE maps, scatter, Taylor diagrams, timeseries,
    boxplots, satellite maps, ...) into one self-contained HTML page --
    SAME grouping/layout/lightbox code as the single-cycle report, just
    labelled for a composite run (composite id instead of a single
    calendar cycle) so a composite run and a single-cycle run can never
    silently disagree on how the report itself is built.

    composite_dir : COMPOSITE_DIR (see 03_composite_validation.ipynb
                     Section 1a) -- every *.png/*.csv directly under it is
                     picked up, exactly like build_html_summary().
    composite_id  : COMPOSITE_ID, shown in the header/filename in place of
                     a single cycle date.
    cycles        : optional list of the composited CYCLES -- appended to
                     extra_title so the report is self-describing about
                     which cycles it pools, without needing
                     validated_cycles.txt open alongside it.
    status        : optional dict (all_pass/criteria), same convention as
                     build_html_summary() -- pass Section 5's summary here
                     for a pass/fail banner; None for no banner.
    extra_title   : optional extra line, appended after the cycles note.

    Returns the path to the written HTML file.
    """
    label = f"composite_{composite_id}" if composite_id else "composite"
    bits = []
    if cycles:
        bits.append(f"{len(cycles)} cycle(s): {', '.join(cycles)}")
    if extra_title:
        bits.append(extra_title)
    combined_extra = " -- ".join(bits) if bits else None
    return val.build_html_summary(composite_dir, cycle=label, config=config,
                                  status=status, extra_title=combined_extra)

# ======================================================================
# NOTEBOOK-SIMPLIFICATION WRAPPERS (03_composite_validation.ipynb)
# ----------------------------------------------------------------------
# Mirrors validation.py's own notebook-simplification wrappers, one lead
# label at a time instead of one calendar day at a time. Every function
# here wraps existing building blocks above (compare_composite_field,
# godae_scorecard_composite, composite_domain_diff, ...) exactly as the
# notebook cells used to call them by hand -- no comparison logic,
# thresholds, or file naming has changed, only where the loop/guard/print
# now lives.
# ======================================================================


def resolve_composite_dir(main_dir, config, cycles, force_id=None, cycle_name="validated_cycles.txt"):
    """Section 0 -- find-or-create the composite output directory for this
    exact set of `cycles`: reuse an existing 'validation_composite_<id>'
    directory if one already lists exactly these cycles, otherwise mint a
    new short uuid. Writes the cycle list to `cycle_name` inside it either
    way. Prints the same messages the inline cell used to.

    Returns (COMPOSITE_ID, COMPOSITE_DIR).
    """
    import glob
    import uuid

    def _find_existing(main_dir, config, cycles):
        wanted = set(cycles)
        base = os.path.join(main_dir, config)
        for d in sorted(glob.glob(os.path.join(base, "validation_composite_*"))):
            f = os.path.join(d, cycle_name)
            if not os.path.isfile(f):
                continue
            with open(f) as fh:
                existing = {line.strip() for line in fh if line.strip()}
            if existing == wanted:
                return d
        return None

    if force_id:
        composite_id = force_id
        composite_dir = os.path.join(main_dir, config, f"validation_composite_{composite_id}")
        print(f"composite ID forced - using composite directory: {composite_dir}")
    else:
        existing_dir = _find_existing(main_dir, config, cycles)
        if existing_dir:
            composite_dir = existing_dir
            composite_id = os.path.basename(composite_dir).replace("validation_composite_", "")
            print(f"Reusing existing composite directory for this exact cycle set: {composite_dir}")
        else:
            composite_id = uuid.uuid4().hex[:8]
            composite_dir = os.path.join(main_dir, config, f"validation_composite_{composite_id}")
            print(f"No existing composite directory matches this cycle set - creating a new one: {composite_dir}")

    os.makedirs(composite_dir, exist_ok=True)
    with open(os.path.join(composite_dir, cycle_name), "w") as f:
        f.write("\n".join(cycles) + "\n")

    print(f"Compositing {len(cycles)} cycle(s): {cycles}")
    print(f"Composite outputs -> {composite_dir}")
    print(f"  (cycle list also recorded in {os.path.join(composite_dir, cycle_name)})")
    return composite_id, composite_dir


def download_and_label_cycles(cycles, config, main_dir, avail, yorig, spinup_days, get_paths_mod):
    """Section 1c -- for every cycle: resolve paths, run the numerical-
    stability check, download the Copernicus Marine Forecast reference (if
    available and not already covering the cycle window) and the satellite
    SST/SSS products, then label that cycle's days by forecast lead time.
    Prints exactly the same progress messages the inline cell used to.

    `get_paths_mod` is the notebook's already-imported `_paths` module.

    Returns (cycles_info, LEADS, FCST_LEADS, stability_ok).
    """
    from datetime import datetime

    cycles_info = []
    stability_ok = True

    for cycle in cycles:
        print(f"\n== cycle {cycle} ==")
        croco_his, reference, _ = get_paths_mod.get_paths(cycle=cycle, config=config, main_dir=main_dir)

        ds_c = pp.open_history(croco_his, Yorig=yorig)
        for v in ("temp", "salt", "zeta"):
            if v in ds_c and (not np.isfinite(ds_c[v].values).all()):
                print(f"  ! NaN/Inf found in {v} - numerical stability check FAILED for this cycle")
                stability_ok = False
        clon_c, clat_c, _ = pp.lonlatmask(ds_c)
        domain = (float(np.nanmin(clon_c)), float(np.nanmax(clon_c)),
                  float(np.nanmin(clat_c)), float(np.nanmax(clat_c)))
        model_times = pd.to_datetime(pp.times(ds_c))
        start_date, end_date = model_times[0].to_pydatetime(), model_times[-1].to_pydatetime()
        cycle_date = datetime.strptime(cycle, "%Y%m%d")
        ds_c.close()

        # ---- (i) Copernicus Marine Forecast (Mercator anfc), combined reference file ----
        if not avail['mercator_forecast']:
            print("  Copernicus Marine Forecast unavailable on the CMEMS platform - "
                  "this cycle's Sections 2/4/5b will contribute no data.")
        elif val.netcdf_covers_time_range(reference, start_date, end_date, max_step_hours=30):
            print(f"  Mercator: already downloaded and covers the full cycle window: {reference}")
        else:
            if os.path.exists(reference):
                os.remove(reference)
            mercator_dir = os.path.dirname(reference)
            fdays = max((end_date.date() - cycle_date.date()).days, 0)
            val.download_mercator_ops(domain, cycle_date, hdays=0, fdays=fdays, outputDir=mercator_dir)
            print(f"  Mercator: downloaded -> {reference}" if os.path.exists(reference)
                  else f"  Mercator: download ran but {reference} wasn't produced - check {mercator_dir}")

        # ---- (ii) Satellite SST: OSTIA & ODYSSEA, one file per day ----
        sat_files = {}
        for product in ("OSTIA", "ODYSSEA"):
            sat_dir = get_paths_mod.satellite_dir(main_dir, config, cycle, product)
            sat_files[product] = val.download_satellite_sst(product, domain, start_date, end_date, sat_dir)

        # ---- (ii-b) Satellite SSS: SMOS L4, one file per day ----
        if not avail['smos_l4_sss']:
            sat_files['SMOS'] = {}
        else:
            smos_dir = get_paths_mod.satellite_dir(main_dir, config, cycle, "SMOS")
            sat_files['SMOS'] = val.download_satellite_sst("SMOS", domain, start_date, end_date, smos_dir)
        print()

        # ---- lead-time labelling for this cycle's days ----
        days, day_label, label_day = label_cycle_days(croco_his, Yorig=yorig, spinup_days=spinup_days)
        print(f"  days: {days}")
        print(f"  lead labels: {day_label}")

        cycles_info.append({
            "cycle": cycle, "croco_his": croco_his, "reference": reference,
            "days": days, "day_label": day_label, "label_day": label_day,
            "sat_files": sat_files,
        })

    LEADS = all_lead_labels(cycles_info)
    FCST_LEADS = all_lead_labels(cycles_info, fcst_only=True)
    print(f"\nAll lead labels in this composite : {LEADS}")
    print(f"Forecast-only lead labels (Section 5): {FCST_LEADS}")
    print(f"Numerical stability across all cycles : {'OK' if stability_ok else 'FAILED - see above'}")
    for lead in LEADS:
        print(f"  {lead:6s}: reached by {cycles_reaching(cycles_info, lead)}")

    return cycles_info, LEADS, FCST_LEADS, stability_ok


def _validate_maps_by_lead(var, cycles_info, leads, yorig, composite_dir, available,
                            prefix, out_var, **extra_kwargs):
    """Shared body of validate_sst_maps/validate_ssh_maps/validate_current_maps/
    validate_sss_maps: guard, loop compare_composite_field() over `leads`,
    print, return {lead: stats}."""
    if not available:
        print("Copernicus Marine Forecast unavailable - skipping this comparison.")
        return {}
    stats_by_lead = {}
    for lead in leads:
        res = compare_composite_field(
            cycles_info, var=var, lead=lead, Yorig=yorig, daily_mean=True,
            out=os.path.join(composite_dir, f'{prefix}_vs_forecast_{lead}.png'),
            **extra_kwargs)
        if res is not None:
            _, stats_by_lead[lead] = res
    print(f"-> {len(stats_by_lead)} figure(s) written, one per lead: {list(stats_by_lead)}")
    return stats_by_lead


def validate_sst_maps(cycles_info, leads, yorig, depth_m, composite_dir, available):
    """Section 2 -- composite SST bias maps, by lead time. Returns {lead: stats}."""
    return _validate_maps_by_lead('temp', cycles_info, leads, yorig, composite_dir,
                                   available, 'sst', 'sst', depth_m=depth_m)


def validate_ssh_maps(cycles_info, leads, yorig, composite_dir, available):
    """Section 2 -- composite SSH bias maps, by lead time. SSH has no depth
    dimension -- always surface, regardless of DEPTH_M. Returns {lead: stats}."""
    return _validate_maps_by_lead('ssh', cycles_info, leads, yorig, composite_dir,
                                   available, 'ssh', 'ssh')


def validate_current_maps(cycles_info, leads, yorig, depth_m, composite_dir, available):
    """Section 2 -- composite surface-current bias maps, by lead time.
    Returns {lead: stats}. Surface velocities pass criterion is QUALITATIVE
    (visual consistency with expected gyre/coastal-jet circulation)."""
    stats_by_lead = _validate_maps_by_lead('speed', cycles_info, leads, yorig, composite_dir,
                                            available, 'currents', 'currents', depth_m=depth_m)
    if stats_by_lead:
        print()
        print("Surface velocities pass criterion is QUALITATIVE (visual consistency with")
        print("expected gyre/coastal-jet circulation) -- inspect the vector maps")
    return stats_by_lead


def validate_sss_maps(cycles_info, leads, yorig, depth_m, composite_dir, available):
    """Section 2 -- composite SSS bias maps, by lead time. Returns {lead: stats}."""
    return _validate_maps_by_lead('salt', cycles_info, leads, yorig, composite_dir,
                                   available, 'sss', 'sss', depth_m=depth_m)


def scatter_maps_by_lead(cycles_info, leads, yorig, composite_dir, available):
    """Section 3 -- composite SST+SSH pointwise scatter, one figure per lead
    (points pooled across every cycle reaching that lead). Returns {lead: out_png}."""
    if not available:
        print("Copernicus Marine Forecast unavailable - skipping this comparison.")
        return {}
    scatter_figs = {}
    for lead in leads:
        out_png = os.path.join(composite_dir, f"scatter_sst_ssh_vs_forecast_{lead}.png")
        res = scatter_composite(cycles_info, lead, Yorig=yorig, daily_mean=True, out=out_png)
        if res is not None:
            scatter_figs[lead] = out_png
    print(f"-> {len(scatter_figs)} figure(s) written, one per lead: {list(scatter_figs)}")
    return scatter_figs


def build_godae_scorecard(cycles_info, leads, yorig, depth_m, available,
                           variables=('temp', 'ssh', 'salt', 'speed')):
    """Section 4 -- composite GODAE scorecard, by lead time (points pooled
    across every cycle reaching that lead). Returns (report, report_by_lead)."""
    if not available:
        print("Copernicus Marine Forecast unavailable - skipping GODAE scorecard "
              "(Sections 4/5 below will be skipped too).")
        empty = pd.DataFrame(columns=['variable', 'lead', 'bias', 'rmsd', 'urmsd', 'corr', 'std_ratio', 'vs', 'layer'])
        return empty, {}

    rows = []
    report_by_lead = {}
    for lead in leads:
        lead_rows = []
        for var in variables:
            s = godae_scorecard_composite(cycles_info, var, lead, depth_m=depth_m, Yorig=yorig)
            if s is not None:
                lead_rows.append({**s, 'vs': 'reference', 'layer': 'all'})
        if lead_rows:
            report_by_lead[lead] = pd.DataFrame(lead_rows)
            rows.extend(lead_rows)
            print(f"-- {lead}  (cycles: {cycles_reaching(cycles_info, lead)}) --")
            val.print_scorecard_table(report_by_lead[lead])

    report = pd.DataFrame(rows)
    return report, report_by_lead


def plot_taylor_diagrams_by_lead(report_by_lead, composite_dir):
    """Section 4 -- one Taylor diagram per lead, all of that lead's
    variables (temp/ssh/salt/speed) plotted together, points pooled across
    cycles. Returns {lead: out_png}."""
    if not report_by_lead:
        print("No GODAE scorecard available (Copernicus Marine Forecast unavailable) - skipping Taylor diagram.")
        return {}

    taylor_figs = {}
    for lead, rep in report_by_lead.items():
        fig = plt.figure(figsize=(7, 7))
        ax = fig.add_subplot(111, polar=True)
        has_neg = bool((rep['corr'] < 0).any())
        thetamax = 180 if has_neg else 90
        corr_ticks = ([-1.0, -0.5, 0, 0.5, 0.8, 0.9, 0.95, 0.99, 1.0] if has_neg
                      else [0, 0.2, 0.4, 0.6, 0.8, 0.9, 0.95, 0.99, 1.0])
        ax.set_thetamin(0); ax.set_thetamax(thetamax)
        ax.set_xticks(np.arccos(corr_ticks)); ax.set_xticklabels([str(c) for c in corr_ticks])
        ax.set_rlabel_position(0)   # keep the std-dev tick labels along the bottom (theta=0) axis
        finite = rep['std_ratio'][np.isfinite(rep['std_ratio'])]
        r_max = max(1.6, finite.max() * 1.2) if len(finite) else 1.6
        ax.set_ylim(0, r_max)
        ax.plot(0, 1, 'k*', ms=16, label='reference')
        for _, row in rep.iterrows():
            theta = np.arccos(np.clip(row['corr'], -1, 1))
            ax.plot(theta, row['std_ratio'], 'o', ms=10,
                    label=f"{row['variable']}  (n_cycles={row['n_cycles']}, RMSD={row['rmsd']:.2f})")

        ax.text(0.5, -0.08, 'Normalised standard deviation (CROCO / reference)',
                transform=ax.transAxes, ha='center', va='top', fontsize=10)
        ax.text(np.radians(thetamax / 6), r_max * 1.25, 'Correlation coefficient',
                ha='center', va='center', fontsize=10,
                rotation=90 - thetamax / 2, rotation_mode='anchor')

        ax.set_title(f'Taylor diagram -- CROCO vs reference ({lead}, composite)', pad=20)
        ax.legend(loc='upper left', bbox_to_anchor=(1.05, 1.0), fontsize=9)
        fig.tight_layout()
        out_png = os.path.join(composite_dir, f"taylor_diagram_{lead}.png")
        fig.savefig(out_png, dpi=150, bbox_inches="tight")
        plt.close(fig)
        taylor_figs[lead] = out_png
    print(f"-> {len(taylor_figs)} figure(s) written, one per lead: {list(taylor_figs)}")
    return taylor_figs


def pass_fail_summary(report, cycles, fcst_leads, stability_ok, spinup_days=2):
    """Section 5 -- automated pass/fail summary, composite, FORECAST LEADS
    ONLY (spin-up excluded via `report['lead'].isin(fcst_leads)` rather than
    "first N calendar days", since spin-up is already a separate set of
    lead labels here). Reads directly from build_godae_scorecard()'s
    `report`, so it can't drift out of sync with the Taylor diagram.

    Returns all_pass, or None if report is empty (nothing to score).
    """
    if report.empty:
        print("Overall V1 status: SKIPPED -- Copernicus Marine Forecast unavailable, "
              "no GODAE scorecard to build the pass/fail summary from.")
        return None

    report_eval = report[report['lead'].isin(fcst_leads)]

    by_var = (report_eval.groupby('variable')[['bias', 'rmsd', 'urmsd', 'corr', 'si', 'si_std', 'std_ratio']]
              .mean().to_dict('index'))

    criteria = [
        ('SST composite-mean domain-avg RMSD < 0.5 degC',       by_var['temp']['rmsd'] < 0.5,  f"{by_var['temp']['rmsd']:.3f} degC"),
        ('SSH composite-mean spatial correlation > 0.90',       by_var['ssh']['corr'] > 0.90, f"{by_var['ssh']['corr']:.3f}"),
        ('Salinity composite-mean domain-avg |bias| < 0.2 PSU', abs(by_var['salt']['bias']) < 0.2, f"{by_var['salt']['bias']:+.3f} PSU"),
        ('Numerical stability (no NaN/Inf, every cycle)',       stability_ok, 'see section 1c'),
    ]

    print(f"Composite scorecard: {len(cycles)} cycle(s) {cycles}, forecast leads only "
          f"(spin-up {[f'sp{i+1}' for i in range(spinup_days)]} excluded): {fcst_leads}")
    print()
    print(f"{'Criterion':50s} {'Result':6s}  Value")
    print('-' * 75)
    all_pass = True
    for name, passed, value in criteria:
        all_pass &= passed
        print(f"{name:50s} {'PASS' if passed else 'FAIL':6s}  {value}")
    print('-' * 75)
    print()

    if len(fcst_leads) > 1:
        temp_rows = report_eval[report_eval['variable'] == 'temp']
        ssh_rows  = report_eval[report_eval['variable'] == 'ssh']
        salt_rows = report_eval[report_eval['variable'] == 'salt']
        worst_temp = temp_rows.loc[temp_rows['rmsd'].idxmax()]
        worst_ssh  = ssh_rows.loc[ssh_rows['corr'].idxmin()]
        worst_salt = salt_rows.loc[salt_rows['bias'].abs().idxmax()]
        print()
        print(f"Worst single lead -- SST RMSD:   {worst_temp['lead']}  ({worst_temp['rmsd']:.3f} degC)")
        print(f"Worst single lead -- SSH corr:   {worst_ssh['lead']}  ({worst_ssh['corr']:.3f})")
        print(f"Worst single lead -- Salt bias:  {worst_salt['lead']}  ({worst_salt['bias']:+.3f} PSU)")

    print()
    print("Note: surface-current skill (by_var['speed']) has no fixed Section 9.3")
    print("threshold -- its criterion is qualitative (inspect the vector maps in")
    print("Section 2), same as 02_validation.ipynb.")
    return all_pass


def stacked_bias_boxplot(cycles_info, leads, yorig, depth_m, composite_dir, cycles, available):
    """Section 5b -- composite domain-wide CROCO-minus-parent bias boxplot,
    one stacked figure with 4 subpanels (SSH, temp, salt, speed), one box
    per lead. Also writes the RMSE-spread companion figure. Matches the
    original cell's behaviour exactly, including neither closing the
    figures nor printing a summary line afterwards.
    """
    if not available:
        print("Copernicus Marine Forecast unavailable - skipping this comparison.")
        return

    dlab = 'surface' if depth_m is None else f'{depth_m:g} m'
    n_cycles = len(cycles)

    ssh_diffs   = [composite_domain_diff(cycles_info, 'ssh',   lead, Yorig=yorig, daily_mean=True) for lead in leads]
    temp_diffs  = [composite_domain_diff(cycles_info, 'temp',  lead, depth_m=depth_m, Yorig=yorig, daily_mean=True) for lead in leads]
    salt_diffs  = [composite_domain_diff(cycles_info, 'salt',  lead, depth_m=depth_m, Yorig=yorig, daily_mean=True) for lead in leads]
    speed_diffs = [composite_domain_diff(cycles_info, 'speed', lead, depth_m=depth_m, Yorig=yorig, daily_mean=True) for lead in leads]

    def _stacked_boxplot(diffs_list, ylabels, titles, out, metric='bias'):
        data = [[np.abs(d) for d in diffs] for diffs in diffs_list] if metric == 'rmse' else diffs_list
        fig, axes = plt.subplots(len(data), 1, figsize=(max(6, 1.1 * len(leads)), 3.3 * len(data)),
                                  sharex=True)
        positions = np.arange(1, len(leads) + 1)
        for ax, diffs, ylabel, title in zip(axes, data, ylabels, titles):
            ax.boxplot(diffs, positions=positions, showfliers=False)
            ax.set_xticks(positions); ax.set_xticklabels(leads)
            if metric == 'bias':
                ax.axhline(0, color='k', ls='--', lw=1)
            ax.set_ylabel(ylabel); ax.set_title(title); ax.grid(alpha=0.3)
        axes[-1].set_xlabel('lead time')
        plt.setp(axes[-1].get_xticklabels(), rotation=45, ha='right')
        fig.suptitle(f"CROCO - parent  (composite of {n_cycles} cycles)" if metric == 'bias'
                     else f"|CROCO - parent|  (composite of {n_cycles} cycles)")
        fig.tight_layout()
        fig.savefig(out, dpi=150, bbox_inches='tight')
        return fig

    diffs_all = [ssh_diffs, temp_diffs, salt_diffs, speed_diffs]
    ylabels = ["SSH' bias (m)", 'temperature bias (degC)', 'salinity bias (PSU)', 'speed bias (m s$^{-1}$)']
    titles = ['SSH anomaly bias', f'temperature bias  ({dlab})', f'salinity bias  ({dlab})',
              f'current speed bias  ({dlab})']

    _stacked_boxplot(diffs_all, ylabels, titles,
                      os.path.join(composite_dir, 'boxplot_bias_vs_forecast_composite.png'), metric='bias')

    rmse_ylabels = [y.replace('bias', '|error|') for y in ylabels]
    rmse_titles = [t.replace('bias', 'RMSE-spread') for t in titles]
    _stacked_boxplot(diffs_all, rmse_ylabels, rmse_titles,
                      os.path.join(composite_dir, 'boxplot_rmse_vs_forecast_composite.png'), metric='rmse')


def satellite_bias_boxplots(cycles_info, avail, leads, yorig, composite_dir, cycles):
    """Section 5c -- composite CROCO-minus-satellite domain-wide bias
    boxplot (per lead), against OSTIA/ODYSSEA (grouped) and SMOS (separate),
    plus the RMSE-spread companions. Reuses val.bias_boxplot_multi()
    directly (same as the original cell did), since
    composite_domain_diff_satellite() already returns plain per-lead
    difference arrays. Returns True if any satellite product was available,
    else False (and prints why, same as before).
    """
    sat_avail_key = {"OSTIA": "ostia_l4", "ODYSSEA": "odyssea_l3s", "SMOS": "smos_l4_sss"}
    any_sat_avail = any(avail.get(sat_avail_key[p], False) and any(c["sat_files"].get(p) for c in cycles_info)
                        for p in sat_avail_key)
    n_cycles = len(cycles)

    if not any_sat_avail:
        print("No satellite product (OSTIA/ODYSSEA/SMOS) available/downloaded for this composite - skipping.")
        return any_sat_avail

    ostia_diffs   = [composite_domain_diff_satellite(cycles_info, 'OSTIA',   lead, Yorig=yorig) for lead in leads]
    odyssea_diffs = [composite_domain_diff_satellite(cycles_info, 'ODYSSEA', lead, Yorig=yorig) for lead in leads]
    smos_diffs    = [composite_domain_diff_satellite(cycles_info, 'SMOS',    lead, Yorig=yorig) for lead in leads]

    sst_groups = {"OSTIA": ostia_diffs, "ODYSSEA": odyssea_diffs}
    sss_groups = {"SMOS": smos_diffs}

    val.bias_boxplot_multi(sst_groups, leads, 'temperature bias (degC)',
                            f'CROCO - satellite SST bias  (composite of {n_cycles} cycles)',
                            colors=['C1', 'C2'],
                            out=os.path.join(composite_dir, 'boxplot_bias_sst_vs_satellite_composite.png'))
    val.bias_boxplot_multi(sss_groups, leads, 'salinity bias (PSU)',
                            f'CROCO - satellite SSS bias  (composite of {n_cycles} cycles)',
                            colors=['C4'],
                            out=os.path.join(composite_dir, 'boxplot_bias_sss_vs_satellite_composite.png'))
    val.bias_boxplot_multi(sst_groups, leads, 'temperature |error| (degC)',
                            f'CROCO - satellite SST RMSE-spread  (composite of {n_cycles} cycles)',
                            colors=['C1', 'C2'], metric='rmse',
                            out=os.path.join(composite_dir, 'boxplot_rmse_sst_vs_satellite_composite.png'))
    val.bias_boxplot_multi(sss_groups, leads, 'salinity |error| (PSU)',
                            f'CROCO - satellite SSS RMSE-spread  (composite of {n_cycles} cycles)',
                            colors=['C4'], metric='rmse',
                            out=os.path.join(composite_dir, 'boxplot_rmse_sss_vs_satellite_composite.png'))
    return any_sat_avail


def plot_composite_timeseries(cycles_info, leads, point_lon, point_lat, depth_m, yorig,
                               composite_dir, cycles, available, any_sat_avail, dlab=None):
    """Sections 5d -- composite mean +/- 1 std timeseries by lead, both
    full-domain and at one point: CROCO vs parent (needs `available`), and
    CROCO vs satellite (OSTIA/ODYSSEA/SMOS, needs `any_sat_avail`, from
    satellite_bias_boxplots()'s return value).

    `dlab` ('surface' or 'N m') is computed from depth_m if not given.
    """
    if dlab is None:
        dlab = 'surface' if depth_m is None else f'{depth_m:g} m'
    n_cycles = len(cycles)

    if not available:
        print("Copernicus Marine Forecast unavailable - skipping this comparison.")
    else:
        vars_specs = [('ssh', "SSH' (m)", 'SSH anomaly', None),
                      ('temp', 'temperature (degC)', f'temperature ({dlab})', depth_m),
                      ('salt', 'salinity (PSU)', f'salinity ({dlab})', depth_m),
                      ('speed', 'speed (m s$^{-1}$)', f'current speed ({dlab})', depth_m)]

        domain_rows = []
        point_rows = []
        for var, ylabel, title, dm in vars_specs:
            cd, p_dom = composite_domain_series(cycles_info, var, leads, depth_m=dm, Yorig=yorig, daily_mean=True)
            domain_rows.append((ylabel, title + ' (full domain)', cd, p_dom))
            cp, p_pt = composite_point_series(cycles_info, var, point_lon, point_lat, leads, depth_m=dm, Yorig=yorig, daily_mean=True)
            point_rows.append((ylabel, title + f' ({point_lon:.2f}, {point_lat:.2f})', cp, p_pt))

        ## full-domain averages over every horizontal level, cumulating localized errors.
        ## The figure shape can then show huge std spread around the mean.
        composite_two_series_timeseries(domain_rows, leads, n_cycles, labels=("CROCO", "parent"),
                                         out=os.path.join(composite_dir, "timeseries_mean_domain_vs_forecast_composite.png"))
        ## single point
        composite_two_series_timeseries(point_rows, leads, n_cycles, labels=("CROCO", "parent"),
                                         out=os.path.join(composite_dir, "timeseries_mean_point_vs_forecast_composite.png"))

    if not any_sat_avail:
        print("No satellite product available/downloaded for this composite - skipping.")
    else:
        sat_domain_rows, sat_point_rows = [], []
        for product, ylabel in (('OSTIA', 'SST bias (degC)'), ('ODYSSEA', 'SST bias (degC)'), ('SMOS', 'SSS (PSU)')):
            cd, sd = composite_satellite_series(cycles_info, product, leads, Yorig=yorig)
            sat_domain_rows.append((ylabel, f'CROCO vs {product} (full domain)', cd, sd))
            cp, sp = composite_satellite_series(cycles_info, product, leads, lon0=point_lon, lat0=point_lat, Yorig=yorig)
            sat_point_rows.append((ylabel, f'CROCO vs {product} ({point_lon:.2f}, {point_lat:.2f})', cp, sp))

        ## full-domain averages over every horizontal level, cumulating localized errors.
        ## The figure shape can then show huge std spread around the mean.
        composite_two_series_timeseries(sat_domain_rows, leads, n_cycles, labels=("CROCO", "satellite"),
                                         out=os.path.join(composite_dir, "timeseries_mean_domain_vs_satellite_composite.png"))
        ## single point
        composite_two_series_timeseries(sat_point_rows, leads, n_cycles, labels=("CROCO", "satellite"),
                                         out=os.path.join(composite_dir, "timeseries_mean_point_vs_satellite_composite.png"))


def composite_profiles(cycles_info, leads, point_lon, point_lat, yorig, composite_dir,
                        variables=('temp', 'salt', 'speed')):
    """Section 5e -- composite point and full-domain vertical profiles, mean
    +/- 1 std across every cycle reaching each lead, for each variable in
    `variables`. One point profile + one full-domain profile figure per
    (lead, variable).
    """
    for lead in leads:
        for var in variables:
            composite_point_profile(cycles_info, var, point_lon, point_lat, lead,
                                     Yorig=yorig, daily_mean=True,
                                     out=os.path.join(composite_dir, f'profile_point_{var}_{lead}_composite.png'))

            # A full-domain mean vertical profile averages over every horizontal point
            # at each depth, cumulating localized errors at every level -- the profile
            # shape then shows a clearer difference between CROCO and reference than
            # the single-point profile above.
            composite_domain_profile(cycles_info, var, lead, Yorig=yorig, daily_mean=True,
                                      out=os.path.join(composite_dir, f'profile_domain_{var}_{lead}_composite.png'))


def validate_satellite(cycles_info, product, avail_flag, leads, yorig, composite_dir):
    """Sections 6/6b -- composite CROCO vs one satellite product (OSTIA,
    ODYSSEA, or SMOS), one figure per lead, plus a CSV of per-lead
    statistics. Guards on `avail_flag` and on whether any cycle in this
    composite downloaded this product at all.

    Returns the stats DataFrame (possibly empty), or None if skipped.
    """
    if not avail_flag:
        print(f"{product} unavailable on the CMEMS platform - skipping.")
        return None
    if not any(c["sat_files"].get(product) for c in cycles_info):
        print(f"{product}: nothing downloaded for any cycle in this composite - skipping.")
        return None

    print(f"\n-- {product} --")
    csv_var = 'sss' if product == 'SMOS' else 'sst'
    figs, rows = {}, []
    for lead in leads:
        out_png = os.path.join(composite_dir, f"{csv_var}_vs_{product.lower()}_{lead}.png")
        res, s = compare_composite_satellite_grid(cycles_info, product, lead, Yorig=yorig, out=out_png)
        if res is not None:
            figs[lead] = res
            rows.append(s)
    stats = pd.DataFrame(rows)
    print(f"  -> {len(figs)} figure(s) written, one per lead: {list(figs)}")
    if not stats.empty:
        stats.to_csv(os.path.join(composite_dir, f"{csv_var}_vs_{product.lower()}_stats_composite.csv"), index=False)
    return stats
