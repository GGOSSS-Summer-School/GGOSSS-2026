"""
glofas_rivers.py -- river discharge for a CROCO hindcast from GloFAS.

Rivers are usually forced with a discharge climatology (e.g. Dai & Trenberth).
This module uses GloFAS instead (Copernicus Emergency Management Service,
`cems-glofas-historical`, LISFLOOD v4 consolidated, daily discharge, 0.05 deg),
so a hindcast gets the real day-to-day discharge of its own period. The rest of
the river chain is the standard croco_pytools one: per-river text files +
river_list.txt -> prepro/make_river_run.py -> croco_runoff.nc + for_croco_in.txt.

Two steps (exposed as ggosss26.py subcommands):

  download_glofas_monthly()  one netcdf per month, glofas_YYYY_MM.nc, over
                             the grid box + pad, from the EWDS data store.
  make_glofas_river_files()  picks the river mouths in the domain and writes
                             <river>.txt ("YYYY/MM/DD hh:mm:ss  Q[m3/s]") +
                             river_list.txt for make_river_run.py.

River mouths: GloFAS is a gridded field, so the MOUTH POSITIONS and NAMES are
taken from the Dai & Trenberth file (lon_mou, lat_mou, riv_name -- positions
only, its discharge is not used), optionally extended by a user list
("name lon lat" per line). For each mouth the GloFAS cell with the largest
mean discharge within `radius` degrees is taken as the river's outlet cell
(the main channel near the coast), and rivers whose mean GloFAS discharge is
below `qmin` are dropped. The position written to river_list.txt is the
mouth itself; make_river_run.py moves it to the nearest wet CROCO cell.

Time stamps: GloFAS dis24 at valid_time t is the mean over [t-24h, t], so it
is written at t-12h (the centre of the averaging window).

Request schema: the post-29-Jul-2026 EWDS form (year/month/day, timespan=
'time_mean', variable=average_river_discharge_in_the_last_24_hours); GLOFAS_VERSION
(default version_4_0, the operational one; version_5_0 = pre-operational v5).

Credentials: EWDS (https://ewds.climate.copernicus.eu) is a different data
store from CDS, same ECMWF account. ~/.ewdsapirc:
    url: https://ewds.climate.copernicus.eu/api
    key: <your personal access token>
and accept the GloFAS historical licence once on the dataset page.
"""
import calendar
import concurrent.futures
import datetime
import os
import threading
import time as _time

import numpy as np

EWDS_RC = os.path.expanduser(os.environ.get('EWDSAPI_RC', '~/.ewdsapirc'))
EWDS_URL = 'https://ewds.climate.copernicus.eu/api'
DATASET = 'cems-glofas-historical'
_lock = threading.Lock()


def _log(msg):
    with _lock:
        print(datetime.datetime.now().strftime('%H:%M:%S') + '  ' + msg, flush=True)


def _ewds_client():
    if not os.path.isfile(EWDS_RC):
        raise SystemExit(
            f"{EWDS_RC} not found. GloFAS is served by EWDS, not CDS. Create it with\n"
            "    url: https://ewds.climate.copernicus.eu/api\n"
            "    key: <your ECMWF personal access token>\n"
            "then `chmod 600 ~/.ewdsapirc` and accept the licence on\n"
            "https://ewds.climate.copernicus.eu/datasets/cems-glofas-historical?tab=download")
    cfg = {}
    with open(EWDS_RC) as f:
        for line in f:
            if ':' in line:
                k, v = line.split(':', 1)
                cfg[k.strip()] = v.strip()
    import cdsapi
    # GloFAS exists only on EWDS: a CDS url in the file (a common copy-paste slip)
    # gives "404 Not Found ... cds.climate.copernicus.eu/.../cems-glofas-historical"
    url = cfg.get('url', '')
    if 'ewds.climate.copernicus.eu' not in url:
        global _URL_WARNED
        if not globals().get('_URL_WARNED'):
            print('  NOTE: %s has url "%s"; GloFAS is on EWDS -- using %s '
                  '(fix the url line in the file)' % (EWDS_RC, url, EWDS_URL))
            _URL_WARNED = True
        url = EWDS_URL
    return cdsapi.Client(url=url, key=cfg['key'], quiet=True, progress=False)


def grid_box(grd_file):
    import netCDF4
    with netCDF4.Dataset(grd_file) as g:
        lon = g['lon_rho'][:]; lat = g['lat_rho'][:]
        return float(lon.min()), float(lon.max()), float(lat.min()), float(lat.max())


def _months(month_start, month_end):
    y, m = map(int, month_start.split('-'))
    ye, me = map(int, month_end.split('-'))
    out = []
    while (y, m) <= (ye, me):
        out.append((y, m))
        m += 1
        if m > 12:
            m, y = 1, y + 1
    return out


_CONS = None
def _constraints():
    """The public EWDS constraints of cems-glofas-historical (cached; None if unreachable)."""
    global _CONS
    if _CONS is None:
        import json, urllib.request
        base = 'https://ewds.climate.copernicus.eu/api/catalogue/v1/collections/' + DATASET
        try:
            coll = json.load(urllib.request.urlopen(base, timeout=60))
            url = [l['href'] for l in coll['links'] if l.get('rel') == 'constraints'][0]
            _CONS = json.load(urllib.request.urlopen(url, timeout=60))
        except Exception as e:
            print('  could not read the EWDS constraints (%s)' % str(e)[:80])
            _CONS = []
    return _CONS or None


def available_days(y, m, version, ptype):
    """Days of month y-m that EWDS offers for this version/product type (sorted),
    or None if the constraints could not be read."""
    cons = _constraints()
    if cons is None:
        return None
    days = set()
    for c in cons:
        if (version in c.get('system_version', []) and ptype in c.get('product_type', [])
                and 'time_mean' in c.get('timespan', []) and str(y) in c.get('year', [])
                and '%02d' % m in c.get('month', [])
                and 'average_river_discharge_in_the_last_24_hours' in c.get('variable', [])):
            days |= set(c.get('day', []))
    return sorted(days)


def choose_product_type(months, version):
    """'consolidated' if EWDS offers it for every month (the last month may be
    partial), else 'intermediate' (recent months exist only as intermediate)."""
    if _constraints() is None:
        return 'consolidated'
    def complete(pt):
        for y, m in months:
            days = available_days(y, m, version, pt)
            if not days:
                return False
            if len(days) < calendar.monthrange(y, m)[1] and (y, m) != months[-1]:
                return False
        return True
    for pt in ('consolidated', 'intermediate'):
        if complete(pt):
            return pt
    print('  WARNING: GloFAS %s does not cover all months in either product type' % version)
    return 'intermediate'


def download_glofas_monthly(grd_file, month_start, month_end, out_dir, pad=0.5, n_parallel=2):
    """GloFAS daily discharge, one file per month: <out_dir>/glofas_YYYY_MM.nc."""
    os.makedirs(out_dir, exist_ok=True)
    version = os.environ.get('GLOFAS_VERSION', 'version_4_0')
    ptype = os.environ.get('GLOFAS_PRODUCT_TYPE', 'auto')
    if ptype == 'auto':
        ptype = choose_product_type(_months(month_start, month_end), version)
    print(' GloFAS %s, product type %s' % (version, ptype))
    lo0, lo1, la0, la1 = grid_box(grd_file)
    # [N, W, S, E], snapped outwards to the 0.05 deg GloFAS grid
    area = [np.ceil((la1 + pad) * 20) / 20, np.floor((lo0 - pad) * 20) / 20,
            np.floor((la0 - pad) * 20) / 20, np.ceil((lo1 + pad) * 20) / 20]
    area = [round(float(a), 2) for a in area]

    jobs = []
    for y, m in _months(month_start, month_end):
        out = os.path.join(out_dir, 'glofas_%04d_%02d.nc' % (y, m))
        if os.path.isfile(out) and os.path.getsize(out) > 0:
            print(' already downloaded, skipping: ' + os.path.basename(out))
            continue
        ndays = calendar.monthrange(y, m)[1]
        # only the days EWDS offers (the current month is partial; asking for
        # missing days makes the whole request fail)
        days = available_days(y, m, version, ptype)
        if days is None:
            days = ['%02d' % d for d in range(1, ndays + 1)]
        if not days:
            raise SystemExit('GloFAS %s %s has no data for %04d-%02d' % (version, ptype, y, m))
        if len(days) < ndays:
            print(' %04d-%02d: only days %s..%s available (partial month)' % (y, m, days[0], days[-1]))
        # EWDS request schema since 29 Jul 2026 (GloFAS v5 release, CEMS page
        # "Changes to the EWDS API Request to download GloFAS Historical"):
        # hyear/hmonth/hday -> year/month/day, mandatory `timespan`, and the
        # variable renamed to average_river_discharge_in_the_last_24_hours.
        # The old keys (still in croco_pytools v2.0.4 download_glofas_river.py)
        # may no longer work. v4 stays the operational version.
        req = {
            'system_version': [version],
            'hydrological_model': ['lisflood'],
            'product_type': [ptype],
            'timespan': ['time_mean'],
            'variable': ['average_river_discharge_in_the_last_24_hours'],
            'year': [str(y)],
            'month': ['%02d' % m],
            'day': days,
            'data_format': 'netcdf',
            'download_format': 'unarchived',
            'area': area,
        }
        jobs.append((req, out))

    def _get(job):
        req, out = job
        name = os.path.basename(out)
        for attempt in range(1, 9):
            try:
                _log('request  ' + name + ('' if attempt == 1 else '  (attempt %d)' % attempt))
                _ewds_client().retrieve(DATASET, req).download(out + '.part')
                os.replace(out + '.part', out)
                _log('done     ' + name)
                return
            except Exception as e:
                _log('FAILED   %s: %s' % (name, str(e).splitlines()[0][:200]))
                if attempt == 8:
                    raise
                _time.sleep(min(60 * attempt, 300))

    print(' %d GloFAS month(s) to do, %d in parallel, area [N,W,S,E] = %s' % (len(jobs), n_parallel, area))
    failed = []
    with concurrent.futures.ThreadPoolExecutor(max_workers=max(1, n_parallel)) as pool:
        futs = {pool.submit(_get, j): j[1] for j in jobs}
        for f in concurrent.futures.as_completed(futs):
            if f.exception() is not None:
                failed.append(os.path.basename(futs[f]))
    if failed:
        raise SystemExit(' GloFAS downloads failed: %s -- re-run (finished months are skipped)'
                         % ', '.join(sorted(failed)))


def open_glofas(glofas_dir, month_start, month_end):
    """(time centre of each daily mean, lat, lon, Q[t,y,x]) over the months."""
    import xarray as xr
    files = [os.path.join(glofas_dir, 'glofas_%04d_%02d.nc' % ym) for ym in _months(month_start, month_end)]
    missing = [f for f in files if not os.path.isfile(f)]
    if missing:
        raise SystemExit('missing GloFAS files: ' + ', '.join(os.path.basename(f) for f in missing))
    parts = []
    for f in files:
        ds = xr.open_dataset(f)
        # discharge variable: 'dis24' (v4 netcdf) or a renamed one after the
        # 2026 MARS restructuring ('avg_dis24', ...): take the one with 'dis'
        dvars = [v for v in ds.data_vars if ds[v].ndim >= 3]
        cand = [v for v in dvars if 'dis' in v.lower()] or dvars
        if len(cand) != 1:
            raise SystemExit('cannot identify the discharge variable in %s: %s' % (f, list(ds.data_vars)))
        da = ds[cand[0]]
        # time axis: the non-spatial dimension; its time is the END of the 24 h
        # window (valid_time) -- use the valid_time coordinate when present
        tdims = [d for d in da.dims if d not in ('latitude', 'longitude')]
        da = da.squeeze([d for d in tdims if da.sizes[d] == 1 and len(tdims) > 1], drop=True)
        tdims = [d for d in da.dims if d not in ('latitude', 'longitude')]
        if len(tdims) != 1 or 'latitude' not in da.dims or 'longitude' not in da.dims:
            raise SystemExit('unexpected GloFAS layout in %s: %s' % (f, da.dims))
        td = tdims[0]
        if 'valid_time' in da.coords and da['valid_time'].dims == (td,):
            da = da.assign_coords({td: da['valid_time'].values})
        da = da.rename({td: 'time'})
        for c in list(da.coords):
            if c not in ('time', 'latitude', 'longitude'):
                da = da.drop_vars(c)
        parts.append(da.transpose('time', 'latitude', 'longitude').load())
        if not parts[1:]:
            print('GloFAS variable %s, first record valid %s (written at -12 h)'
                  % (cand[0], str(da.time.values[0])[:16]))
    da = xr.concat(parts, dim='time').sortby('time')
    _, idx = np.unique(da.time.values, return_index=True)
    da = da.isel(time=np.sort(idx))
    t_centre = da.time.values - np.timedelta64(12, 'h')        # dis24 = mean over [t-24h, t]
    return t_centre, da.latitude.values, da.longitude.values, da.values


def _dai_mouths(dai_file):
    import netCDF4
    with netCDF4.Dataset(dai_file) as d:
        raw = d['riv_name'][:]
        names = [row.tobytes().decode('latin-1', errors='replace').replace('\x00', '').strip()
                 for row in raw]
        return names, np.asarray(d['lon_mou'][:], float), np.asarray(d['lat_mou'][:], float)


def _safe(s):
    import unicodedata
    s = unicodedata.normalize('NFKD', s.split('(')[0].strip())      # Ogooué -> Ogooue
    s = s.replace(' ', '_').replace("'", '').replace('/', '_')
    return s.encode('ascii', 'ignore').decode() or 'river'


def make_glofas_river_files(grd_file, glofas_dir, out_dir, month_start, month_end, dai_file,
                            qmin=100.0, margin=0.0, radius=0.25, extra_list=None):
    """Write <river>.txt + river_list.txt (input of make_river_run.py) from GloFAS."""
    os.makedirs(out_dir, exist_ok=True)
    lo0, lo1, la0, la1 = grid_box(grd_file)
    t, glat, glon, q = open_glofas(glofas_dir, month_start, month_end)
    qmean = np.nanmean(np.where(np.isfinite(q), q, np.nan), axis=0)
    LON, LAT = np.meshgrid(glon, glat)

    cands = []
    names, lonm, latm = _dai_mouths(dai_file)
    for nm, lo, la in zip(names, lonm, latm):
        if lo0 - margin < lo < lo1 + margin and la0 - margin < la < la1 + margin:
            cands.append((nm, float(lo), float(la), 'Dai'))
    if extra_list:
        for line in open(extra_list):
            line = line.split('#')[0].split()
            if len(line) >= 3:
                cands.append((line[0], float(line[1]), float(line[2]), 'user'))

    print('Grid box %.2f..%.2f E  %.2f..%.2f N ; %d candidate mouth(s) ; Qmin %.0f m3/s, radius %.2f deg'
          % (lo0, lo1, la0, la1, len(cands), qmin, radius))
    # outlet cell of each mouth = largest mean GloFAS discharge within `radius`
    found = []
    for nm, lo, la, src in cands:
        near = ((LON - lo) ** 2 + (LAT - la) ** 2 <= radius ** 2) & np.isfinite(qmean)
        if not near.any():
            continue
        jj, ii = np.where(near)
        k = np.argmax(qmean[jj, ii]); j, i = jj[k], ii[k]
        if qmean[j, i] >= qmin:
            found.append((nm, lo, la, src, j, i, float(qmean[j, i])))
    # strongest first; a mouth whose outlet is on/next to (<= 2 cells, 0.1 deg) one
    # already taken is the same channel seen twice -> dropped (no double count)
    found.sort(key=lambda c: -c[-1])
    chosen, used = [], []
    for c in found:
        j, i = c[4], c[5]
        if any(abs(j - u) <= 2 and abs(i - v) <= 2 for u, v in used):
            print('  skip %-20s: same GloFAS channel as a stronger mouth' % c[0][:20])
            continue
        used.append((j, i))
        chosen.append(c)

    chosen.sort(key=lambda c: -c[-1])
    print('%-22s %8s %7s %8s %8s %10s %s' % ('river', 'lon_mou', 'lat_mou', 'glofas_x', 'glofas_y', 'mean m3/s', 'src'))
    lines = ['#  input_file     |    Lon    |  Lat  |   (GloFAS %s, %s..%s)' % (DATASET, month_start, month_end)]
    seen = {}
    for nm, lo, la, src, j, i, qm in chosen:
        base = _safe(nm); seen[base] = seen.get(base, 0) + 1
        fname = base + ('' if seen[base] == 1 else '_%d' % seen[base]) + '.txt'
        series = q[:, j, i]
        with open(os.path.join(out_dir, fname), 'w') as f:
            for tt, qq in zip(t, series):
                if np.isfinite(qq):
                    d = datetime.datetime(1970, 1, 1) + datetime.timedelta(
                        seconds=float((tt - np.datetime64('1970-01-01T00:00:00')) / np.timedelta64(1, 's')))
                    f.write('%s %.3f\n' % (d.strftime('%Y/%m/%d %H:%M:%S'), qq))
        lines.append('%-24s %9.3f %9.3f' % (fname, lo, la))
        print('%-22s %8.2f %7.2f %8.2f %8.2f %10.1f %s' % (nm[:22], lo, la, glon[i], glat[j], qm, src))
    with open(os.path.join(out_dir, 'river_list.txt'), 'w') as f:
        f.write('\n'.join(lines) + '\n')
    print('Wrote %d river file(s) + river_list.txt to %s' % (len(chosen), out_dir))
    return len(chosen)
