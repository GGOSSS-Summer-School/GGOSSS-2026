#!/usr/bin/env python3
"""Post-run verification of a cycled hindcast.

    python3 notebooks/verify_config.py NAME [BUILD]
        BUILD = plain (default) | plain_tides | plain_rivers | plain_tides_rivers

For every cycle folder hindcast/model-runs/NAME/<YYYYMMDD>[_<BUILD>]/ it checks that
the hindcast phase ended with MAIN: DONE, prints its time span (from Yorig) and the
surface T/S range, and checks that consecutive cycles tile the period without gap or
overlap. The last record is compared with the GLORYS month it falls in.
"""
import glob
import os
import re
import sys

import numpy as np
import xarray as xr

if len(sys.argv) < 2:
    sys.exit(__doc__)
NAME = sys.argv[1]
BUILD = sys.argv[2] if len(sys.argv) > 2 else "plain"
ROOT = os.environ.get("CROCO_ROOT", os.path.expanduser("~/croco_work"))
cfg = open(os.path.join(ROOT, "hindcast", "configs", NAME, "domain.cfg")).read()
m = re.search(r"^YORIG=(\d+)", cfg, re.M)
YORIG = int(m.group(1)) if m else 1993
OUT = os.path.join(ROOT, "hindcast", "model-runs", NAME)
GLORYS = os.path.join(ROOT, "hindcast", "scratch", NAME, "downloaded_data", "GLORYS")
suffix = "" if BUILD == "plain" else "_" + BUILD
cycles = sorted(d for d in glob.glob(os.path.join(OUT, "[0-9]" * 8 + suffix)) if os.path.isdir(d))
print(f"{NAME} / {BUILD}: {len(cycles)} cycle folder(s) under {OUT}  (Yorig {YORIG})")

origin = np.datetime64(f"{YORIG}-01-01")
last, prev_end, ok = None, None, 0
for d in cycles:
    tag = os.path.basename(d)
    his = os.path.join(d, "hcast", "CROCO_FILES", "croco_his.nc")
    log = os.path.join(d, "hcast", "croco_hcast.out")
    done = os.path.exists(log) and "MAIN: DONE" in open(log, errors="ignore").read()
    if not (os.path.exists(his) and done):
        print(f"  {tag}: FAILED -- check {log}")
        continue
    with xr.open_dataset(his, decode_times=False) as ds:
        t = origin + (ds.scrum_time.values * 1e9).astype("timedelta64[ns]")
        mk = ds.mask_rho.values == 1
        T = ds.temp.isel(time=-1, s_rho=-1).values[mk]; S = ds.salt.isel(time=-1, s_rho=-1).values[mk]
        gap = "" if prev_end is None else ("" if abs((t[0] - prev_end) / np.timedelta64(1, "h")) < 6.01
                                           else "  <-- GAP/OVERLAP with the previous cycle")
        print(f"  {tag}: OK  {str(t[0])[:16]} -> {str(t[-1])[:16]} ({ds.sizes['time']} rec)"
              f"  SST {T.min():.2f}..{T.max():.2f}  SSS {S.min():.2f}..{S.max():.2f}{gap}")
        prev_end, ok, last = t[-1], ok + 1, his
print(f"{ok}/{len(cycles)} cycle(s) OK")

if last:
    with xr.open_dataset(last, decode_times=False) as ds:
        tend = origin + (ds.scrum_time.values[-1] * 1e9).astype("timedelta64[ns]")
        sst = float(ds.temp.isel(time=-1, s_rho=-1).where(ds.mask_rho == 1).mean())
    f = os.path.join(GLORYS, str(tend)[:7].replace("-", "_") + ".nc")
    if os.path.exists(f):
        with xr.open_dataset(f) as ref:
            r = float(ref.thetao.sel(time=tend, method="nearest").isel(depth=0).mean())
        print(f"last record {str(tend)[:10]}: domain-mean SST CROCO {sst:.2f} degC vs GLORYS {r:.2f} degC (GLORYS box is wider)")
