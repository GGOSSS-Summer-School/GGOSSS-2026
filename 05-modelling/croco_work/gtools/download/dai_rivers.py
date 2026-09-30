import netCDF4 as nc
import numpy as np
import datetime
import os
import sys

'''
dai_rivers.py
=============
Extract the rivers relevant to a CROCO domain from the Dai & Trenberth (2002)
global runoff climatology, formatted for make_rivers.py.

Reads the domain bounds AUTOMATICALLY from the CROCO grid file, so it works
for ANY region -- just point it at that region's croco_grd.nc.

Usage:
  python dai_rivers.py <croco_grd.nc> [output_dir] [Qmin] [margin]

Outputs (in output_dir, default ./DAI_RIVERS/):
  - one <river>.txt per river : YYYY/MM/DD hh:mm:ss  Qsrc  (daily, 1-yr cycle)
  - river_list.txt            : input for make_rivers.py

Only flow is handled (Dai has no T/S); set river T/S in croco.in.
'''

DAI = os.path.expanduser(
    '~/ggosss26/data/DATASETS_CROCOTOOLS/RUNOFF_DAI/Dai_Trenberth_runoff_global_clim.nc')
REF_YEAR = 2000

# ---- args ----
if len(sys.argv) < 2:
    print("Usage: python dai_rivers.py <croco_grd.nc> [output_dir] [Qmin] [margin]")
    sys.exit(1)
grd_file  = sys.argv[1]
OUTDIR    = sys.argv[2] if len(sys.argv) > 2 else './DAI_RIVERS/'
QMIN      = float(sys.argv[3]) if len(sys.argv) > 3 else 300.0
MARGIN    = float(sys.argv[4]) if len(sys.argv) > 4 else 0.6

# ---- read domain bounds from the CROCO grid ----
g = nc.Dataset(grd_file)
glon = g['lon_rho'][:]; glat = g['lat_rho'][:]
LON_MIN, LON_MAX = float(glon.min()), float(glon.max())
LAT_MIN, LAT_MAX = float(glat.min()), float(glat.max())
g.close()
print(f"Grid bounds from {grd_file}:")
print(f"  lon {LON_MIN:.2f} .. {LON_MAX:.2f}   lat {LAT_MIN:.2f} .. {LAT_MAX:.2f}")
print(f"  (margin {MARGIN} deg, Qmin {QMIN} m3/s)")

os.makedirs(OUTDIR, exist_ok=True)
d = nc.Dataset(DAI)
raw = d['riv_name'][:]
names = [row.tobytes().decode('latin-1', errors='replace').replace('\x00','').strip()
         for row in raw]
lonm = d['lon_mou'][:]; latm = d['lat_mou'][:]
flow_clm = d['FLOW_clm'][:]
month_mid_doy = [15,45,74,105,135,166,196,227,258,288,319,349]

def safe_name(s):
    return (s.split('(')[0].strip().replace(' ','_').replace("'","")
             .encode('ascii','ignore').decode() or 'river')

selected = []
for i, nm in enumerate(names):
    lo, la = float(lonm[i]), float(latm[i])
    if not (LON_MIN-MARGIN < lo < LON_MAX+MARGIN and LAT_MIN-MARGIN < la < LAT_MAX+MARGIN):
        continue
    q12 = np.ma.filled(flow_clm[:, i], np.nan).astype(float)
    if np.all(np.isnan(q12)):
        continue
    ann = np.nanmean(q12)
    if ann < QMIN:
        continue
    selected.append((i, nm, lo, la, ann, q12))

selected.sort(key=lambda x: -x[4])
print(f"\nSelected {len(selected)} rivers:")
print(f"{'name':18s} {'lon':>7s} {'lat':>7s} {'ann_m3/s':>10s}")

river_list_lines = ["#  input_file     |    Lon    |  Lat/Qmin  |"]
for i, nm, lo, la, ann, q12 in selected:
    sn = safe_name(nm)
    print(f"{nm[:18]:18s} {lo:7.2f} {la:7.2f} {ann:10.1f}")
    doy = np.arange(1, 366)
    xp = np.array([month_mid_doy[-1]-365] + month_mid_doy + [month_mid_doy[0]+365])
    fp = np.array([q12[-1]] + list(q12) + [q12[0]])
    q_daily = np.interp(doy, xp, fp)
    fname = f"{sn}.txt"
    with open(os.path.join(OUTDIR, fname), 'w') as f:
        base = datetime.date(REF_YEAR, 1, 1)
        for k, dd in enumerate(doy):
            date = base + datetime.timedelta(days=int(dd)-1)
            f.write(f"{date.strftime('%Y/%m/%d')} 12:00:00 {q_daily[k]:.3f}\n")
    river_list_lines.append(f"{fname:20s} {lo:8.2f} {la:9.2f}")

with open(os.path.join(OUTDIR, 'river_list.txt'), 'w') as f:
    f.write('\n'.join(river_list_lines) + '\n')

print(f"\nWrote {len(selected)} river .txt files + river_list.txt to {OUTDIR}")
