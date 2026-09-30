#!/bin/bash
# =============================================================================
# Step 02 -- read the grid and WRITE the derived configuration files
#   hindcast/steps/02_check_grid.sh NAME
#
# From croco_grd.nc it decides / checks, and writes into hindcast/configs/NAME/:
#   obc.cfg                    open boundaries from the land mask (an edge is OPEN
#                              when > 50 % of it is ocean), unless OBC is set in domain.cfg
#   crocotools_param.py        GLORYS ini/bry parameters (inputdata 'mercator', sigma
#                              params and obc_dict from domain.cfg / obc.cfg)
#   crocotools_param_tides.py  TPXO7 parameters (inputdata 'tpxo7_croco'), for step 09
#   CARD.md                    the region card (box, grid, boundaries, dt, CFL, ...)
# and prints: size -> LLm0/MMm0, mask strips, isolated water bodies, rx0, CFL.
# =============================================================================
set -e -o pipefail
source "$(dirname "$0")/../lib/common.sh" "$1"
[ -f "${CF}/croco_grd.nc" ] || die "no grid -- run hindcast/steps/01_make_grid.sh ${CONFIG_NAME}"

python3 - <<'PYEOF'
import os, numpy as np, xarray as xr
from scipy import ndimage
E = os.environ
g = xr.open_dataset(os.path.join(E["CF"], "croco_grd.nc"))
m = g.mask_rho.values
xi, eta = g.sizes["xi_rho"], g.sizes["eta_rho"]
print("=== grid size ===")
print("xi_rho=%d eta_rho=%d  ->  param.h LLm0=%d MMm0=%d N=%s" % (xi, eta, xi - 2, eta - 2, E["N_LEVELS"]))

print("\n=== boundaries from the mask (O ocean, . land) ===")
edges = {"south": m[0, :], "west": m[:, 0], "east": m[:, -1], "north": m[-1, :]}
auto = {k: int(v.mean() > 0.5) for k, v in edges.items()}
for k, v in edges.items():
    s = "".join("O" if x == 1 else "." for x in v)
    if len(s) > 110: s = s[:108] + ".."
    print("%-5s %4d/%-4d ocean (%3.0f%%) -> %s\n      %s" % (k, v.sum(), v.size, 100 * v.mean(),
          "OPEN" if auto[k] else "closed", s))
if E.get("OBC", "").strip():
    vals = [int(x) for x in E["OBC"].split()]
    obc = dict(zip(["south", "west", "east", "north"], vals)); src = "override (OBC in domain.cfg)"
    for k in obc:
        if obc[k] != auto[k]:
            print("  NOTE: %s forced to %d, the mask suggests %d" % (k, obc[k], auto[k]))
else:
    obc = auto; src = "land mask"
print("-> obc (%s): south=%d west=%d east=%d north=%d" % (src, obc["south"], obc["west"], obc["east"], obc["north"]))
if sum(obc.values()) == 0:
    print("  WARNING: all boundaries closed -- no ocean enters from GLORYS; check the box")

print("\n=== isolated water bodies ===")
lab, n = ndimage.label(m)
sizes = ndimage.sum(m, lab, range(1, n + 1)) if n else []
print("ocean %.1f%% of the grid, depth %.0f..%.0f m, %d water body(ies)" % (m.mean() * 100, float(g.h.min()), float(g.h.max()), n))
for i in list(np.argsort(sizes)[::-1][:5]):
    j, k = np.argwhere(lab == i + 1)[0]
    print("  %7d cells at lon %.2f lat %.2f" % (sizes[i], g.lon_rho.values[j, k], g.lat_rho.values[j, k]))

print("\n=== smoothing ===")
hm = np.where(m > 0, g.h.values, np.nan)
rx = np.nanmax(np.abs(np.diff(hm, axis=1)) / (hm[:, 1:] + hm[:, :-1]))
ry = np.nanmax(np.abs(np.diff(hm, axis=0)) / (hm[1:, :] + hm[:-1, :]))
print("rx0 = %.3f  ry0 = %.3f  (target 0.2)" % (rx, ry))

print("\n=== time step / CFL ===")
dt, ndt = float(E["DT"]), int(E["NDTFAST"])
dx = 1.0 / np.maximum(g.pm.values, g.pn.values)
cb = np.nanmax(np.sqrt(9.81 * g.h.values) * (dt / ndt) / dx)
cc = np.nanmax(2.0 * dt / dx)
print("dt = %g s, NDTFAST = %d, min dx = %.0f m" % (dt, ndt, np.nanmin(dx)))
print("barotropic Courant %.2f   baroclinic (2 m/s) %.2f   (keep below ~0.7)" % (cb, cc))
if cb > 0.7 or cc > 0.7:
    print("  WARNING: Courant > 0.7 -- set a smaller DT in domain.cfg")

# ---- write the derived files ---------------------------------------------------
cd = E["CONFIG_DIR"]
with open(os.path.join(cd, "obc.cfg"), "w") as f:
    f.write("# written by steps/02_check_grid.sh from the land mask (%s) -- do not edit;\n"
            "# to force boundaries set OBC=\"south west east north\" in domain.cfg\n" % src)
    for k in ["south", "west", "east", "north"]:
        f.write("OBC_%s=%d\n" % (k.upper(), obc[k]))

hdr = "# written by steps/02_check_grid.sh from domain.cfg + obc.cfg -- do not edit, re-run step 02\n"
with open(os.path.join(cd, "crocotools_param.py"), "w") as f:
    f.write(hdr + "# GLORYS and the Mercator analysis both read through the 'mercator' reader\n"
            "# (zos/thetao/so/uo/vo); there is no 'glorys' key. Ocean source: %s\n"
            "inputdata    = 'mercator'\nNzgoodmin    = 4\nmulti_files  = False\n"
            "tracers      = ['temp', 'salt']\ncroco_grd    = 'croco_grd.nc'\n"
            "sigma_params = dict(theta_s=%s, theta_b=%s, N=%s, hc=%s)\n"
            "ini_prefix   = 'croco_ini_%s'\nbry_prefix   = 'croco_bry_%s'\n"
            "obc_dict     = dict(south=%d, west=%d, east=%d, north=%d)\ncycle_bry    = 0\n"
            % (E["OCEAN_SOURCE"], E["THETA_S"], E["THETA_B"], E["N_LEVELS"], E["SIGMA_HC"],
               E["OCEAN_TAG"], E["OCEAN_TAG"],
               obc["south"], obc["west"], obc["east"], obc["north"]))
with open(os.path.join(cd, "crocotools_param_tides.py"), "w") as f:
    f.write(hdr + "# TPXO7 single file (DATASETS_CROCOTOOLS/TPXO7/TPXO7.nc); reader key 'tpxo7_croco'\n"
            "# ('tpxo7' alone is not a key of gtools/croco_pytools/prepro/Readers/tides_reader.py).\n"
            "# For TPXO10 atlas files: inputdata='tpxo10', multi_files=True, waves_separated=True,\n"
            "# elev_file='h_<tides>_tpxo10_atlas_30_v2.nc', u_file=v_file='u_<tides>_tpxo10_atlas_30_v2.nc'.\n"
            "inputdata       = 'tpxo7_croco'\ninput_file      = 'TPXO7.nc'\ninput_type      = 'Re_Im'\n"
            "multi_files     = False\nwaves_separated = False\nelev_file = ''\nu_file = ''\nv_file = ''\n"
            "croco_grd       = 'croco_grd.nc'\n"
            "tides           = ['M2', 'S2', 'N2', 'K2', 'K1', 'O1', 'P1', 'Q1', 'Mf', 'Mm']\n"
            "cur             = True\npot             = True\nCorrection_ssh  = True\nCorrection_uv   = True\n")

lon, lat = g.lon_rho.values, g.lat_rho.values
opn = [k for k in obc if obc[k]]; cls = [k for k in obc if not obc[k]]
with open(os.path.join(cd, "CARD.md"), "w") as f:
    f.write("""# %s -- region card (generated by steps/02_check_grid.sh)

| | |
|---|---|
| **Box** | lon %s .. %s, lat %s .. %s (grid lon %.2f..%.2f, lat %.2f..%.2f) |
| **Resolution** | 1/%s deg, %s sigma levels (theta_s %s, theta_b %s, hc %s m) |
| **Grid** | %d x %d rho points (LLm0=%d, MMm0=%d) |
| **Boundaries** | open: %s; closed: %s (%s) |
| **Depth** | %.0f .. %.0f m, ocean %.1f %% of the grid, rx0 %.3f |
| **Time step** | dt %s s, NDTFAST %s, barotropic Courant %.2f |
| **Period** | %s -> %s (%s cycles of %s d spin-up + %s d), Yorig %s |
| **Data** | ocean box %s, months %s..%s; ERA5 months %s..%s |
| **Physics** | tides=%s rivers=%s |
| **Build** | `hindcast/new_config.sh %s %s %s %s %s --res %s` |
""" % (E["CONFIG_NAME"], E["LON_MIN"], E["LON_MAX"], E["LAT_MIN"], E["LAT_MAX"],
       lon.min(), lon.max(), lat.min(), lat.max(), E["RES_INV"], E["N_LEVELS"], E["THETA_S"],
       E["THETA_B"], E["SIGMA_HC"], xi, eta, xi - 2, eta - 2, ", ".join(opn) or "none",
       ", ".join(cls) or "none", src, float(g.h.min()), float(g.h.max()), m.mean() * 100, max(rx, ry),
       E["DT"], E["NDTFAST"], cb, E["HC_START"], E["HC_END"], E["NCYCLES"], E["SPINUP_DAYS"],
       E["HCAST_DAYS"], E["YORIG"], E["EXTENTS"], E["GLORYS_MONTH_START"], E["GLORYS_MONTH_END"],
       E["ERA5_MONTH_START"], E["ERA5_MONTH_END"], E["TIDES"], E["RIVERS"], E["CONFIG_NAME"],
       E["LON_MIN"], E["LON_MAX"], E["LAT_MIN"], E["LAT_MAX"], E["RES_INV"]))
print("\nwrote %s/{obc.cfg, crocotools_param.py, crocotools_param_tides.py, CARD.md}" % cd)
PYEOF
cp "${CONFIG_DIR}/crocotools_param.py" "${CF}/crocotools_param.py"
cp "${CONFIG_DIR}/crocotools_param_tides.py" "${CF}/crocotools_param_tides.py"
say "next: hindcast/steps/03_download_data.sh ${CONFIG_NAME}"
