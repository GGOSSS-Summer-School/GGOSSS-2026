#!/usr/bin/env python3
"""
ggosss26.py — the CROCO command-line interface.

A single entry point for CROCO workflows, modelled on the somisana cli.py
pattern. Each subcommand produces an output (downloads, forcing, ini/bry, ...).

Exposes the FORECAST pipeline (Mercator + GFS) and the HINDCAST pipeline
(GLORYS + ERA5, the *_hindcast subcommands), so that `ggosss26.py` remains the
one front door for the whole project.

Usage examples
--------------
    python ggosss26.py download_ocean \
        --domain 3.5,13.5,-6.5,6.5 --run_date "2026-07-01 00:00:00" \
        --hdays 2 --fdays 5 --outputDir <dir>

    python ggosss26.py make_ini \
        --input_file <MERCATOR.nc> --output_dir <CROCO_FILES> \
        --run_date "2026-07-01 00:00:00" --hdays 2 --Yorig 2000

Attribution
-----------
The underlying forecast functions are adapted, with permission, from the SAEON
somisana project (somisana-download, somisana-croco). See
code/ggosss26/forecast/*.py headers and ACKNOWLEDGEMENTS.md.
"""
import argparse
import sys
import os
from datetime import datetime, timedelta

# --- make the CROCO forecast package importable -----------------------
# forecast functions live alongside this file, under ./forecast/
_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _HERE)

from download.cmems import download_mercator_ops, download_cmems_monthly, download_mercator_monthly
from download.gfs import download_gfs_atm
from preprocess import reformat_gfs_atm, make_ini, make_bry, make_tides


# --- small argument parsers -------------------------------------------------
def parse_datetime(value):
    try:
        return datetime.strptime(value, '%Y-%m-%d %H:%M:%S')
    except ValueError:
        raise argparse.ArgumentTypeError(
            "Invalid datetime. Use 'YYYY-MM-DD HH:MM:SS'.")

def parse_list(value):
    return [float(x) for x in value.split(',')]


def main():
    parser = argparse.ArgumentParser(
        description='CROCO command-line interface')
    sub = parser.add_subparsers(dest='function',
                                help='Select the CROCO function to run')

    # =====================================================================
    # FORECAST pipeline
    # =====================================================================

    # ---- download_ocean : Mercator anfc --------------------------------
    p = sub.add_parser('download_ocean',
        help='Download the global ocean forecast (Mercator anfc) and merge it')
    p.add_argument('--domain', required=True, type=parse_list,
                   help='lon_min,lon_max,lat_min,lat_max (pad ~1.5 deg beyond grid)')
    p.add_argument('--run_date', required=True, type=parse_datetime,
                   help='forecast T0, "YYYY-MM-DD HH:MM:SS"')
    p.add_argument('--hdays', required=True, type=int, help='spin-up days before T0')
    p.add_argument('--fdays', required=True, type=int, help='forecast days after T0')
    p.add_argument('--outputDir', required=True, help='where to save the ocean data')
    p.add_argument('--usrname', required=False, default=None,
                   help='CMEMS username (optional; interactive login if omitted)')
    p.add_argument('--passwd', required=False, default=None,
                   help='CMEMS password (optional)')
    def _download_ocean(a):
        download_mercator_ops(a.domain, a.run_date, a.hdays, a.fdays,
                              a.outputDir, usrname=a.usrname, passwd=a.passwd)
    p.set_defaults(func=_download_ocean)

    # ---- download_atmosphere : GFS -------------------------------------
    p = sub.add_parser('download_atmosphere',
        help='Download the global weather forecast (NOAA GFS) as GRIB')
    p.add_argument('--domain', required=True, type=parse_list,
                   help='lon_min,lon_max,lat_min,lat_max')
    p.add_argument('--run_date', required=True, type=parse_datetime)
    p.add_argument('--hdays', required=True, type=int)
    p.add_argument('--fdays', required=True, type=int)
    p.add_argument('--outputDir', required=True, help='where to save the GRIB files')
    def _download_atmosphere(a):
        download_gfs_atm(a.domain, a.run_date, a.hdays, a.fdays, a.outputDir)
    p.set_defaults(func=_download_atmosphere)

    # ---- make_forcing : GFS GRIB -> CROCO online files -----------------
    p = sub.add_parser('make_forcing',
        help='Reshape downloaded GFS into CROCO online (surface forcing) files')
    p.add_argument('--gfsDir', required=True, help='directory of downloaded GFS GRIB')
    p.add_argument('--outputDir', required=True, help='where to write the online files')
    p.add_argument('--Yorig', required=True, type=int,
                   help='time origin year (CROCO time = days/seconds since Yorig-01-01)')
    def _make_forcing(a):
        os.makedirs(a.outputDir, exist_ok=True)
        reformat_gfs_atm(a.gfsDir, a.outputDir, a.Yorig)
    p.set_defaults(func=_make_forcing)

    # ---- make_ini : initial condition ----------------------------------
    # make_ini() itself takes an explicit ini_date + output filename, so this
    # handler reads crocotools_param.py for the prefix, computes the date, and
    # builds the filename (same convention as the somisana *_fcst wrappers).
    p = sub.add_parser('make_ini',
        help='Build the CROCO initial condition from the ocean forecast')
    p.add_argument('--input_file', required=True, help='merged Mercator .nc file')
    p.add_argument('--output_dir', required=True,
                   help='CROCO_FILES dir (must contain crocotools_param.py)')
    p.add_argument('--run_date', required=True, type=parse_datetime)
    p.add_argument('--hdays', required=True, type=int)
    p.add_argument('--Yorig', required=True, type=int)
    def _make_ini(a):
        sys.path.append(a.output_dir)
        import crocotools_param as params
        fname_out = params.ini_prefix + a.run_date.strftime('_%Y%m%d_%H.nc')
        ini_date = a.run_date - timedelta(days=a.hdays)
        make_ini(a.input_file, a.output_dir, ini_date, a.Yorig, fname_out)
    p.set_defaults(func=_make_ini)

    # ---- make_bry : boundary conditions --------------------------------
    p = sub.add_parser('make_bry',
        help='Build the CROCO boundary conditions from the ocean forecast')
    p.add_argument('--input_file', required=True, help='merged Mercator .nc file')
    p.add_argument('--output_dir', required=True,
                   help='CROCO_FILES dir (must contain crocotools_param.py)')
    p.add_argument('--run_date', required=True, type=parse_datetime)
    p.add_argument('--hdays', required=True, type=int)
    p.add_argument('--fdays', required=True, type=int)
    p.add_argument('--Yorig', required=True, type=int)
    def _make_bry(a):
        sys.path.append(a.output_dir)
        import crocotools_param as params
        # 2-day buffer either side (make_bry uses nearest available times)
        hdays = a.hdays + 2
        fdays = a.fdays + 2
        fname_out = params.bry_prefix + a.run_date.strftime('_%Y%m%d_%H.nc')
        start_date = a.run_date - timedelta(days=hdays)
        end_date = a.run_date + timedelta(days=fdays)
        make_bry(a.input_file, a.output_dir, start_date, end_date, a.Yorig, fname_out)
    p.set_defaults(func=_make_bry)
    
    # ---- make_tides : tidal forcing from TPXO -------------------------
    # Regenerated PER CYCLE. With USE_CALENDAR off, the tidal phase epoch is
    # baked into croco_frc.nc at generation time, keyed to run_date -- so the
    # tide file is NOT a build-once asset like the grid; each cycle needs its
    # own, at that cycle's start date.
    #
    # make_tides() reads its parameters from a file named crocotools_param.py
    # in output_dir. Its 'inputdata' is a TPXO tag ('tpxo7_croco'), which would
    # clash with the 'mercator' one make_ini/make_bry read -- so give make_tides
    # its OWN gen dir holding only the tide params (same pattern as child IC).
    p = sub.add_parser('make_tides',
        help='Build the CROCO tidal forcing (croco_frc.nc) from TPXO')
    p.add_argument('--input_dir', required=True,
                   help='dir with the raw TPXO file (e.g. .../DATASETS_CROCOTOOLS/TPXO7/)')
    p.add_argument('--output_dir', required=True,
                   help='gen dir with crocotools_param.py (tide params) + croco_grd.nc')
    p.add_argument('--run_date', required=True, type=parse_datetime,
                   help='cycle start date -- sets the tidal phase reference')
    p.add_argument('--Yorig', required=True, type=int)
    p.add_argument('--fname_out', default='croco_frc.nc')
    def _make_tides(a):
        make_tides(a.input_dir, a.output_dir, a.run_date, a.Yorig, a.fname_out)
    p.set_defaults(func=_make_tides)


    # =====================================================================
    # HINDCAST pipeline (GLORYS + ERA5, Yorig=1993)
    # ---------------------------------------------------------------------
    #   download_ocean_hindcast      -> GLORYS daily, one YYYY_MM.nc per month
    #   download_atmosphere_hindcast -> ERA5 request + convert -> for_croco/
    #   make_ini_hindcast            -> croco_ini_GLORYS_YyyyyMmmDdd.nc
    #   make_bry_hindcast            -> croco_bry_GLORYS_Y..._to_Y....nc
    # =====================================================================
  
    # ---- download_ocean_hindcast : GLORYS monthly reanalysis -----------
    p = sub.add_parser('download_ocean_hindcast',
        help='Download GLORYS monthly ocean reanalysis (CMEMS) for a hindcast')
    p.add_argument('--domain', required=True, type=parse_list,
                   help='lon_min,lon_max,lat_min,lat_max (pad ~1.5 deg beyond grid)')
    p.add_argument('--month_start', required=True, help='first month, "YYYY-MM"')
    p.add_argument('--month_end', required=True, help='last month, "YYYY-MM"')
    p.add_argument('--outputDir', required=True, help='where to save the GLORYS monthly files')
    # Default is the DAILY reanalysis: a monthly-mean file (P1M-m) has one
    # record per month and make_ini_hindcast/make_bry_hindcast cannot bracket
    # a date with it (IndexError: index 1 is out of bounds ...).
    p.add_argument('--product_id', required=False,
                   default='cmems_mod_glo_phy_my_0.083deg_P1D-m',
                   help='CMEMS GLORYS product id (default: DAILY reanalysis P1D-m)')
    p.add_argument('--usrname', required=False, default=None,
                   help='CMEMS username (optional; interactive login if omitted)')
    p.add_argument('--passwd', required=False, default=None,
                   help='CMEMS password (optional)')
    p.add_argument('--source', required=False, default='glorys', choices=['glorys', 'mercator'],
                   help='glorys = daily reanalysis (default); mercator = Mercator global '
                        'analysis (anfc), for periods GLORYS does not cover yet')
    def _download_ocean_hindcast(a):
        start_date = datetime.strptime(a.month_start, '%Y-%m')
        end_date   = datetime.strptime(a.month_end,   '%Y-%m')
        varlist = ['zos', 'uo', 'vo', 'thetao', 'so']
        depths  = [0.493, 5727.918]
        if a.source == 'mercator':
            download_mercator_monthly(a.domain, start_date, end_date, depths, a.outputDir,
                                      usrname=a.usrname, passwd=a.passwd)
            return
        download_cmems_monthly(a.product_id, a.domain, start_date, end_date,
                               varlist, depths, a.outputDir,
                               usrname=a.usrname, passwd=a.passwd)
    p.set_defaults(func=_download_ocean_hindcast)
    # ---- download_atmosphere_hindcast : ERA5 -------------------------
    p = sub.add_parser('download_atmosphere_hindcast',
        help='Download ERA5 atmosphere (CDS) and convert to CROCO online forcing')
    p.add_argument('--domain', required=True, type=parse_list,
                   help='lon_min,lon_max,lat_min,lat_max (grid box; a 2 deg margin is added)')
    p.add_argument('--month_start', required=True, help='"YYYY-MM"')
    p.add_argument('--month_end',   required=True, help='"YYYY-MM"')
    p.add_argument('--outputDir', required=True,
                   help='ERA5 base dir (raw/ + for_croco/ created inside)')
    p.add_argument('--Yorig', required=False, type=int, default=1993)
    p.add_argument('--legacy', action='store_true',
                   help='use the old ERA5_request.py + ERA5_convert.py path (q from the '
                        'slow 1000 hPa pressure-levels product) instead of era5_fast.py')
    def _download_atmosphere_hindcast(a):
        import subprocess
        if not a.legacy:
            sys.path.insert(0, os.path.join(_HERE, 'download', 'ERA5'))
            if os.environ.get('ERA5_SOURCE', 'arco') == 'arco':
                # default: CDS ARCO Zarr store, read directly (no request queue)
                import era5_arco
                era5_arco.run(a.domain, a.month_start, a.month_end, a.outputDir, a.Yorig)
            else:
                # ERA5_SOURCE=cds: CDS API, single-level requests only, q from d2m + msl
                import era5_fast
                era5_fast.run(a.domain, a.month_start, a.month_end, a.outputDir, a.Yorig)
            return
        lon_min, lon_max, lat_min, lat_max = a.domain
        ys, ms = a.month_start.split('-'); ye, me = a.month_end.split('-')
        raw  = os.path.join(a.outputDir, 'raw')
        proc = os.path.join(a.outputDir, 'for_croco')
        os.makedirs(raw, exist_ok=True); os.makedirs(proc, exist_ok=True)
        era5_dir = os.path.join(_HERE, 'download', 'ERA5')
        hours = "[" + ",".join(f"'{h:02d}:00'" for h in range(24)) + "]"
        with open(os.path.join(era5_dir, 'era5_crocotools_param.py'), 'w') as f:
            f.write(
"era5_dir_raw = '" + raw + "'\n"
"era5_dir_processed = '" + proc + "'\n"
"wave_extract = False\n"
"pressure_extract = True\n"
f"year_start = {int(ys)}\n"
f"month_start = {int(ms)}\n"
f"year_end = {int(ye)}\n"
f"month_end = {int(me)}\n"
f"Yorig = {a.Yorig}\n"
"n_overlap = 0\n"
f"times = {hours}\n"
f"lonmin = {lon_min}\n"
f"lonmax = {lon_max}\n"
f"latmin = {lat_min}\n"
f"latmax = {lat_max}\n"
"cff_tp = 1000./3600.\n"
"cff_heat = 1./3600.\n"
"variables = ['lsm','sst','tp','strd','ssr','t2m','q','u10','v10']\n"
"conv_cff  = [1.,1.,cff_tp,cff_heat,cff_heat,1.,1.,1.,1.]\n"
"units     = ['(0-1)','K','kg m-2 s-1','W m-2','W m-2','K','kg kg-1','m s-1','m s-1']\n"
"if pressure_extract:\n"
"    variables.append('msl'); conv_cff.append(1.); units.append('Pa')\n"
            )
        # sys.executable = the interpreter of the active (ggosss26) env
        subprocess.run([sys.executable, 'ERA5_request.py'], cwd=era5_dir, check=True)
        subprocess.run([sys.executable, 'ERA5_convert.py'], cwd=era5_dir, check=True)
    p.set_defaults(func=_download_atmosphere_hindcast)
    # ---- helper: monthly GLORYS files spanning a date window -----------
    def _glorys_months(input_dir, d0, d1):
        files, seen = [], set()
        y, m = d0.year, d0.month
        while (y, m) <= (d1.year, d1.month):
            key = '%04d_%02d' % (y, m)
            if key not in seen:
                f = os.path.join(input_dir, key + '.nc')
                if os.path.exists(f):
                    files.append(f); seen.add(key)
            m += 1
            if m > 12: m = 1; y += 1
        return files

    # ---- make_ini_hindcast : GLORYS initial condition (date-based) -----
    p = sub.add_parser('make_ini_hindcast',
        help='Build CROCO initial condition from monthly GLORYS files (date-based)')
    p.add_argument('--input_dir', required=True,
                   help='directory of monthly GLORYS files (YYYY_MM.nc)')
    p.add_argument('--output_dir', required=True,
                   help='CROCO_FILES dir (crocotools_param.py + croco_grd.nc)')
    p.add_argument('--date', required=True, help='initial-condition date "YYYY-MM-DD"')
    p.add_argument('--Yorig', required=True, type=int)
    def _make_ini_hindcast(a):
        sys.path.append(a.output_dir)
        import crocotools_param as params
        ini_date = datetime.strptime(a.date, '%Y-%m-%d')
        # +/- 1 day so an ini on the 1st (or last) of a month also reads the
        # neighbour month when it is on disk (records on both sides of the date)
        files = _glorys_months(a.input_dir,
                               ini_date - timedelta(days=1),
                               ini_date + timedelta(days=1))
        if not files:
            raise SystemExit('No GLORYS file for %s in %s' % (a.date, a.input_dir))
        fname_in = files if len(files) > 1 else files[0]
        fname_out = params.ini_prefix + ini_date.strftime('_Y%YM%mD%d.nc')
        make_ini(fname_in, a.output_dir, ini_date, a.Yorig, fname_out)
    p.set_defaults(func=_make_ini_hindcast)

    # ---- make_bry_hindcast : GLORYS boundaries (date-based, cross-month) -
    p = sub.add_parser('make_bry_hindcast',
        help='Build CROCO boundaries from monthly GLORYS files (date window)')
    p.add_argument('--input_dir', required=True,
                   help='directory of monthly GLORYS files (YYYY_MM.nc)')
    p.add_argument('--output_dir', required=True,
                   help='CROCO_FILES dir (crocotools_param.py + croco_grd.nc)')
    p.add_argument('--start_date', required=True, help='"YYYY-MM-DD"')
    p.add_argument('--end_date', required=True, help='"YYYY-MM-DD"')
    p.add_argument('--Yorig', required=True, type=int)
    def _make_bry_hindcast(a):
        sys.path.append(a.output_dir)
        import crocotools_param as params
        start_date = datetime.strptime(a.start_date, '%Y-%m-%d')
        end_date   = datetime.strptime(a.end_date,   '%Y-%m-%d')
        files = _glorys_months(a.input_dir,
                               start_date - timedelta(days=1),
                               end_date   + timedelta(days=1))
        if not files:
            raise SystemExit('No GLORYS files for %s..%s' % (a.start_date, a.end_date))
        fname_in = files if len(files) > 1 else files[0]
        tag = start_date.strftime('_Y%YM%mD%d') + end_date.strftime('_to_Y%YM%mD%d.nc')
        fname_out = params.bry_prefix + tag
        make_bry(fname_in, a.output_dir, start_date, end_date, a.Yorig, fname_out)
    p.set_defaults(func=_make_bry_hindcast)
    # =====================================================================
    # RIVERS from GloFAS (hindcast) -- replaces a river discharge
    # climatology by the real daily discharge of the period
    # =====================================================================
    p = sub.add_parser('download_rivers_hindcast',
        help='Download GloFAS daily river discharge (EWDS) -> glofas_YYYY_MM.nc')
    p.add_argument('--grid', required=True, help='croco_grd.nc (sets the download box)')
    p.add_argument('--month_start', required=True, help='"YYYY-MM"')
    p.add_argument('--month_end', required=True, help='"YYYY-MM"')
    p.add_argument('--outputDir', required=True, help='where glofas_YYYY_MM.nc are written')
    p.add_argument('--pad', type=float, default=0.5, help='deg added around the grid box')
    def _download_rivers_hindcast(a):
        from download.glofas_rivers import download_glofas_monthly
        download_glofas_monthly(a.grid, a.month_start, a.month_end, a.outputDir, pad=a.pad,
                                n_parallel=int(os.environ.get('GLOFAS_N_PARALLEL', '2')))
    p.set_defaults(func=_download_rivers_hindcast)

    p = sub.add_parser('make_rivers_hindcast',
        help='GloFAS -> <river>.txt + river_list.txt (input of make_river_run.py)')
    p.add_argument('--grid', required=True, help='croco_grd.nc')
    p.add_argument('--glofas_dir', required=True, help='dir of glofas_YYYY_MM.nc')
    p.add_argument('--outputDir', required=True, help='dir for the river .txt files')
    p.add_argument('--month_start', required=True, help='"YYYY-MM"')
    p.add_argument('--month_end', required=True, help='"YYYY-MM"')
    p.add_argument('--dai_file', default=os.path.join(
        os.environ.get('CROCO_DATA_ROOT', ''), 'DATASETS_CROCOTOOLS', 'RUNOFF_DAI',
        'Dai_Trenberth_runoff_global_clim.nc'),
        help='river-mouth positions/names (Dai & Trenberth; its discharge is not used)')
    p.add_argument('--qmin', type=float, default=100.0, help='min mean GloFAS discharge (m3/s)')
    p.add_argument('--margin', type=float, default=0.0,
                   help='accept mouths up to this far outside the grid (deg); default 0 = inside only')
    p.add_argument('--radius', type=float, default=0.25, help='search radius around a mouth (deg)')
    p.add_argument('--extra_list', default=None, help='optional file: "name lon lat" per line')
    def _make_rivers_hindcast(a):
        from download.glofas_rivers import make_glofas_river_files
        n = make_glofas_river_files(a.grid, a.glofas_dir, a.outputDir, a.month_start, a.month_end,
                                    a.dai_file, qmin=a.qmin, margin=a.margin, radius=a.radius,
                                    extra_list=a.extra_list)
        if n == 0:
            raise SystemExit('no river selected -- lower --qmin or add --extra_list')
    p.set_defaults(func=_make_rivers_hindcast)

    args = parser.parse_args()
    if hasattr(args, 'func'):
        args.func(args)
    else:
        parser.print_help()


if __name__ == '__main__':
    main()
