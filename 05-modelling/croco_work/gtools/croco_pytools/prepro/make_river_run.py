import argparse, os, sys
import numpy as np
import pylab as plt
import cftime
from dateutil.relativedelta import relativedelta
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as py
import netCDF4 as netcdf
sys.path.append("./Modules/")
import croco_class as Croco
import tools_make_river as riv_tools

ap = argparse.ArgumentParser(description='Build CROCO river runoff (non-interactive)')
ap.add_argument('--input_dir',  required=True)
ap.add_argument('--input_file', default='river_list.txt')
ap.add_argument('--croco_dir',  required=True)
ap.add_argument('--croco_grd',  default='croco_grd.nc')
ap.add_argument('--river_filename', default='croco_runoff.nc')
ap.add_argument('--Yorig',  type=int, required=True)
ap.add_argument('--Ystart', type=int, required=True)
ap.add_argument('--Mstart', type=int, required=True)
ap.add_argument('--Yend',   type=int, required=True)
ap.add_argument('--Mend',   type=int, required=True)
ap.add_argument('--rivers_cyl', type=float, default=365.25)
ap.add_argument('--file_format', default='FULL', choices=['MONTHLY','YEARLY','FULL'])
ap.add_argument('--output_frequency', default='DAILY')
a = ap.parse_args()

Yorig = a.Yorig
Ystart, Mstart = a.Ystart, a.Mstart
Yend, Mend = a.Yend, a.Mend
input_dir = a.input_dir if a.input_dir.endswith('/') else a.input_dir + '/'
input_file = a.input_file
croco_dir = a.croco_dir if a.croco_dir.endswith('/') else a.croco_dir + '/'
croco_grd = a.croco_grd
river_filename = a.river_filename
river_file_format = a.file_format
rivers_cyl = a.rivers_cyl
rivers_output_frequency = a.output_frequency

add_ts = False
time_units = 'days since %i-01-01' % Yorig
rstr = plt.datetime.datetime(Ystart, Mstart, 1, 0, 0, 0)
rend = plt.datetime.datetime(Yend, Mend, 1, 12, 0, 0) + relativedelta(months=1, days=-1)

crocogrd = Croco.CROCO_grd(''.join((croco_dir, croco_grd)))

data = np.genfromtxt(input_dir + input_file, dtype=str, comments='#')
if data.ndim == 1:
    data = data.reshape((1, 3))
list_river_files = np.array([x[0] for x in data])
lon_river = np.array([]); lat_river = np.array([])
for x in data:
    try: lon_river = np.append(lon_river, float(x[1]))
    except: lon_river = np.append(lon_river, np.nan)
    try: lat_river = np.append(lat_river, float(x[2]))
    except: lat_river = np.append(lat_river, np.nan)

list_river_files = np.array([input_dir + f for f in list_river_files])

river_obs = riv_tools.read_river(list_river_files, lon_river, lat_river, rstr, rend, time_units)
river_obs = riv_tools.get_river_index(river_obs, crocogrd)
river_ext = riv_tools.fill_period(river_obs, rstr, rend, time_units, rivers_output_frequency)
river_ext = riv_tools.locateji_croco(river_ext, grd=crocogrd, graph=False)

startloc = plt.datetime.datetime(Ystart, Mstart, 1)
endloc = plt.num2date(plt.date2num(rend)).replace(tzinfo=None)

loc_time = river_ext[list(river_ext.keys())[0]]['time']
ind = np.where((loc_time >= cftime.date2num(startloc, time_units)) &
               (loc_time <= cftime.date2num(endloc, time_units)))
dtmin, dtmax = np.min(ind), np.max(ind)

river_outname = croco_dir + river_filename
Croco.CROCO.create_river_nc(None, river_outname, crocogrd, len(river_ext.keys()), add_ts)

nc = netcdf.Dataset(river_outname, 'a')
nc.variables['qbar_time'][:] = river_ext[list(river_ext.keys())[0]]['time'][dtmin:dtmax+1]
nc.variables['qbar_time'].units = time_units
if rivers_cyl > 0:
    nc.variables['qbar_time'].cycle_length = rivers_cyl
for ir, r in enumerate(river_ext):
    nc.variables['runoff_name'][ir] = r
    nc['Qbar'][ir] = river_ext[r]['flow'][dtmin:dtmax+1]
nc.close()

riv_tools.write_croco_in(river_ext, croco_dir, river_filename)

print('make_river_run: wrote %s (cycle=%.2f, %d rivers) + for_croco_in.txt'
      % (river_outname, rivers_cyl, len(river_ext.keys())))
