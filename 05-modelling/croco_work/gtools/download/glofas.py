import calendar
import cdsapi
import datetime
import os
ENV_CP = os.path.dirname(__file__)
ENV_MOD = os.path.join(ENV_CP, '..', 'croco_pytools', 'prepro', 'Modules')   # vendored legacy croco_pytools
import sys
sys.path.append(ENV_MOD)
import croco_class as Croco

'''
glofas.py  (legacy stand-alone GloFAS download; for the hindcast use
  'ggosss26.py download_rivers_hindcast' / glofas_rivers.py instead)
========================

Download monthly grib files of daily river discharges from CMEMS-GLOFAS
https://cds.climate.copernicus.eu/cdsapp#!/dataset/cems-glofas-historical?tab=overview
Tranform each grib file into a netcdf file

Require ECMWF CDS API for data request
https://cds.climate.copernicus.eu/api-how-to#
Require CDO for grib to netcdf conversions

Once the files are dowloaded a mean file should be computed with ncra

Pierrick Penven 17/11/2022
'''


#
################################################################################
################################################################################
#
#
#
# Main Program
# 
# 
#
################################################################################
###########################  USER DATA  ########################################
################################################################################
#
Ystart,Mstart = 2023, 1   # Starting month
Yend,Mend  = 2023, 12      # Ending month 

product = 'cems-glofas-historical'

croco_dir = './'
croco_grd = 'croco_grd.nc'

###
convert2netcdf = True
output_dir = './DATA_RIVER/'
output_name = 'cems_glofas'

#
################################################################################
###########################  END OF USER DATA  #################################
################################################################################
#

if not os.path.exists(output_dir):
    os.makedirs(output_dir)

grd = Croco.CROCO_grd(f"{croco_dir}/{croco_grd}", dict(theta_s=7, theta_b=2, N=50, hc=200))
area = [grd.latmax()+0.1, grd.lonmin()-0.1, grd.latmin()-0.1, grd.lonmax()+0.1]

months=['january','february','march','april',\
        'may','june','july','august',\
	'september','october','november','december']

#
# Start the CDS API
#

_ewds = {}
with open(os.path.expanduser('~/.ewdsapirc')) as _f:
    for _line in _f:
        if ':' in _line:
            _k, _v = _line.split(':', 1)
            _ewds[_k.strip()] = _v.strip()
c = cdsapi.Client(url=_ewds['url'], key=_ewds['key'])

#
# Loop on the years and months
#

for Y in range(Ystart,Yend+1):

  if Y==Ystart: 
    mo_min=Mstart
  else:
    mo_min=1

  if Y==Yend:
    mo_max=Mend
  else:
    mo_max=12

  for M in range(mo_min,mo_max+1):

    outfile_grib= f"{output_dir}/{output_name}_Y{Y}_M{M:02d}.grb"
    outfile_nc  = f"{output_dir}/{output_name}_Y{Y}_M{M:02d}.nc"

#
# Define the request
#
    options = {
        'variable': 'average_river_discharge_in_the_last_24_hours',
        'timespan': 'time_mean',            # mandatory since the Jul-2026 EWDS change
        'data_format': 'grib2',
        'download_format': 'unarchived',
        'system_version': 'version_4_0',
        'hydrological_model': 'lisflood',
        'product_type': 'consolidated',
        'year': str(Y),
        'month': '%02d' % M,
        'day': ['%02d' % d for d in range(1, calendar.monthrange(Y, M)[1] + 1)],
        'data_format': 'grib2',
        'download_format': 'unarchived',
        'area': area,}

#
# Printing message on screen
#

    info_time_clock = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    print('                                                           ')
    print('-----------------------------------------------------------')
    print(''+info_time_clock)
    print(' GLOFAS data request, please wait...                       ')
    print('-----------------------------------------------------------')
    print('Request options: ')
    print(options)

#
# Send the request
#

    c.retrieve(product,options,outfile_grib)

#
# Convert grib files into Netcdf files
#
    if convert2netcdf:
        cmd = 'cdo -f nc copy ' + outfile_grib + ' ' + outfile_nc
        os.system(cmd)
#
# Done
#
if convert2netcdf:
    merge_cmd = \
        f"ncrcat {output_dir}/{output_name}*.nc {output_dir}/{output_name}.nc"
    os.system(merge_cmd)
