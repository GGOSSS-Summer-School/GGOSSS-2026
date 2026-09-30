#!/usr/bin/env python

# Script to download ECMWF ERA5 reanalysis datasets from the Climate Data
#  Store (CDS) of Copernicus https://cds.climate.copernicus.eu
#
#  This script use the CDS Phyton API[*] to connect and download specific ERA5 
#  variables, for a chosen area and monthly date interval, required by CROCO to 
#  perform simulations with atmospheric forcing. Furthermore, this script use 
#  ERA5 parameter names and not parameter IDs as these did not result in stable 
#  downloads. 
#
#  Tested using Python 3.8.6 and Python 3.9.1. This script need the following
#  python libraries pre-installed: "calendar", "datetime", "json" and "os".
#
#  [*] https://cds.climate.copernicus.eu/how-to-api
#
#  Copyright (c) DDONOSO February 2021
#  e-mail:ddonoso@dgeo.udec.cl  
#

#  You may see all available ERA5 variables at the following website
#  https://confluence.ecmwf.int/display/CKB/ERA5%3A+data+documentation#ERA5:datadocumentation-Parameterlistings

# -------------------------------------------------
# Getting libraries and utilities
# -------------------------------------------------
import cdsapi
from ERA5_utilities import *
import calendar
import datetime
import json
import os

# -------------------------------------------------
# Import my crocotools_param_python file
from era5_crocotools_param import *
print('year_start is '+str(year_start))

# -------------------------------------------------
dl=2

lonmin = str(float(lonmin)-dl)
lonmax = str(float(lonmax)+dl)
latmin = str(float(latmin)-dl)
latmax = str(float(latmax)+dl)
print ('lonmin-dl = ', lonmin)
print ('lonmax+dl =', lonmax)
print ('latmin-dl =', latmin)
print ('latmax+dl =', latmax)
# -------------------------------------------------

area = [latmax, lonmin, latmin, lonmax]

# -------------------------------------------------
# Setting raw output directory
# -------------------------------------------------
# Get the current directory
os.makedirs(era5_dir_raw,exist_ok=True)

# -------------------------------------------------
# Loading ERA5 variables's information as 
# python Dictionary from JSON file
# -------------------------------------------------
with open('ERA5_variables.json', 'r') as jf:
    era5 = json.load(jf)

# -------------------------------------------------
# Build the list of requests (one per month x variable)
# -------------------------------------------------
# CDS queues every request for 1-3 min before the (small) transfer, so the
# requests are SUBMITTED IN PARALLEL: the waits overlap instead of adding up.
# Request contents and file names are unchanged (the converter needs them).
import concurrent.futures
import threading
import time as _time

# CDS caps queued requests per user and dataset (a job over the cap is rejected
# with 'Number queued requests for this dataset is temporarily limited'), so keep
# this small; rejected jobs are retried with a growing wait.
N_PARALLEL = int(os.environ.get('ERA5_N_PARALLEL', '2'))   # concurrent CDS requests
MAX_TRIES = 8

monthly_date_start = datetime.datetime(year_start,month_start,1)
monthly_date_end = datetime.datetime(year_end,month_end,1)
len_monthly_dates = (monthly_date_end.year - monthly_date_start.year) * 12 + \
                    (monthly_date_end.month - monthly_date_start.month) + 1

tasks = []
monthly_date = monthly_date_start
for j in range(len_monthly_dates):
    year = monthly_date.year
    month = monthly_date.month
    days_in_month = calendar.monthrange(year,month)[1]
    # ERA5 is released ~5 days behind real time: in the current month only ask
    # for the days that exist (a request for missing days fails as a whole)
    last_avail = datetime.date.today() - datetime.timedelta(days=int(os.environ.get('ERA5_LAG_DAYS', '6')))
    days = [f"{day:02}" for day in range(1, days_in_month + 1)
            if datetime.date(year, month, day) <= last_avail]
    if not days:
        print(' %04d-%02d: no ERA5 data released yet -- skipped' % (year, month))
        monthly_date = addmonths4date(monthly_date,1)
        continue
    if len(days) < days_in_month:
        print(' %04d-%02d: ERA5 only up to day %s (release delay) -- partial month' % (year, month, days[-1]))

    for vname in variables:
        vlong = era5[vname][0]
        options = {
             'product_type': ['reanalysis'],
             'variable': [vlong],
             'year': [str(year)],
             'month': [str(month)],
             'day': days,
             'data_format': 'netcdf',
             'download_format': 'unarchived',
             'area': area,
             }
        # variables without "diurnal variations": one record per day
        if vlong in ('sea_surface_temperature', 'land_sea_mask'):
            options['time'] = ['00:00']
        else:
            options['time'] = times
        # specific / relative humidity come from the pressure-levels product (1000 hPa)
        if vlong in ('specific_humidity', 'relative_humidity'):
            options['pressure_level'] = ['1000']
            product = 'reanalysis-era5-pressure-levels'
        else:
            product = 'reanalysis-era5-single-levels'

        fname = 'ERA5_ecmwf_' + vname.upper() + '_Y' + str(year) + 'M' + str(month).zfill(2) + '.nc'
        output = era5_dir_raw + '/' + fname
        # re-running is safe: skip raw files already downloaded
        if os.path.isfile(output) and os.path.getsize(output) > 0:
            print(' already downloaded, skipping: ' + fname)
            continue
        tasks.append((product, options, output, fname))

    monthly_date = addmonths4date(monthly_date,1)

# -------------------------------------------------
# Run the requests
# -------------------------------------------------
_print_lock = threading.Lock()
def _log(msg):
    with _print_lock:
        print(datetime.datetime.now().strftime('%H:%M:%S') + '  ' + msg, flush=True)

def _retrieve(task):
    product, options, output, fname = task
    for attempt in range(1, MAX_TRIES + 1):
        try:
            _log('request  ' + fname + ('' if attempt == 1 else '  (attempt %d)' % attempt))
            c = cdsapi.Client(quiet=True, progress=False)
            # download to a temp name so an interrupted transfer is never
            # mistaken for a complete file on the next run
            c.retrieve(product, options).download(output + '.part')
            os.replace(output + '.part', output)
            _log('done     ' + fname)
            return fname
        except Exception as e:
            _log('FAILED   ' + fname + ': ' + str(e).splitlines()[0][:200])
            if attempt == MAX_TRIES:
                raise
            _time.sleep(min(60 * attempt, 300))

print(' %d ERA5 request(s) to do, %d in parallel, area %s' % (len(tasks), N_PARALLEL, area))
failed = []
with concurrent.futures.ThreadPoolExecutor(max_workers=max(1, N_PARALLEL)) as pool:
    futures = {pool.submit(_retrieve, t): t[3] for t in tasks}
    for fut in concurrent.futures.as_completed(futures):
        if fut.exception() is not None:
            failed.append(futures[fut])

if failed:
    raise SystemExit(' ERA5 requests failed: ' + ', '.join(sorted(failed)) +
                     ' -- re-run the same command (finished files are skipped)')

print('                                               ')
print(' ERA5 data request has been done successfully! ')
print('                                               ')
