import calendar
import cdsapi
import datetime
import os
'''
download_glofas_africa_clim.py
ONE-TIME reusable Africa-wide GloFAS river-discharge climatology.
Monthly (resumable), reads EWDS creds from ~/.ewdsapirc.
'''
# ---- USER DATA ----
Ystart, Mstart = 2010, 1
Yend,   Mend   = 2025, 12
area = [40.0, -40.0, -45.0, 40.0]
product = 'cems-glofas-historical'
output_dir = os.path.join(os.environ.get('CROCO_DATA_ROOT', os.path.expanduser('~/croco_work/data')),
                          'DATASETS_CROCOTOOLS', 'GLOFAS_AFRICA')
output_name = 'glofas_africa'
# ---- END USER DATA ----

if not os.path.exists(output_dir):
    os.makedirs(output_dir)

_ewds = {}
with open(os.path.expanduser('~/.ewdsapirc')) as _f:
    for _line in _f:
        if ':' in _line:
            _k, _v = _line.split(':', 1)
            _ewds[_k.strip()] = _v.strip()
c = cdsapi.Client(url=_ewds['url'], key=_ewds['key'])

for Y in range(Ystart, Yend+1):
    mo_min = Mstart if Y == Ystart else 1
    mo_max = Mend if Y == Yend else 12
    for M in range(mo_min, mo_max+1):
        outfile_nc = f"{output_dir}/{output_name}_Y{Y}_M{M:02d}.nc"
        if os.path.exists(outfile_nc):
            print(f"skip (exists): {outfile_nc}")
            continue
        request = {
            'system_version': 'version_4_0',
            'hydrological_model': 'lisflood',
            'product_type': 'consolidated',
            'variable': 'average_river_discharge_in_the_last_24_hours',
            'timespan': 'time_mean',        # mandatory since the Jul-2026 EWDS change
            'year': str(Y),
            'month': '%02d' % M,
            'day': ['%02d' % d for d in range(1, calendar.monthrange(Y, M)[1] + 1)],
            'area': area,
            'data_format': 'netcdf',
            'download_format': 'unarchived',
        }
        print('-----------------------------------------------------------')
        print(datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"), f' GLOFAS {Y}-{M:02d}')
        print('-----------------------------------------------------------')
        c.retrieve(product, request, outfile_nc)
        print(f'  -> {outfile_nc}')

print('')
print('ALL MONTHS DOWNLOADED. Next: run build_africa_clim.sh')
