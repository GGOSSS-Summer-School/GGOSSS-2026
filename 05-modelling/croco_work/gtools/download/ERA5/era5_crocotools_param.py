era5_dir_raw = '/home/${USER}/ggosss26/hindcast/scratch/Canary_12/downloaded_data/ERA5/raw'
era5_dir_processed = '/home/${USER}/ggosss26/hindcast/scratch/Canary_12/downloaded_data/ERA5/for_croco'
wave_extract = False
pressure_extract = True
year_start = 2026
month_start = 1
year_end = 2026
month_end = 1
Yorig = 1993
n_overlap = 0
times = ['00:00','01:00','02:00','03:00','04:00','05:00','06:00','07:00','08:00','09:00','10:00','11:00','12:00','13:00','14:00','15:00','16:00','17:00','18:00','19:00','20:00','21:00','22:00','23:00']
lonmin = -22.0
lonmax = -15.5
latmin = 14.0
latmax = 24.0
cff_tp = 1000./3600.
cff_heat = 1./3600.
variables = ['lsm','sst','tp','strd','ssr','t2m','q','u10','v10']
conv_cff  = [1.,1.,cff_tp,cff_heat,cff_heat,1.,1.,1.,1.]
units     = ['(0-1)','K','kg m-2 s-1','W m-2','W m-2','K','kg kg-1','m s-1','m s-1']
if pressure_extract:
    variables.append('msl'); conv_cff.append(1.); units.append('Pa')
