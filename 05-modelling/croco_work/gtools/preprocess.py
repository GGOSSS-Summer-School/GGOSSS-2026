"""
CROCO forecast preprocessing.
Adapted, with permission, from the SAEON somisana project
(somisana-croco / crocotools_py/preprocess.py). Trimmed to: reformat_gfs_atm,
make_ini, make_bry, make_tides. See ACKNOWLEDGEMENTS.md.
"""
import xarray as xr
from datetime import datetime, timedelta
import calendar
import cftime
import numpy as np
import pandas as pd
import os, sys, glob
from os import path
from scipy.interpolate import griddata, interp1d
# croco_pytools is vendored in this repo, right next to this file
_CPT = os.path.join(os.path.dirname(__file__), "croco_pytools")
sys.path.append(os.path.join(_CPT, "prepro"))            # so `Modules.xxx` resolves as a package
sys.path.append(os.path.join(_CPT, "prepro", "Modules"))
sys.path.append(os.path.join(_CPT, "prepro", "Readers"))
import Cgrid_transformation_tools as grd_tools
import interp_tools
import sigmagrid_tools as sig_tools
import croco_class as Croco
import ibc_class as Inp
import tides_class as Inp_tides
from matplotlib.dates import date2num, num2date
import netCDF4 as netcdf



# ============================================================
# make_tides
# ============================================================
def make_tides(input_dir,output_dir,run_ini_date,Yorig,fname_out):
    '''
    This function is mostly taken from croco_pytools/prepro/make_tides.py
    We just adjust how the inputs are handled
    Most of the inputs are contained in a crocotools_param.py file, while
    a few others are direct inputs to this function. The direct inputs are for 
    things which we may want to be update as part of an operational/inter-annual workflow

    input_dir - Path to directory containing the raw tidal files from e.g. TPXO
    output_dir - Path to where the forcing file will be saved. This directory also needs a crocotools_param.py file
    run_ini_date - datetime.datetime object representing the initial time for the run
    Yorig - the Yorig value used in setting up the CROCO model
    fname_out - output filename (only the filename, not the full path)
    '''
    
    # file names below are built as input_dir+name: make sure of the trailing '/'
    input_dir = os.path.join(input_dir, '')

    # read in variables from the crocotools_param.py file
    sys.path.append(output_dir)
    import crocotools_param as params
    inputdata=params.inputdata
    input_type=params.input_type
    multi_files=params.multi_files
    waves_separated=params.waves_separated
    elev_file=params.elev_file
    u_file=params.u_file
    v_file=params.v_file
    croco_grd = os.path.join(output_dir, params.croco_grd)
    tides=params.tides
    cur=params.cur
    pot=params.pot
    Correction_ssh =params.Correction_ssh
    Correction_uv = params.Correction_uv
    
    fname_out = os.path.join(output_dir,fname_out)
    
    #read tides and periods
    tides_txt_file=os.path.join(_CPT, "prepro", "Modules", "tides.txt")
    tides_param=np.loadtxt(tides_txt_file,skiprows=3,comments='^',usecols=(0,4),dtype=str)
    tides_names=tides_param[:,0].tolist()
    tides_periods=tides_param[:,1].tolist()
    tides_names=[x.lower() for x in tides_names]
    tides_names=np.array(tides_names)

    # Load croco_grd
    sigma_params = dict(theta_s=0, theta_b=0, N=1, hc=1) # not relevant as 2d
    crocogrd = Croco.CROCO_grd(croco_grd, sigma_params)

    # --- Load input (restricted to croco_grd) ----------------------------
    if multi_files:
        if waves_separated:
            input_file_ssh=[]
            input_file_u=[]
            input_file_v=[]
            for inp in tides:
                if path.isfile(input_dir+elev_file.replace('<tides>',inp)):
                    input_file_ssh+=[input_dir+elev_file.replace('<tides>',inp)]
                elif path.isfile(input_dir+elev_file.replace('<tides>',inp.lower())):
                    input_file_ssh+=[input_dir+elev_file.replace('<tides>',inp.lower())]
                else:
                    sys.exit('Elevation file %s for wave %s is missing' % (input_dir+elev_file.replace('<tides>',inp), inp))
               
                if cur:
                    if path.isfile(input_dir+u_file.replace('<tides>',inp)):
                        input_file_u+=[input_dir+u_file.replace('<tides>',inp)]
                    elif path.isfile(input_dir+u_file.replace('<tides>',inp.lower())):
                        input_file_u+=[input_dir+u_file.replace('<tides>',inp.lower())]
                    else:
                        sys.exit('Eastward current file for wave %s is missing' % inp)

                    if path.isfile(input_dir+v_file.replace('<tides>',inp)):
                        input_file_v+=[input_dir+v_file.replace('<tides>',inp)]
                    elif path.isfile(input_dir+v_file.replace('<tides>',inp.lower())):
                        input_file_v+=[input_dir+v_file.replace('<tides>',inp.lower())]
                    else:
                        sys.exit('Northward current file for wave %s is missing' % inp)
            
            input_file_ssh = list(input_file_ssh)
            if cur:
                input_file_u = list(input_file_u)
                input_file_v = list(input_file_v)
            else: 
                input_file_u = None
                input_file_v = None
        else:
            input_file_ssh=list(input_dir+elev_file)
            if cur:
                input_file_u = list(input_dir+u_file)
                input_file_v = list(input_dir+v_file)
            else:
                input_file_u = None
                input_file_v = None
    else:
        input_file_ssh=list([input_dir+params.input_file])
        if cur:
            input_file_u=list([input_dir+params.input_file])
            input_file_v=list([input_dir+params.input_file])
        else:
            input_file_u=None
            input_file_v=None

    inpdat=Inp_tides.getdata(inputdata,input_file_ssh,crocogrd,input_type,tides,input_file_u,input_file_v)

    # --- Create the initial file -----------------------------------------

    Croco.CROCO.create_tide_nc(None,fname_out,crocogrd,cur=cur,pot=pot)

    if Correction_ssh or Correction_uv:
        date=cftime.datetime(run_ini_date.year,run_ini_date.month,run_ini_date.day)
        date_orig=cftime.datetime(Yorig,1,1)

    todo=['H']

    if cur:
        todo+=['cur']
    if pot:
        todo+=['pot']
        coslat2=np.cos(np.deg2rad(crocogrd.lat))**2
        sin2lat=np.sin(2.*np.deg2rad(crocogrd.lat))

    # --- Start loop on waves --------------------------------------------
    nc=netcdf.Dataset(fname_out, 'a')

    if Correction_ssh or Correction_uv:
        nc.Nodal_Correction=''.join(('Origin time is ',str(date_orig)))
    else:
        nc.Nodal_Correction='No nodal correction'

    for i,tide in enumerate(tides) :
        print('\nProcessing *%s* wave' %(tide))
        print('-----------------------')
        index=np.argwhere(tides_names==tide.lower())
        if (len(index)>0):
            print("  tides %s is in the list"%(tide))
            # get period
            index=index[0][0]
            period=float(tides_periods[index])
            print("  Period of the wave %s is %f"%(tide,period))
            
            nc.variables['tide_period'][i]=period
        if multi_files:
            # In this case waves had been concatenated in the order they appear in the list
            tndx=i
        else:
            # read ntime/periods dimension and find the closest wave
            # relative tolerance: atlases store periods in float32 with few
            # digits (TPXO7: Mf=327.8599 h, Mm=661.31 h vs tides.txt
            # 327.858980 / 661.309208), which an absolute 1e-4 h test rejects
            tndx=np.argwhere(np.isclose(np.asarray(inpdat.ntides, dtype=float), period,
                                        rtol=1e-5, atol=1e-4))
            if len(tndx)==0:
                sys.exit('  Did not find wave %s in input file' % tide)
            else:
                tndx=tndx[0]

        [pf,pu,mkB]=inpdat.egbert_correction(tide,date)
        # For phase shift time should be in seconds relatively Jan 1 1992
        # As mkB is the wave phase at this date
        t0 = cftime.date2num(date_orig,'seconds since 1992-01-01:00:00:00')
        if Correction_ssh or Correction_uv:      
            correc_amp   = pf
            correc_phase = mkB+np.deg2rad(t0/(period*10)) +         pu
            #              |--- phase at origin time ---|  |nodal cor at ini time|
        else:
            correc_amp=1
            correc_phase= mkB+np.deg2rad(t0/(period*10))

        # --- Start loop on var ------------------------------------------
        for vars in todo:
            # get data
            if vars == 'H':
                print('\n  Processing tidal elevation')
                print('  -------------------------')
                (tide_complex,NzGood) = interp_tools.interp_tides(inpdat,vars,-1,crocogrd,tndx,tndx,input_type)
                if Correction_ssh:
                    tide_amp=np.ma.abs(tide_complex)*correc_amp
                    if 'tpxo' in inputdata :
                        tide_phase=np.mod(np.ma.angle(tide_complex)*-180/np.pi-correc_phase*180/np.pi,360)
                    else:
                        tide_phase=np.mod(np.ma.angle(tide_complex)*180./np.pi-correc_phase*180/np.pi,360)
                else:
                    tide_amp=np.ma.abs(tide_complex)
                    if 'tpxo' in inputdata :
                        tide_phase=np.mod(np.ma.angle(tide_complex)*-180/np.pi,360)
                    else:
                        tide_phase=np.mod(np.ma.angle(tide_complex)*180./np.pi,360)  

                nc.variables['tide_Ephase'][i,:]=tide_phase*crocogrd.maskr
                nc.variables['tide_Eamp'][i,:]=tide_amp*crocogrd.maskr

            #########################
            elif vars == 'cur':
                print('\n  Processing tidal currents')
                print('  -------------------------')
                (u_tide_complex,v_tide_complex,NzGood) = interp_tools.interp_tides(inpdat,vars,-1,crocogrd,tndx,tndx,input_type)
                
                if Correction_uv:
                    u_tide_amp=np.ma.abs(u_tide_complex)*correc_amp
                    v_tide_amp=np.ma.abs(v_tide_complex)*correc_amp
     
                    if 'tpxo' in inputdata:
                        u_tide_phase=np.mod(np.ma.angle(u_tide_complex)*-180/np.pi-correc_phase*180/np.pi,360)
                        v_tide_phase=np.mod(np.ma.angle(v_tide_complex)*-180/np.pi-correc_phase*180/np.pi,360)
                    else:
                        u_tide_phase=np.mod(np.ma.angle(u_tide_complex)*180./np.pi-correc_phase*180/np.pi,360)
                        v_tide_phase=np.mod(np.ma.angle(v_tide_complex)*180./np.pi-correc_phase*180/np.pi,360)
                else:
                    u_tide_amp=np.ma.abs(u_tide_complex)
                    v_tide_amp=np.ma.abs(v_tide_complex)

                    if 'tpxo' in inputdata:
                        u_tide_phase=np.mod(np.ma.angle(u_tide_complex)*-180/np.pi,360)
                        v_tide_phase=np.mod(np.ma.angle(v_tide_complex)*-180/np.pi,360)
                    else:
                        u_tide_phase=np.mod(np.ma.angle(u_tide_complex)*180./np.pi,360)
                        v_tide_phase=np.mod(np.ma.angle(v_tide_complex)*180./np.pi,360)


                major,eccentricity,inclination,phase=inpdat.ap2ep(u_tide_amp,u_tide_phase,v_tide_amp,v_tide_phase)
      
                nc.variables['tide_Cmin'][i,:,:]=major[:,:]*eccentricity[:,:]*crocogrd.maskr
                nc.variables['tide_Cmax'][i,:,:]=major[:,:]*crocogrd.maskr
                nc.variables['tide_Cangle'][i,:,:]=inclination[:,:]*crocogrd.maskr
                nc.variables['tide_Cphase'][i,:,:]=phase[:,:]*crocogrd.maskr
            #########################
            elif vars == 'pot':
                print('\n  Processing equilibrium tidal potential')
                print('  --------------------------------------')
                try:
                    coef=eval(''.join(('inpdat.pot_tide.',tide.lower())))
                except:
                    try:
                        # some waves start with a number (ex: 2N2) and python do not like it
                        coef=eval(''.join(('inpdat.pot_tide._',tide.lower())))
                    except:
                        print('No potential prameter defined for wave %s' % tide)
                        coef=[1 ,0]

                if period<13:  # semidiurnal
                    Pamp=correc_amp*coef[0]*coef[1]*coslat2
                    Ppha=np.mod(-2*crocogrd.lon-correc_phase*180/np.pi,360)
                elif period<26: # diurnal
                    Pamp=correc_amp*coef[0]*coef[1]*sin2lat;
                    Ppha=np.mod(-crocogrd.lon-correc_phase*180/np.pi,360)
                else: # long-term
                    Pamp=correc_amp*coef[0]*coef[1]*(1-1.5*coslat2);
                    Ppha=np.mod(-correc_phase*180/np.pi,360.0)
     
                nc.variables['tide_Pamp'][i,:,:]   = Pamp*crocogrd.maskr
                nc.variables['tide_Pphase'][i,:,:] = Ppha*crocogrd.maskr

    nc.close()


# ============================================================
# reformat_gfs_atm
# ============================================================
def reformat_gfs_atm(gfs_dir,out_dir,Yorig):
    '''
    Convert the GFS atmospheric forecast grb files downloaded by the cli.py function download_gfs_atm
    and convert the data into nc files in a format which can be ingested by CROCO
    using the ONLINE cpp key for online interpolation of the surface forcing
    (we will use the default 'CFSR' file format)
    '''
    
    # Path to the directory containing your GRIB files
    gfs_files = os.path.join(gfs_dir,"*.grb")
    
    # List all GRIB files
    file_paths = sorted(glob.glob(gfs_files))
    
    def open_grib_file(file_path,var_dict):
        return xr.open_dataarray(
            file_path,
            engine='cfgrib',
            filter_by_keys={'shortName': var_dict['shortName'],
                            'stepType': var_dict['stepType'],
                            })

    # (you can see the vars in the files using e.g.
    # grib_ls -P shortName,typeOfLevel,level 2024080106_f001.grb)
    
    variables = {
            "Temperature_height_above_ground": {
                "shortName": "2t",
                "stepType": "instant",
            },
            "Specific_humidity": {
                "shortName": "2sh",
                "stepType": "instant",
            },
            "Precipitation_rate": {
                "shortName": "prate",
                "stepType": "instant",
            },
            "Downward_Short-Wave_Rad_Flux_surface": {
                "shortName": "sdswrf",
                "stepType": "avg",
            },
            "Upward_Short-Wave_Rad_Flux_surface": {
                "shortName": "suswrf",
                "stepType": "avg",
            },
            "Downward_Long-Wave_Rad_Flux": {
                "shortName": "sdlwrf",
                "stepType": "avg",
            },
            "Upward_Long-Wave_Rad_Flux_surface": {
                "shortName": "sulwrf",
                "stepType": "avg",
            },
            "U-component_of_wind": {
                "shortName": "10u",
                "stepType": "instant",
            },
            "V-component_of_wind": {
                "shortName": "10v",
                "stepType": "instant",
            },
            "patm": {
                "shortName": "prmsl",
                "stepType": "instant",
            },
        }
    
    for var in variables:
        
        # get an xarray dataset for this variable
        print('working on '+var)
        var_dict = variables[var]
        datasets = [open_grib_file(fp,var_dict) for fp in file_paths]
        da = xr.concat(datasets, dim='valid_time')
        
        if var_dict['stepType']=='avg':
            # A bunch of variables are (rather annoyingly) written out as the
            # accumulated average over each 6 hour forecast period.
            # We want to convert these accumulated averages into individual one-hour averages
            #
            # see FAQ "How can the individual one-hour averages be computed" from https://rda.ucar.edu/datasets/ds093.0/#docs/FAQs_6hrly.html
            # Excerpt from there:
            # You can compute the one-hour average (X) ending at hour N by using the N-hour average (a) and the (N-1)-hour average (b) as follows:
            # X = N*a - (N-1)*b
            # So if you want the 1-hour Average for the period initial+3 to initial+4 (X), you would use the 4-hour Average (initial+0 to initial+4) as (a) and the 3-hour Average (initial+0 to initial+3) as (b) as follows:
            # X = 4*a - 3*b
                
            # start by extracting the forecast hour for each time-step (can be derived from the 'step' variable)
            # frcst = xr.DataArray(da.step / np.timedelta64(1, 'h'), dims='valid_time')
            frcst = da.step.values / np.timedelta64(1, 'h')
            
            # the averaging period for these variables is (again rather annoyingly) reset every 6 hours, 
            # even if we are looking at forecast hours greater than 6
            # so here we compute the hour within in the 6 hour forecast cycle 
            # (frcst_ave will range from 1-6 by definition)
            frcst_ave = np.mod(frcst - 1, 6) + 1
            
            # get arrays for doing the calcs and updating the output array
            data = da.values
            data_output = data.copy()
            
            for i in range(len(frcst)):
                if frcst_ave[i]>1:
                    # only do the conversion for forecast hours greater than 1
                    if frcst[i]<=120:
                        data_output[i,::] = data[i,::]*frcst_ave[i]-data[i-1,::]*frcst_ave[i-1]
                    else:
                        # after 120 hrs the forecasts are provided at three hourly intervals
                        # so each of these represents the accumulated hourly averages over 3 and 6 hours
                        # so we have to handle this separately to get back to hourly averages
                        # 
                        if frcst_ave[i]==3:
                            data_output[i,::] = data[i,::]/3
                        else: # frcst_ave[i]==6
                            data_output[i,::] = (data[i,::]*frcst_ave[i]-data[i-1,::]*frcst_ave[i-1])/3
            # Create a new DataArray with the same dimensions and coordinates as the original
            # but with the updated output
            da = xr.DataArray(
                data_output,
                dims=da.dims,
                coords=da.coords,
                attrs=da.attrs
            )
        
        # rename the dimensions to match what is expected by CROCO ONLINE option
        da = da.drop_vars(['time','step'])
        da = da.rename({'longitude': 'lon', 'latitude': 'lat', 'valid_time': 'time'})

        # Guard against duplicate timestamps before interpolating: if gfs_dir
        # contains .grb files from more than one GFS cycle/run (e.g. stale
        # files left over from a previous day's download that weren't
        # cleaned out, or two forecast files that happen to share a valid
        # time), the concat above produces a 'time' axis with repeated
        # values. da.interp() then fails with "Reindexing only valid with
        # uniquely valued Index objects" -- sort and drop exact duplicates
        # (keep the first occurrence) so a stale/overlapping file degrades
        # gracefully instead of crashing the whole reformat.
        da = da.sortby('time')
        _, unique_idx = np.unique(da['time'].values, return_index=True)
        if len(unique_idx) != da.sizes['time']:
            n_dupes = da.sizes['time'] - len(unique_idx)
            print(f"  warning: {n_dupes} duplicate 'time' value(s) found for {var} "
                 f"(check {gfs_dir} for .grb files from more than one GFS cycle) "
                 f"-- keeping the first occurrence of each, dropping the rest.")
            da = da.isel(time=np.sort(unique_idx))

        # handle any 3 hourly time-steps at the end of the data 
        # by interpolating onto an hourly time axis (I'm not totally sure this is needed but no harm done)
        time_equidistant = pd.date_range(start=da.time.min().values, end=da.time.max().values, freq='h')
        da = da.interp(time=time_equidistant)
        
        # make a dataset from the dataarray
        ds = xr.Dataset({var: da})
        
        # we need to convert time to days since Yorig!
        # Reference date
        reference_date = np.datetime64(str(Yorig)+'-01-01T00:00:00')
        # time_in_ns = ds['time'].astype('datetime64[ns]')
        # Convert the time dimension to days since the reference date
        ds['time'] = (ds['time'].astype('datetime64[ns]') - reference_date) / np.timedelta64(1, 'D')
        # Set the units attribute for the time coordinate
        ds['time'].attrs['units'] = 'days since 1-Jan-'+str(Yorig)+' 00:00:00'
        
        # write the nc file
        # the ONLINE cppkey is designed for use with monly interannual simulations
        # where the year and month of the file name is appended to the end of the file
        # since we are using this option with forecasts we'll just put dummy values
        # for the year and month, and put these values in the *.in file making it 
        # something we don't have to handle separately
        fname_out = os.path.join(out_dir,var.upper()+"_Y9999M01.nc")
        ds.to_netcdf(fname_out)
        ds.close()
        
        # Delete all temporary files created during reading the .grb files (not sure what these are but they're not needed)
        for gbx9_file in glob.glob(os.path.join(gfs_dir, '*.gbx9')):
            os.remove(gbx9_file)
        for ncx_file in glob.glob(os.path.join(gfs_dir, '*.ncx')):
            os.remove(ncx_file)
        for idx_file in glob.glob(os.path.join(gfs_dir, '*.idx')):
            os.remove(idx_file)


# ============================================================
# make_ini
# ============================================================
def make_ini(input_file,output_dir,ini_date,Yorig,fname_out):
    '''
    Make CROCO initial conditions file from an OGCM file
    
    This function is mostly taken from croco_pytools/prepro/make_ini.py
    We just adjust how the inputs are handled
    Most of the inputs are contained in a crocotools_param.py file, while
    a few others are direct inputs to this function. The direct inputs are for 
    things which we may want to be update as part of an operational/inter-annual workflow
    
    We extract the OGCM data corresponding to the two nearest time-steps to ini_date
    and perform linear interpolation in time before writing the ini file
    
    input_file  - path and filename for the initial file.
    output_dir  - where the output file gets written. This directory must also contain the crocotools_param.py file which contains configurable parameters not already provided as direct inputs.
    ini_date    - the model initialisation time, as a datetime.datetime object. 
    Yorig       - the Yorig value used in setting up the CROCO model
    fname_out   - output ini filename (only the filename, not the full path)
    
    '''

    # imports
    sys.path.append(output_dir)
    import crocotools_param as params

    # Load croco_grd
    # (assumes params.croco_grd is a relative path from output_dir)
    croco_grd = os.path.join(output_dir, params.croco_grd)
    crocogrd = Croco.CROCO_grd(croco_grd, params.sigma_params)

    #--- Load input (restricted to croco_grd) ----------------------------  
    multi_files=params.multi_files
    print(params.inputdata)
    print(input_file)
    print(crocogrd)
    inpdat=Inp.getdata(params.inputdata,input_file,crocogrd,multi_files,params.tracers)
    print(inpdat)

    # --- Create the initial file -----------------------------------------
    fname_out = os.path.join(output_dir,fname_out)
    Croco.CROCO.create_ini_nc(None,fname_out,crocogrd,
                              tracers=params.tracers)

    # --- Handle initial time ---------------------------------------------
    # get Yorig-01-01 in days since 1970-01-01
    ref_datenum = date2num(datetime(Yorig,1, 1))
    
    # convert ini_date to days since Yorig-01-01
    ini_datenum = date2num(ini_date) - ref_datenum
    
    # Read the time from the input file, and convert to days since Yorig-01-01
    input_datenums = date2num(inpdat.ncglo['ssh'].time.values) - ref_datenum
    
    # find the two nearest indices corresponding to the requested ini_date
    sorted_indices = np.argsort(abs(input_datenums - ini_datenum))
    closest_indices = np.sort(sorted_indices[0:2])
    
    # start and end time in days (not sure why we need this in the file?)
    tstart = ini_datenum
    tend   = ini_datenum
    
    # scrumt and oceant, in seconds since Yorig-01-01
    # this is what is used to define the initial time in the model
    scrumt = ini_datenum * 86400
    oceant = ini_datenum * 86400
 
    #  --- Compute and save variables on CROCO grid ---------------

    for vars in ['ssh','tracers','velocity']:
        print('\nProcessing *%s*' %vars)
        nc=netcdf.Dataset(fname_out, 'a')
        if vars == 'ssh' :
            (zeta,NzGood) = interp_tools.interp_tracers(inpdat,vars,-1,crocogrd,\
                                                        closest_indices[0],closest_indices[1]
                                                        )
            
            # interpolate in time to ini_datenum    
            interp_func = interp1d(input_datenums[closest_indices], zeta, axis=0, kind='linear')
            zeta        = interp_func(ini_datenum)
            
            # write to the nc file
            nc.variables['zeta'][0,:,:] = zeta*crocogrd.maskr
            nc.Input_data_type=params.inputdata
            nc.variables['ocean_time'][:] = oceant
            nc.variables['scrum_time'][:] = scrumt
            nc.variables['scrum_time'].units='seconds since %s-01-01 00:00:00'%Yorig
            nc.variables['tstart'][:] = tstart
            nc.variables['tend'][:] = tend
            z_rho = crocogrd.scoord2z_r(zeta=zeta)
            z_w   = crocogrd.scoord2z_w(zeta=zeta)
            
        elif vars == 'tracers':
            for tra in params.tracers:
                print(f'\nIn tracers processing {tra}')
                trac_3d=[]
                # apparently interp_tools.interp can't extract multiple time-steps for 3D 
                # variables so we have to extract the two nearest time-steps in a loop
                for i in range(2):
                    trac_3d_at_i= interp_tools.interp(inpdat,tra,params.Nzgoodmin,z_rho,crocogrd,\
                                                      closest_indices[i],closest_indices[i]
                                                      )
                    trac_3d.append(trac_3d_at_i.squeeze(axis=0)) 
                trac_3d      = np.stack(trac_3d, axis=0)
                
                # interpolate in time to ini_datenum   
                interp_func = interp1d(input_datenums[closest_indices], trac_3d, axis=0, kind='linear')
                trac_3d     = interp_func(ini_datenum)
                
                # write to the nc file
                nc.variables[tra][0,:,:,:] = trac_3d*crocogrd.mask3d()

        elif vars == 'velocity':

            cosa=np.cos(crocogrd.angle)
            sina=np.sin(crocogrd.angle)

            # apparently interp_tools.interp can't extract multiple time-steps for 3D 
            # variables so we have to extract the two nearest time-steps in a loop
            u,v,ubar,vbar=[],[],[],[]
            for i in range(2):
                [u_i,v_i,ubar_i,vbar_i]=interp_tools.interp_uv(inpdat,params.Nzgoodmin,z_rho,cosa,sina,crocogrd,\
                                                               closest_indices[i],closest_indices[i]
                                                               )
                u.append(u_i.squeeze(axis=0)),v.append(v_i.squeeze(axis=0)),ubar.append(ubar_i.squeeze(axis=0)),vbar.append(vbar_i.squeeze(axis=0))
                
            # interpolate in time to ini_datenum   
            u           = np.stack(u, axis=0)
            interp_func = interp1d(input_datenums[closest_indices], u, axis=0, kind='linear')
            u           = interp_func(ini_datenum)
            
            v      = np.stack(v, axis=0)
            interp_func = interp1d(input_datenums[closest_indices], v, axis=0, kind='linear')
            v           = interp_func(ini_datenum)
            
            ubar        = np.stack(ubar, axis=0)
            interp_func = interp1d(input_datenums[closest_indices], ubar, axis=0, kind='linear')
            ubar        = interp_func(ini_datenum)
            
            vbar        = np.stack(vbar, axis=0)
            interp_func = interp1d(input_datenums[closest_indices], vbar, axis=0, kind='linear')
            vbar        = interp_func(ini_datenum)
            
            conserv=1  # Correct the horizontal transport i.e. remove the intergrated tranport and add the OGCM transport          
            if conserv == 1:
                (ubar_croco,h0)=sig_tools.vintegr(u,grd_tools.rho2u(z_w),grd_tools.rho2u(z_rho),np.nan,np.nan)/grd_tools.rho2u(crocogrd.h)
                (vbar_croco,h0)=sig_tools.vintegr(v,grd_tools.rho2v(z_w),grd_tools.rho2v(z_rho),np.nan,np.nan)/grd_tools.rho2v(crocogrd.h)

                u = u - ubar_croco ; u = u + np.tile(ubar,(z_rho.shape[0],1,1))
                v = v - vbar_croco ; v = v + np.tile(vbar,(z_rho.shape[0],1,1))
                
            # write to the nc file
            nc.variables['u'][0,:,:,:] = u *crocogrd.umask3d()
            nc.variables['v'][0,:,:,:] = v * crocogrd.vmask3d()
            nc.variables['ubar'][0,:,:] = ubar *crocogrd.umask
            nc.variables['vbar'][0,:,:] = vbar * crocogrd.vmask

        nc.close()

    print('')
    print(' Initial file created ')
    print(' Path to file is ', fname_out)
    print('')


# ============================================================
# make_bry
# ============================================================
def make_bry(input_file,output_dir,start_date,end_date,Yorig,fname_out):
    '''
    Make CROCO boundary file from an OGCM file
    
    This function is mostly taken from croco_pytools/prepro/make_bry.py
    We just adjust how the inputs are handled
    Most of the inputs are contained in a crocotools_param.py file, while
    a few others are direct inputs to this function. The direct inputs are for 
    things which we may want to be update as part of an operational/inter-annual workflow
    
    input_file  - path and filename for the OGCM file. TODO: I think this should also work with a list of files (to be checked)
    output_dir  - where the output file gets written. This directory must also contain the crocotools_param.py file which contains configurable parameters not already provided as direct inputs.
    start_date  - start of the bry file, as a datetime.datetime object. 
    end_date    - end of the bry file, as a datetime.datetime object. 
    Yorig       - the Yorig value used in setting up the CROCO model
    fname_out   - output bry filename (only the filename, not the full path)
    
    '''  

    # imports
    sys.path.append(output_dir)
    import crocotools_param as params
    
    # --- Load croco_grd --------------------------------------------------
    # (assumes params.croco_grd is a relative path from output_dir)
    croco_grd = os.path.join(output_dir, params.croco_grd)
    crocogrd = Croco.CROCO_grd(croco_grd, params.sigma_params)

    # --- Initialize boundary vars ----------------------------------------
    crocogrd.WEST_grid()
    crocogrd.EAST_grid()
    crocogrd.SOUTH_grid()
    crocogrd.NORTH_grid()

    # --- Initialize input data class -------------------------------------
    multi_files=params.multi_files
    inpdat = Inp.getdata(params.inputdata,input_file,crocogrd,multi_files,
                         params.tracers,
                         bdy=[params.obc_dict,params.cycle_bry])

    fname_out = os.path.join(output_dir,fname_out)
    Croco.CROCO.create_bry_nc(None,fname_out,crocogrd,params.obc_dict,params.cycle_bry,tracers=params.tracers)

    # --- Handle bry_time --------------------------------------------
    
    # get Yorig-01-01 in days since 1970-01-01
    ref_datenum = date2num(datetime(Yorig,1,1))
    # Load full time dataset in days since Yorig
    input_datenums = date2num(inpdat.ncglo['time'].values) - ref_datenum
    # get the start and end times in days since Yorig
    start_datenum = date2num(start_date) - ref_datenum
    end_datenum = date2num(end_date) - ref_datenum    
    # find nearest indices in input_datenums for start_date and end_date
    [dtmin,dtmax] = np.argmin(abs(input_datenums-start_datenum)),np.argmin(abs(input_datenums-end_datenum))
    # get the bry_time as a subset of input_datenums
    bry_time=input_datenums[dtmin:dtmax+1]
    
    print('bry_time: ', bry_time)

    # write the dates to the file
    nc=netcdf.Dataset(fname_out, 'a')
    nc.Input_data_type=params.inputdata
    nc.variables['bry_time'].cycle=params.cycle_bry
    nc.variables['bry_time'][:]=bry_time
    
    if params.cycle_bry==0:
        nc.variables['bry_time'].units='days since %s-01-01 00:00:00' %(Yorig)

    # --- Loop on boundaries ------------------------------------------
    prev=0
    nxt=0

    if len(params.tracers) == 0:
        var_loop = ['ssh','velocity']
    else:
        var_loop = ['ssh','tracers','velocity']

    for boundary, is_open in zip(params.obc_dict.keys(), params.obc_dict.values()):
        if is_open:
            for vars in var_loop:
                print('\n     Processing *%s* for %sern boundary' %(vars, boundary))
                print('     ------------------------------------------')
                if vars == 'ssh':
                    (zeta,NzGood) = interp_tools.interp_tracers(inpdat,vars,-1,crocogrd,dtmin,dtmax,prev,nxt,boundary[0].upper())
                    z_rho = crocogrd.scoord2z_r(zeta=zeta,bdy="_"+boundary)
                    z_w   = crocogrd.scoord2z_w(zeta=zeta,bdy="_"+boundary)
                elif vars == 'tracers':
                    trac_dict = dict()
                    for trc in params.tracers:
                        print(f'\nIn tracers processing {trc}')
                        trac_dict[trc] = interp_tools.interp(inpdat,trc,params.Nzgoodmin,z_rho,crocogrd,dtmin,dtmax,prev,nxt,bdy=boundary[0].upper())

                elif vars == 'velocity':
                    cosa=np.cos(eval(''.join(('crocogrd.angle_',boundary))) )
                    sina=np.sin(eval(''.join(('crocogrd.angle_',boundary))) )

                    [u,v,ubar,vbar]=interp_tools.interp_uv(inpdat,params.Nzgoodmin,z_rho,cosa,sina,crocogrd,dtmin,dtmax,prev,nxt,bdy=boundary[0].upper())

                    conserv=1  # Correct the horizontal transport i.e. remove the intergrated tranport and add the OGCM transport

                    if conserv == 1:
                        ubar_croco=sig_tools.vintegr4D(u,grd_tools.rho2u(z_w),grd_tools.rho2u(z_rho),np.nan,np.nan)[0]/grd_tools.rho2u(eval(''.join(('crocogrd.h_'+boundary))))
                        vbar_croco=sig_tools.vintegr4D(v,grd_tools.rho2v(z_w),grd_tools.rho2v(z_rho),np.nan,np.nan)[0]/grd_tools.rho2v(eval(''.join(('crocogrd.h_'+boundary))))

                        u = u - np.tile(ubar_croco[:,np.newaxis,:,:],(1,z_rho.shape[1],1,1))
                        u = u + np.tile(ubar[:,np.newaxis,:,:],(1,z_rho.shape[1],1,1))

                        v = v - np.tile(vbar_croco[:,np.newaxis,:,:],(1,z_rho.shape[1],1,1))
                        v = v + np.tile(vbar[:,np.newaxis,:,:],(1,z_rho.shape[1],1,1))

    # --- Saving in netcdf ------------------------------------------------

            print('\nSaving %sern boundary in Netcdf' % boundary)
            print('----------------------------------')

             # handle indices (as 2 points where taken next to bdy)
            if str(boundary) == 'west' and is_open:
                indices3D="[:,:,:,0]" # T,N,J,i=0
                indices2D="[:,:,0]"   # T,J,i=0
            elif str(boundary) == 'east' and is_open:
                indices3D="[:,:,:,-1]" # T,N,J,i=last
                indices2D="[:,:,-1]"   # T,J,i=last
            elif str(boundary) == 'south' and is_open:
                indices3D="[:,:,0,:]" # T,N,j=0,I
                indices2D="[:,0,:]"   # T,j=0,I
            elif str(boundary) == 'north' and is_open:
                indices3D="[:,:,-1,:]" # T,N,j=last,I
                indices2D="[:,-1,:]"   # T,j=last,I

            mask_zet = np.tile(eval(''.join(('crocogrd.maskr_',boundary))),[zeta.shape[0],1,1])
            if "velocity" in var_loop:
                mask_u   = np.tile(eval(''.join(('crocogrd.umask_',boundary))),[u.shape[0],u.shape[1],1,1])
                mask_v   = np.tile(eval(''.join(('crocogrd.vmask_',boundary))),[u.shape[0],u.shape[1],1,1])
                mask_ubar   = np.tile(eval(''.join(('crocogrd.umask_',boundary))),[u.shape[0],1,1])
                mask_vbar   = np.tile(eval(''.join(('crocogrd.vmask_',boundary))),[v.shape[0],1,1])

            nc.variables['zeta_'+str(boundary)][:]=eval(''.join(('zeta',indices2D)))*eval(''.join(('mask_zet',indices2D)))
            nc.variables['u_'+str(boundary)][:]   =eval(''.join(('u',indices3D)))*eval(''.join(('mask_u',indices3D)))
            nc.variables['v_'+str(boundary)][:]   =eval(''.join(('v',indices3D)))*eval(''.join(('mask_v',indices3D)))
            nc.variables['ubar_'+str(boundary)][:]=eval(''.join(('ubar',indices2D)))*eval(''.join(('mask_ubar',indices2D)))
            nc.variables['vbar_'+str(boundary)][:]=eval(''.join(('vbar',indices2D)))*eval(''.join(('mask_vbar',indices2D)))

            if 'tracers' in var_loop:
                for varname, value in zip(trac_dict.keys(), trac_dict.values()):
                    mask_tra = np.tile(eval(''.join(('crocogrd.maskr_',boundary))),[value.shape[0],value.shape[1],1,1])
                    nc.variables[f"{varname}_{boundary}"][:] = eval(f'value{indices3D}')*eval(''.join(('mask_tra',indices3D)))

    nc.close()

    print('')
    print(' Boundary file created ')
    print(' Path to file is ', fname_out)
    print('')
