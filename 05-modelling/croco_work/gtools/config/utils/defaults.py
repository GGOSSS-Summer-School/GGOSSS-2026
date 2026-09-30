"""
Set of default configurations for croco in croco_work.
"""
import configparser

croco_data_suffix={
    'topo_file': 'DATASETS_CROCOTOOLS/Topo/etopo2.nc',
    'shp_file': 'DATASETS_CROCOTOOLS/gshhs/GSHHS_shp/i/GSHHS_i_L1.shp',
    'ibc_dir': 'MERCATOR_GLOB_2013/',
    'tide_dir': 'DATASETS_CROCOTOOLS/TPXO7/',
    'rivers_dir': 'DATA_RIVERS/',
}

import configparser
default_config = configparser.ConfigParser()
default_config.optionxform = str  # Keep keys case-sensitive

default_config['Croco_Files'] = {
    'croco_files_dir': 'None',
    'croco_grd_prefix': 'croco_grd',
    'croco_tide_prefix': 'croco_frc',
    'croco_ini_prefix': 'croco_ini',
    'croco_bry_prefix': 'croco_bry',
    'croco_bry_format': 'MONTHLY',
    'croco_rivers_prefix': 'croco_runoff'
}

default_config['Zoom_Options'] = {
    'is_zoom': 'False',
    'is_agrif': 'False',
    'agrif_level': '0',
    'parent_grid': 'None' # path to parent croco_grd.nc'
}

default_config['Grid_Smoothing_Params'] = {
    'hmin': '50.0',
    'hmax': '6000.0',
    'interp_rad': '2',
    'rfact': '0.2',
    'smooth_meth': 'lsmooth'
}

default_config['Grid_Isolated_Waterbodies'] = {
    'mask_isolated_waterbodies': 'False',
    'main_water_body_x_idx': 'None',
    'main_water_body_y_idx': 'None'
}

default_config['Grid_Input_Files'] = {
    'topo_file_reader': 'etopo2',
    'topo_file': 'None', # path to 'DATASETS_CROCOTOOLS/Topo/etopo2.nc',
    'shp_file': 'None', # path to 'DATASETS_CROCOTOOLS/gshhs/GSHHS_shp/i/GSHHS_i_L1.shp'
}

default_config['Grid_Zoom_Params'] = {
    'north_obc': 'True',
    'south_obc': 'True',
    'west_obc': 'True',
    'east_obc': 'True',
    'merging_area': '5'
}

default_config['Grid_Zoom_Agrif'] = {
    'coef': 'None',
    'imin': 'None',
    'imax': 'None',
    'jmin': 'None',
    'jmax': 'None'
}

default_config['Times'] = {
    'Ystart': '2013',
    'Mstart': '1',
    'Dstart': '1',
    'Hstart': '12',
    'Yend': '2013',
    'Mend': '1',
    'Dend': '31',
    'Hend': '12',
    'Yorig': '1970',
    'Morig': '1',
    'Dorig': '1',
    'Horig': '0',
    'use_calendar': 'False'
}

default_config['Sigma_Params'] = {
    'theta_s': '7',
    'theta_b': '2',
    'N': '32',
    'hc': '200'
}

default_config['IBC_Options'] = {
    'obc_dict': "{'south':1, 'north':1,'west':1, 'east':1}",
    'tracers': '["temp", "salt"]',
    'uv_conserv': '1',
    'min_nb_valid_data': '4'
}

default_config['IBC_Input_Files'] = {
    'ibc_reader': 'mercator_croco',
    'ibc_dir': 'None', # path to 'DATA_IBC/MERCATOR_GLOB_2013/',
    'ibc_prefix': 'mercator',
    'ibc_extension': '.cdf',
    'ibc_freq': '1M',
    'ibc_multi_files': 'False',
    'ibc_file_ssh': 'ETAN',
    'ibc_files_tracers': '["THETA", "SALT"]',
    'ibc_file_u': 'EVEL',
    'ibc_file_v': 'NVEL',
    'ini_filedate': 'Y2013M01',
    'ini_idx': '0'
}

default_config['Tide_Options'] = {
    'tide_waves': 'M2,S2,N2,K2,K1,O1,P1,Q1,Mf,Mm',
    'is_tide_current': 'True',
    'is_tide_potential': 'True',
    'is_correction_ssh': 'True',
    'is_correction_uv': 'True'
}

default_config['Tide_Input_Files'] = {
    'tide_reader': 'tpxo7_croco',
    'tide_type': 'Re_Im',
    'tide_dir': 'None', # path to 'DATASETS_CROCOTOOLS/TPXO7/',
    'tide_single_file': 'TPXO7.nc',
    'tide_multi_files': 'False',
    'tide_multi_files_waves_separated': 'True',
    'tide_multi_files_elev_file': 'h_<tide_wave>tpxo9_atlas_30_v5.nc',
    'tide_multi_files_u_file': 'u<tide_wave>tpxo9_atlas_30_v5.nc',
    'tide_multi_files_v_file': 'u<tide_wave>_tpxo9_atlas_30_v5.nc'
}

default_config['Rivers_Options'] = {
    'nb_iter_mask': '3',
    'rivers_freq': 'DAILY'
}

default_config['Rivers_Input_Files'] = {
    'rivers_file_list': 'None' # path to 'Examples/rivers_list.txt'
}

default_config['Download_Options'] = {
    'use_grd_extent': 'True',
    'grd_extent_padding': '1.0',
    'custom_extent': 'None'
}

default_config['Mercator_Download'] = {
    'dataset': 'cmems_mod_glo_phy_my_0.083deg_P1D-m',
    'variables': "['thetao', 'so', 'uo', 'vo', 'zos']",
    'depths': '[0, 50]'
}

default_config['ERA5_Download'] = {
    'era5_dir': 'None', # path to save 'DATA_METEO/ERA5/',
    'n_overlap': '0',
    'time_frames': '00/06/12/18',
    'variables': "['lsm', 'sst', 'tp', 'strd', 'ssr', 't2m', 'q', 'u10', 'v10']",
    'conv_cff': '[1., 1., 1000./3600., 1./3600., 1./3600., 1., 1., 1., 1.]',
    'units': "['(0-1)', 'K', 'kg m-2 s-1', 'W m-2', 'W m-2', 'K', 'kg kg-1', 'm s-1', 'm s-1']",
    'wave_extract': 'False',
    'wave_var': "['swh', 'mwd', 'pp1d', 'cdww']",
    'wave_conv_cff': '[1.0, 1.0, 1.0, 1.0]',
    'wave_units': "['m', 'Degrees true', 's', 'dimensionless']"
}

default_config['Glofas_Download'] = {
    'dataset': 'cems-glofas-historical',
    'rivers_dir': 'None', # path to 'DATA_RIVERS/',
    'rivers_prefix': 'cems_glofas'
}

default_config['Hycom_Download'] = {
    'dataset': 'http://tds.hycom.org/thredds/dodsC/GLBy0.08/expt_93.0'
}
