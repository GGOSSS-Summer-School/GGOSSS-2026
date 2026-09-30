import os
import math
import configparser
import copy

from .sf_utils import get_croco_env
from .defaults import croco_data_suffix, default_config

def create_croco_grid_config(configname, lon_min, lon_max, lat_min, lat_max, dlon, dlat):
    """
    Generates a CROCO grid.ini based on resolution (dlon, dlat).
    Saves to ${CROCO_CONFIGS_ROOT}/configname/grid.ini.
    """
    # 1. Resolve environmental paths
    # CROCO_CONFIGS_ROOT / CROCO_RUNS_ROOT come from hindcast/track.sh,
    # CROCO_DATA_ROOT from env.sh; CROCO_INPUTS_ROOT is not used here.
    env = get_croco_env(
        need_data_root=False,
        need_configs_root=True,
        need_inputs_root=False,
        need_runs_root=True
    )

    # 2. Geometric Calculations
    mean_lat = (lat_min + lat_max) / 2.0
    R = 6371.0  # Earth Radius in km

    dlon_deg = (lon_max - lon_min)
    dlat_deg = (lat_max - lat_min)

    # Calculate number of points based on dlon/dlat
    # N points = (Total Range / Step) + 1
    npoints_x = int(round(dlon_deg / dlon)) + 1
    npoints_y = int(round(dlat_deg / dlat)) + 1

    # size = degree_diff * R * (pi/180). size_x includes cosine correction.
    size_x_km = dlon_deg * math.cos(math.radians(mean_lat)) * R * (math.pi / 180)
    size_y_km = dlat_deg * R * (math.pi / 180)

    # 3. Build Configuration
    config = configparser.ConfigParser()
    config.optionxform = str  # Keep keys case-sensitive

    config['Croco_Files'] = copy.deepcopy(default_config['Croco_Files'])
    config['Croco_Files']['croco_files_dir'] = (
        os.path.join(env['CROCO_RUNS_ROOT'], configname, 'CROCO_FILES/')
        )

    config['Zoom_Options'] = copy.deepcopy(default_config['Zoom_Options'])

    config['Grid_Position'] = {
        'grid_type': 'curvilinear',
        'central_lon': str((lon_min + lon_max) / 2.0),
        'central_lat': str(mean_lat),
        'size_x_km': str(round(size_x_km)),
        'size_y_km': str(round(size_y_km)),
        'npoints_x': str(npoints_x),
        'npoints_y': str(npoints_y),
        'grid_angle': '0.0',
        'lon_min': str(lon_min),
        'lon_max': str(lon_max),
        'lat_min': str(lat_min),
        'lat_max': str(lat_max),
        'dlon': f'{dlon:.6f}', #str(dlon),
        'dlat': f'{dlat:.6f}', #str(dlat)
    }

    config['Grid_Smoothing_Params'] = copy.deepcopy(default_config['Grid_Smoothing_Params'])

    config['Grid_Isolated_Waterbodies'] = {
        'mask_isolated_waterbodies': 'False',
        'main_water_body_x_idx': str(npoints_x // 2), # Default to center
        'main_water_body_y_idx': str(npoints_y // 2)
    }

    config['Grid_Input_Files'] = {
        'topo_file_reader': 'etopo2',
        'topo_file': os.path.join(env['CROCO_DATA_ROOT'], croco_data_suffix['topo_file']),
        'shp_file': os.path.join(env['CROCO_DATA_ROOT'], croco_data_suffix['shp_file'])
    }

    # 4. Save to ${CROCO_CONFIGS_ROOT}/configname/grid.ini
    target_dir = os.path.join(env['CROCO_CONFIGS_ROOT'], configname)
    os.makedirs(target_dir, exist_ok=True)
    save_path = os.path.join(target_dir, 'grid.ini')

    with open(save_path, 'w') as configfile:
        config.write(configfile)

    print(f"Config saved to: {save_path} ({npoints_x}x{npoints_y} points)")
    return save_path
