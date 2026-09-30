import configparser
import os
from copy import deepcopy
from .defaults import default_config
from .sf_utils import get_croco_env

def write_era5_config(config_name, start_date, end_date):
    """
    Writes era5.ini using deepcopy from defaults and paths from sf_utils.
    start_date and end_date are datetime objects.

    Args:
        config_name (str): Name of the simulation (e.g., 'IGOG_12'). Used to
            construct the directory structure in the scratch space.
        start_date (datetime): The beginning of the simulation period; used to
            set model time and select the initial condition data file.
        end_date (datetime): The end of the simulation period.

    Returns:
        None: Writes 'era5.ini' directly to the config directory.

    Note:
        - Paths are resolved via 'get_croco_env' which prioritizes environment
          variables over the 'CROCO_ROOT' fallback.
        - Inherits all non-specified configurations (e.g., Grid Smoothing,
          Tide Options) from 'defaults.py' to ensure consistency.
    """
    # Get environment paths
    # Request only what we need for this specific task
    try:
        env = get_croco_env(
            need_model_dir=False,
            need_pytools_dir=False,
            need_data_root=False,
            need_configs_root=True,
            need_inputs_root=False,
            need_runs_root=True
        )
    except (EnvironmentError, FileNotFoundError) as e:
        print(f"Path resolution error: {e}")
        return
    #
    runs_root = env['CROCO_RUNS_ROOT']
    configs_root = env['CROCO_CONFIGS_ROOT']

    # Define project-specific paths
    croco_files_dir = os.path.join(runs_root, config_name, 'CROCO_FILES/')
    era5_dir = os.path.join(runs_root, config_name, 'DATA_METEO/ERA5/')
    configfile = os.path.join(configs_root, config_name, 'era5.ini')

    # Create the config object
    era5_config = configparser.ConfigParser()
    era5_config.optionxform = str  # Keep keys case-sensitive

    # Deepcopy relevant sections
    relevant_sections = ['Times', 'Download_Options', 'Croco_Files', 'ERA5_Download']
    for section in relevant_sections:
        era5_config[section] = deepcopy(dict(default_config[section]))

    # Update [Times]
    era5_config['Times'].update({
        'Ystart': str(start_date.year),
        'Mstart': str(start_date.month),
        'Dstart': str(start_date.day),
        'Yend': str(end_date.year),
        'Mend': str(end_date.month),
        'Dend': str(end_date.day),
        'use_calendar': 'True'
    })

    # Update [Croco_Files] - Keep only relevant keys
    keep_keys = ['croco_grd_prefix']
    filtered_croco = {k: era5_config['Croco_Files'][k] for k in keep_keys}
    filtered_croco['croco_files_dir'] = croco_files_dir
    era5_config['Croco_Files'] = filtered_croco

    # Update [ERA5_Download]
    era5_config['ERA5_Download'].update({
        'era5_dir': era5_dir,
        'units': "['(0-1)', 'K', 'kg m-2 s-1', 'j m-2', 'j m-2', 'K', 'kg kg-1', 'm s-1', 'm s-1']",
        'conv_cff': "[1., 1. , 1000./3600., 1./3600., 1./3600., 1., 1., 1., 1.]"
    })

    # Ensure the config directory exists
    os.makedirs(os.path.dirname(configfile), exist_ok=True)

    with open(configfile, 'w') as f:
        era5_config.write(f)

    print(f"Config successfully written to: {configfile}")
