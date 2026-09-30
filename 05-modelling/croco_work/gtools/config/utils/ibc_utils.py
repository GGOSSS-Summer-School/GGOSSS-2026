import os
import configparser
from datetime import datetime

from copy import deepcopy

from .sf_utils import get_croco_env
from .defaults import default_config, croco_data_suffix

def write_ibc_config(config_name, n_levels, start_date, end_date):
    """
    Generates a localized 'ibc.ini' configuration file for a specific CROCO simulation.

    This function synchronizes system environment variables, project-wide defaults,
    and simulation-specific parameters into a single INI file used by the model
    preprocessing tools.

    Args:
        config_name (str): Name of the simulation (e.g., 'IGOG_12'). Used to
            construct the directory structure in the scratch space.
        n_levels (int): Number of vertical terrain-following (sigma) layers.
        start_date (datetime): The beginning of the simulation period; used to
            set model time and select the initial condition data file.
        end_date (datetime): The end of the simulation period.

    Returns:
        None: Writes 'ibc.ini' directly to the current working directory.

    Note:
        - Paths are resolved via 'get_croco_env' which prioritizes environment
          variables over the 'CROCO_ROOT' fallback.
        - Inherits all non-specified configurations (e.g., Grid Smoothing,
          Tide Options) from 'defaults.py' to ensure consistency.
    """
    # 1. Retrieve validated system paths (Model, Data, and Run roots)
    env = get_croco_env(
        need_model_dir=False,
        need_pytools_dir=False,
        need_data_root=True,
        need_configs_root=True,
        need_inputs_root=False,
        need_runs_root=True
    )
    runs_root = env['CROCO_RUNS_ROOT']
    data_root = env['CROCO_DATA_ROOT']

    # Build Configuration
    config = configparser.ConfigParser()
    config.optionxform = str  # Keep keys case-sensitive
    #
    # 2. Copy sections from the defaults and update the keys that need updates.
    # This ensures that if a default is updated in defaults.py, it propagates here.

    config['Croco_Files'] = deepcopy(default_config['Croco_Files'])
    # Update paths: Define where preprocessed files will be stored
    config['Croco_Files']['croco_files_dir'] = os.path.join(runs_root, config_name, 'CROCO_FILES')
    #
    config['Zoom_Options'] = deepcopy(default_config['Zoom_Options'])
    #
    config['Times'] = deepcopy(default_config['Times'])
    # Update Time: Map Python datetime objects to INI keys
    config['Times'].update({
        'Ystart': str(start_date.year),
        'Mstart': str(start_date.month),
        'Dstart': str(start_date.day),
        'Yend': str(end_date.year),
        'Mend': str(end_date.month),
        'Dend': str(end_date.day),
        'use_calendar': 'True'
    })
    #
    config['Sigma_Params'] = deepcopy(default_config['Sigma_Params'])
    # Update Physics: Set vertical discretization levels
    config['Sigma_Params']['N'] = str(n_levels)
    #
    config['IBC_Options'] = deepcopy(default_config['IBC_Options'])
    #
    config['IBC_Input_Files'] = deepcopy(default_config['IBC_Input_Files'])
    # Update Initial/Boundary Conditions: Resolve data paths and file labels
    config['IBC_Input_Files'].update({
        'ibc_dir': os.path.join(data_root, croco_data_suffix['ibc_dir']),
        'ibc_file_ssh': 'PHY',
        'ibc_files_tracers': '["PHY", "PHY", "BIO", "BIO", "BIO", "BIO", "BIO"]',
        'ibc_file_u': 'PHY',
        'ibc_file_v': 'PHY',
        'ini_filedate': start_date.strftime('Y%YM%m'),
    })

    # Save to ${CROCO_CONFIGS_ROOT}/config_name/ibc.ini
    target_dir = os.path.join(env['CROCO_CONFIGS_ROOT'], config_name)
    os.makedirs(target_dir, exist_ok=True)
    output_filename = os.path.join(target_dir, 'ibc.ini')
    with open(output_filename, 'w') as configfile:
        config.write(configfile)

    print(f"Successfully wrote {output_filename} for {config_name}")
