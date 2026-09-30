"""
Set of tools used for croco configuration in croco_work.
"""

import os
from datetime import datetime

def get_croco_env0(
    need_model_dir=True,
    need_pytools_dir=True,
    need_data_root=True,
    need_configs_root=True,
    need_inputs_root=True,
    need_runs_root=True
):
    """
    Resolves and validates directory paths for the CROCO modeling workflow.

    The function scans the environment for six specific variables. If a specific variable
    is not set, it attempts to compute the path using the 'CROCO_ROOT' base directory.

    Paths resolved:
    - CROCO_MODEL_DIR: Source code for the CROCO model.
    - CROCO_PYTOOLS_DIR: Source code for the CROCO python tools.
    - CROCO_DATA_ROOT: Sample or reference data.
    - CROCO_CONFIGS_ROOT: Root for region/resolution configurations.
    - CROCO_INPUTS_ROOT: Root for NetCDF inputs (grid, forcing, IBC).
    - CROCO_RUNS_ROOT: Root for model simulation outputs.

    Logic:
    1. Check for specific environment variables.
    2. If missing, use CROCO_ROOT + predefined subdirectories:
       - model: code/croco
       - pytools: code/croco_pytools
       - data: data
       - configs: hindcast/configs
       - inputs: hindcast/scratch
       - runs: hindcast/scratch
    3. Validation: If 'need_xxx' is True, raises EnvironmentError if the path cannot be
       determined. For model, pytools, and data, also verifies that the directory exists.

    Args:
        need_model_dir (bool): If True, validates CROCO_MODEL_DIR. Defaults to True.
        need_pytools_dir (bool): If True, validates CROCO_PYTOOLS_DIR. Defaults to True.
        need_data_root (bool): If True, validates CROCO_DATA_ROOT. Defaults to True.
        need_configs_root (bool): If True, validates CROCO_CONFIGS_ROOT. Defaults to True.
        need_inputs_root (bool): If True, validates CROCO_INPUTS_ROOT. Defaults to True.
        need_runs_root (bool): If True, validates CROCO_RUNS_ROOT. Defaults to True.

    Returns:
        dict: A dictionary mapping environment variable names to their resolved paths.

    Raises:
        EnvironmentError: If a required variable is missing and CROCO_ROOT is not set.
        FileNotFoundError: If model, pytools, or data directories do not exist on disk.
    """
    croco_root = os.environ.get("CROCO_ROOT")

    # Mapping of argument to (Env Var Name, Default Subpath)
    config_map = {
        'model_dir':   (need_model_dir,   "CROCO_MODEL_DIR",   "code/croco"),
        'pytools_dir': (need_pytools_dir, "CROCO_PYTOOLS_DIR", "code/croco_pytools"),
        'data_root':   (need_data_root,   "CROCO_DATA_ROOT",   "data"),
        'configs_root':(need_configs_root,"CROCO_CONFIGS_ROOT","hindcast/configs"),
        'inputs_root': (need_inputs_root, "CROCO_INPUTS_ROOT", "hindcast/scratch"),
        'runs_root':   (need_runs_root,   "CROCO_RUNS_ROOT",   "hindcast/scratch"),
    }

    resolved_paths = {}

    for key, (is_needed, env_var, subpath) in config_map.items():
        val = os.environ.get(env_var)

        # Fallback to CROCO_ROOT if not set
        if not val:
            if croco_root:
                val = os.path.join(croco_root, subpath)
            elif is_needed:
                raise EnvironmentError(
                    f"Required variable {env_var} is not set and CROCO_ROOT is missing."
                )

        # Physical existence validation for source/data directories
        if is_needed and key in ['model_dir', 'pytools_dir', 'data_root']:
            if val and not os.path.isdir(val):
                raise FileNotFoundError(f"Required directory for {env_var} does not exist: {val}")

        resolved_paths[env_var] = val

    return resolved_paths


import logging

logger = logging.getLogger(__name__)

def get_croco_env(
    need_model_dir=True,
    need_pytools_dir=True,
    need_data_root=True,
    need_configs_root=True,
    need_inputs_root=True,
    need_runs_root=True
):
    """
    Scans environment for CROCO paths, falling back to CROCO_ROOT.
    Collects and reports all missing or invalid directories simultaneously.

    The function scans the environment for six specific variables. If a specific variable
    is not set, it attempts to compute the path using the 'CROCO_ROOT' base directory.

    Paths resolved:
    - CROCO_MODEL_DIR: Source code for the CROCO model.
    - CROCO_PYTOOLS_DIR: Source code for the CROCO python tools.
    - CROCO_DATA_ROOT: Sample or reference data.
    - CROCO_CONFIGS_ROOT: Root for region/resolution configurations.
    - CROCO_INPUTS_ROOT: Root for NetCDF inputs (grid, forcing, IBC).
    - CROCO_RUNS_ROOT: Root for model simulation outputs.

    Logic:
    1. Check for specific environment variables.
    2. If missing, use CROCO_ROOT + predefined subdirectories:
       - model: code/croco
       - pytools: code/croco_pytools
       - data: data
       - configs: hindcast/configs
       - inputs: hindcast/scratch
       - runs: hindcast/scratch
    3. Validation: If 'need_xxx' is True, raises EnvironmentError if the path cannot be
       determined. For model, pytools, and data, also verifies that the directory exists.

    Args:
        need_model_dir (bool): If True, validates CROCO_MODEL_DIR. Defaults to True.
        need_pytools_dir (bool): If True, validates CROCO_PYTOOLS_DIR. Defaults to True.
        need_data_root (bool): If True, validates CROCO_DATA_ROOT. Defaults to True.
        need_configs_root (bool): If True, validates CROCO_CONFIGS_ROOT. Defaults to True.
        need_inputs_root (bool): If True, validates CROCO_INPUTS_ROOT. Defaults to True.
        need_runs_root (bool): If True, validates CROCO_RUNS_ROOT. Defaults to True.

    Returns:
        dict: A dictionary mapping environment variable names to their resolved paths.

    Raises:
        EnvironmentError: If a required variable is missing and CROCO_ROOT is not set.
        FileNotFoundError: If model, pytools, or data directories do not exist on disk.
    """
    croco_root = os.environ.get("CROCO_ROOT")

    config_map = {
        'model_dir':   (need_model_dir,   "CROCO_MODEL_DIR",   "code/croco"),
        'pytools_dir': (need_pytools_dir, "CROCO_PYTOOLS_DIR", "code/croco_pytools"),
        'data_root':   (need_data_root,   "CROCO_DATA_ROOT",   "data"),
        'configs_root':(need_configs_root,"CROCO_CONFIGS_ROOT","hindcast/configs"),
        'inputs_root': (need_inputs_root, "CROCO_INPUTS_ROOT", "hindcast/scratch"),
        'runs_root':   (need_runs_root,   "CROCO_RUNS_ROOT",   "hindcast/scratch"),
    }

    resolved_paths = {}
    errors = []

    for key, (is_needed, env_var, subpath) in config_map.items():
        val = os.environ.get(env_var)

        # Fallback to CROCO_ROOT
        if not val:
            if croco_root:
                val = os.path.join(croco_root, subpath)
            elif is_needed:
                errors.append(f"Variable '{env_var}' is not set and CROCO_ROOT is missing.")
                continue # Skip further checks for this specific variable

        # Physical directory validation
        if is_needed and key in ['model_dir', 'pytools_dir', 'data_root']:
            if not os.path.isdir(val):
                errors.append(f"Directory for '{env_var}' does not exist: {val}")

        resolved_paths[env_var] = val
        if is_needed:
            logger.info(f"Resolved {env_var} to {val}")
        else:
            logger.info(f"{key} ({env_var}) not needed")
            pass

    # If any errors were collected, raise them all at once
    if errors:
        error_msg = "Environment validation failed with the following errors:\n" + "\n".join(f"  - {e}" for e in errors)
        raise EnvironmentError(error_msg)

    return resolved_paths
