"""
CROCO Workflow Tools
--------------------
A collection of utilities for managing Coastal and Regional Ocean
Community (CROCO) model configurations, environment paths,
and pre/post-processing tasks.

This package provides automated environment validation and
configuration management to ensure consistent directory structures
across different simulation regions and resolutions.
"""
import logging
import sys
from .sf_utils import * # get_croco_env
from .grid_utils import *
from .defaults import *

# 1. Setup Package-Level Logging
sf_logger = logging.getLogger(__name__)
sf_logger.addHandler(logging.NullHandler())  # Prevents "No handler found" warnings

def setup_sf_logging(level=logging.INFO):
    """
    Optional helper to configure a basic console logger for the package.
    """
    #handler = logging.StreamHandler()
    handler = logging.StreamHandler(sys.stdout) # use the same stream as print
    formatter = logging.Formatter('[%(levelname)s] - %(message)s')
    # Add the name of the package
    # formatter = logging.Formatter('%(name)s - %(levelname)s - %(message)s')
    handler.setFormatter(formatter)
    sf_logger.addHandler(handler)
    sf_logger.setLevel(level)


__all__ = [
    "get_croco_env",
    "setup_sf_logging",
    "sf_logger",
    # grid utils
    "create_croco_grid_config",
    # defaults
    "default_config",
    "croco_data_suffix"
]

__version__ = "0.1.0"
