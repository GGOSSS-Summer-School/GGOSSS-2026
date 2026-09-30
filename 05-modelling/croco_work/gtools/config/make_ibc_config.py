#!/usr/bin/env python3
"""
CLI Script to generate CROCO IBC configurations.

Usage:
    python make_ibc_config.py IGOG_12 50 20130101 20130131
"""

# python3 make_ibc_config.py    "IGOG_12" 50 20130101 20130131
# python3 make_ibc_config.py  "Canary_12" 50 20130101 20130131
# python3 make_ibc_config.py "Agulhas_12" 50 20130101 20130131
# #
# python3 make_ibc_config.py    "IGOG_25" 100 20130101 20130131
# python3 make_ibc_config.py  "Canary_25" 100 20130101 20130131
# python3 make_ibc_config.py "Agulhas_25" 100 20130101 20130131
# #
# python3 make_ibc_config.py    "IGOG_50" 150 20130101 20130131
# python3 make_ibc_config.py  "Canary_50" 150 20130101 20130131
# python3 make_ibc_config.py "Agulhas_50" 150 20130101 20130131

import argparse
import sys
from datetime import datetime
from utils.ibc_utils import write_ibc_config

def valid_date(s):
    """Parses a string in YYYYMMDD format into a datetime object."""
    try:
        return datetime.strptime(s, "%Y%m%d")
    except ValueError:
        msg = f"Not a valid date: '{s}'. Expected format: YYYYMMDD."
        raise argparse.ArgumentTypeError(msg)

def main():
    parser = argparse.ArgumentParser(
        description="Generate an ibc.ini file for CROCO simulations based on environment variables and defaults."
    )

    # Positional Arguments
    parser.add_argument("config_name", type=str, help="Name of the simulation configuration (e.g., IGOG_12)")
    parser.add_argument("n_levels", type=int,help="Number of vertical levels (N)")
    parser.add_argument("start_date", type=valid_date, help="Simulation start date in YYYYMMDD format")
    parser.add_argument("end_date", type=valid_date, help="Simulation end date in YYYYMMDD format")
    #
    args = parser.parse_args()
    #
    try:
        # Call the core configuration function
        write_ibc_config(
            config_name=args.config_name,
            n_levels=args.n_levels,
            start_date=args.start_date,
            end_date=args.end_date
        )
    except EnvironmentError as e:
        print(f"Error: {e}")
        sys.exit(1)
    except Exception as e:
        print(f"An unexpected error occurred: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()
