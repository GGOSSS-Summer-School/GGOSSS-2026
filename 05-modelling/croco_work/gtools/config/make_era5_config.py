"""
CLI Script to generate CROCO era5 configurations.

Usage:
    python make_era5_config.py IGOG_12 20130101 20130131
"""

# python3 make_era5_config.py    "IGOG_12" 20130101 20130131
# python3 make_era5_config.py  "Canary_12" 20130101 20130131
# python3 make_era5_config.py "Agulhas_12" 20130101 20130131
# #
# python3 make_era5_config.py    "IGOG_25" 20130101 20130131
# python3 make_era5_config.py  "Canary_25" 20130101 20130131
# python3 make_era5_config.py "Agulhas_25" 20130101 20130131
# #
# python3 make_era5_config.py    "IGOG_50" 20130101 20130131
# python3 make_era5_config.py  "Canary_50" 20130101 20130131
# python3 make_era5_config.py "Agulhas_50" 20130101 20130131

import argparse
from datetime import datetime
from utils.era5_utils import write_era5_config

def valid_date(s):
    try:
        return datetime.strptime(s, "%Y%m%d")
    except ValueError:
        msg = f"Not a valid date: '{s}'. Expected format: YYYYMMDD."
        raise argparse.ArgumentTypeError(msg)

def main():
    parser = argparse.ArgumentParser(description="Generate era5.ini for CROCO simulation.")
    parser.add_argument("config_name", type=str,
                        help="Name of the configuration (e.g., IGOG_12)")
    parser.add_argument("start_date", type=valid_date,
                        help="Start date in YYYYMMDD format")
    parser.add_argument("end_date", type=valid_date,
                        help="End date in YYYYMMDD format")

    args = parser.parse_args()

    if args.end_date < args.start_date:
        parser.error("end_date must be after start_date")

    write_era5_config(args.config_name, args.start_date, args.end_date)

if __name__ == "__main__":
    main()
