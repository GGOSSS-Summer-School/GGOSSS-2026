"""Generate a CROCO grid.ini configuration file based on geographical extents and resolution.

Usage:
python3 make_grid_config.py configname lon_min lon_max lat_min lat_max dlon dlat
python3 make_grid_config.py "BENGUELA_12" 6.5 23.5 -38 -26 0.083333 0.083333
"""


# # 1/12 degree resolution
# RES_12=$(echo "1/12" | bc -l)
# RES_25=$(echo "1/25" | bc -l)
# RES_50=$(echo "1/50" | bc -l)
# python3 make_grid_config.py    "IGOG_12"   5.3  12.3  -5   5 $RES_12 $RES_12
# python3 make_grid_config.py  "Canary_12" -22   -15.5  14  24 $RES_12 $RES_12
# python3 make_grid_config.py "Agulhas_12"  17    30   -38 -32 $RES_12 $RES_12
# #
# python3 make_grid_config.py    "IGOG_25"   5.3  12.3  -5   5 $RES_25 $RES_25
# python3 make_grid_config.py  "Canary_25" -22   -15.5  14  24 $RES_25 $RES_25
# python3 make_grid_config.py "Agulhas_25"  17    30   -38 -32 $RES_25 $RES_25
# #
# python3 make_grid_config.py    "IGOG_50"   5.3  12.3  -5   5 $RES_50 $RES_50
# python3 make_grid_config.py  "Canary_50" -22   -15.5  14  24 $RES_50 $RES_50
# python3 make_grid_config.py "Agulhas_50"  17    30   -38 -32 $RES_50 $RES_50


import argparse
import logging
from utils import create_croco_grid_config

# Example Usage:
# create_croco_grid_config("BENGUELA_LR", 6.5, 23.5, -38, -26, 1./12., 1./12.)
#
# create_croco_grid_config(   "IGOG_12",   5.3,  12.3,  -5,   5, 1./12., 1./12.)
# create_croco_grid_config( "Canary_12", -22  , -15.5,  14,  24, 1./12., 1./12.)
# create_croco_grid_config("Agulhas_12",  17  ,  30  , -38, -32, 1./12., 1./12.)
#
# create_croco_grid_config(   "IGOG_25",   5.3,  12.3,  -5,   5, 1./25., 1./25.)
# create_croco_grid_config( "Canary_25", -22  , -15.5,  14,  24, 1./25., 1./25.)
# create_croco_grid_config("Agulhas_25",  17  ,  30  , -38, -32, 1./25., 1./25.)
#
# create_croco_grid_config(   "IGOG_50",   5.3,  12.3,  -5,   5, 1./50., 1./50.)
# create_croco_grid_config( "Canary_50", -22  , -15.5,  14,  24, 1./50., 1./50.)
# create_croco_grid_config("Agulhas_50",  17  ,  30  , -38, -32, 1./50., 1./50.)

def main():
    parser = argparse.ArgumentParser(
        description="Generate a CROCO grid.ini configuration file based on geographical extents and resolution.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    # Required Positional Arguments
    parser.add_argument("configname", help="Name of the configuration (subdirectory name)")
    parser.add_argument("lon_min", type=float, help="Minimum longitude of the domain")
    parser.add_argument("lon_max", type=float, help="Maximum longitude of the domain")
    parser.add_argument("lat_min", type=float, help="Minimum latitude of the domain")
    parser.add_argument("lat_max", type=float, help="Maximum latitude of the domain")
    parser.add_argument("dlon", type=float, help="Grid step in longitude (degrees)")
    parser.add_argument("dlat", type=float, help="Grid step in latitude (degrees)")

    args = parser.parse_args()

    try:
        create_croco_grid_config(
            args.configname, args.lon_min, args.lon_max,
            args.lat_min, args.lat_max, args.dlon, args.dlat
        )
    except Exception as e:
        print(f"Error: {e}")
        raise SystemExit(1)   # non-zero exit so scripts using `set -e` stop here

if __name__ == "__main__":
    main()
