"""Create a regional FES2022 subset for CoastSat/PyFES.

This script assumes the original FES2022b NetCDF files are already available
locally from AVISO/CNES. It subsets each NetCDF file to a small longitude/latitude
box and writes a copied YAML file whose NetCDF paths point to the subset files.

It is intended for training workflows where the full global FES2022 archive is
too large to move around. Do not commit the generated NetCDF files to GitHub.

Examples:
    conda activate coastsat

    # Small Benin/Lake Nokoue/Cotonou Channel box.
    python notebooks/prepare_fes2022_regional_subset.py \\
        --fes-root /path/to/fes2022b \\
        --out-root /path/to/fes2022b_benin_subset \\
        --region benin_nokoue

    # Broader Gulf of Guinea box.
    python notebooks/prepare_fes2022_regional_subset.py \\
        --fes-root /path/to/fes2022b \\
        --out-root /path/to/fes2022b_gulf_of_guinea_subset \\
        --region gulf_of_guinea
"""

from __future__ import annotations

import argparse
from pathlib import Path
import re
import shutil

REGIONS = {
    "benin_nokoue": {"lon_min": 1.5, "lon_max": 3.0, "lat_min": 5.8, "lat_max": 7.0},
    "gulf_of_guinea": {"lon_min": -8.0, "lon_max": 14.5, "lat_min": -6.0, "lat_max": 8.0},
    "cameroon_pak": {"lon_min": 8.8, "lon_max": 10.4, "lat_min": 1.8, "lat_max": 4.5},
}


def find_coord_name(ds, candidates: list[str]) -> str:
    for name in candidates:
        if name in ds.coords or name in ds.variables:
            return name
    raise KeyError(f"Coordinate not found. Tried: {candidates}")


def subset_slice(values, low: float, high: float):
    first = float(values[0])
    last = float(values[-1])
    if first <= last:
        return slice(low, high)
    return slice(high, low)


def subset_one_file(src: Path, dst: Path, bounds: dict[str, float]) -> None:
    try:
        import xarray as xr
    except ModuleNotFoundError as exc:
        raise ModuleNotFoundError(
            "xarray is required to subset FES2022 NetCDF files. "
            "Install it in the CoastSat environment with: conda install -n coastsat -c conda-forge xarray netcdf4"
        ) from exc

    dst.parent.mkdir(parents=True, exist_ok=True)
    with xr.open_dataset(src) as ds:
        lon_name = find_coord_name(ds, ["lon", "longitude", "x"])
        lat_name = find_coord_name(ds, ["lat", "latitude", "y"])

        lon = ds[lon_name]
        lat = ds[lat_name]
        lon_min = bounds["lon_min"]
        lon_max = bounds["lon_max"]

        if float(lon.min()) >= 0 and lon_min < 0:
            lon_min = lon_min % 360
            lon_max = lon_max % 360

        if lon_min <= lon_max:
            ds_sub = ds.sel(
                {
                    lon_name: subset_slice(lon.values, lon_min, lon_max),
                    lat_name: subset_slice(lat.values, bounds["lat_min"], bounds["lat_max"]),
                }
            )
        else:
            west = ds.sel({lon_name: subset_slice(lon.values, lon_min, float(lon.max()))})
            east = ds.sel({lon_name: subset_slice(lon.values, float(lon.min()), lon_max)})
            ds_sub = xr.concat([west, east], dim=lon_name)
            ds_sub = ds_sub.sel({lat_name: subset_slice(lat.values, bounds["lat_min"], bounds["lat_max"])})

        encoding = {
            var: {"zlib": True, "complevel": 4}
            for var in ds_sub.data_vars
            if ds_sub[var].ndim > 0
        }
        ds_sub.to_netcdf(dst, encoding=encoding)


def update_yaml_paths(yaml_src: Path, yaml_dst: Path, subset_root: Path) -> None:
    text = yaml_src.read_text(encoding="utf-8")
    subset_files = {path.name.lower(): path.resolve() for path in subset_root.rglob("*.nc")}

    def replace_path(match: re.Match[str]) -> str:
        prefix = match.group(1)
        original_path = match.group(2).strip().strip("'\"")
        basename = Path(original_path).name.lower()
        if basename in subset_files:
            return f"{prefix}{subset_files[basename]}"
        return match.group(0)

    text = re.sub(r"(\s*[A-Za-z0-9_]+:\s*)(.+?\.nc)\s*$", replace_path, text, flags=re.MULTILINE)
    yaml_dst.write_text(text, encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--fes-root", required=True, type=Path, help="Folder containing FES2022b files and fes2022.yaml")
    parser.add_argument("--out-root", required=True, type=Path, help="Output folder for the regional subset")
    parser.add_argument("--region", choices=sorted(REGIONS), default="gulf_of_guinea")
    parser.add_argument("--lon-min", type=float)
    parser.add_argument("--lon-max", type=float)
    parser.add_argument("--lat-min", type=float)
    parser.add_argument("--lat-max", type=float)
    args = parser.parse_args()

    bounds = dict(REGIONS[args.region])
    for key in ["lon_min", "lon_max", "lat_min", "lat_max"]:
        cli_value = getattr(args, key, None)
        if cli_value is not None:
            bounds[key] = cli_value

    fes_root = args.fes_root.expanduser().resolve()
    out_root = args.out_root.expanduser().resolve()
    out_root.mkdir(parents=True, exist_ok=True)

    nc_files = sorted(fes_root.rglob("*.nc"))
    if not nc_files:
        raise FileNotFoundError(f"No NetCDF files found in {fes_root}")

    print("region", args.region)
    print("bounds", bounds)
    print("netcdf_files", len(nc_files))

    for src in nc_files:
        rel = src.relative_to(fes_root)
        dst = out_root / rel
        print("subset", rel)
        subset_one_file(src, dst, bounds)

    yaml_candidates = sorted(fes_root.rglob("fes2022*.yaml"))
    if yaml_candidates:
        update_yaml_paths(yaml_candidates[0], out_root / yaml_candidates[0].name, out_root)
        print("yaml_written", out_root / yaml_candidates[0].name)
    else:
        print("warning", "No fes2022 YAML found; copy and edit the YAML manually before using pyfes.")

    readme = out_root / "README_subset.txt"
    readme.write_text(
        "Regional FES2022 subset generated for GGOSSS CoastSat training.\n"
        f"Region: {args.region}\n"
        f"Bounds: {bounds}\n"
        "Generated NetCDF files are training/runtime data and should not be committed to GitHub.\n",
        encoding="utf-8",
    )
    print("done", out_root)


if __name__ == "__main__":
    main()
