# Install CROCO

One-time setup: build the NetCDF stack CROCO needs, get CROCO, and prepare the
Python environment. Run the scripts **in order** — each depends on the previous
one. This builds the **sequential** (single-process) stack, which is the right
choice for a laptop.

> **Done once per machine.** After this you never repeat it — every working
> session is just:
>
> ```bash
> source ~/croco_work/env.sh
> source ~/croco_work/hindcast/track.sh
> conda activate ggosss26
> ```

## 0. System packages (Ubuntu / WSL2)

```bash
sudo apt update
sudo apt install -y build-essential gfortran m4 curl wget git bc \
                    libcurl4-openssl-dev zlib1g-dev nco gettext-base
```

## 1. Build the libraries and get CROCO

Set the repo root (the scripts install into `opt_seq/` and `code/` beneath it):

```bash
export CROCO_ROOT=${HOME}/croco_work   # where the repo lives
```

`NJOBS` (parallel compile jobs) is **chosen automatically** by the build scripts —
they use `min(cores, RAM_GB / 2)`, since each job needs ~2 GB of RAM. You can
override it: `export NJOBS=6`.

```bash
bash install/00_download_libraries.sh   # download HDF5 / NetCDF-C / NetCDF-Fortran source
bash install/01_build_hdf5.sh           # compile + install HDF5             -> opt_seq/
bash install/02_build_netcdf_c.sh       # compile + install NetCDF-C         -> opt_seq/ (needs HDF5)
bash install/03_build_netcdf_fortran.sh # compile + install NetCDF-Fortran   -> opt_seq/ (needs NetCDF-C)
bash install/04_get_croco.sh            # download CROCO v2.1.3 + croco_pytools v2.0.4 -> code/
```

- **`04`** also compiles the croco_pytools Fortran helpers (`toolsf*.so`, needed
  by `make_grid.py`). Run it with `source env.sh` done **and** the `ggosss26`
  env active (it needs `f2py` from the env and `FFLAGS` from `env.sh`).

✅ **Check the stack** — after `03`:

```bash
source ${CROCO_ROOT}/env.sh
for f in libhdf5.so libnetcdf.so libnetcdff.so; do
  [ -f ${CROCO_ROOT}/opt_seq/lib/$f ] && echo "ok $f" || echo "MISSING $f"; done
nf-config --prefix        # -> .../croco_work/opt_seq
nf-config --flibs         # -> -L.../opt_seq/lib -lnetcdff -lnetcdf ...
```

## 2. Python environment + the package

```bash
cd ${CROCO_ROOT}
conda env create -f environment.yml    # creates the 'ggosss26' env (once)
conda activate ggosss26
pip install --no-deps -e .             # optional: adds the 'ggosss26' command
```

The CLI can always be run as `python ${GTOOLS_DIR}/ggosss26.py <subcommand>`.
The ini/bry/tide builders (`make_ini_hindcast`, `make_bry_hindcast`, `make_tides`) use the
legacy croco_pytools API that ships **vendored** in `gtools/croco_pytools/`; the separate
`code/croco_pytools` (v2.0.4) is what builds the grid.

## 3. CROCO datasets (bathymetry, coastline, TPXO7 tides, river mouths)

One download (several GB; `-c` resumes an interrupted transfer, so it is safe to leave running
overnight):

```bash
mkdir -p ${CROCO_ROOT}/data && cd ${CROCO_ROOT}/data
wget -c https://data-croco.ifremer.fr/DATASETS/DATASETS_CROCOTOOLS.tar.gz
tar -xzf DATASETS_CROCOTOOLS.tar.gz
ls ${CROCO_ROOT}/data/DATASETS_CROCOTOOLS/Topo/etopo2.nc \
   ${CROCO_ROOT}/data/DATASETS_CROCOTOOLS/gshhs/GSHHS_shp/i/GSHHS_i_L1.shp \
   ${CROCO_ROOT}/data/DATASETS_CROCOTOOLS/TPXO7/TPXO7.nc \
   ${CROCO_ROOT}/data/DATASETS_CROCOTOOLS/RUNOFF_DAI/Dai_Trenberth_runoff_global_clim.nc
```
Used for: `Topo/etopo2.nc` + `gshhs/` → the grid; `TPXO7/` → tides (optional);
`RUNOFF_DAI/` → river-mouth positions and names for GloFAS rivers (optional).
Once extracted, you can delete `DATASETS_CROCOTOOLS.tar.gz`.

## Where things end up

| Path | What |
|---|---|
| `${CROCO_ROOT}/install/` | build scripts (the downloaded sources can be deleted after building) |
| `${CROCO_ROOT}/opt_seq/` | installed NetCDF/HDF5 stack — **what CROCO links against** |
| `${CROCO_ROOT}/code/`    | CROCO + croco_pytools |
| `${CROCO_ROOT}/data/`    | DATASETS_CROCOTOOLS (etopo2 + GSHHS) |

## Troubleshooting

- **`conda env create` fails with a YAML error** — the first line of
  `environment.yml` must be `name: ggosss26` (with a space after the colon).
- **`nf-config` shows a conda or `/usr` path, not `opt_seq`** — `conda deactivate`,
  re-`source env.sh`; it must put `opt_seq/bin` first on `PATH`.
- **A build step fails** — run them in order (each needs the previous); ensure
  `gfortran`/`gcc` exist (`which gfortran gcc`). Fix, then re-run the failed step.
- **`02` fails on curl or m4** — `sudo apt install libcurl4-openssl-dev m4`, re-run.
- **Out of memory / machine freezes during compile** — `export NJOBS=2` and re-run.
- **`04` didn't compile the Fortran tools** — `source env.sh; conda activate ggosss26`
  and run `make clean && make` in `code/croco_pytools/prepro/Modules/tools_fort_routines/`.

## Next

Once the checks above pass, build a region: see **`hindcast/README.md`** (quick start:
`hindcast/new_config.sh`, `hindcast/build_config.sh`, `hindcast/run_hindcast_cycle.sh`) or follow
`notebooks/GoG12_hindcast_user_guide.ipynb` step by step.

Optional accounts, needed only by the optional physics or data:
- **Copernicus Marine** (GLORYS, always): `copernicusmarine login` once.
- **CDS** (ERA5, always): `~/.cdsapirc` with `url: https://cds.climate.copernicus.eu/api` and your token; accept the ERA5 single-levels licence (the pressure-levels licence is only needed for `--legacy`).
- **EWDS** (GloFAS, only if `RIVERS=1`): `~/.ewdsapirc` with `url: https://ewds.climate.copernicus.eu/api` and your token; accept the GloFAS historical licence.
- **TPXO7** (tides, only if `TIDES=1`): already in `data/DATASETS_CROCOTOOLS/TPXO7/TPXO7.nc`.
