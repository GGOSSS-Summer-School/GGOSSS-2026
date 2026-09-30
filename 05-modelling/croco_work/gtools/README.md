# gtools — the CROCO Python package

Run the CLI as `python ${GTOOLS_DIR}/ggosss26.py <subcommand>` (or, after
`pip install --no-deps -e .` from the repo root, as `ggosss26 <subcommand>`).
The hindcast scripts in `hindcast/` call it for you.

```
gtools/
├── ggosss26.py       ← the CLI: forecast subcommands + *_hindcast subcommands
│                       (download_ocean_hindcast, download_atmosphere_hindcast, make_ini_hindcast,
│                        make_bry_hindcast, make_tides, download_rivers_hindcast, make_rivers_hindcast)
├── config/           ← make_grid_config.py (writes <configs>/<CONFIG>/grid.ini)
├── download/
│   ├── cmems.py         ← ocean (Mercator forecast / GLORYS reanalysis), map-oriented ARCO store, months in parallel
│   ├── gfs.py           ← atmosphere (GFS, for forecasts)
│   ├── ERA5/            ← atmosphere (ERA5, for hindcasts): request (parallel) + convert, ERA5_variables.json
│   └── glofas_rivers.py ← river discharge (GloFAS, EWDS) → per-river files for make_river_run.py
├── preprocess.py     ← builds ini / bry / tides / GFS forcing
├── croco_pytools/    ← VENDORED legacy croco_pytools (Modules/ + Readers/) used by preprocess.py
├── plotting.py, postprocess.py, validation*.py, animation.py   (optional analysis tools)
```

## Two croco_pytools, two jobs

- `gtools/croco_pytools/` (vendored, legacy API: `ibc_class`, `tides_class`, `interp_tools`,
  `Readers/ibc_reader.py`, `Readers/tides_reader.py`, `make_river_run.py`), originally from SAEON
  somisana-croco. It is imported by `preprocess.py` for `make_ini*`, `make_bry*` and `make_tides`,
  and it contains the fixes for domains that cross 0° (see `CHANGES.md`). Do not delete it.
- `code/croco_pytools/` (v2.0.4, installed by `install/04`) is used only for the grid
  (`prepro/make_grid.py` + its compiled Fortran `toolsf*.so`).

GLORYS reads through the `'mercator'` reader key (zos/thetao/so/uo/vo); there is no `'glorys'` key.
TPXO7 reads through `'tpxo7_croco'`.
