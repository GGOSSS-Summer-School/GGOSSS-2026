# croco_work — CROCO regional ocean hindcasts

A self-contained CROCO ocean-model system that reconstructs the past ocean of **any region** you choose:
- ocean initial and boundary conditions from the **GLORYS** reanalysis (CMEMS);
- atmospheric forcing from **ERA5** (CDS);
- optional **tides** (TPXO) and optional **rivers** (daily **GloFAS** discharge);
- **offline nesting**: a finer child configuration driven by a parent's own hindcast output.

The period is run as one or several cycles, each a short spin-up followed by the hindcast.
Each region is one small settings file; everything else is derived and written by the scripts.

| where | what |
|---|---|
| `install/` + `install/README.md` | one-time setup: NetCDF/HDF5 stack (`opt_seq/`), CROCO v2.1.3 + croco_pytools v2.0.4 (`code/`), conda env `ggosss26`, DATASETS_CROCOTOOLS (`data/`), accounts |
| `env.sh` | shared paths + compilers (sourced by every script) |
| `gtools/` | the Python tools; `gtools/ggosss26.py` is the command-line interface |
| **`hindcast/`** + **`hindcast/README.md`** | **the hindcast system**: `new_config.sh`, `build_config.sh`, `run_hindcast_cycle.sh`, `steps/`, `configs/<NAME>/domain.cfg` |
| `notebooks/GoG12_hindcast_user_guide.ipynb` | step-by-step guide from a bare machine to a finished hindcast, with expected outputs |
| `notebooks/verify_config.py` | checks the cycles of a finished hindcast |
| `CHANGES.md` | everything corrected or added relative to the original system, and the verification record |

```bash
source ~/croco_work/env.sh
~/croco_work/hindcast/new_config.sh MyRegion 4 12.5 -6 5.5 --res 12 [--tides] [--rivers]
~/croco_work/hindcast/build_config.sh MyRegion
nohup ~/croco_work/hindcast/run_hindcast_cycle.sh MyRegion > ~/MyRegion.log 2>&1 &
# a finer child nested in MyRegion (after MyRegion's hindcast has run)
~/croco_work/hindcast/new_config.sh MyChild 6 10 -4 2 --res 36 --parent MyRegion
~/croco_work/hindcast/build_config.sh MyChild --no-test && ~/croco_work/hindcast/run_hindcast_cycle.sh MyChild
```
