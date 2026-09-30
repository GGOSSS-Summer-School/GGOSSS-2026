# croco_work — corrections and additions

The system is self-contained under `~/croco_work`: root variable `CROCO_ROOT`, Python tools in
`gtools/` (variable `GTOOLS_DIR`, CLI `gtools/ggosss26.py`), conda env `ggosss26`. It depends on
nothing outside that folder except the conda env, the online data services (CMEMS, CDS, EWDS) and
the CROCO datasets downloaded into `data/` (install/README.md §3).

## Version 2 — any region, built automatically (see `hindcast/README.md`)

A configuration is now **one hand-written file**, `hindcast/configs/<NAME>/domain.cfg`, created by
`hindcast/new_config.sh NAME lon_min lon_max lat_min lat_max [--res N] [--start D] [--end D] [--tides] [--rivers]`.
`hindcast/lib/common.sh` derives everything else:
- the cpp key, folders, `RES`, the GLORYS and ERA5 boxes, and the time step (300 s at 1/12°, scaled with resolution and latitude);
- the GLORYS and ERA5 months, the number of cycles and the proof-run window;
- the open boundaries, from the land mask (`obc.cfg`, step 02).

Step 02 also writes both parameter files and the region card `CARD.md`. Step 05 reads `LLm0/MMm0`
from the grid. The steps are generic (`hindcast/steps/NN_*.sh NAME`). `hindcast/build_config.sh NAME`
runs them all, and `hindcast/run_hindcast_cycle.sh NAME` runs the cycled hindcast. Compiling removes conda
from the environment by itself, and the Python steps activate `ggosss26` by themselves.

**Tides and rivers are optional**: `TIDES=0/1` and `RIVERS=0/1` in `domain.cfg` (default 0). Each
combination is its own binary (`croco_plain`, `croco_plain_tides`, `croco_plain_rivers`,
`croco_plain_tides_rivers`). The plain build is always built and tested first, and the driver's
`--tides/--no-tides/--rivers/--no-rivers` flags override `domain.cfg`.

The v1 per-configuration scripts named below (`hindcast/configs/GoG_12/00_vars.sh` … `10_rivers_glofas.sh`,
`run_gog12.sh`, `notebooks/verify_gog12.py`) were replaced by `lib/common.sh`, `steps/01_make_grid.sh`
… `steps/10_rivers_glofas.sh` and `notebooks/verify_config.py NAME [BUILD]`. The fixes listed below are
all carried over.

## Offline nesting (new) — a finer child driven by a parent's own hindcast

- **`PARENT=<config>`** in `domain.cfg` (or `new_config.sh … --parent <config>`) makes a configuration a
  nested child. `lib/common.sh` checks that the child box lies inside the parent's, sets the ocean source to
  the parent (`OCEAN_TAG=NEST`, folder `downloaded_data/PARENT_<parent>`), and points `ERA5_DIR` to the
  parent's ERA5.
- **`new_config.sh --parent`** inherits the parent's start, tides, rivers and threads. It sets
  `HC_END` = parent end − 1 day, because the parent's daily means are stamped at 12:00 and must bracket the
  child's last step. It also sets a 1-day spin-up and a sponge about 7 child cells wide (`"20000. 400."`
  at 1/36°).
- **Step 03 of a child** converts instead of downloading. It takes the daily means of every parent cycle
  overlapping the period (spin-up + hindcast `croco_avg.nc`, only for runs that finished with `MAIN: DONE`),
  converts them with `gtools/nesting.py` (`croco_to_mercator`, new `parent_to_monthly`), and writes
  GLORYS-layout monthly files. After that, step 04 and the cycle driver work unchanged. Daily means are used,
  not hourly history, so the parent's tide is not aliased into the boundaries; the child makes its own tide.
- **The driver** (`run_hindcast_cycle.sh`) never downloads GLORYS/ERA5 for a child. It stops with a clear
  message if the converted parent data are missing.
- **Worked example:** `IGOG_36` (4–12.5°E, 5.5°S–5.5°N, 1/36°, 309×399 points, dt 100 s), nested in `GoG_12`. The southern edge was moved from 6°S to 5.5°S, 6 parent cells inside GoG_12's open southern boundary, so the child takes its boundary values from GoG_12's solution rather than from GoG_12's own GLORYS edge data. Both run
  1–15 June 2026, and GoG_12 runs one extra day (to 17 June 00:00) so the child's boundaries cover its last
  day. The GoG_12 hindcast uses GLORYS only. GoG_12's `HCAST_DAYS=16` makes the period one continuous run
  after a 2-day spin-up.
- The notebook was rewritten as a course. A primer comes first, and every code cell is preceded by a
  **What / Why / What to look for** block. It also has a new Part E (OpenMP check) and Part G (nesting).

## Robust GLORYS download: pieces + watchdog (`gtools/download/cmems.py`)

- **Measured:** GLORYS's time-series store uses chunks of 2081 days × 1 level × 16×16 points, and its map store
  uses 1 day × 1 level × 512×2048. So a 2-month request for the GoG_12 box transfers ~15× more data than it
  keeps, and takes 4–5 hours on a saturated 2.3 MB/s link.
- **Two failures seen:**
  1. A connection error restarted the single 1-hour request from zero.
  2. Later, a network drop left the toolbox waiting forever on a dead connection (0 B/s for over 7 minutes,
     no error, so no retry).
- **Pieces:** `download_cmems_monthly` now downloads 21 pieces: one variable × one depth band, plus `zos`.
  This costs no extra data, since each chunk holds one level of one variable. Finished pieces are kept, so a
  re-run resumes; they are merged and split into months at the end.
- **Watchdog:** each piece runs in a child process and is restarted after `CMEMS_STALL_MIN` (10) minutes
  without receiving any byte (`/proc/<pid>/io` `rchar`). The output file is not a usable signal: the toolbox
  writes it in bursts, sometimes only after 20+ minutes of transfer. Up to `CMEMS_MAX_RETRIES` (5) attempts,
  and `download_cmems`' own retries also went from 3 to 5, with a growing wait.

## Faster downloaders (adapted from the user's validation.py and era5_for_exercise.py)

- **CMEMS** (`gtools/download/cmems.py`): `download_cmems` now calls the `copernicusmarine.subset()`
  Python API in-process instead of launching the CLI through `os.system`. No ARCO store is forced
  (`CMEMS_SERVICE` still can force one). An existing file is kept only if it covers the requested dates at
  daily steps (`netcdf_covers_time_range`); a partial or monthly-mean leftover is re-downloaded. Each
  request gets 3 attempts, and the last failure raises, so the pipeline stops instead of continuing
  without data. Months are fetched one at a time by default (`CMEMS_N_PARALLEL=1`), because several
  in-process subsets writing NetCDF at once can crash netCDF4/HDF5. This supersedes the
  `arco-geo-series` default described under *Download speed* below.
- **ERA5 from the CDS ARCO Zarr store (new default, `era5_arco.py`)**, following ECMWF's
  `dss-notebooks` example `reanalysis-era5-single-levels/arco-access`.
  - *Why:* CDS API requests sat in the queue for over an hour each (Sep 2026). The ARCO store is read
    directly: t2m for the GoG_12 box (+2°), May → 22 Sep hourly, took 340 s. It uses the geo-chunked
    store (chunks ≈ 67 000 h × 4×4 points), which suits long series over a box.
  - *Missing fields:* the store has no `ssr` (net solar), no `lsm` and no `q`, so they are derived.
    CROCO uses the SSR file as the net shortwave flux: `online_bulk_var.F` labels it "NET SHORT
    WAVE", and `srflx = radsw`. So SSR = `fdir·(1−α_dir(μ)) + (ssrd−fdir)·(1−0.06)`, with the
    Taylor et al. (1996) open-water albedo `α_dir = 0.037/(1.1 μ^1.4 + 0.15)` and μ taken at
    mid-hour. Q comes from `d2m` + `msl` (Bolton 1980). LSM = 1 where SST is missing.
  - *Dependency:* needs `aiohttp>=3.9`. The env had aiohttp 0.21, which breaks fsspec's HTTP access;
    it is added to `environment.yml`. `ERA5_SOURCE=cds` switches back to the CDS API downloader.
- **Mercator: one request per dataset for the whole period** (`download_mercator_monthly`). The real
  cause of the slow downloads was the ARCO chunk layout, measured on `thetao` anfc:
  - `arco-geo-series` chunks are 1 day × 1 level × 512×2048 points. The GoG_12 box is 205×307 points,
    so the link received ~900 kB/s while only ~30 kB/s reached the file (2 days of `thetao` took 13.5 min).
  - `arco-time-series` chunks are 134 days × 2 levels × 32×64 points. Requesting the whole period
    (e.g. May–Sep, 153 days) in one go wastes very little, and the link runs at full speed (~2 MB/s).

  Each of the 4 datasets is now fetched once for all missing months through `arco-time-series`
  (`CMEMS_SERVICE` overrides it). The result is then split into the usual monthly `YYYY_MM.nc` files.
- **ERA5** (`gtools/download/ERA5/era5_fast.py`, the default of `download_atmosphere_hindcast`): one CDS
  request per variable and month, **single-levels product only**. Specific humidity `Q` is derived from
  2 m dewpoint + mean sea-level pressure (Bolton 1980). The legacy path asked for q at 1000 hPa on the
  pressure-levels product, a separate and much slower CDS queue. SST/LSM are requested once a day. The
  output is what CROCO's `online_get_bulk.F` reads: tags T2M Q U10M V10M TP SSR STRD SST LSM `msl`,
  time in **days** since Yorig, increasing latitude, missing value 9999, TP ×1000/3600, radiation /3600.
  The adaptation changed three things in the exercise script: it wrote `Q2M`, time in hours, and
  decreasing latitude. The raw file names are the legacy ones, so earlier raw downloads are reused.
  Files from a partial month are refreshed automatically. The old path remains available as
  `download_atmosphere_hindcast --legacy`.

## Additions for recent periods and speed (Sep 2026)

- **`OCEAN_SOURCE=mercator`** (`new_config.sh --ocean mercator`): GLORYS ended on 2026-06-23 (checked on
  2026-09-27), so the ocean can come from the **Mercator global analysis** (CMEMS anfc, 1/12°, daily,
  same variable names). Its 4 datasets per month (thetao, so, uo/vo, zos) are merged into `YYYY_MM.nc`
  (`downloaded_data/MERCATOR/`, files `croco_ini_MERCATOR_*`). The `'mercator'` reader is unchanged. CLI:
  `download_ocean_hindcast --source mercator`.
- **OpenMP** (`NTHREADS=8`, `--threads 8`): step 07 defines `OPENMP` + `SPLITTING_X × SPLITTING_ETA` in the
  config block (jobcomp adds `-fopenmp`), always rebuilds from a clean `Compile/`, and the runs export
  `OMP_NUM_THREADS`, `OMP_STACKSIZE=512M` and `ulimit -s unlimited`. `07_compile.sh NAME BUILD 1` gives a
  sequential binary for comparison.
- **GloFAS product type** `auto`: consolidated if it covers the whole period, else intermediate. For
  Jun–Sep 2026 only v4 *intermediate* is complete. Only the days EWDS offers are requested (the current month is partial).
- **ERA5 current month**: requested only up to 6 days ago (`ERA5_LAG_DAYS`); the converter skips months not
  released yet.
- **`--hcast-days`**: choose a cycle length that divides the period (Jun 1 → Sep 1 = 92 d = 23 × 4 d).

## Blocking errors fixed (each one stopped the workflow)

| File | Problem | Fix |
|---|---|---|
| `environment.yml` | `name:ggosss26` (no space) is invalid YAML → `conda env create` fails | `name: ggosss26` |
| `gtools/preprocess.py` (unchanged) needs `gtools/croco_pytools/prepro/{Modules,Readers}` | folder missing → `make_ini_hindcast`/`make_bry_hindcast` crash on `import ibc_class`; `code/croco_pytools` v2.0.4 has a different API (no `ibc_class`, `tides_class`, `Readers/`) | vendored the legacy croco_pytools (Modules + Readers, originally from SAEON somisana-croco) into `gtools/croco_pytools/` |
| `gtools/download/ERA5/` | `ERA5_variables.json` missing → `ERA5_request.py` / `ERA5_convert.py` crash | added (the croco_tools ERA5 variable table) |
| `gtools/config/utils/grid_utils.py` + `sf_utils.py` | `get_croco_env` required a `CROCO_INPUTS_ROOT` / a root variable that was not set → `make_grid_config.py` printed an error and wrote no `grid.ini` (exit code 0) | the fallback root is `CROCO_ROOT` with the croco_work sub-folders (`code/croco`, `hindcast/configs`, `hindcast/scratch`); `grid_utils` no longer asks for the unused inputs root; `make_grid_config.py` exits 1 on error |
| `GoG_12/01_make_grid.sh` | `envsubst` of a `grid.ini` template that is not in the config; hand-written grid.ini used `rectangular`, quoted values, `etopo5` reader for ETOPO2022, `shp_file = None` | `make_grid_config.py` (curvilinear, etopo2 + GSHHS from DATASETS_CROCOTOOLS) → `make_grid.py` |
| `GoG_12/03_download_data.sh`, `04_make_ini_bry.sh` | `cd "${CROCO}` — unterminated quote (syntax error) and undefined variable | `cd "${GTOOLS_DIR}"` |
| `GoG_12/05_patch_source_files.sh` | `${CROCO_ROOT_ROOT}` typo; `require_one_match "LLm0="` always aborts (param.h has many); global seds rewrite every config's LLm0/MMm0/N; **OBC switches never edited** (all 4 boundaries open) ; LLm0/MMm0 hard-coded 270/168 | edits anchored to the block under the config name; `# elif defined GOG_12` inserted above `# else`; LLm0/MMm0 read from `croco_grd.nc` (= **271/169**, grid 273×171) |
| `GoG_12/06_patch_croco_in.sh` | replaced header lines and left the old value lines (corrupted sections); NRREC=0; file names `croco_ini_GLORYS.nc`/`croco_bry_GLORYS.nc`/`croco_blk.nc` that are never produced; ERA5 path on the numbers line; sponge `XXX` left; NDTFAST 30 | rewrite the line(s) below each header; real ini/bry names; NRREC=1; `online` numbers line + path line; sponge `0. 0.`; 2016 300 60 1 |
| `GoG_12/run_gog12.sh` | uses `screen` (not installed); compiles with conda active; ran in `CROCO_FILES` | replaced by `07_compile.sh` (refuses to run inside conda, checks `nf-config`) and `08_run_test.sh` (checks `MAIN: DONE`) |
| `hindcast/run_hindcast_cycle.sh` | missing (no cycled hindcast) | added: spin-up + hindcast cycling driver with automatic data download and `MAIN: DONE` checks (18 cycles from 2018-06-01 for GoG_12) |

## Consistency errors fixed

- `GoG_12/00_vars.sh` set `CROCO_RUNS_ROOT=hindcast/model-runs` and hard-coded roots, overriding `track.sh` (`hindcast/scratch`); `verify_gog12.py` read `hindcast/scratch` → the build and the checks pointed at different trees. Now `00_vars.sh` uses `track.sh`; build in `hindcast/scratch/GoG_12`, cycle outputs in `hindcast/model-runs/GoG_12/<YYYYMMDD>/`.
- Run period started on the 1st of a month with only Jun–Aug downloaded and the bry ended exactly on the run end (no +1 day pad). Now GLORYS 2018-05…2018-09 and ERA5 2018-05…2018-08, bry windows padded ±1 day.
- `ggosss26.py download_ocean_hindcast` defaulted to the monthly-mean product (P1M-m); default is now the daily P1D-m (the scripts still pass it explicitly).
- `make_ini_hindcast` read only the ini date's month; it now also reads the neighbour month (±1 day) when present.
- ERA5 converter/request launched with bare `python`; now `sys.executable` (the active env). `ERA5_request.py` skips raw files already downloaded (safe re-runs) and writes to a `.part` file first.
- `download/cmems.py` retry loop re-raised on the first failure (never retried) — fixed.
- `pyproject.toml` entry point `gtools.cli:main` (no such module) → `gtools.ggosss26:main`.
- `install/00_download_libraries.sh` defaulted to `~/ggosss26` instead of `~/croco_work`.
- `install/README.md` referred to a `croco_py` env and a non-existent `notes/README_parallel.md`.
- `hindcast/track.sh` now refuses to run if `env.sh` was not sourced.
- `02_check_grid.py` derives the open boundaries from the mask (instead of printing a fixed expectation) and adds the CFL check.
- `crocotools_param.py` comment pointed to a non-existent reader file.

## Download speed

Measured on this machine (~4.4 MB/s link):

- **ERA5**: `ERA5_request.py` sent the 10 requests per month one after another. Each one waits 1–3 min
  in the CDS queue before a ~12 MB transfer, so the wait was mostly queue time (~30 min per month).
  The requests are now submitted **in parallel** (`ERA5_N_PARALLEL`, default 2). CDS enforces a per-user,
  per-dataset limit on queued requests; above it, jobs are rejected with *"Number queued requests for this
  dataset is temporarily limited"*. 6 in parallel was rejected, 2 works, so the gain is about 2×.
  Rejected or failed requests are retried (up to 8 times, with a growing wait). Finished files are skipped on re-run.
- **GLORYS**: one month = one large subset (50 levels × 5 variables × 30/31 days over the 25.5°×17° box).
  This is bandwidth-bound. Months are now fetched in parallel (`CMEMS_N_PARALLEL`, default 2), and
  finished months are skipped.
- **GLORYS storage format**: left alone, copernicusmarine picked the *time-chunked* ARCO store (made for
  point time series). For a full-depth regional map it pulls far more data. Same request (2 days, full
  depth, GoG_12 box): time-chunked reached 2 % in 6 min (estimated hours); `--service arco-geo-series`
  finished in 19 min (49 MB) on a slower link. `download_cmems` now forces `arco-geo-series`
  (`CMEMS_SERVICE` overrides it; empty = toolbox choice).
- **Why earlier downloads felt faster**: the CLI default used to be the *monthly-mean* product (1 record per
  month, ~30× less data). That is unusable for a hindcast (make_ini/bry need daily records), so the daily
  `P1D-m` product is now the default.
- `gtools/download/era5_for_exercise.py` is **not** used: it is sequential, says in its header that it is an
  unverified reconstruction, and writes `Q2M`, which CROCO v2.1.3 does not read (it reads `Q`).
- Connection drops or machine sleep break long downloads. Re-run the same command; it resumes at the
  first missing file or month.

## Added: Tides (TPXO) and Rivers from GloFAS

| File | What |
|---|---|
| `crocotools_param_tides.py` (written by `steps/02_check_grid.sh`) | TPXO7 parameters (`inputdata='tpxo7_croco'`; `'tpxo7'` alone is **not** a key of the vendored `tides_reader.py`) |
| `hindcast/steps/09_tides.sh` | make_tides in its own gen dir, checks (M2 amplitude, 0 fill values); `TIDES` is set per build by step 07, `forcing:` → croco_frc.nc by step 06, hourly his / daily avg by step 08 and the driver |
| `gtools/download/glofas_rivers.py` + CLI `download_rivers_hindcast`, `make_rivers_hindcast` | GloFAS v4 daily discharge from EWDS → per-river text files + `river_list.txt` |
| `hindcast/steps/10_rivers_glofas.sh` | GloFAS → `make_river_run.py` (real dates, Yorig, no cycle) → `croco_runoff.nc`; the `psource_ncfile` block (T = `RIVER_TEMP`, S = `RIVER_SALT`); `PSOURCE`/`PSOURCE_NCFILE` set per build by step 07, `_TS` off |
| `steps/07_compile.sh NAME <build>` / `steps/08_run_test.sh NAME <build>` | builds and runs `croco_plain`, `croco_plain_tides`, `croco_plain_rivers`, `croco_plain_tides_rivers`; test outputs kept apart in `test_<build>/` |
| `hindcast/run_hindcast_cycle.sh NAME [--tides] [--rivers]` | per-cycle tide file at the spin-up start, staged into both phases; runoff staged; output retimed; folders `<date>_<build>` |

**GloFAS request schema.** Since 29 Jul 2026 (GloFAS v5 release) EWDS requires `year/month/day`
(not `hyear/hmonth/hday`), a mandatory `timespan` (`'time_mean'` for discharge) and the variable
`average_river_discharge_in_the_last_24_hours`. Checked against the live EWDS form and constraints:
2018-05…09 are complete for v4 consolidated. `glofas_rivers.py`, and your older `glofas.py` and
`download_glofas_africa_clim.py`, now use this schema. The old keys are still in croco_pytools
v2.0.4 `download_glofas_river.py`, so don't use that one. The other fixes in `glofas.py`:
`~/ggosss26` paths, the wrong merge command, days per month.

### Bugs found and fixed in the vendored croco_pytools (triggered because GoG_12 crosses 0°)

1. `tides_class.var_periodicity` / `handle_periodicity` (same code in `ibc_class`): for a domain crossing
   0° in a 0–360 atlas (TPXO), the eastern part copied `imax` columns instead of `imax+1`. The last
   column was left at lon = 0 with data = 0, so the longitude axis was non-monotonic and a false zero-tide
   strip was treated as valid data. Result: 669 ocean cells with no tidal currents, for every wave.
2. Same branch, single wave: squeeze dropped the period axis → `IndexError: too many indices`.
3. `tides_class`: TPXO transport ÷ depth gives ±inf where the atlas depth is 0 next to the coast; these
   are now NaN, so they are gap-filled like land.
4. `preprocess.make_tides`: the period match `abs(ΔT) < 1e-4 h` rejected Mf and Mm (TPXO7 stores
   327.8599 / 661.31 h vs tides.txt 327.858980 / 661.309208). It now uses a relative tolerance of 1e-5.
   `input_dir` also gets its trailing `/`.

## Verification record (WSL Ubuntu, CROCO v2.1.3, croco_pytools v2.0.4, opt_seq NetCDF)

Done in an isolated copy (`/tmp/cwv`), reusing your built `opt_seq` and CROCO source read-only:

- ✅ `bash -n` on every script; `py_compile` on every changed Python file.
- ✅ **Version 2 automation, three regions built from `new_config.sh` alone** (steps 01, 02, 05, 07; plus 06/09/10 for GoG_12):

  | region | derived dt | grid → LLm0 × MMm0 | boundaries from the mask | compile |
  |---|---|---|---|---|
  | GoG_12 (−10…12.5°E, −6…8°N, 1/12°) | 300 s | 273×171 → 271×169 | S, W open; E, N closed | plain + plain_tides_rivers OK |
  | IGOG_25 (4…12.5°E, −6…5.5°N, 1/25°) | 144 s | 215×291 → 213×289 | S, W open; E, N closed | plain OK |
  | Canary_12 (−22…−15.5°E, 14…24°N, 1/12°) | 300 s | 81×123 → 79×121 | S, W, N open; E closed | plain OK |

  The guards give clear messages for: a missing domain.cfg, a bad name, an inverted box, an existing config, an unknown option, a missing binary.
- ✅ Step 1 grid: xi_rho × eta_rho = 273 × 171 → LLm0 = 271, MMm0 = 169. Mask: south 272/273 and
  west 146/171 ocean → open; north 1/273 and east 1/171 → closed. rx0 = ry0 = 0.200; barotropic CFL 0.13.
- ✅ Step 8 patches = exactly the documented edits, in the GOG_12 block only; `cpp` gives LLm0=271, MMm0=169, N=50.
- ✅ Compile: `CROCO is OK` for `croco_plain`, `croco_plain_tides`, `croco_plain_rivers`, `croco_plain_tides_rivers`.
- ✅ Tides: 10 waves, 8 variables, **0** fill values on ocean cells. The generated M2 amplitude divided by
  raw TPXO7 is 1.0237 at every test point, west and east of 0° (= the M2 nodal factor for June 2018).
  M2 amplitude 0.445 m on average, 0.759 m at most.
- ✅ Rivers plumbing (selection → text files → `make_river_run.py` → `croco_runoff.nc` → `psource_ncfile`
  block), tested with **synthetic** discharge in the GloFAS file layout: 5 rivers placed on wet cells,
  qbar_time 2018-05-01 12:00 → 2018-09-30 12:00 in days since 1993, real dates.
- ✅ Downloads: ERA5 parallel requests (2 at a time) completed May 2018 and 9/10 files of June 2018
  (the last one was pending when the test was stopped). The GLORYS `arco-geo-series` test subset completed.
- ⚠️ **Not yet run** (stopped at your request while the link was down to ~0.1 MB/s): a full GLORYS month,
  the ERA5 conversion, `make_ini_hindcast` / `make_bry_hindcast` on real GLORYS, the 7-day `08_run_test.sh`
  runs, a real GloFAS download (needs `~/.ewdsapirc`), and a driver cycle. Run them with the notebook;
  each step prints the checks to compare against.
