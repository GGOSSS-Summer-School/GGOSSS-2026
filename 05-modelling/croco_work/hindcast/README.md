# The hindcast track — build and run any region with (almost) no hand-editing

This folder builds and runs a **CROCO regional hindcast** (GLORYS ocean + ERA5 atmosphere, Yorig 1993),
with optional **tides** (TPXO) and **rivers** (discharge from **GloFAS**), for any
region you choose, and finer **nested children** inside it (offline nesting, section 7). You write **one small file per region** (`domain.cfg`, generated for you). Every other
number is derived from it, and every file is written by a script.

> Step-by-step tutorial with explanations and expected outputs: `notebooks/GoG12_hindcast_user_guide.ipynb`.
> One-time machine setup (NetCDF, CROCO, conda env, datasets): `install/README.md`.
> Every change relative to the original system: `CHANGES.md`.

---

## 1. Quick start

```bash
source ~/croco_work/env.sh                     # (the scripts also do this themselves)

# 1. describe the region: name, box (deg), options
~/croco_work/hindcast/new_config.sh GoG_12 -10 12.5 -6 8 --start 2026-06-01 --end 2026-06-17 \
        --hcast-days 16 --tides --rivers --threads 8                              # 1/12 deg, one 16-day run
~/croco_work/hindcast/new_config.sh Canary_12 -22 -15.5 14 24 --start 2025-12-02 --end 2026-01-09
# a finer CHILD nested in GoG_12 (period, physics, threads from the parent; see section 7)
~/croco_work/hindcast/new_config.sh IGOG_36 4 12.5 -5.5 5.5 --res 36 --parent GoG_12

# 2. build + prove it (grid, boundaries, data, ini/bry, source files, compile, 7-day test run)
~/croco_work/hindcast/build_config.sh GoG_12

# 3. run the whole cycled hindcast (hours: run it in the background)
nohup ~/croco_work/hindcast/run_hindcast_cycle.sh GoG_12 > ~/GoG_12_hindcast.log 2>&1 &

# 4. check it
python3 ~/croco_work/notebooks/verify_config.py GoG_12
```

You don't need `conda activate` or `conda deactivate`. The scripts activate the `ggosss26` env for
the Python tools and remove conda from the environment for the compile (conda's NetCDF would break the link).

---

## 2. What you decide, and what is derived

### You decide — `hindcast/configs/<NAME>/domain.cfg`

`new_config.sh` writes it. Every line is commented, and you can edit it at any time.

| setting | meaning | default |
|---|---|---|
| `CONFIG_NAME` | name of the region (folder names; upper-cased it becomes the cpp key, e.g. `GOG_12`) | argument |
| `LON_MIN LON_MAX LAT_MIN LAT_MAX` | the grid box in degrees (longitude −180..180) | arguments |
| `RES_INV` | resolution = 1/RES_INV degree (12 → 1/12°, ≈ 9 km) | 12 (`--res`) |
| `HC_START HC_END` | start of the period / end of the period (00:00 of that day; `2026-09-01` covers all of August) | 2018-06-01 / 2018-08-30 (`--start`, `--end`) |
| `SPINUP_DAYS HCAST_DAYS TEST_DAYS` | spin-up and hindcast length per cycle (choose `HCAST_DAYS` dividing the period, e.g. 92 d = 23 × 4); proof-run length | 2, 5 (`--hcast-days`), 7 |
| `THETA_S THETA_B N_LEVELS SIGMA_HC` | vertical grid | 7, 2, 50, 200 |
| `TIDES` | 1 = TPXO tides | 0 (`--tides`) |
| `RIVERS` | 1 = GloFAS rivers | 0 (`--rivers`) |
| `RIVER_TEMP RIVER_SALT RIVER_QMIN RIVER_RADIUS RIVER_MARGIN RIVER_EXTRA` | river water T/S, selection threshold, search radius, extra mouths | 27, 0, 100 m³/s, 0.25°, 0°, none |
| `OCEAN_SOURCE` | `glorys` = GLORYS12 daily reanalysis; `mercator` = Mercator global analysis (anfc, 1/12°, same variables), for periods GLORYS does not cover yet (it ends a few months before the present, e.g. 2026-06-23 as of Sep 2026) | glorys (`--ocean`) |
| `NTHREADS` | > 1 = build and run CROCO with OpenMP on that many threads (≤ `nproc`); tiles `SPLITTING_X × SPLITTING_ETA` chosen automatically | 1 (`--threads`) |
| `GLOFAS_PRODUCT_TYPE` | `auto` = consolidated if EWDS offers it for the whole period, else intermediate (recent months); or force one | auto |
| `PARENT` | offline nesting: the parent configuration whose hindcast output replaces GLORYS as this configuration's ocean (section 7) | empty (`--parent`) |
| `PARENT_RUN` | *optional*: one parent cycle folder to nest in | empty = every parent cycle overlapping the period |
| `OBC` | *optional* override of the open boundaries, `"south west east north"` as 1/0 | empty = from the mask |
| `DT NDTFAST SPONGE` | *optional* overrides | automatic / 60 / `0. 0.` |

### Derived automatically (`lib/common.sh`, printed at the start of every build/run)

| derived | rule | example GoG_12 |
|---|---|---|
| cpp key `CONFIG_CPP` | name upper-cased, non-alphanumerics → `_` | `GOG_12` |
| grid spacing `RES` | `1/RES_INV` at full precision (`bc -l`) | 0.08333… |
| GLORYS download box `EXTENTS` | grid box ± 1.5° | −11.5, 14, −7.5, 9.5 |
| ERA5 box `ERA5_BOX` | the grid box (the downloader adds 2°) | −10, 12.5, −6, 8 |
| time step `DT` | 300 s × 12/RES_INV, × cos(lat)/cos(35°) poleward of 35°, rounded down to a divisor of 3600 s | 300 s (1/25° → 144 s) |
| GLORYS months | from (start − spin-up − 2 d) to (end + 2 d) | 2026-05 … 2026-06 |
| ERA5 months | from (start − spin-up) to end | 2026-05 … 2026-06 |
| number of cycles | ⌈ period / HCAST_DAYS ⌉ | 1 (1–17 June, HCAST_DAYS 16) |
| proof-run window | day 2 of the period, TEST_DAYS long | 2026-06-02 → 06-09 |
| **open boundaries** | step 02: an edge is open if > 50 % of it is ocean in the land mask → `obc.cfg` | S W open, E N closed |
| **grid size** `LLm0 MMm0` | step 05: `xi_rho − 2`, `eta_rho − 2` read from `croco_grd.nc` | 271, 169 |
| `crocotools_param.py` | step 02: `obc_dict` + `sigma_params` written for you | — |
| `crocotools_param_tides.py` | step 02: TPXO7 parameters | — |
| `cppdefs.h`, `param.h`, `jobcomp` | step 05: the documented edits | — |
| `croco.in` | step 06: title, time step, S-coord, files, sponge, ERA5 `online` block, rivers block | — |
| binary | step 07: `croco_<build>`, build = `plain[_tides][_rivers]` | `croco_plain` |
| region card | step 02: `configs/<NAME>/CARD.md` (box, grid, boundaries, dt, CFL, data, period) | — |

The consistency rules of a CROCO configuration therefore hold **by construction**. The config name is the same in
`cppdefs.h` and `param.h`, the grid size comes from the file, the vertical grid comes from one place, and
the boundaries come from one place: the mask → `obc.cfg` → `obc_dict`, `OBC_*`, and the boundary file contents.

---

## 3. Folder layout

```
hindcast/
├── new_config.sh            create configs/<NAME>/domain.cfg
├── build_config.sh          run steps 01-10 for a configuration (logs in scratch/<NAME>/build_logs/)
├── run_hindcast_cycle.sh    the cycled hindcast for a configuration
├── track.sh                 CROCO_CONFIGS_ROOT / CROCO_RUNS_ROOT of the hindcast track
├── lib/common.sh            loads domain.cfg and derives everything (read it: it is short)
├── steps/                   the individual steps, each: steps/NN_xxx.sh <NAME> [...]
│   ├── 01_make_grid.sh          grid.ini + croco_grd.nc
│   ├── 02_check_grid.sh         mask -> obc.cfg, params, CARD.md, CFL check
│   ├── 03_download_data.sh      GLORYS + ERA5 for the whole period
│   ├── 04_make_ini_bry.sh       ini + bry for the proof run
│   ├── 05_patch_source_files.sh cppdefs.h / param.h / jobcomp
│   ├── 06_patch_croco_in.sh     croco.in
│   ├── 07_compile.sh            croco_<build>, conda removed automatically
│   ├── 08_run_test.sh           proof run -> scratch/<NAME>/test_<build>/
│   ├── 09_tides.sh              OPTIONAL  croco_frc.nc from TPXO
│   └── 10_rivers_glofas.sh      OPTIONAL  croco_runoff.nc from GloFAS
├── configs/<NAME>/          the RECIPE: domain.cfg (yours) + generated grid.ini, obc.cfg,
│                            crocotools_param*.py, cppdefs.h, param.h, croco.in, jobcomp, CARD.md
├── scratch/<NAME>/          the WORKBENCH: CROCO_FILES/ (grid, ini/bry, frc, runoff), downloaded_data/
│                            (GLORYS, ERA5, GLOFAS), croco_<build> binaries, test_<build>/ proof runs
└── model-runs/<NAME>/       the RESULTS: <YYYYMMDD>[_<build>]/{spinup,hcast,gen_*}/ per cycle
```

---

## 4. The steps in detail

`build_config.sh NAME` runs them in this order and stops at the first failure. After fixing it, continue
with `--from N`. `--only N` runs one step. `--no-download` skips 03, and the driver then downloads what it
needs. `--no-test` skips the proof runs.

1. **Grid (01).** `make_grid_config.py` writes `grid.ini`: a curvilinear grid with `npoints = span × RES_INV + 1`,
   ETOPO2 bathymetry and GSHHS coastline from `DATASETS_CROCOTOOLS`. `make_grid.py` then builds the mask
   and smooths the topography to rx0 ≤ 0.2. The grid gets 2 extra rows, so 1/12° over 22.5° × 14° gives 273 × 171.
2. **Boundaries and parameters (02).** Prints the mask along each edge, isolated water bodies, rx0 and the
   barotropic and baroclinic Courant numbers. It writes `obc.cfg`, `crocotools_param.py`,
   `crocotools_param_tides.py` and `CARD.md`. **Look at its output once for every new region**: if an edge
   is borderline (40–60 % ocean) or the Courant number is above 0.7, set `OBC` or `DT` in `domain.cfg` and
   re-run from 02 (or 05 for DT).
3. **Data (03).** GLORYS daily reanalysis (`cmems_mod_glo_phy_my_0.083deg_P1D-m`; the monthly-mean product
   cannot be used). It is fetched with the `copernicusmarine.subset()` Python API from the time-chunked
   ARCO store (`arco-time-series`; `CMEMS_SERVICE` overrides it), for all missing months at once, as
   **21 pieces**: one variable × one depth band (0–50, 50–250, 250–1100, 1100–3400, 3400–6000 m), plus
   sea level. Every chunk of that store holds one level of one variable, so the pieces transfer no more than
   one big request would. A network failure repeats one piece (minutes), not the whole request (hours), and
   a re-run skips the pieces already on disk. A **watchdog** runs each piece in a child process and restarts
   it after `CMEMS_STALL_MIN` (10) minutes without receiving data. This matters because after a network drop
   the toolbox can hang forever on a dead connection. Pieces get up to `CMEMS_MAX_RETRIES` (5) attempts.
   The pieces are then merged and split into `YYYY_MM.nc` (a month at the end of GLORYS holds the days that
   exist). GLORYS chunks are 2081 days long, so even this layout moves far more data than it keeps: expect
   hours for a month or two over a GoG_12-size box. By default ERA5 is read directly from the **CDS ARCO Zarr store** by
   `gtools/download/ERA5/era5_arco.py` (`ERA5_SOURCE=arco`). There is no request queue: one read per
   variable for the whole period, from the geo-chunked store, 3 variables at a time. The token comes
   from `~/.cdsapirc`. Fields the store lacks are derived:
   - `Q` from dewpoint + mean sea-level pressure (Bolton 1980);
   - `SSR`, the NET solar flux that CROCO uses as `srflx`, from downward `ssrd` and direct `fdir`
     with the Taylor et al. (1996) open-water albedo;
   - `LSM` from where SST is missing.

   The store lags real time by about 5 days; the last complete day is used. `ERA5_SOURCE=cds` uses
   `era5_fast.py` instead: CDS API requests, single-levels product only, 2 in parallel. The legacy
   path (`--legacy`) also asks the slow pressure-levels queue for humidity. The files are converted to CROCO online files
   `<VAR>_YyyyyMmm.nc` (T2M, Q, U10M, V10M, TP, SSR, STRD, SST, LSM, `msl`; time in days since Yorig,
   latitude increasing, missing value 9999). Re-running is safe: finished files are skipped. This step can take hours.
   With `OCEAN_SOURCE=mercator` the ocean comes from the Mercator global analysis. It uses 4 CMEMS
   datasets (thetao, so, uo/vo, zos), **each requested once for the whole period** from the time-chunked
   ARCO store (`arco-time-series`, chunks of 134 days × 32×64 points: small waste for a long request,
   full link speed). The map store wastes 17–60× on a regional box. The data are split and merged into
   the same monthly `YYYY_MM.nc` files (folder `downloaded_data/MERCATOR`,
   files `croco_ini_MERCATOR_*`). For the current month ERA5 is requested only up to 6 days ago (its release
   delay, `ERA5_LAG_DAYS`), and GloFAS only for the days EWDS offers.
4. **Initial and boundary conditions (04)** for the proof run. The initial condition is at day 2. The
   boundary window is the run ± 1 day, because CROCO needs a record bracketing every step. Only open edges get data.
5. **Source files (05).** `cppdefs.h`: config name, the OBC block right under it (never the other
   configurations' blocks), `ONLINE` + `ERA_ECMWF` on, `AROME` off. `param.h`: an `# elif defined <KEY>` line
   with the sizes from the grid. `jobcomp`: an absolute `SOURCE1`.
6. **croco.in (06).** Everything from `domain.cfg`. The ERA5 `online` block uses real dates and spans
   months (`bmonthend` = last month). Sponge is off by default. If a runoff file exists, the rivers block is added.
7. **Compile (07)** a build: `plain`, `plain_tides`, `plain_rivers` or `plain_tides_rivers`. The TIDES /
   PSOURCE switches are set for that build, so one recipe gives every variant.
8. **Proof run (08).** 7 days. It must end with `MAIN: DONE`, and the first `time[DAYS]` must equal the days
   since `YORIG`-01-01. Outputs go to `test_<build>/`. With tides it prints the sea-level period (≈ 12.4 h).
9. **Tides (09, optional).** See section 5.
10. **Rivers (10, optional).** See section 6.

**The cycled hindcast** (`run_hindcast_cycle.sh NAME`). Each cycle T is a spin-up (T − 2 d → T, initial
condition from GLORYS) followed by a hindcast (T → T + 5 d, initial condition = the spin-up restart). The
cycles tile the period with no gap or overlap, and each one starts again from GLORYS, so the run stays
anchored to the reanalysis. Missing data months are downloaded on the fly. Options: `--tides/--no-tides`,
`--rivers/--no-rivers` (default: `domain.cfg`), `--start`, `--ncycles`. A cycle that already finished is
skipped when you re-run.

---

## 5. Tides (optional — `TIDES=1`)

* **What:** TPXO7 harmonics (10 waves: M2 S2 N2 K2 K1 O1 P1 Q1 Mf Mm), for elevation, currents and potential,
  with nodal corrections. They are written to `croco_frc.nc`. CROCO switches: `TIDES` (plus `SSH_TIDES`,
  `UV_TIDES`, `POT_TIDES` and `TIDERAMP`, already set in the template block).
* **Data:** `DATASETS_CROCOTOOLS/TPXO7/TPXO7.nc`. It comes with the CROCO datasets, so there's nothing to download.
* **The date matters:** phases are referenced to Yorig and the nodal corrections are computed for the run's
  start date, so every run or cycle needs its own tide file. `08` builds one for the proof run; the driver
  builds one per cycle at its spin-up start.
* **Output:** hourly history (to resolve M2) and daily averages (a full tidal day averages the tide out).
  Compare daily means only with GLORYS, which has no tide.
* **Without tides:** `TIDES=0`. Nothing tide-related is built or read, and the binary is `croco_plain`.
* Domains that cross 0° (like GoG_12) triggered three bugs in the vendored tide reader and one in the period
  matching (Mf, Mm). They are fixed here; step 09 checks that there are 0 fill values on ocean cells.

## 6. Rivers from GloFAS (optional — `RIVERS=1`)

* **What:** daily river discharge from **GloFAS v4** (`cems-glofas-historical`, consolidated) for the real
  dates of your period, instead of a climatology such as Dai & Trenberth. CROCO switches: `PSOURCE` + `PSOURCE_NCFILE`.
  `PSOURCE_NCFILE_TS` stays off (the file carries flow only; T/S come from `croco.in`).
* **Once:** GloFAS is served by **EWDS**, not CDS. Same ECMWF account, different store and licence:
  ```bash
  printf "url: https://ewds.climate.copernicus.eu/api\nkey: <your ECMWF token>\n" > ~/.ewdsapirc
  chmod 600 ~/.ewdsapirc
  # accept the licence: https://ewds.climate.copernicus.eu/datasets/cems-glofas-historical?tab=download
  ```
  Request schema since 29 Jul 2026: `year/month/day`, `timespan='time_mean'`,
  `variable='average_river_discharge_in_the_last_24_hours'`. The old `hyear/hmonth/hday` requests may fail.
  `GLOFAS_VERSION=version_5_0` selects the pre-operational v5.
* **How rivers are chosen:** the mouth positions and names come from the Dai & Trenberth file (its discharge is
  not used). For each mouth inside the grid, the GloFAS cell with the highest mean discharge within
  `RIVER_RADIUS` is its outlet. Rivers below `RIVER_QMIN` are dropped, and two mouths that land on the same
  channel are counted once. `make_river_run.py` then puts each river on the nearest wet CROCO cell and writes
  `croco_runoff.nc` with real dates in days since Yorig. Add a missing river with `RIVER_EXTRA` (a file with
  `name lon lat` lines).
* **Water properties:** `RIVER_TEMP` (27 °C, tropical) and `RIVER_SALT` (0; the builder's default of 5 PSU would
  inject brackish water).
* **Without rivers:** `RIVERS=0`. Nothing is downloaded or built, and the binary has no `PSOURCE`.
* **Rivers at the edge of the domain:** a mouth just outside the grid is dropped (`RIVER_MARGIN=0`). Its
  water already enters through an open boundary with GLORYS. For GoG_12 this is the Congo, on the southern
  edge; add it only deliberately.

---

## 7. Offline nesting — a finer child inside a parent

**Idea.** GoG_12 downscales GLORYS. A nested child downscales *GoG_12*: its initial state and open boundaries
come from the parent's output instead of GLORYS. Everything else is an ordinary configuration. It is *offline*:
the parent runs first, then the child reads the parent's saved output (the models do not exchange data while
running, unlike AGRIF).

```bash
# 1. the parent, as usual (its period must cover the child's + 1 day)
~/croco_work/hindcast/build_config.sh GoG_12 && ~/croco_work/hindcast/run_hindcast_cycle.sh GoG_12
# 2. the child: one command describes it ...
~/croco_work/hindcast/new_config.sh IGOG_36 4 12.5 -5.5 5.5 --res 36 --parent GoG_12
# 3. ... the grid (step 01), then the usual steps; step 03 CONVERTS the parent instead of downloading
~/croco_work/hindcast/build_config.sh IGOG_36 --no-test
~/croco_work/hindcast/run_hindcast_cycle.sh IGOG_36
```

**What `--parent` sets** (all in the child's `domain.cfg`, all editable):

| setting | value | reason |
|---|---|---|
| `PARENT` | `GoG_12` | switches the ocean source to the parent (`OCEAN_SOURCE` is then ignored) |
| box | checked to lie inside the parent's | the parent must cover every child point |
| `HC_START`, `HC_END` | the parent's start; the parent's end − 1 day | the parent's daily means are stamped at 12:00, so its last one (end − 12 h) must come after the child's last step |
| `HCAST_DAYS` | the parent's, capped at the child's period | one cycle per parent cycle |
| `SPINUP_DAYS` | 1 | the interpolated parent state has no fine structure; the child grows its own in about a day |
| `SPONGE` | ≈ 7 child cells wide, 400 m² s⁻¹ (1/36° → `"20000. 400."`) | absorbs the parent's sharp structures at the open edges instead of reflecting them |
| `TIDES`, `RIVERS`, `NTHREADS` | the parent's | the child makes its own tide file (TPXO) and river file on its grid |
| `DT` | derived as usual (1/36° → 100 s) | CFL: a 3× finer grid needs a 3× shorter step |

**Step 03 of a child** (`steps/03_download_data.sh`) downloads nothing:
* it finds the parent cycles overlapping the child period (`hindcast/model-runs/<PARENT>/<date>_<build>/`),
  checks they finished (`MAIN: DONE`), and converts their **daily means** (`spinup/` and `hcast/`
  `croco_avg.nc`) with `gtools/nesting.py`: sigma levels → the 50 standard GLORYS depths, currents rotated to
  east/north, names `thetao so uo vo zos`. The result is written as monthly `YYYY_MM.nc` files in
  `scratch/<CHILD>/downloaded_data/PARENT_<PARENT>/`, the same layout as the GLORYS folder. Step 04 and the
  cycle driver then build the child's ini/bry exactly as from GLORYS (`crocotools_param.py` keeps
  `inputdata='mercator'`; the files carry the NEST tag: `croco_ini_NEST_…`, `croco_bry_NEST_…`);
* daily means, not the hourly history: averaging removes the parent's tide. The child adds the tide itself
  (TPXO) at its boundaries, so hourly snapshots taken once a day would **alias** the tide into the boundary data;
* ERA5: the parent's `for_croco` files are used as they are (`ERA5_DIR` points to the parent's). CROCO
  interpolates ERA5 onto the running grid (`ONLINE`), so no conversion is needed; nesting refines the ocean,
  not the atmosphere.

**Rules** (checked or set by the scripts): child box inside the parent; the parent's Yorig (1993) for the
conversion; child `DT` ≈ parent `DT` × (child cell / parent cell); child period inside the parent period, ending
at least one day earlier. Keep the child a few parent cells away from the parent's **open** edges when you can:
there the parent's solution is mostly its own boundary data. That is why IGOG_36 stops at 5.5°S: its southern
edge lies 6 parent cells (~55 km) inside GoG_12's open southern boundary (6°S), so the child's boundary values
come from GoG_12's own solution, not from the GLORYS data GoG_12 receives at its edge. Its eastern edge
(12.5°E) coincides with the parent's, but both are land there (closed).

---

## 8. Making a new region — checklist

1. `new_config.sh NAME lon_min lon_max lat_min lat_max [--res N] [--start/--end] [--tides] [--rivers]`
2. `build_config.sh NAME --only 1` and `--only 2`. Read step 02's output: are the boundaries physically
   sensible, is the Courant number below 0.7, are there odd isolated water bodies? Adjust `domain.cfg` if needed.
3. `build_config.sh NAME --from 3` (downloads, then build and prove).
4. `run_hindcast_cycle.sh NAME` → `verify_config.py NAME [build]`.
5. A finer child: `new_config.sh CHILD … --res N --parent NAME`, then `build_config.sh CHILD --no-test` and
   `run_hindcast_cycle.sh CHILD` (section 7).

## 9. Troubleshooting

| message | meaning / fix |
|---|---|
| `domain.cfg not found` | create it with `new_config.sh` |
| `open boundaries unknown` | run step 02 (or set `OBC` in domain.cfg) |
| `missing …/etopo2.nc` / `GSHHS` | install DATASETS_CROCOTOOLS (install/README.md §3) |
| Courant > 0.7 warning | set a smaller `DT` in domain.cfg, re-run from step 05 |
| `IndexError: index 1 is out of bounds` | monthly-mean GLORYS: re-download the daily `P1D-m` |
| `ERROR in get_bry … 'bry_time'` | a boundary file doesn't cover the run ± 1 day: missing neighbour month |
| `Abnormal termination: netCDF INPUT` after *Open Meteo file* | ERA5 file missing, or wrong path or month in `online` |
| `nf-config is …, expected …/opt_seq/bin/nf-config` | re-source env.sh; rebuild the NetCDF stack (install/) |
| `binary not found: …/croco_<build>` | `steps/07_compile.sh NAME <build>` |
| CDS `Number queued requests … temporarily limited` | normal under load; retried automatically (keep `ERA5_N_PARALLEL=2`) |
| `.ewdsapirc not found` / GloFAS 401 | set up EWDS (section 6) |
| NaN / `BLOW UP` | boundaries vs mask (step 02), DT, or try `SPONGE="50000. 400."` |
| `no <PARENT> run … covers …` / `parent run not finished` | run (or finish) the parent's hindcast over the child's period first |
| `the child box must lie inside the parent` | shrink the child box (or enlarge the parent) |
| child: `ERROR in get_bry … 'bry_time'` at the end | the child ends too late: at most the parent's end − 1 day |
| ERA5 `TimeoutError` / connection reset | temporary network trouble: retried 5 times; re-run to resume (finished variables are cached in `ERA5/raw/`) |
