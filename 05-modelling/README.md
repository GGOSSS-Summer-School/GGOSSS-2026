# croco_work — CROCO regional ocean hindcasts and Lagrangian drift

A self-contained CROCO ocean-model system that reconstructs the past ocean of **any region** you choose:
- ocean initial and boundary conditions from the **GLORYS** reanalysis (CMEMS);
- atmospheric forcing from **ERA5** (CDS);
- optional **tides** (TPXO) and optional **rivers** (daily **GloFAS** discharge);
- **offline nesting**: a finer child configuration driven by a parent's own hindcast output;
- **Lagrangian trajectory modelling with OpenDrift**: particles (oil, people and objects at sea, plastics) advected by
  the CROCO currents or by GLORYS, with ERA5 winds.

The period is run as one or several cycles, each a short spin-up followed by the hindcast.
Each region is one small settings file; everything else is derived and written by the scripts.

| where | what |
|---|---|
| `install/` + `install/README.md` | one-time setup: NetCDF/HDF5 stack (`opt_seq/`), CROCO v2.1.3 + croco_pytools v2.0.4 (`code/`), conda env `ggosss2026`, DATASETS_CROCOTOOLS (`data/`), accounts |
| `env.sh` | shared paths + compilers (sourced by every script) |
| `gtools/` | the Python tools; `gtools/ggosss26.py` is the command-line interface |
| **`hindcast/`** + **`hindcast/README.md`** | **the hindcast system**: `new_config.sh`, `build_config.sh`, `run_hindcast_cycle.sh`, `steps/`, `configs/<NAME>/domain.cfg` |
| `notebooks/GoG12_hindcast_user_guide.ipynb` | step-by-step guide from a bare machine to a finished hindcast, with expected outputs |
| `notebooks/verify_config.py` | checks the cycles of a finished hindcast |
| OpenDrift (in the `ggosss2026` env) | Lagrangian drift driven by the hindcast outputs, see [Lagrangian modelling with OpenDrift](#lagrangian-modelling-with-opendrift) |
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

---

## Lagrangian modelling with OpenDrift

The hindcast gives **Eulerian** fields: currents, temperature and salinity known at grid points every output step.
**OpenDrift** (Dagestad et al. 2018, *GMD* 11, 1405–1420; <https://opendrift.github.io>) adds the **Lagrangian** view.
It follows particles ("elements") through those fields, which answers questions such as *where does an oil slick
go*, *where should we search for a person in the water*, and *where do river plastics end up on the coast*.

OpenDrift 1.14.11 is installed in the `ggosss2026` environment (`conda activate ggosss2026`).

### How a particle moves
Each element is advected by the interpolated forcing, plus a random walk for motions the grid does not resolve:

  dx/dt = u_current(x, t) + α·U₁₀ + u_Stokes + u_leeway      x_{n+1} = x_n + v·Δt + √(2KΔt)·ξ,  ξ ~ N(0, 1)

- u_current comes from the ocean model; α·U₁₀ is the windage (≈ 2–3 % of the 10 m wind for floating material).
- u_Stokes is the wave-induced drift (from a wave product, or estimated from the wind).
- u_leeway is the object-specific drift of Search & Rescue targets.
- Integration is Euler or Runge–Kutta 4. In 3-D, vertical turbulent mixing (Kz) and buoyancy (droplet rise)
  act as well.

### Forcing: what this system provides

| forcing | source in `croco_work` | OpenDrift reader |
|---|---|---|
| ocean currents u, v (+ T, S, bathymetry) | the CROCO hindcast: `hindcast/model-runs/<NAME>/<date>_<build>/hcast/CROCO_FILES/croco_his.nc` (hourly with tides) — e.g. GoG_12 at 1/12°, IGOG_36 at 1/36° | `reader_ROMS_native` (reads CROCO/ROMS sigma-grid output) |
| ocean currents from the reanalysis | GLORYS12: `hindcast/scratch/<NAME>/downloaded_data/GLORYS/YYYY_MM.nc` | `reader_netCDF_CF_generic` |
| 10 m wind | ERA5: `…/downloaded_data/ERA5/for_croco/U10M_Y*.nc`, `V10M_Y*.nc` | `reader_netCDF_CF_generic` (see the note below) |
| coastline | GSHHS (OpenDrift's built-in `reader_global_landmask`) | automatic |

Running the same experiment with GLORYS and with the CROCO parent and child isolates the effect of resolution
(1/12° vs 1/36°) and of the time sampling (daily means vs hourly, including tides) on the trajectories.

### The models (same engine, different element types)

| application | OpenDrift model | specific physics | typical product |
|---|---|---|---|
| oil spill | `OpenOil` | evaporation, emulsification, natural dispersion, 2-D or 3-D | oil budget, beached mass |
| search & rescue | `Leeway` | object-specific leeway coefficients (person in water, life raft, fishing vessel) | search area at T + 24 h |
| marine plastics / passive tracers | `OceanDrift` (or `PlastDrift`) | windage, stranding and resuspension at the coast; river sources such as the Sanaga, Nyong and Wouri | coastal accumulation sectors |

### Workflow
1. **Readers**: attach the forcing (CROCO or GLORYS currents, ERA5 or GFS wind).
2. **Seed** elements: position, time, number, radius; or along a river mouth, continuously.
3. **Configure**: 2-D or 3-D, weathering, `general:coastline_action` (`stranding` or `previous`), diffusivity.
4. **Run**: time step, duration, output time step.
5. **Analyse**: `o.result` is an xarray Dataset of the trajectories.
6. **Compare** forcings: resolution vs temporal sampling, parent vs child.

```python
import os
from datetime import timedelta
from opendrift.models.oceandrift import OceanDrift
from opendrift.readers import reader_ROMS_native

o = OceanDrift(loglevel=20)
croco = reader_ROMS_native.Reader(os.path.expanduser(       # the reader does not expand "~" itself
    "~/croco_work/hindcast/model-runs/GoG_12/20260601_plain_tides_rivers/hcast/CROCO_FILES/croco_his.nc"))
o.add_reader(croco)
o.set_config("general:coastline_action", "stranding")
o.set_config("environment:constant:horizontal_diffusivity", 10)   # m2/s: the random walk for unresolved motions
o.seed_elements(lon=9.4, lat=3.8, radius=2000, number=1000, time=croco.start_time,   # e.g. off the Sanaga mouth
                wind_drift_factor=0.02)                          # windage α = 2 % (used once a wind reader is added)
o.run(duration=timedelta(days=10), time_step=timedelta(minutes=30),
      time_step_output=timedelta(hours=1), outfile="drift_GoG12.nc")
o.plot(fast=True)
print(o.result)                                              # xarray: lon, lat, status per element and time
```

**Diagnostics** commonly derived from the output:
- time to the first coastal impact, and the stranded fraction;
- area of the particle cloud (convex hull) and separation of the centroids between two forcings;
- probability (density) maps from ensembles of seedings.

### Notes
- **Time origin of CROCO files.** Our CROCO files count time in seconds from **1993-01-01**. If OpenDrift cannot
  date them (the `time` units must read `seconds since 1993-01-01 00:00:00`), set the attribute once with
  `ncatted -a units,time,o,c,"seconds since 1993-01-01 00:00:00" croco_his.nc`.
- **Wind from ERA5.** The `for_croco` files use CROCO's names (`U10M`, `V10M`, time in days since 1993). OpenDrift
  identifies the wind through the CF `standard_name`. Either add `standard_name = x_wind` / `y_wind` (and CF time
  units), or read ERA5 through OpenDrift's own readers from the original files.
- **Forcing time step.** OpenDrift interpolates in time between records. Hourly CROCO history resolves the
  tidal currents, while daily GLORYS means do not: that difference is part of what a forcing comparison shows.