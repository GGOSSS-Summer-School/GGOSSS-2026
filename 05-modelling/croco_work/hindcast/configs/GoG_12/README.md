# GoG_12 — Gulf of Guinea, 1/12°

The only hand-written file here is **`domain.cfg`**. It was created with

```bash
hindcast/new_config.sh GoG_12 -10 12.5 -6 8
```

i.e. 10°W–12.5°E, 6°S–8°N, 1/12°, hindcast 2018-06-01 → 2018-08-30, with tides and rivers **off**.
To switch them on, set `TIDES=1` and/or `RIVERS=1` in `domain.cfg` (rivers need `~/.ewdsapirc`, see
`hindcast/README.md` §6). Then:

```bash
hindcast/build_config.sh GoG_12
nohup hindcast/run_hindcast_cycle.sh GoG_12 > ~/GoG_12.log 2>&1 &
```

The build writes the rest of this folder: `grid.ini`, `obc.cfg`, `crocotools_param.py`,
`crocotools_param_tides.py`, `cppdefs.h`, `param.h`, `croco.in`, `jobcomp`, and the region card `CARD.md`.

What the build produced when this was verified (WSL, CROCO v2.1.3):

| | |
|---|---|
| Grid | 273 × 171 (LLm0 = 271, MMm0 = 169), depth 50–5641 m, rx0 0.200 |
| Boundaries (from the mask) | south 272/273 and west 146/171 ocean → **open**; east 1/171 and north 1/273 → **closed** |
| Time step | dt 300 s, NDTFAST 60, barotropic Courant 0.13 |
| Period | 18 cycles × (2 d spin-up + 5 d); GLORYS 2018-05…09, ERA5 2018-05…08 |
| Tides | M2 amplitude 0.45 m mean / 0.76 m max, 0 fill values |
| Rivers (inside the grid) | Niger, Ogooué, Sanaga, Volta, Cross, … (depending on GloFAS discharge ≥ 100 m³/s); the Congo mouth is just outside the southern edge |
