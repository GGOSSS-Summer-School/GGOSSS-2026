# Installing the `ggosss26` conda environment

This guide installs **the Python/conda environment** defined in
`ggosss2026_conda_env.yaml`. That environment is the one used by

- **croco_work**: the CROCO pre-processing tools (`gtools/`, `croco_pytools`), downloads (CMEMS, ERA5), plotting and validation; and
- the **OpenDrift exercises** (`open_drift/*.ipynb`).

> **Done once per machine.** After the install, every session only needs `conda activate ggosss26`.

---

## 0. Pick your platform

| Platform | What you use |
|---|---|
| **Linux** (Ubuntu 20.04+, etc.) | A normal terminal. |
| **macOS** 12+ (Intel or Apple Silicon) | The Terminal app. Every package in the YAML has macOS builds (`osx-64` / `osx-arm64`) on conda-forge. |
| **Windows 10/11** | **WSL2 + Ubuntu**. CROCO and its tools are built for Linux, so don't use a native Windows conda. |

### Windows only: install WSL2

1. Open **PowerShell as Administrator** and run:
   ```powershell
   wsl --install -d Ubuntu
   ```
   You can also install **Ubuntu** from the Microsoft Store.
2. Restart if asked. Open **Ubuntu** from the Start menu and create a UNIX username and password. The password isn't shown while you type, which is expected.
3. Check that you're in Linux:
   ```bash
   uname -a     # should mention "microsoft-standard-WSL2"
   ```

From here on, **every command runs in the Ubuntu terminal**.

Put your files in the Linux home directory (`~/`), not under `/mnt/c/...`. Working from the Windows drive is
many times slower. To copy the YAML over from Windows downloads (replace `<WinUser>`
with your Windows user name):

```bash
cp /mnt/c/Users/<WinUser>/Downloads/ggosss2026_conda_env.yaml ~/
```

To see the Linux files from Windows Explorer, go to `\\wsl$\Ubuntu\home\<linux-user>` in your Windows Explorer zdress bar.
If your  `linux-user` is for example `student`, then `\\wsl$\Ubuntu\home\student`

---

## 1. Install conda (Miniconda)

Install Miniconda into the default location, **`~/miniconda3`**. The croco_work scripts

Download the installer for your machine:

| Machine | Installer |
|---|---|
| Linux / WSL2, x86_64 | `Miniconda3-latest-Linux-x86_64.sh` |
| Linux, ARM (aarch64) | `Miniconda3-latest-Linux-aarch64.sh` |
| macOS, Apple Silicon (M1–M4) | `Miniconda3-latest-MacOSX-arm64.sh` |
| macOS, Intel | `Miniconda3-latest-MacOSX-x86_64.sh` |

```bash
cd ~
# Linux / WSL2 (x86_64). On macOS, replace the file name using the table above;
# curl is used because macOS has no wget by default.
curl -LO https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
```

- Accept the licence.
- Keep the default location (`~/miniconda3`).
- Answer **yes** when asked whether to initialise conda.

Close and reopen the terminal, or run `source ~/.bashrc` (on macOS, `source ~/.zshrc`). The prompt
should now start with `(base)`.

```bash
conda --version
```

### Use conda-forge only (recommended)

The YAML only uses `conda-forge`. Taking Anaconda's `defaults` channel out of the configuration
prevents mixed-channel conflicts. It also avoids the “Terms of Service” prompt that recent Miniconda
versions show for `defaults`:

```bash
conda config --add channels conda-forge
conda config --remove channels defaults      # ignore "not found" if it isn't there
conda config --set channel_priority strict
```

If conda still asks you to accept the Terms of Service for `repo.anaconda.com` channels, either accept
them (`conda tos accept`) or run the `--remove channels defaults` line again.

---

## 2. Create the environment

Go to the folder that holds the YAML (`~` if you copied it there):

```bash
cd ~
conda env create -f ggosss2026_conda_env.yaml
```

- The environment name comes from the first line of the file, `name: ggosss26`. The croco_work scripts activate exactly that name, so don't rename it.
- Creating the environment takes about **10–30 minutes**: it solves roughly 50 conda packages, then pip installs `regionmask`, `pydap`, `copernicusmarine` and `opendrift` (OpenDrift pulls in `cmocean`, `roaring-landmask`, `adios_db` and more).
- Don't interrupt it. If it fails part-way.

Activate it:

```bash
conda activate ggosss26
```

The prompt now shows `(ggosss26)`.


---

## 3. Verify the environment

### 3a. Imports (croco_work + OpenDrift)

```bash
conda activate ggosss26
python - <<'EOF'
import xarray, netCDF4, numpy, scipy, pandas, dask, h5netcdf, zarr, cftime
import cartopy, geopandas, shapely, pyproj, pyinterp, regionmask
import copernicusmarine, cdsapi, cfgrib, pydap
import xgcm, xrft, cf_xarray, numba, pyamg
import opendrift, cmocean
from opendrift.models.oceandrift import OceanDrift
from opendrift.models.leeway import Leeway
from opendrift.models.openoil import OpenOil
from opendrift.readers import reader_netCDF_CF_generic, reader_ROMS_native
print("numpy", numpy.__version__, "| xarray", xarray.__version__, "| opendrift", opendrift.__version__)
print("ggosss26 env OK")
EOF
```

You should see **`ggosss26 env OK`**.

### 3a. OpenDrift smoke test (no data needed)

```bash
python - <<'EOF'
from datetime import datetime, timedelta
from opendrift.models.oceandrift import OceanDrift
o = OceanDrift(loglevel=50)
o.set_config('environment:fallback:x_sea_water_velocity', 0.2)
o.set_config('environment:fallback:y_sea_water_velocity', 0.1)
o.seed_elements(lon=3.0, lat=4.0, number=10, time=datetime(2026, 9, 18))
o.run(duration=timedelta(hours=6), time_step=3600)
print("OpenDrift OK:", o.elements.lon.mean().round(3), o.elements.lat.mean().round(3))
EOF
```

The first run downloads nothing. The `roaring-landmask` coastline ships with the package.

---


## 5. Troubleshooting

| Symptom | Fix |
|---|---|
| `conda: command not found` | Reopen the terminal, or run `source ~/miniconda3/etc/profile.d/conda.sh`. |
| YAML parse error on `conda env create` | The first line must be exactly `name: ggosss26`, with a space after the colon. Don't edit the file in Word or Notepad with “smart” formatting. |
| Terms of Service error for `repo.anaconda.com` | See “Use conda-forge only” in step 1, or run `conda tos accept`. |
| Conda part OK, **pip part failed** (network or proxy) | Finish it by hand: `conda activate ggosss26 && pip install regionmask pydap copernicusmarine opendrift` |
| Solve takes forever or conflicts | Update conda (`conda update -n base conda`) and make sure `channel_priority` is `strict`. Then start over: `conda env remove -n ggosss26` and create it again. |
| Environment exists but is broken or old | `conda env remove -n ggosss26`, then `conda env create -f ggosss2026_conda_env.yaml` |
| WSL2: install very slow, or files “missing” | Work in `~` (Linux filesystem), not in `/mnt/c/...`. |
| WSL2: runs out of memory | Create `C:\Users\<WinUser>\.wslconfig` with `[wsl2]` and `memory=8GB` on the next line. Then run `wsl --shutdown` in PowerShell. |
| `import wx` fails on a server with no display | Only `croco_pyvisu`'s GUI needs wxPython. The rest of the environment is unaffected. |

---

## Summary
# example for Linux / Windows (WSL2)

```bash
# once per machine
cd ~
curl -LO https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-linux-x86_64.sh              # -> ~/miniconda3
conda config --add channels conda-forge && conda config --set channel_priority strict

cp /mnt/c/Users/<WinUser>/Downloads/ggosss2026_conda_env.yaml ~/
conda env create -f ggosss2026_conda_env.yaml      # -> env "ggosss26"
conda activate ggosss26
python -m ipykernel install --user --name ggosss26 --display-name "Python (ggosss26)"

# every session
conda activate ggosss26
```
# Installing the `ggosss26` conda environment

This guide installs **the Python/conda environment** defined in
`ggosss2026_conda_env.yaml`. That environment is the one used by

- **croco_work**: the CROCO pre-processing tools (`gtools/`, `croco_pytools`), downloads (CMEMS, ERA5), plotting and validation; and
- the **OpenDrift exercises** (`open_drift/*.ipynb`).

> **Done once per machine.** After the install, every session only needs `conda activate ggosss26`.

---

## 0. Pick your platform

| Platform | What you use |
|---|---|
| **Linux** (Ubuntu 20.04+, etc.) | A normal terminal. |
| **macOS** 12+ (Intel or Apple Silicon) | The Terminal app. Every package in the YAML has macOS builds (`osx-64` / `osx-arm64`) on conda-forge. |
| **Windows 10/11** | **WSL2 + Ubuntu**. CROCO and its tools are built for Linux, so don't use a native Windows conda. |

### Windows only: install WSL2

1. Open **PowerShell as Administrator** and run:
   ```powershell
   wsl --install -d Ubuntu
   ```
   You can also install **Ubuntu** from the Microsoft Store.
2. Restart if asked. Open **Ubuntu** from the Start menu and create a UNIX username and password. The password isn't shown while you type, which is expected.
3. Check that you're in Linux:
   ```bash
   uname -a     # should mention "microsoft-standard-WSL2"
   ```

From here on, **every command runs in the Ubuntu terminal**.

Put your files in the Linux home directory (`~/`), not under `/mnt/c/...`. Working from the Windows drive is
many times slower. To copy the YAML over from Windows downloads (replace `<WinUser>`
with your Windows user name):

```bash
cp /mnt/c/Users/<WinUser>/Downloads/ggosss2026_conda_env.yaml ~/
```

To see the Linux files from Windows Explorer, go to `\\wsl$\Ubuntu\home\<linux-user>` in your Windows Explorer zdress bar.
If your  `linux-user` is for example `student`, then `\\wsl$\Ubuntu\home\student`

---

## 1. Install conda (Miniconda)

Install Miniconda into the default location, **`~/miniconda3`**. The croco_work scripts

Download the installer for your machine:

| Machine | Installer |
|---|---|
| Linux / WSL2, x86_64 | `Miniconda3-latest-Linux-x86_64.sh` |
| Linux, ARM (aarch64) | `Miniconda3-latest-Linux-aarch64.sh` |
| macOS, Apple Silicon (M1–M4) | `Miniconda3-latest-MacOSX-arm64.sh` |
| macOS, Intel | `Miniconda3-latest-MacOSX-x86_64.sh` |

```bash
cd ~
# Linux / WSL2 (x86_64). On macOS, replace the file name using the table above;
# curl is used because macOS has no wget by default.
curl -LO https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
```

- Accept the licence.
- Keep the default location (`~/miniconda3`).
- Answer **yes** when asked whether to initialise conda.

Close and reopen the terminal, or run `source ~/.bashrc` (on macOS, `source ~/.zshrc`). The prompt
should now start with `(base)`.

```bash
conda --version
```

### Use conda-forge only (recommended)

The YAML only uses `conda-forge`. Taking Anaconda's `defaults` channel out of the configuration
prevents mixed-channel conflicts. It also avoids the “Terms of Service” prompt that recent Miniconda
versions show for `defaults`:

```bash
conda config --add channels conda-forge
conda config --remove channels defaults      # ignore "not found" if it isn't there
conda config --set channel_priority strict
```

If conda still asks you to accept the Terms of Service for `repo.anaconda.com` channels, either accept
them (`conda tos accept`) or run the `--remove channels defaults` line again.

---

## 2. Create the environment

Go to the folder that holds the YAML (`~` if you copied it there):

```bash
cd ~
conda env create -f ggosss2026_conda_env.yaml
```

- The environment name comes from the first line of the file, `name: ggosss26`. The croco_work scripts activate exactly that name, so don't rename it.
- Creating the environment takes about **10–30 minutes**: it solves roughly 50 conda packages, then pip installs `regionmask`, `pydap` and `opendrift` (OpenDrift pulls in `roaring-landmask`, `adios_db` and more).
- Don't interrupt it. If it fails part-way.

Activate it:

```bash
conda activate ggosss26
```

The prompt now shows `(ggosss26)`.


---

## 3. Verify the environment

### 3a. Imports (croco_work + OpenDrift)

```bash
conda activate ggosss26
python - <<'EOF'
import xarray, netCDF4, numpy, scipy, pandas, dask, h5netcdf, zarr, cftime
import cartopy, geopandas, shapely, pyproj, pyinterp, regionmask
import copernicusmarine, cdsapi, cfgrib, pydap
import xgcm, xrft, cf_xarray, numba, pyamg
import opendrift, cmocean, gsw, sklearn
from opendrift.models.oceandrift import OceanDrift
from opendrift.models.leeway import Leeway
from opendrift.models.openoil import OpenOil
from opendrift.readers import reader_netCDF_CF_generic, reader_ROMS_native
print("numpy", numpy.__version__, "| xarray", xarray.__version__, "| opendrift", opendrift.__version__)
print("ggosss26 env OK")
EOF
```

You should see **`ggosss26 env OK`**.

### 3a. OpenDrift smoke test (no data needed)

```bash
python - <<'EOF'
from datetime import datetime, timedelta
from opendrift.models.oceandrift import OceanDrift
o = OceanDrift(loglevel=50)
o.set_config('environment:fallback:x_sea_water_velocity', 0.2)
o.set_config('environment:fallback:y_sea_water_velocity', 0.1)
o.seed_elements(lon=3.0, lat=4.0, number=10, time=datetime(2026, 9, 18))
o.run(duration=timedelta(hours=6), time_step=3600)
print("OpenDrift OK:", o.elements.lon.mean().round(3), o.elements.lat.mean().round(3))
EOF
```

The first run downloads nothing. The `roaring-landmask` coastline ships with the package.

---


## 4. Troubleshooting

| Symptom | Fix |
|---|---|
| `conda: command not found` | Reopen the terminal, or run `source ~/miniconda3/etc/profile.d/conda.sh`. |
| YAML parse error on `conda env create` | The first line must be exactly `name: ggosss26`, with a space after the colon. Don't edit the file in Word or Notepad with “smart” formatting. |
| Terms of Service error for `repo.anaconda.com` | See “Use conda-forge only” in step 1, or run `conda tos accept`. |
| Conda part OK, **pip part failed** (network or proxy) | Finish it by hand: `conda activate ggosss26 && pip install regionmask pydap opendrift` |
| Solve takes forever or conflicts | Update conda (`conda update -n base conda`) and make sure `channel_priority` is `strict`. Then start over: `conda env remove -n ggosss26` and create it again. |
| Environment exists but is broken or old | `conda env remove -n ggosss26`, then `conda env create -f ggosss2026_conda_env.yaml` |
| WSL2: install very slow, or files “missing” | Work in `~` (Linux filesystem), not in `/mnt/c/...`. |
| WSL2: runs out of memory | Create `C:\Users\<WinUser>\.wslconfig` with `[wsl2]` and `memory=8GB` on the next line. Then run `wsl --shutdown` in PowerShell. |
| `import wx` fails on a server with no display | Only `croco_pyvisu`'s GUI needs wxPython. The rest of the environment is unaffected. |

---

## Summary
# example for Linux / Windows (WSL2)

```bash
# once per machine
cd ~
curl -LO https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-linux-x86_64.sh              # -> ~/miniconda3
conda config --add channels conda-forge && conda config --set channel_priority strict

cp /mnt/c/Users/<WinUser>/Downloads/ggosss2026_conda_env.yaml ~/
conda env create -f ggosss2026_conda_env.yaml      # -> env "ggosss26"
conda activate ggosss26
python -m ipykernel install --user --name ggosss26 --display-name "Python (ggosss26)"

# every session
conda activate ggosss26
```
