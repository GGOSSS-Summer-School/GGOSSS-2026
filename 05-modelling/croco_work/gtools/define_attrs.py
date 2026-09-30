"""
gtools/define_attrs.py — CF metadata + display defaults for CROCO fields.

Each variable declares, in one place:
  - CF metadata:  long_name, units, standard_name
  - display:      cmap, vmin, vmax (None = auto from data), diverging

so plots label AND colour themselves from the variable being drawn. Velocity,
anomaly and difference fields are 'diverging' -> a symmetric range about zero is
chosen from the data when vmin/vmax are None.
"""
from __future__ import annotations
from dataclasses import dataclass
from typing import Optional
import numpy as np


@dataclass
class VarAttrs:
    long_name: str
    units: str
    standard_name: str
    cmap: str = "viridis"
    vmin: Optional[float] = None      # None -> auto from data
    vmax: Optional[float] = None
    diverging: bool = False           # symmetric range about 0 when auto


# ----------------------------------------------------------------------
# Scalar and grid variables
# ----------------------------------------------------------------------
ATTRS = {
    # --- grid / coordinates ---
    'xi_rho':  VarAttrs('x-dimension of the grid', '1', 'x_grid_index'),
    'eta_rho': VarAttrs('y-dimension of the grid', '1', 'y_grid_index'),
    'lon_rho': VarAttrs('Longitude', 'degrees_east', 'longitude'),
    'lat_rho': VarAttrs('Latitude', 'degrees_north', 'latitude'),
    'h':       VarAttrs('Depth of the sea floor', 'm', 'sea_floor_depth',
                        cmap='Blues', vmin=0.0, vmax=None),
    'mask':    VarAttrs('Land-sea mask (1=water, 0=land)', '1', 'land_binary_mask',
                        cmap='gray', vmin=0.0, vmax=1.0),
    'depth':   VarAttrs('Depth', 'm', 'depth'),
    'f':       VarAttrs('Coriolis parameter', '$s^{-1}$', 'coriolis_parameter',
                        cmap='RdBu_r', diverging=True),

    # --- ocean state ---
    'zeta':    VarAttrs('Sea Surface Elevation', 'm', 'sea_surface_elevation',
                        cmap='RdBu_r', diverging=True),
    'ssh':     VarAttrs('Sea Surface Height', 'm', 'sea_surface_height',
                        cmap='RdBu_r', diverging=True),
    'temp':    VarAttrs('Sea Water Temperature', 'degC', 'sea_water_temperature',
                        cmap='RdYlBu_r', vmin=None, vmax=None),
    'sst':     VarAttrs('Sea Surface Temperature', 'degC', 'sea_surface_temperature',
                        cmap='RdYlBu_r'),
    'salt':    VarAttrs('Sea Water Salinity', 'PSU', 'sea_water_salinity',
                        cmap='viridis', vmin=None, vmax=None),
    'rho':     VarAttrs('Sea Water Density', '$kg.m^{-3}$', 'sea_water_density',
                        cmap='cmo.dense', vmin=None, vmax=None),
    'sigma_t': VarAttrs('Sea Water Sigma-t (density - 1000)', 'kg.m-3',
                        'sea_water_sigma_t', cmap='cmo.dense'),
    'w':       VarAttrs('Upward seawater velocity', '$m.s^{-1}$', 'upward_sea_water_velocity',
                        cmap='RdBu_r', diverging=True),
    'speed':   VarAttrs('Current speed', '$m.s^{-1}$', 'sea_water_speed',
                        cmap='viridis', vmin=0.0, vmax=None),

    # --- derived ---
    'vort':    VarAttrs('Relative vorticity', '$\zeta$ $(s^{-1})$', 'ocean_relative_vorticity',
                        cmap='RdBu_r', diverging=True),
    'vort_f':  VarAttrs('Relative vorticity ($\zeta$)/ f', '', 'normalized_relative_vorticity',
                        cmap='RdBu_r', vmin=-1.0, vmax=1.0, diverging=True),
    'eke':     VarAttrs('Eddy kinetic energy', '$m^{2}.s^{-2}$', 'eddy_kinetic_energy',
                        cmap='hot_r', vmin=0.0, vmax=None),
    'mld':     VarAttrs('Mixed layer depth', 'm', 'ocean_mixed_layer_thickness',
                        cmap='viridis_r', vmin=0.0, vmax=None),
    'n2':      VarAttrs('Brunt-Vaisala frequency squared', '$s^{-2}$',
                        'square_of_brunt_vaisala_frequency_in_sea_water',
                        cmap='viridis', vmin=0.0, vmax=None),

    # --- atmosphere / fluxes ---
    'wind_speed': VarAttrs('10 m wind speed', '$m.s^{-1}$', 'wind_speed',
                           cmap='viridis', vmin=0.0, vmax=None),
    'shflux':  VarAttrs('Surface net heat flux', '$W.m^{-2}$', 'surface_downward_heat_flux_in_sea_water',
                        cmap='RdBu_r', diverging=True),
    'swflux':  VarAttrs('Surface freshwater flux (E-P)', 'cm.day^{-1}',
                        'water_flux_out_of_sea_water', cmap='BrBG', diverging=True),
    'swrad':   VarAttrs('Shortwave radiation', '$W.m^{-2}$', 'surface_net_downward_shortwave_flux',
                        cmap='hot_r', vmin=0.0, vmax=None),

    # --- differences / anomalies (generic) ---
    'diff':    VarAttrs('Difference (model - reference)', '', 'difference',
                        cmap='RdBu_r', diverging=True),
    'anom':    VarAttrs('Anomaly', '', 'anomaly', cmap='RdBu_r', diverging=True),
}


# ----------------------------------------------------------------------
# Vector variables: (grid-aligned, east/north, section) triples
# ----------------------------------------------------------------------
def _vel(lnx, sn):
    return VarAttrs(lnx, '$m.s^{-1}$', sn, cmap='RdBu_r', diverging=True)


VECTOR_ATTRS = {
    'u':   (_vel('Sea water velocity in x direction', 'sea_water_x_velocity'),
            _vel('Eastward sea water velocity', 'eastward_sea_water_velocity'),
            _vel('Across-section sea water velocity', 'across_section_sea_water_velocity')),
    'v':   (_vel('Sea water velocity in y direction', 'sea_water_y_velocity'),
            _vel('Northward sea water velocity', 'northward_sea_water_velocity'),
            _vel('Along-section sea water velocity', 'along_section_sea_water_velocity')),
    'ubar':(_vel('Barotropic velocity in x direction', 'barotropic_sea_water_x_velocity'),
            _vel('Eastward barotropic velocity', 'barotropic_eastward_sea_water_velocity'),
            _vel('Across-section barotropic velocity', 'barotropic_across_section_sea_water_velocity')),
    'vbar':(_vel('Barotropic velocity in y direction', 'barotropic_sea_water_y_velocity'),
            _vel('Northward barotropic velocity', 'barotropic_northward_sea_water_velocity'),
            _vel('Along-section barotropic velocity', 'barotropic_along_section_sea_water_velocity')),
    'u10': (VarAttrs('10 m wind in x direction', '$m.s^{-1}$', 'x_wind', cmap='RdBu_r', diverging=True),
            VarAttrs('Eastward 10 m wind', '$m.s^{-1}$', 'eastward_wind', cmap='RdBu_r', diverging=True),
            VarAttrs('Across-section 10 m wind', '$m.s^{-1}$', 'across_section_wind', cmap='RdBu_r', diverging=True)),
    'v10': (VarAttrs('10 m wind in y direction', '$m.s^{-1}$', 'y_wind', cmap='RdBu_r', diverging=True),
            VarAttrs('Northward 10 m wind', '$m.s^{-1}$', 'northward_wind', cmap='RdBu_r', diverging=True),
            VarAttrs('Along-section 10 m wind', '$m.s^{-1}$', 'along_section_wind', cmap='RdBu_r', diverging=True)),
    'sustr':(VarAttrs('Wind stress in x direction', '$N.m^{-2}$', 'surface_downward_x_stress', cmap='RdBu_r', diverging=True),
             VarAttrs('Eastward surface stress', '$N.m^{-2}$', 'surface_eastward_stress', cmap='RdBu_r', diverging=True),
             VarAttrs('Across-section surface stress', '$N.m^{-2}$', 'surface_across_section_stress', cmap='RdBu_r', diverging=True)),
    'svstr':(VarAttrs('Wind stress in y direction', '$N.m^{-2}$', 'surface_downward_y_stress', cmap='RdBu_r', diverging=True),
             VarAttrs('Northward surface stress', '$N.m^{-2}$', 'surface_northward_stress', cmap='RdBu_r', diverging=True),
             VarAttrs('Along-section surface stress', '$N.m^{-2}$', 'surface_along_section_stress', cmap='RdBu_r', diverging=True)),
}


# ----------------------------------------------------------------------
# apply + resolve
# ----------------------------------------------------------------------
def _meta_for(var_str, rotated=False, section=False):
    if var_str in ATTRS:
        return ATTRS[var_str]
    if var_str in VECTOR_ATTRS:
        if section:
            return VECTOR_ATTRS[var_str][2]
        if rotated:
            return VECTOR_ATTRS[var_str][1]
        return VECTOR_ATTRS[var_str][0]
    return None


def apply_attrs(da, var_str, rotated=False, section=False):
    """Attach CF + display attributes to a DataArray (in place)."""
    meta = _meta_for(var_str, rotated=rotated, section=section)
    if meta is None:
        return da
    da.attrs['long_name'] = meta.long_name
    da.attrs['units'] = meta.units
    da.attrs['standard_name'] = meta.standard_name
    da.attrs['cmap'] = meta.cmap
    if meta.vmin is not None:
        da.attrs['vmin'] = meta.vmin
    if meta.vmax is not None:
        da.attrs['vmax'] = meta.vmax
    da.attrs['diverging'] = int(meta.diverging)
    return da


def display_for(var):
    """Return {cmap, vmin, vmax, diverging} for a variable name or DataArray."""
    if hasattr(var, "attrs") and not isinstance(var, str):
        a = var.attrs
        return {"cmap": a.get("cmap", "viridis"),
                "vmin": a.get("vmin", None),
                "vmax": a.get("vmax", None),
                "diverging": bool(a.get("diverging", 0))}
    meta = _meta_for(var)
    if meta is None:
        return {"cmap": "viridis", "vmin": None, "vmax": None, "diverging": False}
    return {"cmap": meta.cmap, "vmin": meta.vmin, "vmax": meta.vmax,
            "diverging": meta.diverging}


def resolve_limits(da, cmap=None, vmin=None, vmax=None):
    """Resolve (cmap, vmin, vmax) for a plot.

    Priority: explicit kwarg > variable attribute > auto from data.
    Diverging variables with no fixed range get a symmetric range about 0.
    """
    disp = display_for(da)
    cmap = cmap if cmap is not None else disp["cmap"]
    vmin = vmin if vmin is not None else disp["vmin"]
    vmax = vmax if vmax is not None else disp["vmax"]

    data = np.asarray(da.values) if hasattr(da, "values") else np.asarray(da)
    finite = data[np.isfinite(data)]
    if finite.size == 0:
        return cmap, vmin, vmax

    if disp["diverging"] and (vmin is None or vmax is None):
        m = float(np.nanmax(np.abs(finite)))
        vmin = -m if vmin is None else vmin
        vmax = m if vmax is None else vmax
    else:
        if vmin is None:
            vmin = float(np.nanmin(finite))
        if vmax is None:
            vmax = float(np.nanmax(finite))
    return cmap, vmin, vmax
