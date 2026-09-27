from pathlib import Path

import pandas as pd


def find_data_dir() -> Path:
    for candidate in [Path("../data"), Path("../Datasets"), Path("data"), Path("Datasets")]:
        if candidate.exists():
            return candidate
    return Path(__file__).resolve().parents[1] / "Datasets"


data_dir = find_data_dir()
rates = pd.read_csv(data_dir / "pak_shoreline_change_rates_2026.csv")
boxes = pd.read_csv(data_dir / "pak_coastsat_box_inventory.csv")
wind = pd.read_csv(data_dir / "pak_wind_zonal_stats.csv")
currents = pd.read_csv(data_dir / "pak_current_points_stats.csv")
factors = pd.read_csv(data_dir / "pak_integrated_coastal_factors.csv")
localities = pd.read_csv(data_dir / "selected_localities_for_maps.csv")

classes = pd.cut(
    rates["endpoint_rate_myr"],
    bins=[float("-inf"), -0.5, 0.5, float("inf")],
    labels=["erosion", "stable", "accretion"],
)

print("nombre_lignes_shoreline", len(rates))
print("transects_en_erosion", int((classes == "erosion").sum()))
print("transects_stables", int((classes == "stable").sum()))
print("transects_en_accretion", int((classes == "accretion").sum()))
print("boites_coastsat", len(boxes))
print("zones_vent", len(wind))
print("points_courant_surface", len(currents))
print("localites_cartographie", len(localities))
print("epr_moyen", round(rates["endpoint_rate_myr"].mean(), 3))
print("lignes_table_cinq_facteurs", len(factors))
print("colonnes_facteurs_ok", all(col in factors.columns for col in [
    "shoreline_mean_epr_myr",
    "swv_mean_current_speed_mps",
    "tide_range_m",
    "wave_swh_p95_m",
    "wind_speed_p95_ms",
]))
print("indice_pression_max", round(factors["coastal_dynamics_pressure_index"].max(), 3))

score_cols = ["shoreline_erosion_score", "swv_score", "tide_score", "wave_score", "wind_score"]
weight_scenarios = {
    "balanced": {"shoreline_erosion_score": 0.20, "swv_score": 0.20, "tide_score": 0.20, "wave_score": 0.20, "wind_score": 0.20},
    "observed_erosion_dominant": {"shoreline_erosion_score": 0.50, "swv_score": 0.15, "tide_score": 0.10, "wave_score": 0.15, "wind_score": 0.10},
    "hydrodynamic_forcing": {"shoreline_erosion_score": 0.20, "swv_score": 0.25, "tide_score": 0.15, "wave_score": 0.25, "wind_score": 0.15},
}

scenario_results = factors[["Point"]].copy()
for scenario, weights in weight_scenarios.items():
    assert round(sum(weights.values()), 6) == 1
    scenario_results[scenario] = sum(factors[col] * weight for col, weight in weights.items())

priority_by_scenario = scenario_results.set_index("Point").idxmax()
print("scenarios_ponderation", len(weight_scenarios))
print("priorite_par_scenario")
print(priority_by_scenario.to_string())

assert factors[score_cols].notna().all().all()
assert len(priority_by_scenario) == 3
assert {"Name", "Lon", "Lat", "Zone"}.issubset(localities.columns)
