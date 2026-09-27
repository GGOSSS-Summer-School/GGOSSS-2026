# Suivi côtier et machine learning

**Jour 2 - GGOSSS 2026**

Instructeurs: Frédéric Bonou, Nourdi Njutapvoui

## Contenu de la session

Cette session introduit des approches pratiques pour suivre les processus côtiers et évaluer la vulnérabilité climatique au Bénin et dans le golfe de Guinée. Elle relie l'évolution du trait de côte, le transport sédimentaire, les inondations côtières, l'exposition aux surcotes, les échanges lagune-océan et les pressions humaines aux décisions que doivent prendre les scientifiques et gestionnaires du littoral.

Les notebooks de ce dossier servent de supports de démonstration et de prolongement pratique. Ils utilisent des jeux de données préparés sur PAK/Kribi et le littoral camerounais afin d'illustrer comment des observations de terrain, des produits satellites, des données historiques de trait de côte et des cas locaux peuvent être combinés pour identifier le changement côtier et communiquer la vulnérabilité.

## Objectifs pédagogiques

A la fin de la session, les participants doivent pouvoir:

1. Décrire les principaux processus qui contrôlent le changement du trait de côte et les échanges lagune-océan dans le golfe de Guinée.
2. Identifier les facteurs climatiques et humains qui augmentent la vulnérabilité côtière: érosion, inondation, exposition aux surcotes et pression d'occupation du sol.
3. Choisir des indicateurs et sources de données adaptés au suivi du trait de côte et de l'exposition aux aléas.
4. Interpréter un cas d'étude simple en distinguant changement physique, exposition, sensibilité et capacité d'adaptation.
5. Proposer un plan de suivi de base reliant observations, produits de données et besoins de décision pour un site côtier vulnérable.

## Organisation

| Dossier | Description |
|---|---|
| `notebooks/` | Notebooks de TP et scripts de test rapide pour la dynamique côtière, la cartographie et le machine learning. |
| `data/` | Petits jeux de données CSV utilisés par les exercices. |

## Notebooks et scripts

Ordre recommandé:

1. `notebooks/GGOSSS2026_Coastal_Dynamics_PAK_Practical.ipynb`
2. `notebooks/GGOSSS2026_ML_Coastal_Vulnerability_BD_Practical.ipynb`
3. `notebooks/GGOSSS2026_Cartographie_Resultats_Spatiaux.ipynb`
4. Module optionnel CoastSat:
   - `notebooks/GGOSSS2026_CoastSat_Optional_Setup_Download_Extraction.ipynb`
   - `notebooks/smoke_test_coastsat_environment.py`
   - `notebooks/prepare_fes2022_regional_subset.py`
5. Scripts optionnels de vérification:
   - `notebooks/smoke_test_pak_coastal_dynamics.py`
   - `notebooks/smoke_test_bd_vulnerability_ml.py`

Les tests rapides peuvent être lancés depuis ce dossier:

```bash
python notebooks/smoke_test_pak_coastal_dynamics.py
python notebooks/smoke_test_bd_vulnerability_ml.py
```

La sortie ML attendue contient trois sections: validation aléatoire avec tous les prédicteurs, sensibilité au leakage après retrait des indices dérivés, et validation spatiale avec le jeu de prédicteurs réduit. Les scores peuvent varier légèrement selon la plateforme, mais les six modèles doivent tourner:

```text
random_validation_all_predictors
ANN_MLP accuracy ...
...

leakage_sensitivity_f1_drop
...

spatial_validation_reduced_predictors
...
```

La sortie PAK attendue contient le nombre de transects en érosion/stabilité/accrétion et le secteur prioritaire selon trois scénarios de pondération.

Pour la session orale d'une heure, ces notebooks ne sont pas obligatoires à exécuter en entier. Ils peuvent être utilisés comme démonstrations courtes, support de discussion, ou base de travail pour les activités de groupe après la session.

Le module CoastSat est séparé de l'environnement `ggosss2026`. Il doit être lancé dans un environnement `coastsat` déjà configuré avec CoastSat, Earth Engine et PyFES:

```bash
conda activate coastsat
python notebooks/smoke_test_coastsat_environment.py
```

Le notebook CoastSat désactive les téléchargements par défaut. Les interrupteurs `RUN_CHECK_AVAILABLE`, `RUN_DOWNLOAD` et `RUN_EXTRACTION` doivent être passés explicitement à `True` par l'instructeur pour éviter de créer des fichiers lourds pendant la session.

## Données

Les notebooks utilisent ces petits fichiers locaux:

- `data/cameroon_coastal_vulnerability_ml.csv`
- `data/cameroon_coastal_vulnerability_ml_sample.csv`
- `data/pak_shoreline_change_rates_2026.csv`
- `data/pak_coastsat_box_inventory.csv`
- `data/pak_coastsat_transect_inventory.csv`
- `data/pak_wind_zonal_stats.csv`
- `data/pak_current_points_stats.csv`
- `data/pak_integrated_coastal_factors.csv`
- `data/pak_priorites_decision_suivi.csv`
- `data/pak_typologie_processus_workshop_par_zone.csv`
- `data/selected_localities_for_maps.csv`
- `data/cameroon_ml_predictions_spatiales.csv`
- `data/cameroon_ml_predictions_spatiales_sample.csv`

Les slides, handouts, cartes interactives HTML et gros fichiers de référence sont distribués séparément, conformément aux consignes du dépôt GGOSSS. Le notebook cartographique utilise Cartopy/Natural Earth pour les cartes statiques; un premier lancement peut nécessiter Internet si les couches Natural Earth ne sont pas encore en cache.

## Avant la session

Les participants doivent créer l'environnement commun GGOSSS 2026 avant l'école, comme demandé dans le README racine du dépôt:

```bash
conda env create -f environment.yml
conda activate ggosss2026
jupyter lab
```

La session orale ne nécessite pas de logiciel en direct. Les notebooks de démonstration utilisent les packages déjà listés dans l'environnement du dépôt: `numpy`, `pandas`, `matplotlib`, `scikit-learn`, `cartopy`, `pyproj` et `jupyterlab`.

Le module optionnel CoastSat nécessite un environnement séparé `coastsat`. Les données FES2022 ne doivent pas être committées dans le dépôt; si elles sont disponibles localement, le script `notebooks/prepare_fes2022_regional_subset.py` permet de préparer un sous-ensemble Bénin/lac Nokoué, PAK/Kribi ou golfe de Guinée. Ce sous-ensemble requiert `xarray` et `netcdf4` dans l'environnement CoastSat:

```bash
conda install -n coastsat -c conda-forge xarray netcdf4
```

## Lancer les notebooks

Depuis la racine du dépôt:

```bash
cd 06-coastal-monitoring
conda activate ggosss2026
jupyter lab
```

Les notebooks utilisent uniquement des chemins relatifs.

## Figures pedagogiques

La figure de typologie par zone est reproductible a partir de `data/pak_typologie_processus_workshop_par_zone.csv` dans le notebook cartographique. Elle reprend le style des figures AHP/CVI existantes, mais elle n'est pas un nouveau calcul AHP/CVI: elle synthetise les cinq facteurs du TP PAK/Kribi.
