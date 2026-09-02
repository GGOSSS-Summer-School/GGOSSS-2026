# GGOSSS 2026 — Gulf of Guinea Ocean Sciences Summer School

**Second edition · Cotonou, Benin · 3–9 October 2026**
Hosted by the ICMPA–UNESCO Chair, Université d'Abomey-Calavi.

This repository holds the notebooks, scripts and code used during GGOSSS 2026.
Slides, handouts and large datasets are distributed separately — see
[www.ggosss.org](https://www.ggosss.org).

---

## Getting started

Clone the repository:

```bash
git clone https://github.com/GGOSSS-Summer-School/GGOSSS-2026.git
cd GGOSSS-2026
```

Create the environment:

```bash
conda env create -f environment.yml
conda activate ggosss2026
jupyter lab
```

If you prefer pip:

```bash
python -m venv .venv
source .venv/bin/activate      # Windows: .venv\Scripts\activate
pip install -r requirements.txt
```

**Please do this before you arrive in Cotonou.** Bandwidth on site is limited and
time spent installing software on Day 1 is time taken from the training. If the
installation fails, write to contact@ggosss.org rather than waiting until October.

You will also need a free **Copernicus Marine Service** account —
register at [marine.copernicus.eu](https://marine.copernicus.eu) and check that
your login works before you travel.

---

## Repository structure

```
GGOSSS-2026/
├── 01-python-tutorial/        Python for ocean data analysis
├── 02-in-situ/                Instrumentation, CTD/ADCP processing, quality control
├── 03-satellite/              Satellite remote sensing, altimetry, SST, SSS
├── 04-cmems/                  Copernicus Marine Service — access and applications
├── 05-modelling/              Numerical modelling with CROCO
├── 06-coastal-monitoring/     Shoreline change, coastal processes, machine learning
├── 07-field-campaign/         Lake Nokoué campaign — processing the data collected
├── 08-group-projects/         Trainee group work
├── environment.yml
└── requirements.txt
```

Each session folder contains its own `README.md` describing the session, the
notebooks in the order they are used, and any data the notebooks expect.

---

## A note on data

Keep this repository light. **Do not commit large files.**

- Datasets under about 10 MB may go in `data/`.
- Anything larger is distributed on USB on arrival, or via a download link given
  in the relevant session README.
- Never commit `.nc` files of model or satellite output.

The `.gitignore` excludes common data formats by default. If you need to add a
dataset, raise it with the organising committee first.

---

## For instructors

Work on a branch and open a pull request rather than pushing to `main`:

```bash
git checkout -b session/satellite-tchonang
# ... add your notebooks ...
git add .
git commit -m "Add satellite remote sensing notebooks"
git push -u origin session/satellite-tchonang
```

Before submitting:

- **Clear all outputs** before committing (`Kernel → Restart and Clear Output`).
  Notebooks with embedded figures bloat the repository and produce unreadable diffs.
- **Use relative paths** — `data/nokoue_ctd.nc`, never `/Users/yourname/Desktop/...`.
- **Check it runs from a clean environment**, not just from yours.
- **Add every dependency** to `environment.yml`.
- Write a short `README.md` in your session folder.

Final materials are due **Friday 25 September 2026**.

---

## Previous edition

The first edition, held at the University of Dschang, Cameroon in August 2025:
[GGOSSS-2025](https://github.com/GGOSSS-Summer-School/GGOSSS-2025)

---

## Contact

Gulf of Guinea Ocean Sciences Summer School
[contact@ggosss.org](mailto:contact@ggosss.org) · [www.ggosss.org](https://www.ggosss.org)

---

## Acknowledgments

GGOSSS 2026 is supported by the Partnership for Observation of the Global Ocean
(POGO), the Institut de Recherche pour le Développement (IRD) and the Scientific
Committee on Oceanic Research (SCOR).

The school is also supported by
[OPERA](https://www.unoceanprediction.org/en/ocean-prediction-enhancement-regions-africa-opera)
— Ocean Prediction Enhancement in Regions of Africa — a project funded by the
European Union and implemented by
[Mercator Ocean International](https://www.mercator-ocean.eu/) through its
coordinating role in the
[OceanPrediction Decade Collaborative Centre](https://www.unoceanprediction.org/en/),
and draws throughout on the
[Copernicus Marine Service](https://marine.copernicus.eu/).

Institutional and technical support is provided by the ICMPA–UNESCO Chair at the
Université d'Abomey-Calavi, IRHOB, COESSING and PRODATA SARL.

*Funded by the European Union. Views and opinions expressed are however those of
the author(s) only and do not necessarily reflect those of the European Union.
Neither the European Union nor the granting authority can be held responsible for
them.*
