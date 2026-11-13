# EGFR Inhibitor

EGFR Inhibitor is an early-stage project for exploring EGFR inhibitor activity data and developing a compound-ranking workflow. The current implementation inspects a local CSV, filters measurements, and calculates pIC50. Fingerprint generation, model training, evaluation, and ranked predictions are not implemented yet.

Read the [project history](PROJECT_HISTORY.md) for the research question, dataset preparation, results, and limitations.

## Quick start

Run these PowerShell commands from the **repository root (`EGFR-Inhibitor`)**:

```powershell
python -m venv .venv
.\.venv\Scripts\python.exe -m pip install pandas numpy rdkit
.\.venv\Scripts\python.exe src/egfrank/preprocess.py
```

The preprocessor uses pandas, NumPy, and RDKit. The repository-root `requirements.txt` contains a broader pinned environment.

## Input data

Place the CSV at `data/raw/egfr_sample.csv`, relative to the repository root.

Each row requires an activity ID, compound ID, and assay ID (`activity_id`, `molecule_chembl_id`, and `assay_chembl_id`), a parseable `canonical_smiles` structure, and a positive, finite `standard_value`. The screen requires `target_chembl_id=CHEMBL203`, `standard_type=IC50`, an exact `standard_relation` of `=`, and `standard_units=nM`.

## Current workflow

The script adds `pIC50 = 9 - log10(IC50 in nM)` and writes passing rows to `data/processed/egfr_sample_screened.csv`. It preserves the input columns and keeps repeated measurements as separate rows. The terminal shows the raw dataset counts, the number of rows excluded, and the first five retained measurements.

## Repository layout

```text
EGFR-Inhibitor/
|-- data/raw/egfr_sample.csv       # Current input dataset
|-- src/egfrank/preprocess.py      # Working preprocessing entry point
|-- requirements.txt              # Broader pinned Python environment
|-- frontend/                     # React and TypeScript app
|-- docs/
|-- models/
|-- notebooks/
`-- tests/
```

The remaining Python modules live alongside the working preprocessor in `src/egfrank/`. Run `.\.venv\Scripts\python.exe src/egfrank/fetch_data.py` from the repository root to download the full ChEMBL EGFR IC50 query. Each successful run creates a dated directory under `data/raw` containing `egfr_ic50.csv` and `manifest.json`. The manifest records the query, retrieval time, page and ID counts, and SHA-256 checksum. The downloader preserves different units and measurement relations for later screening. Fingerprinting, training, and evaluation are placeholders.

## Planned work


- Add molecular fingerprints.
- Define compound-level aggregation and dataset splitting.
- Implement and evaluate a baseline compound-ranking model.
- Add tests for data processing and model evaluation.
