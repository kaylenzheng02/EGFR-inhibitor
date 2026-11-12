# EGFR Inhibitor

EGFR Inhibitor is an early-stage project for exploring EGFR inhibitor activity data and developing a compound-ranking workflow. The current implementation inspects a local CSV, filters measurements, and calculates pIC50. Fingerprint generation, model training, evaluation, and ranked predictions are not implemented yet.

## Quick start

Run these PowerShell commands from the **repository root (`EGFR-Inhibitor`)**, one directory above this README:

```powershell
python -m venv .venv
.\.venv\Scripts\python.exe -m pip install pandas numpy
.\.venv\Scripts\python.exe src/egfrank/preprocess.py
```

Only pandas and NumPy are required for the current script. The repository-root `requirements.txt` contains a broader pinned environment; `EGFRank/requirements.txt` is currently empty.

## Input data

Place the CSV at `data/raw/egfr_sample.csv`, relative to the repository root. If that file is absent, the script also checks `EGFRank/data/raw/egfr_sample.csv`.

Required columns:

| Column | Description |
| --- | --- |
| `molecule_chembl_id` | Compound identifier used to count unique compounds. |
| `canonical_smiles` | Molecular structure represented as a SMILES string. |
| `standard_value` | Activity value, interpreted as IC50 in nM. |

When present, `standard_type` and `standard_units` are used to retain only `IC50` measurements in `nM`. If these columns are absent, the script assumes the supplied values have that type and unit.

## Current workflow

The preprocessing script:

1. Checks that the input file and required columns exist.
2. Prints total measurements, unique compounds, and missing IC50 and SMILES counts before filtering.
3. Converts activity values to numeric values and removes missing SMILES and missing, nonpositive, or nonfinite activity values.
4. Filters measurement type and units when the corresponding columns are available.
5. Calculates `pIC50 = 9 - log10(IC50 in nM)` and prints the first five retained measurements.

Results are printed to the terminal. The script does not save processed data, validate chemical structures, aggregate repeated compounds, or filter by target or measurement relation.

## Repository layout

```text
EGFR-Inhibitor/
|-- data/raw/egfr_sample.csv       # Current input dataset
|-- src/egfrank/preprocess.py      # Working preprocessing entry point
|-- requirements.txt              # Broader pinned Python environment
|-- frontend/                     # Separate frontend scaffold
`-- EGFRank/
    |-- README.md
    |-- requirements.txt          # Empty placeholder
    |-- data/
    |-- docs/
    |-- models/
    |-- notebooks/
    |-- tests/
    `-- src/egfrank/
        |-- __init__.py
        |-- fetch_data.py         # Initial local CSV inspection script
        |-- preprocess.py         # Placeholder
        |-- fingerprints.py       # Placeholder
        |-- train.py              # Placeholder
        `-- evaluate.py           # Placeholder
```

Despite its name, `EGFRank/src/egfrank/fetch_data.py` reads a local CSV from the nested data directory; it does not download data. Use the repository-root preprocessing entry point shown above for the current workflow.

## Planned work

- Consolidate the source and dependency layout.
- Implement reproducible data acquisition and processed-data export.
- Add chemical structure validation and molecular fingerprints.
- Define compound-level aggregation and dataset splitting.
- Implement and evaluate a baseline compound-ranking model.
- Add tests for data processing and model evaluation.
