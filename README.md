# EGFRank

EGFRank uses chemical structures to rank compounds by their likelihood of inhibiting human EGFR in laboratory assays. Non-small cell lung cancer (NSCLC) motivates the choice of target. The model learns from published EGFR measurements in ChEMBL, and its benchmark evaluates laboratory activity. Predicting patient response would require clinical data.

## Demo

Download or open [`demo/EGFRank-Demo.html`](demo/EGFRank-Demo.html) in a browser. This standalone demo includes its data and needs no installation. See the [demo instructions](demo/README.md).

## Dataset and model

The [project history](PROJECT_HISTORY.md) documents how 26,600 ChEMBL records yielded 1,600 compounds with consistent activity labels. A compound qualified as active when all its eligible IC50 measurements were at or below 1,000 nM. Of the labeled compounds, 1,396 were active and 204 inactive.

RDKit represents each structure as a 2,048-bit Morgan fingerprint with radius 2 and chirality enabled. A 100-tree Random Forest learns from 1,280 compounds. The remaining 320 form the test set. The split keeps each chemical scaffold (the core molecular framework) entirely on one side so the benchmark tests compounds with scaffolds absent from training.

## First benchmark

| Method | Correct predictions | Accuracy | ROC-AUC | Average precision (active class) |
| --- | ---: | ---: | ---: | ---: |
| Always predict active | 279 / 320 | 87.2% | 0.500 | 0.872 |
| Random Forest | 306 / 320 | 95.6% | 0.929 | 0.981 |

| Experimental label | Predicted inactive | Predicted active |
| --- | ---: | ---: |
| Inactive | 29 | 12 |
| Active | 2 | 277 |

The forest identified 29 of 41 inactive compounds while retaining 277 of 279 active compounds. For the active class, precision was 0.96, recall was 0.99, and F1 was 0.98. ROC-AUC and average precision summarize the ordering of the full test set by score. The number of experimentally active compounds among the top 12 remains to be calculated.

## Run the benchmark in PowerShell

Run these commands from the repository root. For a new environment:

```powershell
python -m venv .venv
.\.venv\Scripts\python.exe -m pip install pandas numpy rdkit scikit-learn joblib
```

Starting with the labeled dataset at `data/processed/egfr_training.csv`, run:

```powershell
.\.venv\Scripts\python.exe src/egfrank/fingerprints.py
.\.venv\Scripts\python.exe src/egfrank/split_data.py
.\.venv\Scripts\python.exe src/egfrank/train.py
.\.venv\Scripts\python.exe src/egfrank/evaluate.py
```

Training saves `models/egfr_random_forest.joblib`. If the saved model, fingerprints and split are already present, run just the evaluation command to print both classification reports and confusion matrices.

To rebuild the labels, run `src/egfrank/review_assays.py` with the project Python interpreter. It expects the original local CSV at `data/raw/egfr_ic50_20261113T184940536954Z/egfr_ic50.csv`. The [snapshot manifest](data/raw/egfr_ic50_20261113T184940536954Z/manifest.json) records its checksum. The [downloader](src/egfrank/fetch_data.py) creates a new dated snapshot on each run. To analyze a new snapshot, update the review script’s `SNAPSHOT` path first. A newer ChEMBL release can yield different counts and scores.

## Explore the pilot data

The React app displays experimental data from the original pilot. From the repository root:

```powershell
cd frontend
npm.cmd install
npm.cmd run dev
```

Open the local URL printed by Vite. To refresh the pilot from `data/raw/egfr_sample.csv`, run `src/egfrank/preprocess.py` with the project Python interpreter from the repository root. It screens measurements and calculates `pIC50 = 9 - log10(IC50 in nM)`.

The next milestone is to export the model’s ranked test compounds and show their scores alongside the experimental labels in the dashboard.
