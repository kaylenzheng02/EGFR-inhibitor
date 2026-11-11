from pathlib import Path
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
path = ROOT / "data" / "raw" / "egfr_sample.csv"

df = pd.read_csv(path)

print("Total measurements:", len(df))
print("Unique compounds:", df["molecule_chembl_id"].nunique())
print("Missing IC50:", df["standard_value"].isna().sum())
print("Missing SMILES:", df["canonical_smiles"].isna().sum())

# Initial screening only
df["standard_value"] = pd.to_numeric(
    df["standard_value"], errors="coerce"
)

df = df.dropna(
    subset=["canonical_smiles", "standard_value"]
)

df = df[df["standard_value"] > 0].copy()

# IC50 is measured in nM
import numpy as np

df["pIC50"] = 9 - np.log10(df["standard_value"])

print(df[[
    "molecule_chembl_id",
    "standard_value",
    "pIC50"
]].head())