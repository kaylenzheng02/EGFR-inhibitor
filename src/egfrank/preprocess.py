from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[2]


def main():
    path = ROOT / "data" / "raw" / "egfr_sample.csv"
    nested_path = ROOT / "EGFRank" / "data" / "raw" / path.name
    if not path.is_file() and nested_path.is_file():
        path = nested_path
    if not path.is_file():
        raise SystemExit(f"Dataset missing. Save the EGFR sample CSV at {path}")

    df = pd.read_csv(path)
    required = {"molecule_chembl_id", "standard_value", "canonical_smiles"}
    missing = required.difference(df.columns)
    if missing:
        raise SystemExit(f"Missing CSV columns: {', '.join(sorted(missing))}")

    print("Total measurements:", len(df))
    print("Unique compounds:", df["molecule_chembl_id"].nunique())
    print("Missing IC50:", df["standard_value"].isna().sum())
    print("Missing SMILES:", df["canonical_smiles"].isna().sum())

    df["standard_value"] = pd.to_numeric(df["standard_value"], errors="coerce")
    df = df.dropna(subset=["canonical_smiles", "standard_value"])
    if "standard_type" in df:
        df = df[df["standard_type"] == "IC50"]
    if "standard_units" in df:
        df = df[df["standard_units"] == "nM"]
    df = df[df["standard_value"].gt(0) & np.isfinite(df["standard_value"])].copy()
    df["pIC50"] = 9 - np.log10(df["standard_value"])
    print(df[["molecule_chembl_id", "standard_value", "pIC50"]].head())


if __name__ == "__main__":
    main()
