from pathlib import Path
from tempfile import NamedTemporaryFile
from time import sleep

import numpy as np
import pandas as pd
from rdkit import Chem


ROOT = Path(__file__).resolve().parents[2]
FILTERS = {
    "target_chembl_id": "CHEMBL203",
    "standard_type": "IC50",
    "standard_units": "nM",
    "standard_relation": "=",
}
IDS = ["activity_id", "molecule_chembl_id", "assay_chembl_id"]


def valid_smiles(value):
    if pd.isna(value) or not value.strip():
        return False
    molecule = Chem.MolFromSmiles(value.strip())
    return molecule is not None and molecule.GetNumAtoms() > 0


def screen(df):
    required = {*FILTERS, *IDS, "canonical_smiles", "standard_value"}
    missing = required.difference(df.columns)
    if missing:
        raise ValueError(f"Missing CSV columns: {', '.join(sorted(missing))}")

    values = pd.to_numeric(df["standard_value"], errors="coerce")
    keep = values.gt(0) & np.isfinite(values)
    for column, expected in FILTERS.items():
        keep &= df[column].eq(expected)
    for column in IDS:
        keep &= df[column].astype("string").str.strip().ne("").fillna(False)

    result = df.loc[keep].copy()
    result = result.loc[result["canonical_smiles"].map(valid_smiles).astype(bool)].copy()
    result["pIC50"] = 9 - np.log10(values.loc[result.index])
    return result


def write_if_changed(path, content):
    data = content.encode("utf-8")
    if path.is_file() and path.read_bytes() == data:
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    with NamedTemporaryFile(dir=path.parent, delete=False) as temp:
        temporary = Path(temp.name)
    try:
        temporary.write_bytes(data)
        for attempt in range(6):
            try:
                temporary.replace(path)
                return
            except PermissionError:
                if attempt == 5:
                    raise SystemExit(f"Windows denied replacement of {path}. Close any app holding the file and rerun.")
                sleep(0.1 * 2**attempt)
    finally:
        temporary.unlink(missing_ok=True)


def main():
    path = ROOT / "data" / "raw" / "egfr_sample.csv"
    output = ROOT / "data" / "processed" / "egfr_sample_screened.csv"
    if not path.is_file():
        raise SystemExit(f"Dataset missing. Save the EGFR sample CSV at {path}")
    df = pd.read_csv(path, dtype="string")
    try:
        result = screen(df)
    except ValueError as exc:
        raise SystemExit(str(exc)) from exc

    web_output = ROOT / "frontend" / "public" / "data" / "egfr_sample.json"
    write_if_changed(output, result.to_csv(index=False))
    write_if_changed(web_output, result.to_json(orient="records", double_precision=15))
    print(f"Website data ready: {web_output}")

    print("Total measurements:", len(df))
    print("Unique compounds:", df["molecule_chembl_id"].nunique())
    print("Missing IC50:", df["standard_value"].isna().sum())
    print("Missing SMILES:", df["canonical_smiles"].isna().sum())
    print(f"Saved {len(result)} measurements to {output} ({len(df) - len(result)} excluded)")
    print(result[["molecule_chembl_id", "standard_value", "pIC50"]].head())


if __name__ == "__main__":
    main()
