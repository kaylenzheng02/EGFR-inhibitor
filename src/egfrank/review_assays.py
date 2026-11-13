"""Preview compound counts for a conservative set of reviewed EGFR assay descriptions."""
from pathlib import Path

import numpy as np
import pandas as pd
from rdkit import Chem

ROOT = Path(__file__).resolve().parents[2]
SNAPSHOT = ROOT / "data/raw/egfr_ic50_20261113T184940536954Z/egfr_ic50.csv"


def main():
    data = pd.read_csv(SNAPSHOT, low_memory=False)
    descriptions = data.assay_description.fillna("")
    reviewed = descriptions.isin([
        "DRUGMATRIX: Protein Tyrosine Kinase, EGF Receptor enzyme inhibition (substrate: Poly(Glu:Tyr))",
        "Inhibition of recombinant human N-terminal GST His6-tagged thrombin cleavage site-fused EGFR (H672 to A1210 residues) expressed in baculovirus using poly-EY as substrate after 210 mins by luciferase coupled chemiluminescence assay",
    ]) | descriptions.str.startswith("High Throughput Biochemical Assay: EGFR activity was measured using KinEASE")
    rows = data.loc[reviewed].copy()
    print("Reviewed description groups:", rows.assay_description.nunique())
    print("Candidate assay IDs:", rows.assay_chembl_id.nunique(), "rows:", len(rows))
    values = pd.to_numeric(rows.standard_value, errors="coerce")
    keep = (rows.standard_relation.eq("=") & rows.standard_units.eq("nM")
            & values.gt(0) & np.isfinite(values) & rows.data_validity_comment.isna()
            & rows.potential_duplicate.eq(0) & rows.molecule_chembl_id.notna())
    rows = rows.loc[keep].copy()
    rows["ic50_nm"] = values.loc[keep]
    valid = {s for s in rows.canonical_smiles.dropna().unique()
             if s.strip() and Chem.MolFromSmiles(s) is not None}
    rows = rows[rows.canonical_smiles.isin(valid)]
    grouped = rows.groupby("molecule_chembl_id").ic50_nm.agg(["min", "max", "count"])
    print("Screened measurements:", len(rows), "compounds:", len(grouped))
    print("Repeated compounds:", (grouped["count"] > 1).sum())
    print("All values <=1000:", (grouped["max"] <= 1000).sum())
    print("All values >1000:", (grouped["min"] > 1000).sum())
    print("Cross cutoff:", ((grouped["min"] <= 1000) & (grouped["max"] > 1000)).sum())


if __name__ == "__main__":
    main()
