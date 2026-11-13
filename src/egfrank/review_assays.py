"""Review EGFR assays and export one activity label per eligible compound."""
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
    # Extend the original three with descriptions explicitly reporting a substrate
    # assay for wild-type EGFR. Mixed variant descriptions require manual review.
    def mentions(pattern):
        return descriptions.str.contains(pattern, case=False, regex=True)

    wild_type = mentions(r"wild[ -]?type|\bWT\b")
    egfr = mentions(r"EGFR|EGF[ -]?R|epidermal growth factor receptor")
    enzyme = mentions(r"substrate|ADP[ -]?Glo|Z[ -]?Lyte")
    ambiguous = mentions(r"mutant|mutation|\b[A-Z]\d{2,4}[A-Z]\b|exon|deletion|insertion|proliferation|viability|binding affinity|tracer|degradation|in.cells|whole.cell|cell.based")
    reviewed |= wild_type & egfr & enzyme & ~ambiguous & data.assay_type.eq("B")
    audit = data[["assay_chembl_id", "assay_description"]].copy()
    audit["candidate"] = reviewed
    audit.drop_duplicates().to_csv(ROOT / "docs/research/phase2-assay-review.csv", index=False)
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
    compounds = rows.groupby("molecule_chembl_id")
    if compounds.canonical_smiles.nunique().gt(1).any():
        raise ValueError("Resolve inconsistent structures before exporting compound labels")
    dataset = grouped.rename(columns={"min": "min_ic50_nm", "max": "max_ic50_nm", "count": "measurements"})
    dataset["canonical_smiles"] = compounds.canonical_smiles.first()
    dataset["activity_ids"] = compounds.activity_id.agg(lambda ids: "|".join(ids.astype(str)))
    dataset["assay_ids"] = compounds.assay_chembl_id.agg(lambda ids: "|".join(sorted(set(ids))))
    active = dataset.max_ic50_nm.le(1000)
    inactive = dataset.min_ic50_nm.gt(1000)
    dataset["label"] = np.select([active, inactive], ["active", "inactive"], default="conflict")
    output = ROOT / "data/processed"
    output.mkdir(parents=True, exist_ok=True)
    dataset.loc[active | inactive].to_csv(output / "egfr_training.csv")
    dataset.loc[~(active | inactive)].to_csv(output / "egfr_conflicts.csv")
    print("Screened measurements:", len(rows), "compounds:", len(grouped))
    print("Repeated compounds:", (grouped["count"] > 1).sum())
    print("All values <=1000:", (grouped["max"] <= 1000).sum())
    print("All values >1000:", (grouped["min"] > 1000).sum())
    print("Cross cutoff:", ((grouped["min"] <= 1000) & (grouped["max"] > 1000)).sum())


if __name__ == "__main__":
    main()
