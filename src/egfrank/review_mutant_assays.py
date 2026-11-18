"""Review EGFR mutant assays and export one label per compound per variant.

NSCLC is driven by specific EGFR variants, not by wild-type receptor: the
sensitizing mutations L858R and exon 19 deletions, and the resistance
mutations T790M and C797S. review_assays.py deliberately screens those out.
This script screens them in, keeping one variant per label so that potency
against different variants is never averaged together.

The snapshot does not carry ChEMBL's variant_sequence field, so the variant is
read from the assay description. A description that names a mutant alongside
wild-type cannot be attributed to one entity from text alone and is excluded.
"""
from pathlib import Path

import numpy as np
import pandas as pd
from rdkit import Chem

ROOT = Path(__file__).resolve().parents[2]
SNAPSHOT = ROOT / "data/raw/egfr_ic50_20261113T184940536954Z/egfr_ic50.csv"
CUTOFF_NM = 1000

# Variants that matter in NSCLC: sensitizing, then acquired resistance.
VARIANTS = {
    "L858R": r"L858R",
    "T790M": r"T790M",
    "C797S": r"C797S",
    "ex19del": r"exon\s*19|delE746|del\s*E746|ex19del|E746[_ -]?A750",
    "L861Q": r"L861Q",
    "G719x": r"G719[SAC]",
}


def main():
    data = pd.read_csv(SNAPSHOT, low_memory=False)
    descriptions = data.assay_description.fillna("")

    def mentions(pattern):
        return descriptions.str.contains(pattern, case=False, regex=True)

    # Variant signature: every listed mutation the description names, joined.
    found = pd.DataFrame({name: mentions(p) for name, p in VARIANTS.items()})
    signature = found.apply(lambda row: "+".join(found.columns[row.values]), axis=1)

    wild_type = mentions(r"wild[ -]?type|\bWT\b")
    cell_based = mentions(
        r"proliferation|viability|in .{0,12}cells|cell[- ]based"
        r"|\bPC-?9\b|\bHCC827\b|H1975|H3255|H1650|HCC4006|A431|A549"
    )
    enzyme = mentions(
        r"substrate|poly\s*\(?\s*(?:Glu|E)|polyE|Omnia|Z.?-?Lyte|ADP[- ]?Glo"
        r"|KinEASE|TR-?FRET|kinase (?:assay|inhibition)"
    )
    # A mutant named beside wild-type in one protocol cannot be attributed from
    # text alone; a cell readout answers a different question than the enzyme.
    candidate = signature.ne("") & ~wild_type & ~cell_based & enzyme & data.assay_type.eq("B")

    audit = data[["assay_chembl_id", "assay_description"]].copy()
    audit["variant"] = signature
    audit["candidate"] = candidate
    audit.drop_duplicates().to_csv(ROOT / "docs/research/mutant-assay-review.csv", index=False)

    rows = data.loc[candidate].copy()
    rows["variant"] = signature.loc[candidate]
    print("Candidate assay IDs:", rows.assay_chembl_id.nunique(), "rows:", len(rows))
    print("Excluded, mutant named beside wild-type:",
          data.loc[signature.ne("") & wild_type, "assay_chembl_id"].nunique(), "assays")

    values = pd.to_numeric(rows.standard_value, errors="coerce")
    keep = (rows.standard_relation.eq("=") & rows.standard_units.eq("nM")
            & values.gt(0) & np.isfinite(values) & rows.data_validity_comment.isna()
            & rows.potential_duplicate.eq(0) & rows.molecule_chembl_id.notna())
    rows = rows.loc[keep].copy()
    rows["ic50_nm"] = values.loc[keep]
    valid = {s for s in rows.canonical_smiles.dropna().unique()
             if s.strip() and Chem.MolFromSmiles(s) is not None}
    rows = rows[rows.canonical_smiles.isin(valid)]

    # One row per compound PER VARIANT: potency against L858R and against
    # L858R/T790M are different measurements of different entities.
    pairs = rows.groupby(["molecule_chembl_id", "variant"])
    grouped = pairs.ic50_nm.agg(["min", "max", "count"])
    if pairs.canonical_smiles.nunique().gt(1).any():
        raise ValueError("Resolve inconsistent structures before exporting compound labels")
    dataset = grouped.rename(columns={"min": "min_ic50_nm", "max": "max_ic50_nm", "count": "measurements"})
    dataset["canonical_smiles"] = pairs.canonical_smiles.first()
    dataset["activity_ids"] = pairs.activity_id.agg(lambda ids: "|".join(ids.astype(str)))
    dataset["assay_ids"] = pairs.assay_chembl_id.agg(lambda ids: "|".join(sorted(set(ids))))
    active = dataset.max_ic50_nm.le(CUTOFF_NM)
    inactive = dataset.min_ic50_nm.gt(CUTOFF_NM)
    dataset["label"] = np.select([active, inactive], ["active", "inactive"], default="conflict")

    output = ROOT / "data/processed"
    output.mkdir(parents=True, exist_ok=True)
    dataset.loc[active | inactive].to_csv(output / "egfr_mutant_training.csv")
    dataset.loc[~(active | inactive)].to_csv(output / "egfr_mutant_conflicts.csv")

    labeled = dataset.loc[active | inactive]
    print("\nScreened measurements:", len(rows))
    print("Compound-variant pairs:", len(dataset), "labeled:", len(labeled),
          "conflicts:", len(dataset) - len(labeled))
    print("Distinct compounds:", labeled.index.get_level_values(0).nunique())
    print("\nBy variant:")
    summary = labeled.groupby(["variant", "label"]).size().unstack(fill_value=0)
    print(summary.to_string())

    # Wild-type comparison: the same compound measured against both.
    wt_path = output / "egfr_training.csv"
    if wt_path.exists():
        wt = pd.read_csv(wt_path)
        shared = set(wt.molecule_chembl_id) & set(labeled.index.get_level_values(0))
        print("\nCompounds in both the wild-type and mutant sets:", len(shared))


if __name__ == "__main__":
    main()
