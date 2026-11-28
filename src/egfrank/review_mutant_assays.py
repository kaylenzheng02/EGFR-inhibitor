"""Review EGFR mutant assays and export compound activity labels.

Running without arguments preserves the original all-variant screen. Running
with ``--variant L858R`` applies the fixed single-L858R target specification,
writes a reviewable assay audit, and exports provisional L858R-only labels.
The L858R exports remain provisional until every automated candidate assay has
an explicit manual include or exclude decision.
"""

import argparse
import json
import re
from pathlib import Path

import numpy as np
import pandas as pd
from rdkit import Chem

ROOT = Path(__file__).resolve().parents[2]
SNAPSHOT = ROOT / "data/raw/egfr_ic50_20261113T184940536954Z/egfr_ic50.csv"
CUTOFF_NM = 1000

# These patterns preserve the original all-variant screen when no CLI target
# is supplied.
LEGACY_VARIANTS = {
    "L858R": r"L858R",
    "T790M": r"T790M",
    "C797S": r"C797S",
    "ex19del": r"exon\s*19|delE746|del\s*E746|ex19del|E746[_ -]?A750",
    "L861Q": r"L861Q",
    "G719x": r"G719[SAC]",
}

# Canonical order makes T790M/L858R and L858R/T790M the same signature.
VARIANT_PATTERNS = {
    "ex19del": (
        r"exon\s*19[^.;,]{0,24}(?:del|delet)|\bex\s*19\s*del\b|"
        r"\bdel\s*19\b|\b19\s*del\b|\bE746[_ -]?A750(?:del)?\b|"
        r"\bdel\s*E?746[_ -]?(?:A?750)?\b|\bd746[_ -]?750\b"
    ),
    "L858R": r"\bL\s*858\s*R\b|\bLeu(?:cine)?\s*858\s*Arg(?:inine)?\b",
    "T790M": r"\bT\s*790\s*M\b|\bThr(?:eonine)?\s*790\s*Met(?:hionine)?\b",
    "C797S": r"\bC\s*797\s*S\b|\bCys(?:teine)?\s*797\s*Ser(?:ine)?\b",
    "L861Q": r"\bL\s*861\s*Q\b|\bLeu(?:cine)?\s*861\s*Gln(?:utamine)?\b",
    "G719X": r"\bG\s*719\s*[SAC]\b",
}

KNOWN_SUBSTITUTIONS = {
    "L858R", "T790M", "C797S", "L861Q", "G719S", "G719A", "G719C"
}
AMINO_ACID_CODES = "ACDEFGHIKLMNPQRSTVWY"
SUBSTITUTION_PATTERN = re.compile(
    rf"\b([{AMINO_ACID_CODES}])\s*(\d{{3,4}})\s*([{AMINO_ACID_CODES}])\b",
    re.I,
)
WILD_TYPE_PATTERN = re.compile(r"wild[ -]?type|\bWT\b", re.I)

CELL_PATTERN = re.compile(
    r"proliferation|viability|cell[- ]based|cellular|in .{0,16}cells|lysate|"
    r"\bPC-?9\b|\bHCC827\b|H1975|H3255|H1650|HCC4006|A431|A549",
    re.I,
)
BIOCHEMICAL_PATTERN = re.compile(
    r"substrate|poly\s*\(?\s*(?:Glu|E)|polyE|Omnia|Z.?-?Lyte|"
    r"ADP[- ]?Glo|KinEASE|Kinase[- ]?Glo|TR-?FRET|HTRF|radiometric|"
    r"scintillation|ATP incorporation|kinase (?:assay|inhibition|activity)",
    re.I,
)
BINDING_PATTERN = re.compile(
    r"binding affinity|surface plasmon|\bSPR\b|tracer|displacement|"
    r"dissociation constant|\bK[dD]\b",
    re.I,
)
DEGRADATION_PATTERN = re.compile(r"degrad|PROTAC", re.I)

AUDIT_PATH = ROOT / "docs/research/l858r-assay-review.csv"
DECISIONS_PATH = ROOT / "docs/research/l858r-assay-decisions.csv"
PROCESSED = ROOT / "data/processed"


class IncompleteReviewError(RuntimeError):
    """Raised when strict output is requested before assay review is complete."""


def parse_variant_signature(description):
    """Return a canonical mutation signature and unrecognized substitutions."""
    text = "" if pd.isna(description) else str(description)
    variants = [
        name for name, pattern in VARIANT_PATTERNS.items()
        if re.search(pattern, text, flags=re.I)
    ]

    substitutions = {
        f"{left.upper()}{position}{right.upper()}"
        for left, position, right in SUBSTITUTION_PATTERN.findall(text)
    }
    unknown = sorted(substitutions - KNOWN_SUBSTITUTIONS)

    # Named insertions and non-exon-19 alterations also make a target compound.
    if re.search(
        r"\binsertion\b|\b[A-Z]\d{2,4}(?:_[A-Z]\d{2,4})?ins[A-Z]{1,8}\b",
        text,
        re.I,
    ):
        unknown.append("named_insertion")
    if re.search(r"\bexon\s*(?!19\b)\d+", text, re.I):
        unknown.append("other_exon_alteration")

    return "+".join(variants), "|".join(sorted(set(unknown)))


def classify_endpoint(description):
    """Classify the biological endpoint conservatively from assay text."""
    text = "" if pd.isna(description) else str(description)
    if DEGRADATION_PATTERN.search(text):
        return "degradation"
    if re.search(r"proliferation", text, re.I):
        return "proliferation"
    if re.search(r"viability", text, re.I):
        return "viability"
    if BINDING_PATTERN.search(text):
        return "binding"
    if CELL_PATTERN.search(text) or re.search(r"autophosphorylation", text, re.I):
        return "cellular_phosphorylation"
    if BIOCHEMICAL_PATTERN.search(text):
        return "biochemical_enzyme"
    return "ambiguous"


def classify_l858r_assay(description, assay_type):
    """Return target metadata and the automated inclusion decision."""
    signature, unknown = parse_variant_signature(description)
    endpoint = classify_endpoint(description)
    mentions_wild_type = bool(WILD_TYPE_PATTERN.search(str(description)))

    if "L858R" not in signature.split("+"):
        reason = "does_not_name_L858R"
    elif signature != "L858R":
        reason = "combination_variant"
    elif unknown:
        reason = "unrecognized_additional_alteration"
    elif mentions_wild_type:
        reason = "mixed_wild_type_and_mutant"
    elif endpoint != "biochemical_enzyme":
        reason = f"endpoint_{endpoint}"
    elif assay_type != "B":
        reason = "assay_type_not_B"
    else:
        reason = ""

    return {
        "parsed_variant": signature,
        "unrecognized_alterations": unknown,
        "endpoint_type": endpoint,
        "mentions_wild_type": mentions_wild_type,
        "automated_candidate": reason == "",
        "exclusion_reason": reason,
    }


def label_from_range(minimum, maximum):
    if maximum <= CUTOFF_NM:
        return "active"
    if minimum > CUTOFF_NM:
        return "inactive"
    return "conflict"


def measurement_exclusion_reason(row):
    """Return all activity-level reasons that make a row ineligible."""
    reasons = []
    value = pd.to_numeric(pd.Series([row.get("standard_value")]), errors="coerce").iloc[0]
    if row.get("standard_type") != "IC50":
        reasons.append("standard_type_not_IC50")
    if row.get("standard_relation") != "=":
        reasons.append("relation_not_exact")
    if row.get("standard_units") != "nM":
        reasons.append("units_not_nM")
    if pd.isna(value) or not np.isfinite(value) or value <= 0:
        reasons.append("value_not_positive_finite")
    validity = row.get("data_validity_comment")
    if pd.notna(validity) and str(validity).strip():
        reasons.append("data_validity_comment")
    duplicate = pd.to_numeric(
        pd.Series([row.get("potential_duplicate")]), errors="coerce"
    ).iloc[0]
    if pd.isna(duplicate) or duplicate != 0:
        reasons.append("potential_duplicate")
    compound_id = row.get("molecule_chembl_id")
    if pd.isna(compound_id) or not str(compound_id).strip():
        reasons.append("missing_compound_id")
    smiles = row.get("canonical_smiles")
    if pd.isna(smiles) or not str(smiles).strip():
        reasons.append("missing_smiles")
    elif Chem.MolFromSmiles(str(smiles)) is None:
        reasons.append("invalid_smiles")
    return "|".join(reasons)


def _load_manual_decisions(audit):
    """Preserve human decisions and notes when the generated audit is rerun."""
    columns = ["assay_chembl_id", "manual_decision", "review_notes"]
    sources = []
    if AUDIT_PATH.exists():
        prior = pd.read_csv(AUDIT_PATH, dtype=str).reindex(columns=columns)
        sources.append(prior[prior.manual_decision.isin(["include", "exclude"])])
    if DECISIONS_PATH.exists():
        decisions = pd.read_csv(DECISIONS_PATH, dtype=str).reindex(columns=columns)
        if decisions.assay_chembl_id.isna().any() or decisions.assay_chembl_id.duplicated().any():
            raise ValueError(f"Missing or duplicate assay IDs in {DECISIONS_PATH.name}")
        sources.append(decisions)

    if sources:
        prior = pd.concat(sources, ignore_index=True).drop_duplicates(
            "assay_chembl_id", keep="last"
        )
        audit = audit.merge(prior, on="assay_chembl_id", how="left")
    else:
        audit["manual_decision"] = pd.NA
        audit["review_notes"] = pd.NA

    defaults = pd.Series(
        np.where(audit.automated_candidate, "pending", "not_applicable"),
        index=audit.index,
    )
    # Only explicit human include/exclude decisions survive a classifier
    # change. Generated pending/not-applicable states follow the new result.
    explicit = audit.manual_decision.isin(["include", "exclude"])
    audit.loc[~explicit, "manual_decision"] = defaults.loc[~explicit]
    audit["review_notes"] = audit.review_notes.fillna("")
    allowed = {"pending", "include", "exclude", "not_applicable"}
    invalid = sorted(set(audit.manual_decision) - allowed)
    if invalid:
        raise ValueError(f"Invalid manual decisions in {AUDIT_PATH.name}: {invalid}")
    return audit


def build_l858r_audit(data):
    """Build one review row per assay that names L858R in any signature."""
    assay_rows = data[["assay_chembl_id", "assay_description", "assay_type"]].drop_duplicates()
    description_counts = assay_rows.groupby("assay_chembl_id").assay_description.nunique(
        dropna=False
    )
    if description_counts.gt(1).any():
        ambiguous = description_counts[description_counts.gt(1)].index.tolist()
        raise ValueError(f"Assay IDs have inconsistent descriptions: {ambiguous[:5]}")
    assay_rows = assay_rows.drop_duplicates("assay_chembl_id").copy()
    classified = assay_rows.apply(
        lambda row: classify_l858r_assay(row.assay_description, row.assay_type),
        axis=1,
        result_type="expand",
    )
    audit = pd.concat(
        [assay_rows.reset_index(drop=True), classified.reset_index(drop=True)], axis=1
    )
    audit["parsed_variant"] = audit.parsed_variant.fillna("")
    audit["unrecognized_alterations"] = audit.unrecognized_alterations.fillna("")
    audit["exclusion_reason"] = audit.exclusion_reason.fillna("")
    audit = audit[
        audit.parsed_variant.str.split("+").apply(lambda values: "L858R" in values)
    ].copy()
    audit = _load_manual_decisions(audit)
    columns = [
        "assay_chembl_id", "assay_description", "assay_type", "parsed_variant",
        "unrecognized_alterations", "endpoint_type", "mentions_wild_type",
        "automated_candidate", "exclusion_reason", "manual_decision", "review_notes",
    ]
    return audit[columns].sort_values("assay_chembl_id")


def _join_unique(values):
    return "|".join(sorted({str(value) for value in values if pd.notna(value)}))


def _build_l858r_dataset(rows):
    if rows.activity_id.duplicated().any():
        duplicates = rows.loc[rows.activity_id.duplicated(), "activity_id"].tolist()
        raise ValueError(f"Duplicate activity IDs in eligible L858R rows: {duplicates[:5]}")

    pairs = rows.groupby("molecule_chembl_id")
    if pairs.canonical_smiles.nunique().gt(1).any():
        raise ValueError("Resolve inconsistent structures before exporting L858R labels")

    grouped = pairs.ic50_nm.agg(["min", "max", "count"])
    dataset = grouped.rename(
        columns={"min": "min_ic50_nm", "max": "max_ic50_nm", "count": "measurements"}
    )
    dataset["variant"] = "L858R"
    dataset["canonical_smiles"] = pairs.canonical_smiles.first()
    dataset["activity_ids"] = pairs.activity_id.agg(_join_unique)
    dataset["assay_ids"] = pairs.assay_chembl_id.agg(_join_unique)
    dataset["document_ids"] = pairs.document_chembl_id.agg(_join_unique)
    assay_counts = pairs.assay_chembl_id.nunique()
    document_counts = pairs.document_chembl_id.nunique()
    dataset["evidence_grade"] = np.select(
        [dataset.measurements.eq(1), assay_counts.gt(1) | document_counts.gt(1)],
        ["C", "A"],
        default="B",
    )
    dataset["label"] = [
        label_from_range(low, high)
        for low, high in zip(dataset.min_ic50_nm, dataset.max_ic50_nm)
    ]
    return dataset.reset_index()[[
        "molecule_chembl_id", "variant", "canonical_smiles", "label",
        "min_ic50_nm", "max_ic50_nm", "measurements", "evidence_grade",
        "activity_ids", "assay_ids", "document_ids",
    ]]


def run_l858r(data, require_reviewed=False):
    """Write the review audit and provisional single-L858R label exports."""
    audit = build_l858r_audit(data)
    AUDIT_PATH.parent.mkdir(parents=True, exist_ok=True)
    audit.to_csv(AUDIT_PATH, index=False)

    pending = audit.loc[audit.automated_candidate & audit.manual_decision.eq("pending")]
    if require_reviewed and not pending.empty:
        raise IncompleteReviewError(
            f"{len(pending)} candidate L858R assays still require manual review in "
            f"{AUDIT_PATH.name}"
        )

    effective_include = audit.manual_decision.eq("include") | (
        audit.automated_candidate & audit.manual_decision.eq("pending")
    )
    included_assays = set(audit.loc[effective_include, "assay_chembl_id"])
    candidate_rows = data[data.assay_chembl_id.isin(included_assays)].copy()
    candidate_rows["exclusion_reason"] = candidate_rows.apply(
        measurement_exclusion_reason, axis=1
    )
    excluded = candidate_rows[candidate_rows.exclusion_reason.ne("")].copy()
    eligible = candidate_rows[candidate_rows.exclusion_reason.eq("")].copy()
    eligible["ic50_nm"] = pd.to_numeric(eligible.standard_value)

    dataset = _build_l858r_dataset(eligible)
    training = dataset[dataset.label.isin(["active", "inactive"])].copy()
    conflicts = dataset[dataset.label.eq("conflict")].copy()

    PROCESSED.mkdir(parents=True, exist_ok=True)
    training.to_csv(PROCESSED / "egfr_L858R_training.csv", index=False)
    conflicts.to_csv(PROCESSED / "egfr_L858R_conflicts.csv", index=False)
    excluded.to_csv(PROCESSED / "egfr_L858R_excluded.csv", index=False)

    label_counts = training.label.value_counts().to_dict()
    evidence_counts = training.evidence_grade.value_counts().sort_index().to_dict()
    summary = {
        "target": "EGFR L858R",
        "endpoint": "biochemical kinase inhibition IC50",
        "cutoff_nm": CUTOFF_NM,
        "source_snapshot": str(SNAPSHOT.relative_to(ROOT)).replace("\\", "/"),
        "review_complete": pending.empty,
        "status": "reviewed" if pending.empty else "provisional_pending_manual_assay_review",
        "assays_naming_L858R": int(len(audit)),
        "automated_candidate_assays": int(audit.automated_candidate.sum()),
        "pending_candidate_assays": int(len(pending)),
        "included_assays": int(len(included_assays)),
        "candidate_measurements": int(len(candidate_rows)),
        "eligible_measurements": int(len(eligible)),
        "excluded_candidate_measurements": int(len(excluded)),
        "labeled_compounds": int(len(training)),
        "active_compounds": int(label_counts.get("active", 0)),
        "inactive_compounds": int(label_counts.get("inactive", 0)),
        "conflicting_compounds": int(len(conflicts)),
        "evidence_grades": {key: int(value) for key, value in evidence_counts.items()},
    }
    (PROCESSED / "egfr_L858R_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8"
    )

    print(f"L858R assay audit: {AUDIT_PATH.relative_to(ROOT)}")
    print("Automated candidate assays:", summary["automated_candidate_assays"])
    print("Pending manual review:", summary["pending_candidate_assays"])
    print("Eligible measurements:", summary["eligible_measurements"])
    print("Labeled compounds:", summary["labeled_compounds"])
    print("Active:", summary["active_compounds"], "Inactive:", summary["inactive_compounds"])
    print("Conflicts:", summary["conflicting_compounds"])
    print("Status:", summary["status"])


def run_all_variants(data):
    """Preserve the original all-variant screen and output filenames."""
    descriptions = data.assay_description.fillna("")

    def mentions(pattern):
        return descriptions.str.contains(pattern, case=False, regex=True)

    found = pd.DataFrame({
        name: mentions(pattern) for name, pattern in LEGACY_VARIANTS.items()
    })
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
    candidate = signature.ne("") & ~wild_type & ~cell_based & enzyme & data.assay_type.eq("B")

    audit = data[["assay_chembl_id", "assay_description"]].copy()
    audit["variant"] = signature
    audit["candidate"] = candidate
    audit.drop_duplicates().to_csv(ROOT / "docs/research/mutant-assay-review.csv", index=False)

    rows = data.loc[candidate].copy()
    rows["variant"] = signature.loc[candidate]
    print("Candidate assay IDs:", rows.assay_chembl_id.nunique(), "rows:", len(rows))
    print(
        "Excluded, mutant named beside wild-type:",
        data.loc[signature.ne("") & wild_type, "assay_chembl_id"].nunique(),
        "assays",
    )

    values = pd.to_numeric(rows.standard_value, errors="coerce")
    keep = (
        rows.standard_relation.eq("=") & rows.standard_units.eq("nM")
        & values.gt(0) & np.isfinite(values) & rows.data_validity_comment.isna()
        & rows.potential_duplicate.eq(0) & rows.molecule_chembl_id.notna()
    )
    rows = rows.loc[keep].copy()
    rows["ic50_nm"] = values.loc[keep]
    valid = {
        smiles for smiles in rows.canonical_smiles.dropna().unique()
        if smiles.strip() and Chem.MolFromSmiles(smiles) is not None
    }
    rows = rows[rows.canonical_smiles.isin(valid)]

    pairs = rows.groupby(["molecule_chembl_id", "variant"])
    grouped = pairs.ic50_nm.agg(["min", "max", "count"])
    if pairs.canonical_smiles.nunique().gt(1).any():
        raise ValueError("Resolve inconsistent structures before exporting compound labels")
    dataset = grouped.rename(
        columns={"min": "min_ic50_nm", "max": "max_ic50_nm", "count": "measurements"}
    )
    dataset["canonical_smiles"] = pairs.canonical_smiles.first()
    dataset["activity_ids"] = pairs.activity_id.agg(lambda values: "|".join(values.astype(str)))
    dataset["assay_ids"] = pairs.assay_chembl_id.agg(_join_unique)
    dataset["label"] = [
        label_from_range(low, high)
        for low, high in zip(dataset.min_ic50_nm, dataset.max_ic50_nm)
    ]

    PROCESSED.mkdir(parents=True, exist_ok=True)
    labeled = dataset[dataset.label.isin(["active", "inactive"])]
    conflicts = dataset[dataset.label.eq("conflict")]
    labeled.to_csv(PROCESSED / "egfr_mutant_training.csv")
    conflicts.to_csv(PROCESSED / "egfr_mutant_conflicts.csv")

    print("\nScreened measurements:", len(rows))
    print(
        "Compound-variant pairs:", len(dataset), "labeled:", len(labeled),
        "conflicts:", len(conflicts),
    )
    print("Distinct compounds:", labeled.index.get_level_values(0).nunique())
    print("\nBy variant:")
    print(labeled.groupby(["variant", "label"]).size().unstack(fill_value=0).to_string())

    wt_path = PROCESSED / "egfr_training.csv"
    if wt_path.exists():
        wt = pd.read_csv(wt_path)
        shared = set(wt.molecule_chembl_id) & set(labeled.index.get_level_values(0))
        print("\nCompounds in both the wild-type and mutant sets:", len(shared))


def parse_args(args=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--variant", choices=["all", "L858R"], default="all",
        help="Use 'L858R' for the fixed single-variant target; default preserves the all-variant screen.",
    )
    parser.add_argument(
        "--require-reviewed", action="store_true",
        help="For L858R, stop unless every automated candidate assay has a manual decision.",
    )
    return parser.parse_args(args)


def main(args=None):
    options = parse_args(args)
    data = pd.read_csv(SNAPSHOT, low_memory=False)
    if options.variant == "L858R":
        try:
            run_l858r(data, require_reviewed=options.require_reviewed)
        except IncompleteReviewError as error:
            raise SystemExit(str(error)) from None
    else:
        if options.require_reviewed:
            raise ValueError("--require-reviewed is only valid with --variant L858R")
        run_all_variants(data)


if __name__ == "__main__":
    main()
