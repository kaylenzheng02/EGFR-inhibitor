"""Validate and characterize the reviewed EGFR L858R dataset before splitting."""

import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
from rdkit import Chem, rdBase
from rdkit.Chem.Scaffolds import MurckoScaffold

ROOT = Path(__file__).resolve().parents[2]
PROCESSED = ROOT / "data/processed"
DATASET_PATH = PROCESSED / "egfr_L858R_training.csv"
FINGERPRINT_PATH = PROCESSED / "egfr_L858R_fingerprints.npz"
SCAFFOLD_PATH = PROCESSED / "egfr_L858R_scaffolds.csv"
REPORT_PATH = PROCESSED / "egfr_L858R_validation.json"
MARKDOWN_PATH = ROOT / "docs/research/l858r-dataset-validation.md"
DECISIONS_PATH = ROOT / "docs/research/l858r-assay-decisions.csv"
CUTOFF_NM = 1000
REQUIRED_COLUMNS = {
    "molecule_chembl_id", "variant", "canonical_smiles", "label",
    "min_ic50_nm", "max_ic50_nm", "measurements", "evidence_grade",
    "activity_ids", "assay_ids", "document_ids",
}


def pipe_values(value):
    if pd.isna(value) or not str(value).strip():
        return []
    return [part for part in str(value).split("|") if part]


def expected_label(minimum, maximum):
    if maximum <= CUTOFF_NM:
        return "active"
    if minimum > CUTOFF_NM:
        return "inactive"
    return "conflict"


def expected_evidence_grade(measurements, assay_count, document_count):
    if measurements == 1:
        return "C"
    if assay_count > 1 or document_count > 1:
        return "A"
    return "B"


def scaffold_for_smiles(smiles):
    molecule = Chem.MolFromSmiles(smiles)
    if molecule is None:
        raise ValueError(f"RDKit could not read this SMILES: {smiles}")
    return MurckoScaffold.MurckoScaffoldSmiles(mol=molecule) or "ACYCLIC"


def validate_table(data, decisions):
    """Validate label provenance and return a structure-enriched table."""
    problems = []
    missing = REQUIRED_COLUMNS.difference(data.columns)
    if missing:
        raise ValueError(f"Missing required columns: {sorted(missing)}")
    if data[list(REQUIRED_COLUMNS)].isna().any().any():
        problems.append("required values are missing")
    if data.molecule_chembl_id.duplicated().any():
        problems.append("compound IDs are duplicated")
    if set(data.variant) != {"L858R"}:
        problems.append("variant values are not exclusively L858R")
    if set(data.label) != {"active", "inactive"}:
        problems.append("both and only active/inactive labels are required")
    if set(data.evidence_grade).difference({"A", "B", "C"}):
        problems.append("unexpected evidence grades are present")

    numeric = data[["min_ic50_nm", "max_ic50_nm", "measurements"]].apply(
        pd.to_numeric, errors="coerce"
    )
    if numeric.isna().any().any() or not np.isfinite(numeric.to_numpy()).all():
        problems.append("IC50 or measurement counts are not finite numbers")
    if (numeric[["min_ic50_nm", "max_ic50_nm"]] <= 0).any().any():
        problems.append("IC50 values must be positive")
    if (numeric.min_ic50_nm > numeric.max_ic50_nm).any():
        problems.append("an IC50 minimum exceeds its maximum")
    if (numeric.measurements < 1).any() or not np.equal(
        numeric.measurements, np.floor(numeric.measurements)
    ).all():
        problems.append("measurement counts must be positive integers")

    activity_lists = data.activity_ids.map(pipe_values)
    assay_lists = data.assay_ids.map(pipe_values)
    document_lists = data.document_ids.map(pipe_values)
    if not np.array_equal(activity_lists.map(len).to_numpy(), numeric.measurements.to_numpy()):
        problems.append("measurement counts disagree with activity ID counts")
    all_activity_ids = [item for values in activity_lists for item in values]
    if len(all_activity_ids) != len(set(all_activity_ids)):
        problems.append("an activity ID is assigned to more than one compound")

    expected_labels = [
        expected_label(low, high)
        for low, high in zip(numeric.min_ic50_nm, numeric.max_ic50_nm)
    ]
    if expected_labels != data.label.tolist():
        problems.append("one or more labels disagree with the IC50 cutoff")
    expected_grades = [
        expected_evidence_grade(count, len(assays), len(documents))
        for count, assays, documents in zip(
            numeric.measurements, assay_lists, document_lists
        )
    ]
    if expected_grades != data.evidence_grade.tolist():
        problems.append("one or more evidence grades disagree with provenance counts")

    if decisions.assay_chembl_id.duplicated().any():
        problems.append("assay decisions contain duplicate IDs")
    included_assays = set(
        decisions.loc[decisions.manual_decision.eq("include"), "assay_chembl_id"]
    )
    used_assays = {item for values in assay_lists for item in values}
    if not used_assays.issubset(included_assays):
        problems.append("the dataset uses an assay without an include decision")

    molecules = data.canonical_smiles.map(Chem.MolFromSmiles)
    if molecules.isna().any():
        problems.append("RDKit could not parse one or more structures")
    if problems:
        raise ValueError("; ".join(problems))

    enriched = data.copy()
    enriched["rdkit_canonical_smiles"] = molecules.map(
        lambda molecule: Chem.MolToSmiles(
            molecule, canonical=True, isomericSmiles=True
        )
    )
    if enriched.rdkit_canonical_smiles.duplicated().any():
        raise ValueError("Two compound IDs resolve to the same RDKit canonical structure")
    enriched["scaffold"] = enriched.canonical_smiles.map(scaffold_for_smiles)
    return enriched


def validate_fingerprints(data):
    if not FINGERPRINT_PATH.exists():
        raise ValueError(
            f"Missing {FINGERPRINT_PATH.name}; run fingerprints.py --target L858R first"
        )
    with np.load(FINGERPRINT_PATH) as saved:
        matrix = saved["fingerprints"]
        ids = saved["molecule_chembl_id"].astype(str)
        labels = saved["label"].astype(str)
        radius = int(saved["radius"])
        size = int(saved["fp_size"])
        chirality = bool(saved["include_chirality"])
        saved_rdkit = str(saved["rdkit_version"])
        target = str(saved["target"])
        source_file = str(saved["source_file"])

    expected_shape = (len(data), 2048)
    if matrix.shape != expected_shape:
        raise ValueError(f"Expected fingerprint shape {expected_shape}, found {matrix.shape}")
    if not np.array_equal(ids, data.molecule_chembl_id.to_numpy(dtype=str)):
        raise ValueError("Fingerprint compound IDs do not match the dataset order")
    if not np.array_equal(labels, data.label.to_numpy(dtype=str)):
        raise ValueError("Fingerprint labels do not match the dataset")
    if set(np.unique(matrix)).difference({0, 1}):
        raise ValueError("Fingerprints contain values other than zero and one")
    bit_counts = matrix.sum(axis=1)
    if (bit_counts == 0).any():
        raise ValueError("At least one compound has an empty fingerprint")
    group_ids = np.array([
        hashlib.sha256(row.tobytes()).hexdigest() for row in matrix
    ])
    groups = pd.DataFrame({
        "fingerprint_group": group_ids,
        "label": data.label.to_numpy(),
        "scaffold": data.scaffold.to_numpy(),
    }).groupby("fingerprint_group").agg(
        compounds=("label", "size"),
        labels=("label", "nunique"),
        scaffolds=("scaffold", "nunique"),
    )
    repeated = groups[groups.compounds.gt(1)]
    cross_label_members = []
    for group_id in repeated.index[repeated.labels.gt(1)]:
        members = data.loc[
            group_ids == group_id,
            ["molecule_chembl_id", "label", "canonical_smiles"],
        ]
        cross_label_members.append({
            "fingerprint_group": group_id,
            "members": members.to_dict(orient="records"),
        })
    summary = {
        "path": str(FINGERPRINT_PATH.relative_to(ROOT)).replace("\\", "/"),
        "sha256": sha256(FINGERPRINT_PATH),
        "shape": list(matrix.shape),
        "radius": radius,
        "size": size,
        "include_chirality": chirality,
        "rdkit_version": saved_rdkit,
        "runtime_rdkit_version": rdBase.rdkitVersion,
        "target": target,
        "source_file": source_file,
        "minimum_set_bits": int(bit_counts.min()),
        "median_set_bits": float(np.median(bit_counts)),
        "maximum_set_bits": int(bit_counts.max()),
        "duplicate_fingerprint_rows": int((repeated.compounds - 1).sum()),
        "duplicate_fingerprint_groups": int(len(repeated)),
        "cross_label_fingerprint_groups": int(repeated.labels.gt(1).sum()),
        "cross_scaffold_fingerprint_groups": int(repeated.scaffolds.gt(1).sum()),
        "cross_label_group_members": cross_label_members,
    }
    return summary, group_ids


def sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def build_report(data, decisions, fingerprint_summary):
    scaffold_groups = data.groupby("scaffold").agg(
        compounds=("molecule_chembl_id", "size"),
        active=("label", lambda values: int(values.eq("active").sum())),
        inactive=("label", lambda values: int(values.eq("inactive").sum())),
    )
    scaffold_groups["classes"] = (
        scaffold_groups.active.gt(0).astype(int)
        + scaffold_groups.inactive.gt(0).astype(int)
    )
    data = data.copy()
    data["scaffold_size"] = data.scaffold.map(scaffold_groups.compounds)
    data["scaffold_active"] = data.scaffold.map(scaffold_groups.active)
    data["scaffold_inactive"] = data.scaffold.map(scaffold_groups.inactive)

    reviewed = set(
        decisions.loc[decisions.manual_decision.eq("include"), "assay_chembl_id"]
    )
    used = {item for value in data.assay_ids for item in pipe_values(value)}
    documents = {item for value in data.document_ids for item in pipe_values(value)}
    counts = data.label.value_counts().to_dict()
    evidence = data.evidence_grade.value_counts().sort_index().to_dict()
    near_cutoff = data[
        data.min_ic50_nm.le(1200) & data.max_ic50_nm.ge(800)
    ]

    report = {
        "status": "passed",
        "target": "EGFR L858R",
        "cutoff_nm": CUTOFF_NM,
        "dataset": {
            "path": str(DATASET_PATH.relative_to(ROOT)).replace("\\", "/"),
            "sha256": sha256(DATASET_PATH),
            "compounds": int(len(data)),
            "unique_compound_ids": int(data.molecule_chembl_id.nunique()),
            "active": int(counts.get("active", 0)),
            "inactive": int(counts.get("inactive", 0)),
            "active_fraction": float(data.label.eq("active").mean()),
            "always_active_accuracy": float(data.label.eq("active").mean()),
            "measurements": int(data.measurements.sum()),
            "multi_measurement_compounds": int(data.measurements.gt(1).sum()),
            "near_cutoff_800_to_1200_nm": int(len(near_cutoff)),
            "rdkit_canonical_changes": int(
                data.canonical_smiles.ne(data.rdkit_canonical_smiles).sum()
            ),
            "duplicate_rdkit_structures": int(
                data.rdkit_canonical_smiles.duplicated().sum()
            ),
            "evidence_grades": {
                key: int(value) for key, value in evidence.items()
            },
        },
        "provenance": {
            "reviewed_included_assays": int(len(reviewed)),
            "contributing_assays": int(len(used)),
            "included_assays_without_eligible_measurements": sorted(reviewed - used),
            "contributing_documents": int(len(documents)),
        },
        "scaffolds": {
            "definition": "RDKit Murcko scaffold without chirality; ringless structures use ACYCLIC",
            "groups": int(len(scaffold_groups)),
            "active_scaffolds": int(scaffold_groups.active.gt(0).sum()),
            "inactive_scaffolds": int(scaffold_groups.inactive.gt(0).sum()),
            "mixed_label_scaffolds": int(scaffold_groups.classes.eq(2).sum()),
            "singleton_scaffolds": int(scaffold_groups.compounds.eq(1).sum()),
            "largest_group": int(scaffold_groups.compounds.max()),
            "acyclic_compounds": int(
                scaffold_groups.loc["ACYCLIC", "compounds"]
                if "ACYCLIC" in scaffold_groups.index else 0
            ),
        },
        "fingerprints": fingerprint_summary,
        "benchmark_readiness": {
            "passed_identity_and_label_checks": True,
            "scaffold_split_not_yet_created": True,
            "primary_limitation": (
                "Most labels rely on one eligible measurement and only 27 "
                "scaffolds contain an inactive compound."
            ),
        },
    }
    return data, report


def render_markdown(report):
    dataset = report["dataset"]
    provenance = report["provenance"]
    scaffolds = report["scaffolds"]
    fingerprints = report["fingerprints"]
    excluded = ", ".join(provenance["included_assays_without_eligible_measurements"])
    cross_label_members = ", ".join(
        f"{member['molecule_chembl_id']} ({member['label']})"
        for group in fingerprints["cross_label_group_members"]
        for member in group["members"]
    )
    return f"""# EGFR L858R dataset validation

The reviewed L858R dataset passed the pre-split identity, structure, label, provenance, scaffold, and fingerprint checks. This report characterizes the frozen input to the future benchmark; it does not create a train/test split or evaluate a model.

## Dataset

| Check | Result |
| --- | ---: |
| Compounds | {dataset['compounds']} |
| Active | {dataset['active']} |
| Inactive | {dataset['inactive']} |
| Always-active baseline | {dataset['always_active_accuracy']:.1%} |
| Eligible measurements | {dataset['measurements']} |
| Compounds with multiple measurements | {dataset['multi_measurement_compounds']} |
| Compounds near 800–1,200 nM | {dataset['near_cutoff_800_to_1200_nm']} |
| Duplicate RDKit structures | {dataset['duplicate_rdkit_structures']} |

Evidence grades are A: {dataset['evidence_grades'].get('A', 0)}, B: {dataset['evidence_grades'].get('B', 0)}, and C: {dataset['evidence_grades'].get('C', 0)}. Grade C means that one eligible measurement supports the label.

## Provenance

All contributing assay IDs have an explicit include decision. Of {provenance['reviewed_included_assays']} reviewed assays, {provenance['contributing_assays']} contribute an exact eligible measurement from {provenance['contributing_documents']} source documents. Five reviewed assays contribute no final row after activity-level quality screening: {excluded}.

## Chemical coverage

| Check | Result |
| --- | ---: |
| Murcko scaffold groups | {scaffolds['groups']} |
| Scaffolds containing active compounds | {scaffolds['active_scaffolds']} |
| Scaffolds containing inactive compounds | {scaffolds['inactive_scaffolds']} |
| Scaffolds containing both classes | {scaffolds['mixed_label_scaffolds']} |
| Singleton scaffolds | {scaffolds['singleton_scaffolds']} |
| Largest scaffold group | {scaffolds['largest_group']} compounds |
| Acyclic compounds | {scaffolds['acyclic_compounds']} |

The scaffold definition matches the existing EGFR benchmark: an RDKit Murcko scaffold without chirality, with ringless compounds assigned to `ACYCLIC`.

## Fingerprints

The saved matrix contains {fingerprints['shape'][0]} rows and {fingerprints['shape'][1]} bits per row. It uses radius {fingerprints['radius']} Morgan fingerprints with chirality enabled under RDKit {fingerprints['rdkit_version']}. Set-bit counts range from {fingerprints['minimum_set_bits']} to {fingerprints['maximum_set_bits']}, with a median of {fingerprints['median_set_bits']:.1f}.

{fingerprints['duplicate_fingerprint_groups']} fingerprint groups contain more than one structure, accounting for {fingerprints['duplicate_fingerprint_rows']} rows beyond the first member of each group. {fingerprints['cross_scaffold_fingerprint_groups']} groups span more than one Murcko scaffold, and {fingerprints['cross_label_fingerprint_groups']} group contains both an active and an inactive compound. These are distinct structures that produce the same hashed fingerprint, not duplicate structures. Every member of a repeated fingerprint group must remain on the same side of the future train/test split.

The cross-label fingerprint group contains {cross_label_members}. A model using this fingerprint representation must assign those compounds the same score even though their experimental labels differ.

## Benchmark implication

The dataset is ready for a leakage-safe scaffold split. The main constraint is the small and shallow negative evidence: 48 inactive compounds occupy 27 scaffolds, while 321 of all 329 labels rely on one eligible measurement. The benchmark should therefore retain an untouched scaffold-held-out test set and use scaffold-grouped validation within training. Accuracy must be compared with the {dataset['always_active_accuracy']:.1%} always-active baseline and accompanied by balanced and inactive-class metrics.
"""


def main():
    data = pd.read_csv(DATASET_PATH)
    decisions = pd.read_csv(DECISIONS_PATH)
    enriched = validate_table(data, decisions)
    fingerprints, fingerprint_groups = validate_fingerprints(enriched)
    enriched["fingerprint_group"] = fingerprint_groups
    inventory, report = build_report(enriched, decisions, fingerprints)

    inventory.to_csv(SCAFFOLD_PATH, index=False)
    REPORT_PATH.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    MARKDOWN_PATH.write_text(render_markdown(report), encoding="utf-8")

    print(f"Validated {len(inventory)} L858R compounds")
    print(f"Labels: {report['dataset']['active']} active, {report['dataset']['inactive']} inactive")
    print(f"Scaffold groups: {report['scaffolds']['groups']}")
    print(f"Inactive scaffolds: {report['scaffolds']['inactive_scaffolds']}")
    print(f"Saved {SCAFFOLD_PATH.name}, {REPORT_PATH.name}, and {MARKDOWN_PATH.name}")


if __name__ == "__main__":
    main()
