"""Create reproducible scaffold-separated EGFRank train/test splits."""

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem.Scaffolds import MurckoScaffold

ROOT = Path(__file__).resolve().parents[2]
DATA_DIR = ROOT / "data/processed"
TEST_FRACTION = 0.20
SEED = 42
SEARCH_ATTEMPTS = 25_000

GENERAL_INPUT = DATA_DIR / "egfr_training.csv"
GENERAL_OUTPUT = DATA_DIR / "egfr_scaffold_split.csv"
L858R_INPUT = DATA_DIR / "egfr_L858R_scaffolds.csv"
L858R_OUTPUT = DATA_DIR / "egfr_L858R_scaffold_split.csv"
L858R_SUMMARY = DATA_DIR / "egfr_L858R_split_summary.json"
L858R_REPORT = ROOT / "docs/research/l858r-split-validation.md"


def get_scaffold(smiles):
    molecule = Chem.MolFromSmiles(smiles)
    if molecule is None:
        raise ValueError(f"RDKit could not read this SMILES: {smiles}")
    scaffold = MurckoScaffold.MurckoScaffoldSmiles(mol=molecule)
    return scaffold or "ACYCLIC"


def validate_base_table(data):
    required = {"molecule_chembl_id", "canonical_smiles", "label"}
    missing = required.difference(data.columns)
    if missing:
        raise ValueError(f"Missing required columns: {sorted(missing)}")
    if data[list(required)].isna().any().any():
        raise ValueError("A compound ID, structure, or label is missing")
    if data.molecule_chembl_id.duplicated().any():
        raise ValueError("The input contains duplicate compound IDs")
    if set(data.label) != {"active", "inactive"}:
        raise ValueError("Expected both and only active/inactive labels")


def choose_test_groups(data, group_column, test_fraction=TEST_FRACTION,
                       seed=SEED, attempts=SEARCH_ATTEMPTS):
    """Choose whole groups near the requested fraction of both classes."""
    groups = pd.crosstab(data[group_column], data.label).reindex(
        columns=["active", "inactive"], fill_value=0
    ).sort_index()
    if len(groups) < 2:
        raise ValueError("At least two split groups are required")
    totals = groups.sum().to_numpy(dtype=int)
    desired = np.rint(totals * test_fraction).astype(int)
    desired_groups = round(len(groups) * test_fraction)
    counts = groups.to_numpy(dtype=int)
    names = groups.index.to_numpy(dtype=str)

    best_key = None
    best_test_groups = None
    for attempt in range(attempts):
        rng = np.random.default_rng(seed + attempt)
        selected = rng.random(len(groups)) < test_fraction
        if not selected.any() or selected.all():
            continue
        test_counts = counts[selected].sum(axis=0)
        train_counts = totals - test_counts
        if (test_counts == 0).any() or (train_counts == 0).any():
            continue
        chosen_names = tuple(sorted(names[selected]))
        key = (
            int(np.abs(test_counts - desired).sum()),
            int(abs(test_counts.sum() - desired.sum())),
            int(abs(selected.sum() - desired_groups)),
            chosen_names,
        )
        if best_key is None or key < best_key:
            best_key = key
            best_test_groups = chosen_names

    if best_test_groups is None:
        raise ValueError("Could not find a grouped split containing both classes")
    return set(best_test_groups)


def assign_split(data, group_column, test_groups):
    output = data.copy()
    output["split"] = np.where(
        output[group_column].astype(str).isin(test_groups), "test", "train"
    )
    if output.groupby(group_column).split.nunique().gt(1).any():
        raise ValueError(f"A {group_column} value crosses the split boundary")
    if set(output.split) != {"train", "test"}:
        raise ValueError("Expected both train and test rows")
    if (pd.crosstab(output.split, output.label) == 0).any().any():
        raise ValueError("Both classes must appear in training and testing")
    return output


def connected_split_groups(data):
    """Connect compounds sharing either a scaffold or a fingerprint row."""
    required = {"molecule_chembl_id", "scaffold", "fingerprint_group"}
    missing = required.difference(data.columns)
    if missing:
        raise ValueError(f"Missing grouping columns: {sorted(missing)}")

    parent = list(range(len(data)))

    def find(index):
        while parent[index] != index:
            parent[index] = parent[parent[index]]
            index = parent[index]
        return index

    def union(left, right):
        left_root = find(left)
        right_root = find(right)
        if left_root != right_root:
            parent[right_root] = left_root

    seen = {}
    for position, row in enumerate(data.itertuples(index=False)):
        keys = (("scaffold", row.scaffold), ("fingerprint", row.fingerprint_group))
        for key in keys:
            if key in seen:
                union(position, seen[key])
            else:
                seen[key] = position

    components = {}
    for position in range(len(data)):
        components.setdefault(find(position), []).append(position)

    group_ids = np.empty(len(data), dtype=object)
    for positions in components.values():
        compound_ids = sorted(data.iloc[positions].molecule_chembl_id.astype(str))
        digest = hashlib.sha256("|".join(compound_ids).encode("utf-8")).hexdigest()
        group_id = digest[:16]
        for position in positions:
            group_ids[position] = group_id
    return pd.Series(group_ids, index=data.index, name="split_group")


def run_general():
    """Preserve the original general-EGFR scaffold split behavior."""
    data = pd.read_csv(GENERAL_INPUT)
    validate_base_table(data)
    data["scaffold"] = data.canonical_smiles.map(get_scaffold)

    groups = pd.crosstab(data.scaffold, data.label)
    totals = groups.sum()
    best_score = float("inf")
    best_test_groups = None
    for attempt in range(100):
        shuffled = groups.sample(frac=1, random_state=SEED + attempt)
        running_counts = shuffled.cumsum()
        valid = (running_counts > 0) & (running_counts < totals)
        candidates = running_counts.loc[valid.all(axis=1)]
        if candidates.empty:
            continue
        fractions = candidates.div(totals)
        scores = (fractions - TEST_FRACTION).abs().mean(axis=1)
        last_group = scores.idxmin()
        if scores.loc[last_group] < best_score:
            best_score = scores.loc[last_group]
            last_position = shuffled.index.get_loc(last_group)
            best_test_groups = set(shuffled.index[:last_position + 1])
    if best_test_groups is None:
        raise ValueError("Could not find a scaffold split with both classes")

    output = assign_split(data, "scaffold", best_test_groups)
    output.to_csv(GENERAL_OUTPUT, index=False)
    print(pd.crosstab(output.split, output.label))
    print(f"Saved scaffold split to {GENERAL_OUTPUT.name}")


def sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def validate_l858r_split(data):
    checks = {
        "compound_ids_crossing_split": int(
            data.groupby("molecule_chembl_id").split.nunique().gt(1).sum()
        ),
        "canonical_structures_crossing_split": int(
            data.groupby("rdkit_canonical_smiles").split.nunique().gt(1).sum()
        ),
        "scaffolds_crossing_split": int(
            data.groupby("scaffold").split.nunique().gt(1).sum()
        ),
        "fingerprint_groups_crossing_split": int(
            data.groupby("fingerprint_group").split.nunique().gt(1).sum()
        ),
        "connected_groups_crossing_split": int(
            data.groupby("split_group").split.nunique().gt(1).sum()
        ),
    }
    if any(checks.values()):
        raise ValueError(f"Leakage checks failed: {checks}")
    return checks


def summarize_partition(data, split_name):
    rows = data[data.split.eq(split_name)]
    return {
        "compounds": int(len(rows)),
        "active": int(rows.label.eq("active").sum()),
        "inactive": int(rows.label.eq("inactive").sum()),
        "active_fraction": float(rows.label.eq("active").mean()),
        "connected_groups": int(rows.split_group.nunique()),
        "scaffolds": int(rows.scaffold.nunique()),
        "fingerprint_groups": int(rows.fingerprint_group.nunique()),
        "evidence_grade_A": int(rows.evidence_grade.eq("A").sum()),
        "evidence_grade_B": int(rows.evidence_grade.eq("B").sum()),
        "evidence_grade_C": int(rows.evidence_grade.eq("C").sum()),
    }


def render_l858r_report(summary):
    train = summary["partitions"]["train"]
    test = summary["partitions"]["test"]
    checks = summary["leakage_checks"]
    return f"""# EGFR L858R split validation

This split freezes the held-out set for the first L858R benchmark. It assigns whole connected groups rather than individual compounds. A connected group combines compounds that share a Murcko scaffold or an identical 2,048-bit Morgan fingerprint, including transitive connections between those relationships.

## Split settings

| Setting | Value |
| --- | ---: |
| Test target | {summary['settings']['test_fraction']:.0%} |
| Random seed | {summary['settings']['seed']} |
| Candidate assignments searched | {summary['settings']['search_attempts']:,} |
| Connected groups | {summary['connected_groups']} |

## Final assignment

| Set | Active | Inactive | Total | Scaffolds | Connected groups |
| --- | ---: | ---: | ---: | ---: | ---: |
| Training | {train['active']} | {train['inactive']} | {train['compounds']} | {train['scaffolds']} | {train['connected_groups']} |
| Test | {test['active']} | {test['inactive']} | {test['compounds']} | {test['scaffolds']} | {test['connected_groups']} |

The test partition contains {summary['test_fractions']['overall']:.1%} of all compounds, {summary['test_fractions']['active']:.1%} of active compounds, and {summary['test_fractions']['inactive']:.1%} of inactive compounds. Its always-active accuracy is {test['active_fraction']:.1%}.

The held-out test set is now frozen. Model selection, feature decisions, hyperparameter choices, and score-threshold decisions must use only the training partition and training-only grouped validation.

## Leakage checks

| Check | Crossing groups |
| --- | ---: |
| Compound IDs | {checks['compound_ids_crossing_split']} |
| Canonical structures | {checks['canonical_structures_crossing_split']} |
| Murcko scaffolds | {checks['scaffolds_crossing_split']} |
| Identical fingerprints | {checks['fingerprint_groups_crossing_split']} |
| Connected split groups | {checks['connected_groups_crossing_split']} |

All leakage checks pass at zero. Both activity classes occur in both partitions. The split tests chemical frameworks absent from training while also preventing identical model inputs from appearing on both sides.

## Interpretation

The test set contains {test['inactive']} inactive compounds. Consequently, each inactive error changes inactive recall by {100 / test['inactive']:.1f} percentage points. The final benchmark must report raw confusion-matrix counts alongside summary metrics and should not use this test set for iterative model development.
"""


def run_l858r():
    data = pd.read_csv(L858R_INPUT)
    validate_base_table(data)
    required = {
        "rdkit_canonical_smiles", "scaffold", "fingerprint_group", "evidence_grade"
    }
    missing = required.difference(data.columns)
    if missing:
        raise ValueError(
            f"Missing validated L858R columns: {sorted(missing)}; run validate_l858r_dataset.py"
        )
    data["split_group"] = connected_split_groups(data)
    test_groups = choose_test_groups(data, "split_group")
    output = assign_split(data, "split_group", test_groups)
    checks = validate_l858r_split(output)
    output.to_csv(L858R_OUTPUT, index=False)

    total_active = int(output.label.eq("active").sum())
    total_inactive = int(output.label.eq("inactive").sum())
    test_rows = output[output.split.eq("test")]

    summary = {
        "status": "passed",
        "target": "EGFR L858R",
        "source_file": str(L858R_INPUT.relative_to(ROOT)).replace("\\", "/"),
        "source_sha256": sha256(L858R_INPUT),
        "output_file": str(L858R_OUTPUT.relative_to(ROOT)).replace("\\", "/"),
        "output_sha256": sha256(L858R_OUTPUT),
        "settings": {
            "test_fraction": TEST_FRACTION,
            "seed": SEED,
            "search_attempts": SEARCH_ATTEMPTS,
            "grouping": "connected components of Murcko scaffold and identical fingerprint",
        },
        "connected_groups": int(output.split_group.nunique()),
        "test_fractions": {
            "overall": float(len(test_rows) / len(output)),
            "active": float(test_rows.label.eq("active").sum() / total_active),
            "inactive": float(test_rows.label.eq("inactive").sum() / total_inactive),
        },
        "partitions": {
            "train": summarize_partition(output, "train"),
            "test": summarize_partition(output, "test"),
        },
        "leakage_checks": checks,
        "test_split_groups": sorted(test_groups),
    }
    L858R_SUMMARY.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    L858R_REPORT.write_text(render_l858r_report(summary), encoding="utf-8")

    print(pd.crosstab(output.split, output.label))
    print(f"Connected groups: {summary['connected_groups']}")
    print(f"Saved {L858R_OUTPUT.name}, {L858R_SUMMARY.name}, and {L858R_REPORT.name}")


def parse_args(args=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--target", choices=["general", "L858R"], default="general",
        help="Dataset target; default preserves the original general-EGFR split.",
    )
    return parser.parse_args(args)


def main(args=None):
    options = parse_args(args)
    if options.target == "L858R":
        run_l858r()
    else:
        run_general()


if __name__ == "__main__":
    main()
