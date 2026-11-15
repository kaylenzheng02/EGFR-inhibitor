from pathlib import Path

import pandas as pd
from rdkit import Chem
from rdkit.Chem.Scaffolds import MurckoScaffold


ROOT = Path(__file__).resolve().parents[2]
DATA_DIR = ROOT / "data" / "processed"
INPUT_FILE = DATA_DIR / "egfr_training.csv"
OUTPUT_FILE = DATA_DIR / "egfr_scaffold_split.csv"

TEST_FRACTION = 0.20
SEED = 42


def get_scaffold(smiles):
    molecule = Chem.MolFromSmiles(smiles)
    if molecule is None:
        raise ValueError(f"RDKit could not read this SMILES: {smiles}")


    # molecules without rings have an empty scaffold
    # Keep them together as one group for splitting.
    scaffold = MurckoScaffold.MurckoScaffoldSmiles(mol=molecule)

    return scaffold or "ACYCLIC"  # Return "ACYCLIC" if the scaffold is empty


def main():
    data = pd.read_csv(INPUT_FILE)
    id_column = "molecule_chembl_id"
    required = [id_column, "canonical_smiles", "label"]

    missing = set(required).difference(data.columns)
    if missing:
        raise ValueError(f"Missing required columns: {sorted(missing)}")
    if data[required].isna().any().any():
        raise ValueError("A compound ID, structure, or label is missing.")
    if data[id_column].duplicated().any():
        raise ValueError("The training CSV contains duplicate compound IDs.")
    if set(data["label"]) != {"active", "inactive"}:
        raise ValueError("Expected both active and inactive labels.")

    data["scaffold"] = data["canonical_smiles"].map(get_scaffold)
    groups = pd.crosstab(data["scaffold"], data["label"])
    if len(groups) < 2:
        raise ValueError("At least two scaffold groups are needed.")

    totals = groups.sum()
    best_score = float("inf")
    best_test_groups = None

    # Try repeatable group orderings. Choose the split closest to
    # holding out 20% of each class without splitting scaffolds.
    for attempt in range(100):
        shuffled = groups.sample(frac=1, random_state=SEED + attempt)
        running_counts = shuffled.cumsum()

        # Both training and testing must contain both classes.
        valid = (running_counts > 0) & (running_counts < totals)
        candidates = running_counts.loc[valid.all(axis=1)]
        if candidates.empty:
            continue

        fractions = candidates.div(totals)
        scores = (fractions - TEST_FRACTION).abs().mean(axis=1)
        last_group = scores.idxmin()
        score = scores.loc[last_group]

        if score < best_score:
            best_score = score
            last_position = shuffled.index.get_loc(last_group)
            best_test_groups = shuffled.index[:last_position + 1]

    if best_test_groups is None:
        raise ValueError("Could not find a scaffold split with both classes.")

    data["split"] = "train"
    data.loc[data["scaffold"].isin(best_test_groups), "split"] = "test"

    if data.groupby("scaffold")["split"].nunique().gt(1).any():
        raise ValueError("A scaffold appears in both training and test sets.")

    data.to_csv(OUTPUT_FILE, index=False)
    print(pd.crosstab(data["split"], data["label"]))
    print(f"Saved scaffold split to {OUTPUT_FILE}")


if __name__ == "__main__":
    main()
