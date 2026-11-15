from pathlib import Path

import numpy as np
import pandas as pd
import joblib
from sklearn.ensemble import RandomForestClassifier


DATA_DIR = Path(__file__).resolve().parents[2] / "data" / "processed"
ID = "molecule_chembl_id"


def load_training_data():
    split = pd.read_csv(DATA_DIR / "egfr_scaffold_split.csv")
    split = split.set_index(ID)
    if not split.index.is_unique:
        raise ValueError("duplicate compound IDs in the scaffold split file")

    with np.load(DATA_DIR / "egfr_fingerprints.npz") as saved:
        fingerprints = saved["fingerprints"]
        ids = pd.Index(saved[ID], name=ID)
        labels = saved["label"]

    if not ids.is_unique:
        raise ValueError("duplicate compound IDs in the fingerprint file")
    if set(ids) != set(split.index):
        raise ValueError("the two files contain different compounds")
    if fingerprints.shape != (len(ids), 2048):
        raise ValueError("expected one 2048 bit fingerprint per compound")

    # put the labels and split assignments in the same order as the fingerprints
    split = split.loc[ids]

    if not np.array_equal(split["label"].to_numpy(), labels):
        raise ValueError("Activity labels disagree between the two files")

    return fingerprints, split


def main():
    fingerprints, split = load_training_data()
    print(f"Matched {len(split)} compounds with {fingerprints.shape[1]} features each")
    print(pd.crosstab(split["label"], split["split"]))

    train_rows = split["split"].eq("train").to_numpy()
    X_train = fingerprints[train_rows]
    y_train = split.loc[train_rows, "label"]

    model = RandomForestClassifier(n_estimators=100, random_state=42)
    model.fit(X_train, y_train)

    model_dir = Path(__file__).resolve().parents[2] / "models"
    model_dir.mkdir(exist_ok=True)
    model_path = model_dir / "egfr_random_forest.joblib"
    joblib.dump(model, model_path)

    print(f"Trained on {len(y_train)} compounds.")
    print(f"Saved model to {model_path}")


if __name__ == "__main__":
    main()
