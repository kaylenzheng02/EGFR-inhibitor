from pathlib import Path

import joblib
import numpy as np
from sklearn.metrics import (
    average_precision_score,
    classification_report,
    confusion_matrix,
    roc_auc_score,
)

from train import load_training_data


MODEL_PATH = Path(__file__).resolve().parents[2] / "models" / "egfr_random_forest.joblib"


def report(name, actual, predicted, scores):
    print(f"\n{name}")
    print(classification_report(
        actual, predicted, labels=["inactive", "active"], zero_division=0
    ))

    is_active = actual == "active"
    print(f"ROC-AUC: {roc_auc_score(is_active, scores):.3f}")
    print(f"Average precision (PR summary): {average_precision_score(is_active, scores):.3f}")

    matrix = confusion_matrix(actual, predicted, labels=["inactive", "active"])
    print("Confusion matrix (rows = actual, columns = predicted)")
    print("Label order: inactive, active")
    print(matrix)


def main():
    fingerprints, split = load_training_data()
    test_rows = split["split"].eq("test").to_numpy()
    X_test = fingerprints[test_rows]
    actual = split.loc[test_rows, "label"].to_numpy()

    model = joblib.load(MODEL_PATH)
    predicted = model.predict(X_test)
    active_column = list(model.classes_).index("active")
    scores = model.predict_proba(X_test)[:, active_column]

    print(f"Evaluating {len(actual)} held-out compounds.")

    # every compound receives the same score and an "active" prediction.
    report("Always-active baseline", actual, np.repeat("active", len(actual)),
           np.ones(len(actual)))
    report("Random Forest", actual, predicted, scores)


if __name__ == "__main__":
    main()