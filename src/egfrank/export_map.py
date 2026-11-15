import json
from pathlib import Path

import joblib
import numpy as np
from sklearn.manifold import TSNE

from train import load_training_data


def main():
    root = Path(__file__).resolve().parents[2]
    features, split = load_training_data()
    test = split["split"].eq("test")
    rows = split.loc[test]
    bits = features[test].astype(np.float64)
    model = joblib.load(root / "models/egfr_random_forest.joblib")
    scores = model.predict_proba(bits)[:, list(model.classes_).index("active")]

    overlap = bits @ bits.T
    total = bits.sum(axis=1)
    distance = 1 - overlap / (total[:, None] + total - overlap)
    coords = TSNE(metric="precomputed", init="random", random_state=42,
                  perplexity=15).fit_transform(distance.copy())
    coords = (coords - coords.min(axis=0)) / np.ptp(coords, axis=0)
    np.fill_diagonal(distance, np.inf)
    compounds = []
    for i, (compound_id, row) in enumerate(rows.iterrows()):
        compounds.append(dict(id=compound_id, score=float(scores[i]), label=row.label,
                              low=row.min_ic50_nm, high=row.max_ic50_nm,
                              x=float(coords[i, 0]), y=float(coords[i, 1]),
                              neighbors=np.argsort(distance[i], kind="stable")[:4].tolist()))
    target = root / "frontend/public/data/compound_map.json"
    target.write_text(json.dumps(compounds, allow_nan=False), encoding="utf-8")
    print(f"Exported {len(compounds)} held-out compounds")


if __name__ == "__main__":
    main()
