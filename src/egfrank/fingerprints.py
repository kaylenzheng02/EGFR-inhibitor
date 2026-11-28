import argparse
from pathlib import Path
from tempfile import NamedTemporaryFile

import numpy as np
import pandas as pd
from rdkit import Chem, rdBase
from rdkit.Chem import rdFingerprintGenerator


TARGET_FILES = {
    "general": ("egfr_training.csv", "egfr_fingerprints.npz"),
    "L858R": ("egfr_L858R_training.csv", "egfr_L858R_fingerprints.npz"),
}


def generate_fingerprints(data):
    """Return the project's chirality-aware Morgan fingerprint matrix."""
    generator = rdFingerprintGenerator.GetMorganGenerator(
        radius=2, fpSize=2048, includeChirality=True
    )
    fingerprints = np.empty((len(data), 2048), dtype=np.uint8)
    for index, row in enumerate(data.itertuples(index=False)):
        molecule = Chem.MolFromSmiles(row.canonical_smiles)
        if molecule is None:
            raise ValueError(f"Invalid SMILES for {row.molecule_chembl_id}")
        fingerprints[index] = generator.GetFingerprintAsNumPy(molecule)
    return fingerprints


def parse_args(args=None):
    parser = argparse.ArgumentParser(
        description="Generate molecular fingerprints for an EGFRank target dataset."
    )
    parser.add_argument(
        "--target", choices=sorted(TARGET_FILES), default="general",
        help="Dataset target; default preserves the original general-EGFR command.",
    )
    return parser.parse_args(args)


def main(args=None):
    options = parse_args(args)
    folder = Path(__file__).resolve().parents[2] / "data" / "processed"
    input_name, output_name = TARGET_FILES[options.target]
    data = pd.read_csv(folder / input_name)
    required = {"molecule_chembl_id", "canonical_smiles", "label"}
    missing = required.difference(data.columns)
    if missing:
        raise ValueError(f"Missing required columns: {sorted(missing)}")
    if data.molecule_chembl_id.duplicated().any():
        raise ValueError("The fingerprint input contains duplicate compound IDs")
    fingerprints = generate_fingerprints(data)

    output = folder / output_name
    with NamedTemporaryFile(dir=folder, delete=False) as temporary:
        path = Path(temporary.name)
    try:
        with path.open("wb") as stream:
            np.savez_compressed(
                stream, fingerprints=fingerprints,
                molecule_chembl_id=data.molecule_chembl_id.to_numpy(dtype=str),
                label=data.label.to_numpy(dtype=str), radius=2,
                fp_size=2048, include_chirality=True, rdkit_version=rdBase.rdkitVersion,
                target=options.target, source_file=input_name,
            )
        path.replace(output)
    finally:
        path.unlink(missing_ok=True)
    print(
        f"Saved {len(data)} {options.target} fingerprints of 2048 bits to {output.name}"
    )


if __name__ == "__main__":
    main()
