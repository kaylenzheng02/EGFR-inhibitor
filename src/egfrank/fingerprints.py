from pathlib import Path
from tempfile import NamedTemporaryFile

import numpy as np
import pandas as pd
from rdkit import Chem, rdBase
from rdkit.Chem import rdFingerprintGenerator


def main():
    folder = Path(__file__).resolve().parents[2] / "data" / "processed"
    data = pd.read_csv(folder / "egfr_training.csv")
    generator = rdFingerprintGenerator.GetMorganGenerator(
        radius=2, fpSize=2048, includeChirality=True
    )
    fingerprints = np.empty((len(data), 2048), dtype=np.uint8)
    for index, row in enumerate(data.itertuples(index=False)):
        molecule = Chem.MolFromSmiles(row.canonical_smiles)
        if molecule is None:
            raise ValueError(f"Invalid SMILES for {row.molecule_chembl_id}")
        fingerprints[index] = generator.GetFingerprintAsNumPy(molecule)

    output = folder / "egfr_fingerprints.npz"
    with NamedTemporaryFile(dir=folder, delete=False) as temporary:
        path = Path(temporary.name)
    try:
        with path.open("wb") as stream:
            np.savez_compressed(
                stream, fingerprints=fingerprints,
                molecule_chembl_id=data.molecule_chembl_id.to_numpy(dtype=str),
                label=data.label.to_numpy(dtype=str), radius=2,
                fp_size=2048, include_chirality=True, rdkit_version=rdBase.rdkitVersion,
            )
        path.replace(output)
    finally:
        path.unlink(missing_ok=True)
    print(f"Saved {len(data)} fingerprints of 2048 bits to {output.name}")


if __name__ == "__main__":
    main()
