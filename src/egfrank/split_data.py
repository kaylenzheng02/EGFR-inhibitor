from pathlib import Path

import pandas as pd
from rdkit import Chem
from rdkit.Chem.Scaffolds import MurckoScaffold


ROOT = Path(__file__).resolve().parents[2]
DATA_DIR = ROOT / "data" / "processed"
INPUT_FILE = DATA_DIR / "egfr_training.csv"
OUTPUT_FILE = DATA_DIR / "egfr_scaffold_split.csv"

TESt_FRACTION = 0.20
SEED = 42


def get_scaffold(smiles):
    molecule = Chem.MolFromSmiles(smiles)
    if molecule is None:
        raise ValueError(f"RDKit could not read this SMILES: {smiles}")

    scaffold = MurckoScaffold.MurckoScaffoldSmiles(mol=molecule)