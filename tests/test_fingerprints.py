import unittest

import pandas as pd

from src.egfrank.fingerprints import generate_fingerprints


class FingerprintTests(unittest.TestCase):
    def test_generates_binary_2048_bit_rows(self):
        data = pd.DataFrame({
            "molecule_chembl_id": ["CHEMBL1", "CHEMBL2"],
            "canonical_smiles": ["c1ccccc1O", "CCN"],
        })
        matrix = generate_fingerprints(data)
        self.assertEqual(matrix.shape, (2, 2048))
        self.assertTrue(set(matrix.ravel()).issubset({0, 1}))
        self.assertTrue((matrix.sum(axis=1) > 0).all())

    def test_invalid_smiles_fails_with_compound_id(self):
        data = pd.DataFrame({
            "molecule_chembl_id": ["CHEMBL_BAD"],
            "canonical_smiles": ["not smiles"],
        })
        with self.assertRaisesRegex(ValueError, "CHEMBL_BAD"):
            generate_fingerprints(data)


if __name__ == "__main__":
    unittest.main()
