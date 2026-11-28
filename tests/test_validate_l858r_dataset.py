import unittest

import pandas as pd

from src.egfrank.validate_l858r_dataset import (
    expected_evidence_grade,
    expected_label,
    scaffold_for_smiles,
    validate_table,
)


def valid_table():
    return pd.DataFrame([
        {
            "molecule_chembl_id": "CHEMBL1",
            "variant": "L858R",
            "canonical_smiles": "c1ccccc1O",
            "label": "active",
            "min_ic50_nm": 100,
            "max_ic50_nm": 1000,
            "measurements": 1,
            "evidence_grade": "C",
            "activity_ids": "1",
            "assay_ids": "ASSAY1",
            "document_ids": "DOC1",
        },
        {
            "molecule_chembl_id": "CHEMBL2",
            "variant": "L858R",
            "canonical_smiles": "c1ccncc1",
            "label": "inactive",
            "min_ic50_nm": 1001,
            "max_ic50_nm": 2000,
            "measurements": 2,
            "evidence_grade": "A",
            "activity_ids": "2|3",
            "assay_ids": "ASSAY1|ASSAY2",
            "document_ids": "DOC1",
        },
    ])


def decisions():
    return pd.DataFrame({
        "assay_chembl_id": ["ASSAY1", "ASSAY2"],
        "manual_decision": ["include", "include"],
        "review_notes": ["reviewed", "reviewed"],
    })


class LabelAndEvidenceTests(unittest.TestCase):
    def test_label_boundary(self):
        self.assertEqual(expected_label(10, 1000), "active")
        self.assertEqual(expected_label(1001, 2000), "inactive")
        self.assertEqual(expected_label(900, 1100), "conflict")

    def test_evidence_grades(self):
        self.assertEqual(expected_evidence_grade(1, 1, 1), "C")
        self.assertEqual(expected_evidence_grade(2, 1, 1), "B")
        self.assertEqual(expected_evidence_grade(2, 2, 1), "A")


class DatasetValidationTests(unittest.TestCase):
    def test_valid_dataset_is_enriched_with_scaffolds(self):
        enriched = validate_table(valid_table(), decisions())
        self.assertEqual(len(enriched), 2)
        self.assertTrue(enriched.scaffold.ne("").all())

    def test_duplicate_compound_id_fails(self):
        data = valid_table()
        data.loc[1, "molecule_chembl_id"] = "CHEMBL1"
        with self.assertRaisesRegex(ValueError, "compound IDs are duplicated"):
            validate_table(data, decisions())

    def test_duplicate_structure_under_two_ids_fails(self):
        data = valid_table()
        data.loc[1, "canonical_smiles"] = data.loc[0, "canonical_smiles"]
        with self.assertRaisesRegex(ValueError, "same RDKit canonical structure"):
            validate_table(data, decisions())

    def test_label_disagreement_fails(self):
        data = valid_table()
        data.loc[0, "label"] = "inactive"
        with self.assertRaisesRegex(ValueError, "labels disagree"):
            validate_table(data, decisions())

    def test_unreviewed_assay_fails(self):
        data = valid_table()
        data.loc[0, "assay_ids"] = "UNREVIEWED"
        with self.assertRaisesRegex(ValueError, "without an include decision"):
            validate_table(data, decisions())

    def test_scaffold_for_ringless_structure(self):
        self.assertEqual(scaffold_for_smiles("CCO"), "ACYCLIC")


if __name__ == "__main__":
    unittest.main()
