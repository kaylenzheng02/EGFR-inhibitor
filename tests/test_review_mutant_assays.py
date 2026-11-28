import unittest

import pandas as pd

from src.egfrank.review_mutant_assays import (
    build_l858r_audit,
    classify_endpoint,
    classify_l858r_assay,
    label_from_range,
    measurement_exclusion_reason,
    parse_variant_signature,
)


class VariantParsingTests(unittest.TestCase):
    def test_single_l858r_aliases(self):
        self.assertEqual(parse_variant_signature("EGFR L858R mutant")[0], "L858R")
        self.assertEqual(parse_variant_signature("EGFR Leu858Arg mutant")[0], "L858R")

    def test_combination_has_canonical_order(self):
        forward = parse_variant_signature("EGFR L858R/T790M/C797S mutant")[0]
        reverse = parse_variant_signature("EGFR C797S/T790M/L858R mutant")[0]
        self.assertEqual(forward, "L858R+T790M+C797S")
        self.assertEqual(reverse, forward)

    def test_exon_19_alias(self):
        signature, unknown = parse_variant_signature("EGFR delE746_A750/T790M")
        self.assertEqual(signature, "ex19del+T790M")
        self.assertEqual(unknown, "")

    def test_unrecognized_substitution_is_flagged(self):
        signature, unknown = parse_variant_signature("EGFR L858R/L861R mutant")
        self.assertEqual(signature, "L858R")
        self.assertEqual(unknown, "L861R")

    def test_residue_range_is_not_a_substitution(self):
        signature, unknown = parse_variant_signature(
            "EGFR L858R mutant (669 to 1210 residues) using peptide substrate"
        )
        self.assertEqual(signature, "L858R")
        self.assertEqual(unknown, "")

    def test_insect_expression_is_not_an_insertion(self):
        signature, unknown = parse_variant_signature(
            "EGFR L858R mutant expressed in Sf9 insect cells using substrate"
        )
        self.assertEqual(signature, "L858R")
        self.assertEqual(unknown, "")

    def test_audit_ignores_assays_without_a_recognized_variant(self):
        data = pd.DataFrame({
            "assay_chembl_id": ["CHEMBL_A", "CHEMBL_B"],
            "assay_description": [
                "Wild-type EGFR kinase assay using peptide substrate",
                "EGFR L858R kinase assay using peptide substrate",
            ],
            "assay_type": ["B", "B"],
        }, index=[10, 20])
        audit = build_l858r_audit(data)
        self.assertEqual(len(audit), 1)
        self.assertEqual(audit.assay_chembl_id.tolist(), ["CHEMBL_B"])


class EndpointTests(unittest.TestCase):
    def test_biochemical_endpoint(self):
        description = "Inhibition of EGFR L858R using poly(Glu:Tyr) substrate by kinase assay"
        self.assertEqual(classify_endpoint(description), "biochemical_enzyme")

    def test_cellular_endpoint(self):
        description = "Inhibition of EGFR L858R phosphorylation in H3255 cells"
        self.assertEqual(classify_endpoint(description), "cellular_phosphorylation")

    def test_combination_is_not_single_l858r_candidate(self):
        result = classify_l858r_assay(
            "Inhibition of EGFR L858R/T790M using peptide substrate", "B"
        )
        self.assertFalse(result["automated_candidate"])
        self.assertEqual(result["exclusion_reason"], "combination_variant")

    def test_wild_type_mixture_is_not_candidate(self):
        result = classify_l858r_assay(
            "Kinase assay of wild-type EGFR and L858R using peptide substrate", "B"
        )
        self.assertFalse(result["automated_candidate"])
        self.assertEqual(result["exclusion_reason"], "mixed_wild_type_and_mutant")


class LabelAndQualityTests(unittest.TestCase):
    def test_cutoff_labels(self):
        self.assertEqual(label_from_range(100, 1000), "active")
        self.assertEqual(label_from_range(1001, 2000), "inactive")
        self.assertEqual(label_from_range(800, 1200), "conflict")

    def test_valid_measurement_has_no_exclusion_reason(self):
        row = pd.Series({
            "standard_type": "IC50",
            "standard_relation": "=",
            "standard_units": "nM",
            "standard_value": 100,
            "data_validity_comment": None,
            "potential_duplicate": 0,
            "molecule_chembl_id": "CHEMBL1",
            "canonical_smiles": "CCO",
        })
        self.assertEqual(measurement_exclusion_reason(row), "")

    def test_invalid_measurement_reports_all_reasons(self):
        row = pd.Series({
            "standard_type": "IC50",
            "standard_relation": ">",
            "standard_units": "uM",
            "standard_value": -1,
            "data_validity_comment": "Outside typical range",
            "potential_duplicate": 1,
            "molecule_chembl_id": None,
            "canonical_smiles": "not smiles",
        })
        reasons = set(measurement_exclusion_reason(row).split("|"))
        self.assertEqual(reasons, {
            "relation_not_exact", "units_not_nM", "value_not_positive_finite",
            "data_validity_comment", "potential_duplicate", "missing_compound_id",
            "invalid_smiles",
        })


if __name__ == "__main__":
    unittest.main()
