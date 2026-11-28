import unittest

import pandas as pd

from src.egfrank.split_data import (
    assign_split,
    choose_test_groups,
    connected_split_groups,
    validate_l858r_split,
)


def example_groups():
    rows = []
    for index in range(20):
        rows.append({
            "molecule_chembl_id": f"CHEMBL{index}",
            "canonical_smiles": f"SMILES{index}",
            "rdkit_canonical_smiles": f"SMILES{index}",
            "label": "inactive" if index % 5 == 0 else "active",
            "scaffold": f"S{index}",
            "fingerprint_group": f"F{index}",
        })
    return pd.DataFrame(rows)


class ConnectedGroupTests(unittest.TestCase):
    def test_transitive_scaffold_and_fingerprint_connections(self):
        data = pd.DataFrame({
            "molecule_chembl_id": ["A", "B", "C", "D"],
            "scaffold": ["S1", "S1", "S2", "S3"],
            "fingerprint_group": ["F1", "F2", "F2", "F3"],
        })
        groups = connected_split_groups(data)
        self.assertEqual(groups.iloc[0], groups.iloc[1])
        self.assertEqual(groups.iloc[1], groups.iloc[2])
        self.assertNotEqual(groups.iloc[2], groups.iloc[3])

    def test_group_ids_do_not_depend_on_row_order(self):
        data = pd.DataFrame({
            "molecule_chembl_id": ["A", "B", "C"],
            "scaffold": ["S1", "S1", "S2"],
            "fingerprint_group": ["F1", "F2", "F3"],
        })
        first = dict(zip(data.molecule_chembl_id, connected_split_groups(data)))
        shuffled = data.sample(frac=1, random_state=7).reset_index(drop=True)
        second = dict(zip(shuffled.molecule_chembl_id, connected_split_groups(shuffled)))
        self.assertEqual(first, second)


class SplitSelectionTests(unittest.TestCase):
    def test_split_selection_is_reproducible(self):
        data = example_groups()
        data["split_group"] = connected_split_groups(data)
        first = choose_test_groups(data, "split_group", attempts=1000)
        second = choose_test_groups(data, "split_group", attempts=1000)
        self.assertEqual(first, second)

    def test_assignment_keeps_all_grouping_keys_together(self):
        data = example_groups()
        data["split_group"] = connected_split_groups(data)
        test_groups = choose_test_groups(data, "split_group", attempts=1000)
        output = assign_split(data, "split_group", test_groups)
        checks = validate_l858r_split(output)
        self.assertTrue(all(value == 0 for value in checks.values()))
        self.assertTrue((pd.crosstab(output.split, output.label) > 0).all().all())


if __name__ == "__main__":
    unittest.main()
