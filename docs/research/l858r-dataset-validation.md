# EGFR L858R dataset validation

The reviewed L858R dataset passed the pre-split identity, structure, label, provenance, scaffold, and fingerprint checks. This report characterizes the frozen input to the future benchmark; it does not create a train/test split or evaluate a model.

## Dataset

| Check | Result |
| --- | ---: |
| Compounds | 329 |
| Active | 281 |
| Inactive | 48 |
| Always-active baseline | 85.4% |
| Eligible measurements | 343 |
| Compounds with multiple measurements | 8 |
| Compounds near 800–1,200 nM | 8 |
| Duplicate RDKit structures | 0 |

Evidence grades are A: 7, B: 1, and C: 321. Grade C means that one eligible measurement supports the label.

## Provenance

All contributing assay IDs have an explicit include decision. Of 47 reviewed assays, 42 contribute an exact eligible measurement from 40 source documents. Five reviewed assays contribute no final row after activity-level quality screening: CHEMBL4342867, CHEMBL4379304, CHEMBL4603464, CHEMBL4771365, CHEMBL5150414.

## Chemical coverage

| Check | Result |
| --- | ---: |
| Murcko scaffold groups | 161 |
| Scaffolds containing active compounds | 142 |
| Scaffolds containing inactive compounds | 27 |
| Scaffolds containing both classes | 8 |
| Singleton scaffolds | 109 |
| Largest scaffold group | 18 compounds |
| Acyclic compounds | 0 |

The scaffold definition matches the existing EGFR benchmark: an RDKit Murcko scaffold without chirality, with ringless compounds assigned to `ACYCLIC`.

## Fingerprints

The saved matrix contains 329 rows and 2048 bits per row. It uses radius 2 Morgan fingerprints with chirality enabled under RDKit 2026.09.1. Set-bit counts range from 35 to 105, with a median of 68.0.

7 fingerprint groups contain more than one structure, accounting for 9 rows beyond the first member of each group. 6 groups span more than one Murcko scaffold, and 1 group contains both an active and an inactive compound. These are distinct structures that produce the same hashed fingerprint, not duplicate structures. Every member of a repeated fingerprint group must remain on the same side of the future train/test split.

The cross-label fingerprint group contains CHEMBL5192038 (active), CHEMBL5205073 (inactive). A model using this fingerprint representation must assign those compounds the same score even though their experimental labels differ.

## Benchmark implication

The dataset is ready for a leakage-safe scaffold split. The main constraint is the small and shallow negative evidence: 48 inactive compounds occupy 27 scaffolds, while 321 of all 329 labels rely on one eligible measurement. The benchmark should therefore retain an untouched scaffold-held-out test set and use scaffold-grouped validation within training. Accuracy must be compared with the 85.4% always-active baseline and accompanied by balanced and inactive-class metrics.
