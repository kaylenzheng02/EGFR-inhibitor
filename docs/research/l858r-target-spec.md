# EGFR L858R prediction target specification

## Status

This specification fixes the first variant-specific prediction target for EGFRank. It defines the scientific question and the evidence that may be used before the dataset, model, or frontend is changed.

## Prediction question

Given a chemical structure, how strongly should it be prioritized for a biochemical experiment measuring inhibition of the human EGFR L858R kinase variant?

The first variant model will rank compounds by whether they meet the project's activity definition for **single-mutant EGFR L858R**. It will not pool L858R with wild-type EGFR or with compound variants such as L858R/T790M.

## Prediction unit and output

- One labeled example is one compound tested against the single L858R variant.
- The compound is represented by its canonical SMILES and a molecular fingerprint.
- The supervised outcome is `active` or `inactive` under the labeling rule below.
- The model output is a ranking score for laboratory prioritization, not a prediction of patient response, clinical efficacy, safety, or dose.

## Eligible evidence

An activity record is eligible only when all of the following are true:

1. The assay can be attributed to human EGFR carrying L858R and no other named EGFR alteration.
2. The endpoint is biochemical kinase or enzyme inhibition. The description must provide evidence of an enzyme-activity experiment; a ChEMBL assay-type classification alone is insufficient.
3. The standardized activity type is IC50, the relation is exact (`=`), and the unit is nM.
4. The IC50 is positive and finite.
5. The ChEMBL data-validity comment is empty and the potential-duplicate flag is zero.
6. The record has a compound ID and a canonical SMILES that RDKit can parse.
7. All eligible records assigned to one compound ID use a consistent structure.

The assay audit must preserve the activity, assay, and source-document identifiers so every label can be traced to its evidence.

## Excluded evidence

The first target excludes:

- Wild-type EGFR.
- L858R combined with T790M, C797S, an exon 19 deletion, or any other named alteration.
- T790M, C797S, exon 19 deletions, L861Q, G719X, and other variants without L858R as the sole named alteration.
- Assays whose descriptions name several EGFR proteins but do not permit the activity row to be assigned to one exact variant.
- Cellular phosphorylation, autophosphorylation, proliferation, viability, binding-affinity, degradation, and other non-biochemical endpoints.
- Approximate or censored measurements such as `<`, `>`, `<=`, or `>=`.
- Measurements that fail the quality, identity, unit, or structure checks above.

Ambiguous records enter the assay-review table with an exclusion reason; they are not silently assigned to L858R.

## Labeling rule

Eligible measurements are grouped by compound ID and exact variant signature.

- `active`: every eligible IC50 is at or below 1,000 nM.
- `inactive`: every eligible IC50 is above 1,000 nM.
- `conflict`: eligible measurements occur on both sides of 1,000 nM.

Exactly 1,000 nM belongs to the active class. Conflicts are preserved for review and excluded from binary model training. The terms active and inactive describe this project cutoff rather than a universal pharmacological boundary.

## Why L858R is first

The original all-variant screen contains 294 consistently labeled single-L858R compound pairs: 246 active and 48 inactive. The target-specific screen expands mutation-alias and biochemical-method recognition and produces 329 pairs: 281 active and 48 inactive. The difference comes from six assays that the two automated rules classify differently, including Leu858Arg notation and HTRF descriptions. All 47 target-specific candidate assays have been reviewed and included, so these are the final counts for the frozen snapshot under the current target and screening rules. L858R has the strongest combination of relevance and two-class coverage among the single variants in the snapshot; the current extraction has too few single-variant exon 19 deletion or C797S examples to support a two-class model.

## Relationship to the twelve-compound goal

This model will support the later selection of twelve compounds by providing a target-specific ranking. The held-out benchmark set is used to evaluate the ranking because its experimental labels are already known. A prospective twelve-compound shortlist must instead come from a separate candidate library whose outcomes were not used for training or evaluation.

## Completion gate for target preparation

The target is ready for modeling only after:

1. Every candidate L858R assay has an explicit reviewed decision and exclusion reason where applicable.
2. The labeled dataset passes identity, structure, unit, cutoff, conflict, and traceability checks.
3. Final active, inactive, conflict, measurement, assay, compound, and scaffold counts are recorded.
4. The split procedure keeps identical compounds, structures, and chemical scaffolds on one side of the train/test boundary.
5. The approved dataset and its source snapshot can be reproduced without changing this target definition.
