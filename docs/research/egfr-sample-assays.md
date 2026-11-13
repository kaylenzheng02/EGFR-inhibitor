# EGFR sample assay review

The 100 screened measurements use the same target identifier, units, and IC50 relation. The six assays measure different responses. Use assay ID as the default comparison group. Keep each activity ID and its assay ID when displaying or exporting measurements.

| ChEMBL assay | Sample rows | What the assay measures | Source document |
| --- | ---: | --- | --- |
| [CHEMBL674637](https://www.ebi.ac.uk/chembl/api/data/assay/CHEMBL674637.json) | 17 | Inhibition of EGFR kinase tyrosine phosphorylation. The description gives no cell context. | CHEMBL1134862 |
| [CHEMBL621151](https://www.ebi.ac.uk/chembl/api/data/assay/CHEMBL621151.json) | 10 | Inhibition of human EGFR autophosphorylation in A431 cells. | CHEMBL1134862 |
| [CHEMBL615325](https://www.ebi.ac.uk/chembl/api/data/assay/CHEMBL615325.json) | 8 | Inhibition of ligand-induced proliferation in 3T3 cells expressing human EGFR. | CHEMBL1134862 |
| [CHEMBL674643](https://www.ebi.ac.uk/chembl/api/data/assay/CHEMBL674643.json) | 48 | Inhibition of EGF-dependent receptor autophosphorylation. The description gives no cell context. | CHEMBL1124610 |
| [CHEMBL677833](https://www.ebi.ac.uk/chembl/api/data/assay/CHEMBL677833.json) | 16 | In vitro inhibition of receptor kinase isolated from A431 cells. A431 identifies the source of the receptor here. | CHEMBL1145114 |
| [CHEMBL679674](https://www.ebi.ac.uk/chembl/api/data/assay/CHEMBL679674.json) | 1 | EGFR inhibition. The record gives too little detail to classify the readout or preparation. | CHEMBL1132739 |

## Worked example

All three measurements of CHEMBL68920 come from source document CHEMBL1134862. The kinase phosphorylation assay reports 41 nM, the A431 receptor autophosphorylation assay reports 300 nM, and the 3T3 proliferation assay reports 7,820 nM. These values describe concentrations associated with three distinct responses. The records establish that the endpoints differ. They do not establish a mechanism for the potency differences. Preserve all three values and label each with its endpoint.

## Comparison decisions

- Compare compounds within a single assay ID first. CHEMBL674643 provides the largest such group with 48 measurements. Preserve activity-level records if a compound occurs more than once.
- Across assay IDs, show values side by side with endpoint labels. Pooling kinase-related records would require checking the source methods for receptor preparation, substrate, ATP conditions, and readout. Shared phosphorylation terminology alone does not establish equivalence.
- Keep the 3T3 proliferation measurements separate from kinase and receptor phosphorylation measurements when defining potency labels. Keep the single underspecified CHEMBL679674 measurement available for exploration and defer its inclusion in an endpoint-specific training set.
- None of the six ChEMBL records supplies a variant sequence. Treat mutation status as unspecified unless the primary publication confirms it.
