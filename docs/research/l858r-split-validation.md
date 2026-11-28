# EGFR L858R split validation

This split freezes the held-out set for the first L858R benchmark. It assigns whole connected groups rather than individual compounds. A connected group combines compounds that share a Murcko scaffold or an identical 2,048-bit Morgan fingerprint, including transitive connections between those relationships.

## Split settings

| Setting | Value |
| --- | ---: |
| Test target | 20% |
| Random seed | 42 |
| Candidate assignments searched | 25,000 |
| Connected groups | 156 |

## Final assignment

| Set | Active | Inactive | Total | Scaffolds | Connected groups |
| --- | ---: | ---: | ---: | ---: | ---: |
| Training | 225 | 38 | 263 | 130 | 125 |
| Test | 56 | 10 | 66 | 31 | 31 |

The test partition contains 20.1% of all compounds, 19.9% of active compounds, and 20.8% of inactive compounds. Its always-active accuracy is 84.8%.

The held-out test set is now frozen. Model selection, feature decisions, hyperparameter choices, and score-threshold decisions must use only the training partition and training-only grouped validation.

## Leakage checks

| Check | Crossing groups |
| --- | ---: |
| Compound IDs | 0 |
| Canonical structures | 0 |
| Murcko scaffolds | 0 |
| Identical fingerprints | 0 |
| Connected split groups | 0 |

All leakage checks pass at zero. Both activity classes occur in both partitions. The split tests chemical frameworks absent from training while also preventing identical model inputs from appearing on both sides.

## Interpretation

The test set contains 10 inactive compounds. Consequently, each inactive error changes inactive recall by 10.0 percentage points. The final benchmark must report raw confusion-matrix counts alongside summary metrics and should not use this test set for iterative model development.
