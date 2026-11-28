# EGFRank project history

**Research question.** Can the chemical structure of a compound help us prioritize it for an experiment measuring inhibition of human EGFR kinase?

**Abstract.** I used previously reported experiments to construct a dataset for answering this question. A pilot of 100 measurements exposed a central problem: the same compound can give different results when scientists measure different biological events. I traced those measurements to their experiments, built a small explorer, and applied an assay selection and quality screen to a larger download. From 26,600 records, I obtained 1,600 compounds with consistent labels under a project cutoff of 1,000 nM. I generated molecular fingerprints for all 1,600 compounds and reserved 320 compounds for testing, keeping each chemical scaffold entirely on one side of the split. A Random Forest trained on the other 1,280 compounds correctly classified 306 of the 320 held-out compounds and achieved a ROC-AUC of 0.929.

## Step 1. Research question and project setup

Can chemical structure help prioritize compounds for testing against human epidermal growth factor receptor (EGFR) kinase? I chose general EGFR inhibition for the first model. I organized a Git project in VS Code with a Python environment for processing data and a React and TypeScript frontend for exploring the results. I consolidated the original nested folders so commands run from one repository root. · [Project setup](README.md)

### Biological background

Epidermal growth factor receptor (EGFR) is a protein that spans the membrane surrounding a cell. A growth factor is a signaling molecule. When it binds the outside portion of EGFR, receptor molecules can pair up and activate the part inside the cell. This part is an **enzyme**, a molecule that speeds up a chemical reaction. EGFR belongs to the **tyrosine kinase** family. A kinase transfers a phosphate group, typically from a molecule called ATP, onto a receiving molecule called a substrate. Tyrosine is an amino acid, one of the building blocks of proteins. EGFR can attach phosphate to particular tyrosines. The transfer is called **phosphorylation**. [UniProt record for human EGFR](https://www.uniprot.org/uniprotkb/P00533/entry)

In an activated receptor pair, one receptor can phosphorylate tyrosines on its partner. This is called **autophosphorylation** because the receptor protein is itself the recipient. Other proteins can recognize these sites and pass the signal onward, influencing processes such as cell growth. This sequence gives scientists several things to measure: the enzyme transferring phosphate to a test substrate, the receptor becoming phosphorylated inside a cell, or the eventual growth of a population of cells. Each is a different observation. [EGFR function and signaling](https://www.uniprot.org/uniprotkb/P00533/entry)

Non-small cell lung cancer (NSCLC) motivates this project’s focus on EGFR. The analysis uses published measurements of EGFR inhibition from a broader set of laboratory experiments.

The intended application is early compound prioritization. A researcher with a limited testing budget could eventually use rankings to choose which structures to study first. The outcome defined in this project is inhibition in selected laboratory experiments. Patient response, safety, and the dose needed to treat a person would require different evidence.

## Step 2. Pilot data and checks

I downloaded 100 real ChEMBL measurements with exact IC50 values in nanomolar units for human EGFR (CHEMBL203). The sample covered 82 compounds in six assays. I checked IDs, units, numerical values, and whether RDKit could read each molecular structure. All 100 passed. The preprocessor preserves the original values and calculates pIC50 as 9 minus the base-10 logarithm of IC50 in nM. For example, 41 nM gives 7.3872. [Pilot query](https://www.ebi.ac.uk/chembl/api/data/activity.json?target_chembl_id=CHEMBL203&standard_type=IC50&standard_units=nM&standard_relation=%3D&limit=100) · [Screening script](src/egfrank/preprocess.py)

### Reading a ChEMBL record

**ChEMBL** is a public database maintained at the European Bioinformatics Institute. It connects chemical structures with measurements of **bioactivity**, meaning effects observed when a substance is tested in a biological experiment. Its records come from sources including scientific publications and deposited datasets. The experiments in EGFRank were conducted by the original researchers. My contribution was collecting, checking, interpreting, and organizing the records for a specific question. [What is ChEMBL?](https://www.ebi.ac.uk/training/online/courses/chembl-quick-tour/what-is-chembl/) · [How ChEMBL curates data](https://www.ebi.ac.uk/training/online/courses/chembl-quick-tour/what-is-chembl/how-is-chembl-data-curated/)

| Field | Question it answers | Why I kept it |
| --- | --- | --- |
| Target ID | Which biological entity did ChEMBL associate with the result? | CHEMBL203 identifies human EGFR. |
| Compound ID | Which molecule did the researchers test? | The same ID connects results across experiments. |
| Assay ID | Which experimental test produced the result? | Looking up the assay reveals the measured endpoint. |
| Activity ID | Which individual result is this? | It lets me trace a training label back to specific observations. |
| Canonical SMILES | How is the compound structured? | SMILES encodes atoms, bonds, branches, and rings as text. RDKit reads this text as a molecular structure. |
| Type, relation, value, units | What quantity was reported, and how precisely? | All four are needed to distinguish an exact IC50 of 50 nM from a lower bound above 50 nM. |

### Interpreting the numbers

**IC50** is the concentration that produces half-maximal inhibition in the specified test. Concentration describes how much substance is present per volume. One nanomolar (nM) is one billionth of a mole per liter. A value of 41 nM means that the assay reached its halfway inhibition point at that concentration. A relation of '=' means that the database reports a numerical estimate directly. Experimental uncertainty still exists. A relation of '>' tells us only that the estimate lies above a boundary.

The logarithmic value **pIC50** makes large differences in concentration easier to read. Starting with 41 nM, convert to molar units by multiplying by one billionth, then take the negative base-10 logarithm. This gives approximately 7.3872. As a second example, 100 nM corresponds to 7 and 1,000 nM corresponds to 6. A difference of one pIC50 unit represents a tenfold concentration difference. Meaningful potency comparisons still depend on the assay context.

I used pandas to read the file as a table, NumPy for numerical calculations, and RDKit to check that each molecular string parsed. Parsing confirms that the software can interpret the structure as written. It cannot confirm that the original laboratory tested the correct material. [RDKit guide](https://www.rdkit.org/docs/GettingStartedInPython.html)

## Step 3. Assay interpretation

I looked up what each pilot experiment measured. CHEMBL68920 had values of 41 nM for kinase phosphorylation, 300 nM for receptor autophosphorylation in A431 cells, and 7,820 nM for proliferation in 3T3 cells. These are distinct experimental outcomes. Keeping the assay ID with each value lets a reader see why averaging all three would obscure useful information. [Six-assay review with links to the original records](docs/research/egfr-sample-assays.md)

| Observation for CHEMBL68920 | What the experiment followed | What can be concluded |
| --- | --- | --- |
| 41 nM, CHEMBL674637 | Tyrosine phosphorylation by EGFR kinase | This concentration describes inhibition of the phosphorylation readout under that assay's conditions. |
| 300 nM, CHEMBL621151 | EGFR autophosphorylation in A431 cells, a laboratory cell line | This is a receptor-level response measured in a cellular environment. |
| 7,820 nM, CHEMBL615325 | Growth-factor-induced proliferation of 3T3 cells expressing human EGFR | Proliferation means an increase in cell number through division. This result describes the response of the cell population. |

These results come from the same source document and span roughly 190-fold from lowest to highest concentration. A whole-cell growth measurement can depend on several processes between adding the compound and observing division. The records tell us which outcome changed. They do not establish which process caused this particular difference. Averaging the three would yield a number with no single experimental endpoint. This case is why the later labeling rule operates only after assay selection. [Kinase assay](https://www.ebi.ac.uk/chembl/api/data/assay/CHEMBL674637.json) · [Receptor assay](https://www.ebi.ac.uk/chembl/api/data/assay/CHEMBL621151.json) · [Cell growth assay](https://www.ebi.ac.uk/chembl/api/data/assay/CHEMBL615325.json)

## Step 4. Interactive demo

I connected the screened sample to a React dashboard with counts calculated from the data and a compound search showing the assay ID, endpoint, IC50, and pIC50. I verified the dashboard counts and the three search results for CHEMBL68920. The frontend passed its lint and production build checks. I also fixed a Windows file replacement issue so rerunning the preprocessor skips identical exports and retries a temporarily blocked update. [React app](frontend/src/App.tsx) · [Data export](src/egfrank/preprocess.py)

The pilot also tested the path from a published measurement to what a user actually sees. The raw CSV is a table with one row per result. Running the Python script produces a screened CSV for inspection and a JSON file, a structured text format the browser can load. React reads the JSON to count distinct compound and assay IDs. Searching a compound filters the measurement rows but preserves each assay label. The dashboard therefore shows 100 measurements even though the search index contains only 82 distinct compounds. [Raw pilot CSV](data/raw/egfr_sample.csv) · [Browser data](frontend/public/data/egfr_sample.json)

The three-result search provided an end-to-end check with a known answer. A working search box alone would not have exposed an accidental merge of repeated compounds. Checking for all three concentrations and their corresponding assays tested precisely the information the explorer was built to preserve.

## Step 5. Full dataset acquisition

On November 13, 2026, the downloader collected 26,600 IC50 measurements for target CHEMBL203 across 54 API pages. The snapshot contains 14,776 distinct compound IDs and 2,425 assay IDs. I verified the reported total against the saved CSV and checked that every activity ID was present and unique. A dated manifest records the exact query, UTC retrieval time, and file checksum. The full query preserves units and measurement bounds for later screening. [Exact full query](https://www.ebi.ac.uk/chembl/api/data/activity.json?target_chembl_id=CHEMBL203&standard_type=IC50&order_by=activity_id&limit=500&only=activity_id%2Cmolecule_chembl_id%2Ccanonical_smiles%2Cstandard_type%2Cstandard_relation%2Cstandard_value%2Cstandard_units%2Cassay_chembl_id%2Ctarget_chembl_id%2Cdocument_chembl_id%2Cassay_type%2Cassay_description%2Ctarget_organism%2Ctarget_pref_name%2Cpchembl_value%2Cdata_validity_comment%2Cpotential_duplicate%2Cactivity_comment%2Ctype%2Crelation%2Cvalue%2Cunits) · [Snapshot manifest](data/raw/egfr_ic50_20261113T184940536954Z/manifest.json) · [Downloader](src/egfrank/fetch_data.py)

An **application programming interface** (API) lets a program ask a database for records using explicit conditions. The larger query selected human EGFR and the standardized measurement type IC50. Unlike the pilot, it retained the reported units and relations for later review. The 26,600 count therefore describes this specific query at the time of retrieval. It includes records that would fail the later training criteria.

The API returned results in pages of up to 500 records. The downloader ordered them by activity ID, followed each next-page link, checked the expected offset, and compared the final number of distinct activity IDs with the reported total. These checks catch problems such as silently downloading only the first page or counting one page twice. A **checksum** is a value computed from the saved file's bytes. Another reader can compare that value against the manifest to verify that a local copy matches this snapshot. A matching checksum establishes file identity; the scientific suitability of the records is a separate question addressed in Step 6.

## Step 6. Training dataset preparation

I first reviewed three enzyme assays, which yielded 950 eligible compounds. I then screened descriptions across all 2,425 assay IDs. The expanded rule kept the three original assays and added binding-type records whose descriptions identify wild-type EGFR and a substrate or kinase activity assay while excluding ambiguous descriptions. This automated selection can miss relevant experiments and cannot independently confirm every assay protocol or variant.

| Stage | Result |
| --- | --- |
| Assay description screen | 74 candidate assays, 3,128 measurements |
| Exact, positive, finite nM values, accepted quality flags, compound IDs and readable structures | 1,693 measurements, 1,602 compound IDs |
| Repeated measurements | 70 compound IDs with multiple eligible measurements |
| Every eligible IC50 at or below 1,000 nM | 1,396 active |
| Every eligible IC50 above 1,000 nM | 204 inactive |
| Measurements on both sides | 2 compounds set aside |

I approved 1,000 nM as the project cutoff. The training CSV has one row per compound, its structure, label, IC50 range, measurement count, and source activity and assay IDs. The terms active and inactive describe this project's cutoff. The active share is 87.25%, which is the full dataset accuracy of always guessing active. [Screening and label code](src/egfrank/review_assays.py) · [Assay decisions](docs/research/phase2-assay-review.csv) · [Local training CSV](data/processed/egfr_training.csv) · [Two conflicts](data/processed/egfr_conflicts.csv)

### Why the eligibility rules matter

The original three assays entered through reviewed descriptions. For additional assays, the script looked for a ChEMBL binding-type classification ('B'), an EGFR name, an explicit wild-type designation, and evidence of an enzyme activity test, such as a substrate or a named kinase assay method. **Wild type** refers here to the reference form of EGFR. Descriptions mentioning mutations, cell proliferation, viability, binding affinity, or other ambiguous endpoints were flagged out by the rule. A binding-type classification alone is insufficient: measuring how tightly a compound binds a protein answers a different question from measuring how much it slows a reaction. The three original assays followed their own inclusion rule, so the full set should not be described as independently verified wild-type assays.

Of the 3,128 candidate measurements, 1,693 passed the combined checks. The screen required exact positive finite values in nM, readable structures and compound IDs, an empty ChEMBL data-validity comment, and a potential-duplicate flag of zero. A flagged duplicate is a possible repeated report of the same evidence. Multiple eligible results for the same compound are still valuable observations. These two types of repetition were handled separately. The difference of 1,435 records is the total removed by all checks combined; it cannot be attributed to any one reason without a separate exclusion analysis.

### Turning measurements into labels

One row in the final table represents one ChEMBL compound ID. I checked for consistent structure strings within each ID, counted its measurements, and recorded the smallest and largest qualifying IC50. If the maximum was at or below 1,000 nM, every qualifying measurement met the project definition of active. If the minimum exceeded 1,000 nM, every result fell on the inactive side. Exactly 1,000 nM belongs to the active group. The cutoff equals one micromolar or pIC50 6. It defines a project decision, not a universal boundary between useful and ineffective compounds.

For example, two hypothetical results of 200 and 800 nM agree under this rule despite differing fourfold. Results of 800 and 1,200 nM conflict even though they are closer numerically. The rule tests agreement about the cutoff, not whether laboratories obtained identical measurements. The export keeps the activity IDs, assay IDs, range, and measurement count so readers can revisit how any label was assigned.

## Step 7. Interpret the results

Phase 2 produced 1,600 examples with structures and labels. The saved files were checked for unique compound IDs, compliance with the cutoff, and separation of the two conflicting compounds. The training and conflict files together account for all 1,693 eligible measurements. These checks establish that the export implements the agreed rule. They do not measure predictive performance.

**Limitations.** Chemical series, assay conditions, and publication choices influence which compounds appear in this collection. Automated description matching may overlook eligible experiments or admit a poorly specified one. Measurements close to the cutoff can change class with small experimental differences. Only 204 labeled compounds fall in the inactive group. Consequently, someone predicting “active” for every row would be correct 1,396 times out of 1,600 while identifying zero inactive compounds. The resulting 87.25% is a property of the dataset, not a score achieved by a trained model.

## Step 8. Represent each compound as a molecular fingerprint

I converted the SMILES structure for each of the 1,600 labeled compounds into a Morgan fingerprint using RDKit. The algorithm examines the neighborhood around each atom, extending up to two bonds away. It records the features it finds in a row of 2,048 bits. A bit switches on when a feature maps to that position. Several features can map to the same bit. I enabled chirality so the calculation can account for three-dimensional arrangements specified in the input structure.

The saved fingerprint matrix has 1,600 rows and 2,048 columns. I verified that each row contained at least one set bit and that the compound IDs and activity labels matched the training CSV in order. The compressed NumPy file stores those IDs and labels alongside the matrix, the fingerprint settings, and the RDKit version. This keeps each numerical representation traceable to its original experimental record. [Fingerprint script](src/egfrank/fingerprints.py) · [Saved fingerprints](data/processed/egfr_fingerprints.npz)

## Step 9. Reserve chemical scaffolds for testing

Two compounds can have different ChEMBL IDs while sharing the same central ring structure. Placing close relatives on opposite sides of a random split can make a test easier because the model has already seen examples from the same chemical family. I used RDKit to calculate a Murcko scaffold, the framework formed by a molecule's rings and the chains connecting them. I kept all compounds with the same scaffold together and assigned molecules without rings to one group called `ACYCLIC`.

The script tried 100 repeatable orderings of the scaffold groups, starting with random seed 42. It chose the assignment closest to reserving 20% of each activity class for testing, with both classes represented on each side. I checked the saved assignments directly.

| Set | Active | Inactive | Total compounds | Scaffold groups |
| --- | ---: | ---: | ---: | ---: |
| Training | 1,117 | 163 | 1,280 | 435 |
| Test | 279 | 41 | 320 | 105 |

All 1,600 compound IDs appear exactly once, and the two sets share zero scaffold groups. The 320 test compounds therefore provide a held-out set for evaluating predictions on other chemical frameworks. These counts describe the dataset split. [Split script](src/egfrank/split_data.py) · [Saved assignments](data/processed/egfr_scaffold_split.csv)

## Step 10. Train the first Random Forest

I matched each fingerprint to its activity label and scaffold assignment using its ChEMBL compound ID. The training script checks for duplicate IDs, missing matches, unexpected fingerprint dimensions, and labels that disagree between files. I then selected the 1,280 training compounds. Each example supplied 2,048 fingerprint bits as inputs and one experimental label, active or inactive, as the answer.

I fitted a Random Forest with 100 decision trees and a random seed of 42. Each tree learns a series of decisions based on structural features. The forest averages the trees’ class probability estimates to score each compound. A higher active-class score expresses greater model confidence that the compound meets the project’s 1,000 nM activity cutoff. The classifier predicts a category rather than an exact IC50 concentration. I saved the fitted model to `models/egfr_random_forest.joblib` so I could evaluate it without retraining. [Training script](src/egfrank/train.py)

## Step 11. Benchmark the model on held-out scaffolds

I evaluated the saved model on 320 compounds from 105 scaffolds held out during training. The reference labels came from the same experimental cutoff used for training. Because 279 of the 320 test compounds were active, a simple baseline could achieve 87.2% accuracy by calling every compound active. I evaluated both methods on the same compounds.

| Method | Correct predictions | Accuracy | ROC-AUC | Average precision for the active class |
| --- | ---: | ---: | ---: | ---: |
| Always predict active | 279 / 320 | 87.2% | 0.500 | 0.872 |
| Random Forest | 306 / 320 | 95.6% | 0.929 | 0.981 |

The confusion matrix shows where the forest made its 14 errors. Rows give the experimental label and columns give the model prediction.

| Experimental label | Predicted inactive | Predicted active |
| --- | ---: | ---: |
| Inactive | 29 | 12 |
| Active | 2 | 277 |

Of the 289 compounds the forest called active, 277 met the experimental cutoff. That fraction is **active precision** (95.8%). Of the 279 experimentally active compounds, it found 277. That fraction is **active recall** (99.3%). The active-class F1 score, which combines precision and recall, was 0.98. The smaller inactive class was harder to identify. The forest recognized 29 of 41 inactive compounds (70.7% recall), compared with zero for the baseline.

For compound prioritization, the order of the scores also matters. ROC-AUC measures how often an active compound ranks above an inactive one, counting ties as half credit. Average precision summarizes how concentrated the active compounds are toward the top as progressively more compounds are selected. Both figures describe the full ranking on this single scaffold split. The experimental endpoint is EGFR inhibition under the selected assay conditions. Assessing response in an NSCLC patient would require clinical evidence. [Evaluation script](src/egfrank/evaluate.py)

## Step 12. Export predictions and build an interactive compound map

I exported the saved model's activity scores for all 320 held-out compounds to a JSON file that the React app can read. Each record contains its ChEMBL ID, experimental label, minimum and maximum eligible IC50, model score, map coordinates, and neighboring compound indices. The export connects the trained model to the interface. [Export script](src/egfrank/export_map.py) · [Map dataset](frontend/public/data/compound_map.json)

To calculate chemical similarity, I compared the compounds' Morgan fingerprints using Tanimoto similarity. This divides the number of shared set bits by the number set in either fingerprint. I converted similarity to distance by subtracting it from one, then used t-SNE to place the compounds on a two-dimensional map. t-SNE compresses a high-dimensional representation while emphasizing local neighborhoods. Its plotted coordinates support exploration, while the fingerprint comparison supplies the actual similarity measurement. The export uses random seed 42 and perplexity 15. Edges connect each compound to its four closest fingerprint neighbors.

The React prototype supports searching by ChEMBL ID, selecting a compound through the map or ranked list, and viewing its measured IC50 range and activity label. A control highlights twelve compounds to illustrate the fictional laboratory's testing budget. Scores are sorted from highest to lowest, with ChEMBL ID used to order equal scores. [Interactive map](frontend/src/CompoundMap.tsx)

### Interpret tied scores

The export contains 164 compounds with the maximum model score of 1.00. Their relative positions in the list follow the ID ordering. That ordering gives a repeatable display choice, while the model assigns the same score to every member of this group. Consequently, the highlighted twelve represent one selection from a much larger tied group. A score of 1.00 describes the forest's output for the active class under this model and training setup. Experimental confirmation remains necessary.

### Prototype verification and availability

The frontend production build passed after the interactive map was added. The local browser connection subsequently failed, so a standalone HTML snapshot was created as a demonstration artifact. The November 18 documentation check looked for that snapshot in the project root and recorded it as absent. It had been moved to `demo/EGFRank-Demo.html` three days earlier, so that entry describes where the check looked rather than a file the project had lost.

I rechecked the prototype on November 25. The snapshot carries both datasets inline as base64 `data:application/json` URIs instead of requesting them, which is why it opens with no server. Decoding them returned 320 compound-map records and 100 pilot measurements, and their SHA-256 digests match `frontend/public/data/compound_map.json` and `frontend/public/data/egfr_sample.json`. The snapshot therefore presents the same data the pipeline exported. The linter reported no findings, and the production build compiled 18 modules. A development server returned the application shell, the module entry point, and both data files, and the served map parsed to 320 records.

These checks establish the production build and the data path from export to browser request. They do not exercise rendering or interaction, so the visual behavior of the map and the ranked list rests on manual use.

## Step 13. Fix the first variant prediction target

I fixed the first variant-specific target as biochemical inhibition of human EGFR carrying L858R as its sole named alteration. One example will represent one compound and its eligible exact IC50 evidence against that protein. L858R/T790M, L858R/C797S, wild-type EGFR, and other signatures remain separate targets because their measurements describe different protein forms. Cellular phosphorylation, proliferation, viability, binding-affinity, and degradation endpoints are outside this first target.

The existing project cutoff remains unchanged. A compound is active when every eligible single-L858R IC50 is at or below 1,000 nM, inactive when every eligible value is above 1,000 nM, and conflicting when measurements fall on both sides. Exactly 1,000 nM is active. The current automated screen yields 246 active and 48 inactive single-L858R compound pairs, but these counts are provisional until each candidate L858R assay has an explicit reviewed decision.

The model's role is to rank compounds for a later laboratory shortlist. Its held-out labeled compounds will measure retrospective performance; they are not a prospective candidate library. The final twelve-compound selection must come from a separate collection whose outcomes were not used for training or evaluation. [L858R target specification](docs/research/l858r-target-spec.md) · [Current mutant screen](src/egfrank/review_mutant_assays.py)

## Step 14. Build the L858R assay-review pipeline

I added a target-specific command that recognizes common L858R mutation notation, keeps mutation combinations separate, classifies assay endpoints, records exclusion reasons, and applies the existing activity-quality checks. Its assay audit preserves explicit manual decisions and review notes between runs. A strict review gate stops the workflow while any automated candidate remains pending. The original no-argument command continues to reproduce the earlier all-variant exports and counts.

The first run found 473 assay descriptions naming L858R, of which 47 were automated candidates for the single-L858R biochemical endpoint. I reviewed all 47 against the fixed target. Each description assigns a biochemical enzyme readout to L858R without naming wild type, another mutation, or a cellular endpoint; the corresponding activity rows consistently identify human EGFR target CHEMBL203, assay type B, and a source document. I included all 47 and recorded the decisions separately from the generated audit. The strict review gate now passes with zero pending candidates.

The reviewed assays contain 470 measurements before activity-level screening. The quality rules retain 343 measurements and produce 329 consistent compound labels: 281 active and 48 inactive, with no cross-cutoff conflicts. Of these labels, 321 rely on one eligible measurement. These are the final counts for the frozen snapshot under the current target and screening rules; their reliance on mostly single measurements remains a limitation for later evaluation. [L858R assay decisions](docs/research/l858r-assay-decisions.csv) · [Generated assay audit](docs/research/l858r-assay-review.csv) · [Screening and labeling code](src/egfrank/review_mutant_assays.py) · [Target specification](docs/research/l858r-target-spec.md)

## Step 15. Validate structures and chemical coverage

I validated the reviewed dataset before creating any benchmark split. All 329 compound IDs are unique, every SMILES parses with RDKit, the stored structures already match RDKit's canonical isomeric SMILES, and no two IDs resolve to the same canonical structure. Activity IDs are unique across compounds, measurement counts match their provenance IDs, every contributing assay has an include decision, evidence grades match their measurement and source counts, and every active or inactive label agrees with the 1,000 nM cutoff. The dataset contains 281 active and 48 inactive compounds, so an always-active baseline would achieve 85.4% accuracy.

The 329 structures occupy 161 Murcko scaffold groups. Of those, 142 contain an active compound, 27 contain an inactive compound, eight contain both classes, and 109 contain only one compound. The largest scaffold group contains 18 compounds. This provides enough scaffold diversity to attempt a held-out split, but the inactive class has limited coverage and must be evaluated with balanced and inactive-class metrics.

I generated 2,048-bit radius-2 Morgan fingerprints with chirality enabled. All rows contain set bits, but seven groups of distinct structures produce identical fingerprint rows. Six such groups span different Murcko scaffolds, and one contains both an active and an inactive compound. The scaffold inventory therefore carries a stable fingerprint-group identifier. A future split must keep connected scaffold and identical-fingerprint groups together so the same model input cannot appear in both training and testing. [Dataset validation report](docs/research/l858r-dataset-validation.md) · [Validation code](src/egfrank/validate_l858r_dataset.py) · [Fingerprint code](src/egfrank/fingerprints.py)

## Step 16. Freeze a leakage-safe L858R benchmark split

I connected compounds whenever they share a Murcko scaffold or an identical Morgan fingerprint. This transitive grouping prevents a fingerprint-equivalent pair from crossing the boundary even when its structures have different scaffolds. The 329 compounds form 156 connected groups. With random seed 42, the split search evaluated 25,000 reproducible whole-group assignments and selected one matching the rounded 20% targets for both activity classes.

The frozen training partition contains 263 compounds: 225 active and 38 inactive across 130 scaffolds. The held-out test partition contains 66 compounds: 56 active and 10 inactive across 31 scaffolds. It therefore holds 20.1% of all compounds, 19.9% of active compounds, and 20.8% of inactive compounds. Compound IDs, canonical structures, scaffolds, identical-fingerprint groups, and connected split groups each have zero overlap between training and testing. The fingerprint-equivalent active/inactive pair remains together in training.

This split is now fixed for the first L858R benchmark. Model and hyperparameter selection must use the training partition with training-only grouped validation. The test set will be evaluated only after those decisions are frozen. Because it contains 10 inactive compounds, one inactive error changes test inactive recall by 10 percentage points; raw counts must accompany every final metric. [Split validation report](docs/research/l858r-split-validation.md) · [Split code](src/egfrank/split_data.py)

## Reproduce the work in PowerShell

From the repository root, run the following with the project environment installed. The first command refreshes the demo. The second rebuilds the labels from the dated snapshot above. The final two start the website. To download a **new** dated snapshot, run the downloader separately with the same Python interpreter. A new download may reflect a later ChEMBL release.

```powershell
.\.venv\Scripts\python.exe src/egfrank/preprocess.py
.\.venv\Scripts\python.exe src/egfrank/review_assays.py
.\.venv\Scripts\python.exe src/egfrank/fingerprints.py
.\.venv\Scripts\python.exe src/egfrank/split_data.py
.\.venv\Scripts\python.exe src/egfrank/train.py
.\.venv\Scripts\python.exe src/egfrank/evaluate.py
.\.venv\Scripts\python.exe src/egfrank/export_map.py
cd frontend
npm.cmd run dev
```

## Resources

- Experimental data and metadata: [ChEMBL EGFR target](https://www.ebi.ac.uk/chembl/explore/target/CHEMBL203) · [ChEMBL API documentation](https://www.ebi.ac.uk/chembl/api/data/docs)
- Table processing, logarithms, and structures: [pandas](https://pandas.pydata.org/docs/) · [NumPy](https://numpy.org/doc/stable/) · [RDKit](https://www.rdkit.org/docs/GettingStartedInPython.html)
- Website tools: [React](https://react.dev/learn) · [TypeScript](https://www.typescriptlang.org/docs/) · [Vite](https://vite.dev/guide/)
