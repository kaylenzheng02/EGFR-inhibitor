# EGFRank project history

**Research question.** Can the chemical structure of a compound help us prioritize it for an experiment measuring inhibition of human EGFR kinase?

**Abstract.** I used previously reported experiments to construct a dataset for answering this question. A pilot of 100 measurements exposed a central problem: the same compound can give different results when scientists measure different biological events. I traced those measurements to their experiments, built a small explorer, and applied an assay selection and quality screen to a larger download. From 26,600 records, I obtained 1,600 compounds with consistent labels under a project cutoff of 1,000 nM. The completed work establishes the dataset and its provenance.

## Step 1. Research question and project setup

Can chemical structure help prioritize compounds for testing against human epidermal growth factor receptor (EGFR) kinase? I chose general EGFR inhibition for the first model. I organized a Git project in VS Code with a Python environment for processing data and a React and TypeScript frontend for exploring the results. I consolidated the original nested folders so commands run from one repository root. · [Project setup](README.md)

### Biological background

Epidermal growth factor receptor (EGFR) is a protein that spans the membrane surrounding a cell. A growth factor is a signaling molecule. When it binds the outside portion of EGFR, receptor molecules can pair up and activate the part inside the cell. This part is an **enzyme**, a molecule that speeds up a chemical reaction. EGFR belongs to the **tyrosine kinase** family. A kinase transfers a phosphate group, typically from a molecule called ATP, onto a receiving molecule called a substrate. Tyrosine is an amino acid, one of the building blocks of proteins. EGFR can attach phosphate to particular tyrosines. The transfer is called **phosphorylation**. [UniProt record for human EGFR](https://www.uniprot.org/uniprotkb/P00533/entry)

In an activated receptor pair, one receptor can phosphorylate tyrosines on its partner. This is called **autophosphorylation** because the receptor protein is itself the recipient. Other proteins can recognize these sites and pass the signal onward, influencing processes such as cell growth. This sequence gives scientists several things to measure: the enzyme transferring phosphate to a test substrate, the receptor becoming phosphorylated inside a cell, or the eventual growth of a population of cells. Each is a different observation. [EGFR function and signaling](https://www.uniprot.org/uniprotkb/P00533/entry)

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

## Reproduce the work in PowerShell

From the repository root, run the following with the project environment installed. The first command refreshes the demo. The second rebuilds the labels from the dated snapshot above. The final two start the website. To download a **new** dated snapshot, run the downloader separately with the same Python interpreter. A new download may reflect a later ChEMBL release.

```powershell
.\.venv\Scripts\python.exe src/egfrank/preprocess.py
.\.venv\Scripts\python.exe src/egfrank/review_assays.py
cd frontend
npm.cmd run dev
```

## Resources

- Experimental data and metadata: [ChEMBL EGFR target](https://www.ebi.ac.uk/chembl/explore/target/CHEMBL203) · [ChEMBL API documentation](https://www.ebi.ac.uk/chembl/api/data/docs)
- Table processing, logarithms, and structures: [pandas](https://pandas.pydata.org/docs/) · [NumPy](https://numpy.org/doc/stable/) · [RDKit](https://www.rdkit.org/docs/GettingStartedInPython.html)
- Website tools: [React](https://react.dev/learn) · [TypeScript](https://www.typescriptlang.org/docs/) · [Vite](https://vite.dev/guide/)
