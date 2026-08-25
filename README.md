# Depression symptoms × metabolomics in the Estonian Biobank

Cleaned analysis code accompanying the manuscript **“Metabolic profile of 14 lifetime depression symptoms and body mass index in the Estonian Biobank.”**

This repository contains the minimal code needed to reproduce the *reported analysis logic* from analysis-ready Estonian Biobank data. It intentionally omits exploratory scripts, abandoned SOM configurations, one-off plotting experiments, and participant-level data.

## Analysis overview

The analysis follows the manuscript in four steps:

1. **Metabolome-wide association study (MetWAS)** across 14 lifetime depressive symptoms and 249 Nightingale NMR metabolites using three progressively adjusted logistic regression models.
2. **BMI attenuation and effect modification**, including Model 2 → Model 3 coefficient attenuation and metabolite × BMI interaction tests among the 691 metabolite–symptom pairs significant in Model 2 or Model 3.
3. **Systemic metabolic profiles**, using a 21 × 21 SOM followed by Ward.D2 hierarchical clustering into six profiles.
4. **Conditional SOM follow-up**, testing metabolite × profile heterogeneity among the 211 pairs with FDR-significant metabolite × BMI interactions while retaining the continuous metabolite × BMI interaction in the model.

The final reported progression is:

- 1,367 Bonferroni-significant pairs in Model 1
- 660 in Model 2
- 136 in Model 3
- 105 of the 660 Model 2 associations retained in Model 3
- 691 unique Model 2/3 candidate pairs for BMI interaction testing
- 211/691 FDR-significant metabolite × BMI interactions
- 3/211 FDR-significant additional metabolite × SOM-profile interactions

## Repository structure

```text
R/
  01_metwas_models.R          Fit the three MetWAS models
  02_manuscript_analysis.R    Final counts, attenuation, BMI interactions,
                              figures, supplementary tables, SOM follow-up
  03_som_training.R           Clean SOM training implementation
  04_som_clustering_qc.R      Ward.D2 clustering, K=6 profiles, QE/TE and
                              cluster-validation metrics
results/
  manuscript_numbers_verified.txt
  significant_counts_by_symptom.csv
  model2_to_model3_retention.csv
  bmi_interactions_fdr_significant.csv
  som_interactions_fdr_significant.csv
data/
  README.md                   Required input structure and data-access note
docs/
  ANALYSIS_MAP.md             Mapping from manuscript analyses to code
CHECK_BEFORE_PUBLICATION.md   Final provenance checks before making repo public
```

## Data availability

Individual-level Estonian Biobank data cannot be distributed through GitHub. The scripts expect an **analysis-ready individual-level file** containing the 249 metabolite measures, 14 symptom phenotypes, BMI, covariates, and an individual identifier. Access to Estonian Biobank data is subject to the Biobank's data-access procedures and applicable data-protection requirements.

The `results/` directory contains only aggregate, non-identifying outputs used to verify that the code corresponds to the manuscript.

## Running the analysis

### 1. Fit MetWAS Models 1–3

`01_metwas_models.R` expects the prepared analysis dataset. It can be run separately for each model:

```bash
Rscript R/01_metwas_models.R --model 1
Rscript R/01_metwas_models.R --model 2
Rscript R/01_metwas_models.R --model 3
```

Set the input/output locations with environment variables, for example:

```bash
export METWAS_DATA=/path/to/analysis_ready_data.csv
export MODEL_OUT_DIR=/path/to/model_outputs
```

### 2. Reproduce manuscript outputs

```bash
export M1_PATH=/path/to/model_outputs/Model1.csv
export M2_PATH=/path/to/model_outputs/Model2.csv
export M3_PATH=/path/to/model_outputs/Model3.csv
export METWAS_INDIVIDUAL_PATH=/path/to/analysis_ready_data.csv
export SOM_ASSIGNMENT_PATH=/path/to/som_assignments.csv
Rscript R/02_manuscript_analysis.R
```

This script reconstructs the final manuscript counts and Tables S2–S6/S9, reruns the 691 BMI interaction tests, applies global BH FDR, performs the corrected 211-pair SOM follow-up, and generates the manuscript figures.

### 3. SOM analysis

The cleaned SOM scripts are separated into training and clustering/QC:

```bash
export METWAS_DATA=/path/to/analysis_ready_data.csv
Rscript R/03_som_training.R

export SOM_MODEL_PATH=/path/to/som_model.rds
Rscript R/04_som_clustering_qc.R
```


## Statistical conventions

- Metabolites: 249 Nightingale NMR measures.
- Outcomes: 14 primary lifetime depressive symptoms.
- Bonferroni threshold: `0.05 / (249 × 14) = 1.4343087e-05`.
- Model 2 → Model 3 attenuation is calculated on the absolute log-odds scale.
- BMI interaction tests are restricted to the union of Model 2- or Model 3-Bonferroni-significant pairs (691 tests) and corrected jointly using Benjamini–Hochberg FDR.
- SOM-profile interaction follow-up is restricted to BMI-interaction FDR-significant pairs (211 tests) and corrected jointly using Benjamini–Hochberg FDR.
- In the SOM follow-up, the continuous metabolite × BMI term is retained in both nested models, so the likelihood-ratio test evaluates additional profile-dependent heterogeneity beyond BMI moderation.

## Software

The manuscript analysis was performed in R. Core packages used by the cleaned scripts include `dplyr`, `tidyr`, `readr`, `purrr`, `broom`, `ggplot2`, `patchwork`, `scales`, `kohonen`, `cluster`, and `clusterCrit`.


## License

Code is released under the MIT License. No license is granted for Estonian Biobank data.
