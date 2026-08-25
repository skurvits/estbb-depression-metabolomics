# Manuscript-to-code map

| Manuscript analysis | Main script | Key output |
|---|---|---|
| MetWAS Models 1–3 | `R/01_metwas_models.R` | Model 1/2/3 coefficient tables |
| Significant-association counts | `R/02_manuscript_analysis.R` | Figure 1A; verified count table |
| Model 2/3 heatmaps | `R/02_manuscript_analysis.R` | Figure 1B–C |
| BMI attenuation | `R/02_manuscript_analysis.R` | Figure 2; Supplementary Table S5 |
| BMI candidate-set construction | `R/02_manuscript_analysis.R` | 691 unique pairs |
| Metabolite × BMI interactions | `R/02_manuscript_analysis.R` | Figure 3; Supplementary Table S6 |
| SOM training | `R/03_som_training.R` | Trained SOM RDS |
| K=6 SOM profiles and QC | `R/04_som_clustering_qc.R` | SOM assignments; validation metrics |
| Conditional metabolite × SOM-profile interactions | `R/02_manuscript_analysis.R` | Figure 4; Supplementary Table S9 |

## Files intentionally omitted

The original analysis was developed iteratively. The public repository does **not** include exploratory scripts that are not part of the final reported workflow.
