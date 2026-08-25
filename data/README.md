# Data inputs

Participant-level Estonian Biobank data are **not included** in this repository.

The cleaned code expects an analysis-ready file with:

- 249 Nightingale NMR metabolites, already quality-controlled and standardized as used in the manuscript;
- 14 binary depressive symptom variables;
- BMI at Mental Health Online Survey completion;
- Model covariates;
- an individual identifier used only to merge SOM assignments back to the analysis data.

## Primary symptom variables

- `deprDepressedEver`
- `deprLossOfInterestEver`
- `deprHopelessness`
- `deprWorthlessness`
- `deprConfidence`
- `deprAttention`
- `deprThoughtsOfDeath`
- `deprTired`
- `deprSleepGain`
- `deprSleepLoss`
- `deprAppetite1`
- `deprGainWeight`
- `deprLostWeight`
- `deprSleep1`

## Covariates

Model 1:
- `age_at_sample`
- `Person.gender.code`
- `deltaTime`

Model 2 additionally includes:
- `PersonPortrait.lastSmokingStatus.name`
- `AlcoholConsumptionGroup`
- `PersonPortrait.educationGroup.code`
- `A10_before_blood`
- `C02_before_blood`
- `C10_before_blood`

Model 3 additionally includes:
- `PersonPortrait.lastBmi`

If column names differ in the controlled EstBB release, adapt the configuration block at the top of the relevant script rather than changing the statistical logic.
