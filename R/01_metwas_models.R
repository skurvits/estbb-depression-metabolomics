#!/usr/bin/env Rscript

# Metabolome-wide association study (MetWAS)
# Final manuscript specification: 14 symptoms × 249 metabolites, Models 1–3.
#
# Input is an analysis-ready EstBB dataset. Nightingale QC/normalization and
# sex-specific standardization are treated as upstream data-preparation steps.

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(purrr)
  library(tibble)
})

args <- commandArgs(trailingOnly = TRUE)
model_arg <- grep("^--model$", args)
if (length(model_arg) != 1L || model_arg == length(args)) {
  stop("Usage: Rscript R/01_metwas_models.R --model {1|2|3}")
}
model_id <- as.integer(args[model_arg + 1L])
if (!model_id %in% 1:3) stop("--model must be 1, 2, or 3")

input_path <- Sys.getenv("METWAS_DATA", unset = "data/private/analysis_ready_data.csv")
out_dir <- Sys.getenv("MODEL_OUT_DIR", unset = "outputs/model_outputs")
delim <- Sys.getenv("METWAS_DELIM", unset = ";")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

primary_symptoms <- c(
  "deprDepressedEver", "deprLossOfInterestEver", "deprHopelessness",
  "deprWorthlessness", "deprConfidence", "deprAttention",
  "deprThoughtsOfDeath", "deprTired", "deprSleepGain", "deprSleepLoss",
  "deprAppetite1", "deprGainWeight", "deprLostWeight", "deprSleep1"
)

cov_model1 <- c("age_at_sample", "Person.gender.code", "deltaTime")
cov_model2 <- c(
  cov_model1,
  "PersonPortrait.lastSmokingStatus.name", "AlcoholConsumptionGroup",
  "PersonPortrait.educationGroup.code", "A10_before_blood",
  "C02_before_blood", "C10_before_blood"
)
cov_model3 <- c(cov_model2, "PersonPortrait.lastBmi")

covariates <- switch(as.character(model_id),
  "1" = cov_model1,
  "2" = cov_model2,
  "3" = cov_model3
)

if (!file.exists(input_path)) stop("Input not found: ", input_path)
dat <- readr::read_delim(input_path, delim = delim, show_col_types = FALSE)

# The analysis-ready file keeps the 249 Nightingale measures contiguously.
start_col <- match("Total_C", names(dat))
end_col <- match("S_HDL_TG_pct", names(dat))
if (is.na(start_col) || is.na(end_col) || start_col > end_col) {
  stop("Could not identify metabolite span Total_C ... S_HDL_TG_pct")
}
metabolites <- names(dat)[start_col:end_col]
if (length(metabolites) != 249L) {
  stop("Expected 249 metabolites; found ", length(metabolites))
}

required <- unique(c(primary_symptoms, metabolites, covariates))
missing <- setdiff(required, names(dat))
if (length(missing) > 0) stop("Missing columns: ", paste(missing, collapse = ", "))

# Preserve categorical coding from the analysis dataset.
dat[["Person.gender.code"]] <- as.factor(dat[["Person.gender.code"]])
for (v in intersect(
  c("PersonPortrait.lastSmokingStatus.name", "AlcoholConsumptionGroup",
    "PersonPortrait.educationGroup.code"), names(dat))) {
  dat[[v]] <- as.factor(dat[[v]])
}

fit_one <- function(outcome, metabolite) {
  vars <- unique(c(outcome, metabolite, covariates))
  d <- dat[, vars, drop = FALSE]
  d <- d[stats::complete.cases(d), , drop = FALSE]

  y <- d[[outcome]]
  if (is.factor(y)) y <- as.character(y)
  y <- suppressWarnings(as.numeric(y))
  if (!all(sort(unique(y[!is.na(y)])) %in% c(0, 1))) {
    stop("Outcome is not coded 0/1: ", outcome)
  }
  d[[outcome]] <- y

  form <- stats::reformulate(c(metabolite, covariates), response = outcome)
  fit <- stats::glm(form, data = d, family = stats::binomial(link = "logit"))
  sm <- summary(fit)$coefficients

  if (!metabolite %in% rownames(sm)) {
    return(tibble(
      variable = outcome, feature = metabolite, p_value = NA_real_,
      odds_ratio = NA_real_, CI_low = NA_real_, CI_high = NA_real_,
      N = nrow(d), N_cases = sum(y == 1), N_controls = sum(y == 0)
    ))
  }

  beta <- stats::coef(fit)[[metabolite]]
  # Match the historical workflow: profile-likelihood CI when available.
  ci_beta <- tryCatch(
    suppressMessages(stats::confint(fit, parm = metabolite)),
    error = function(e) beta + c(-1, 1) * 1.96 * sm[metabolite, "Std. Error"]
  )

  tibble(
    variable = outcome,
    feature = metabolite,
    p_value = sm[metabolite, "Pr(>|z|)"],
    odds_ratio = exp(beta),
    CI_low = exp(ci_beta[1]),
    CI_high = exp(ci_beta[2]),
    N = nrow(d),
    N_cases = sum(y == 1),
    N_controls = sum(y == 0)
  )
}

message("Running Model ", model_id, " with ", length(covariates), " covariates")
message(length(primary_symptoms), " symptoms × ", length(metabolites), " metabolites")

results <- purrr::map_dfr(primary_symptoms, function(symptom) {
  message("  ", symptom)
  purrr::map_dfr(metabolites, ~ fit_one(symptom, .x))
})

bonf_alpha <- 0.05 / (length(primary_symptoms) * length(metabolites))
results <- results %>%
  mutate(
    adjusted_p_value = p.adjust(p_value, method = "bonferroni"),
    sig_bonf_249x14 = !is.na(p_value) & p_value < bonf_alpha
  ) %>%
  arrange(variable, p_value)

out_file <- file.path(out_dir, paste0("Model", model_id, ".csv"))
readr::write_csv(results, out_file)
message("Wrote: ", out_file)
message("Bonferroni-significant pairs: ", sum(results$sig_bonf_249x14, na.rm = TRUE))
