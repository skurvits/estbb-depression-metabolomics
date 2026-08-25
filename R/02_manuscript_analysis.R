#!/usr/bin/env Rscript

# Final downstream analysis corresponding to the submitted manuscript.
# Builds Figures 1–4 and Supplementary Tables S2–S6/S9 from the final
# MetWAS outputs, analysis-ready individual data, and six-profile SOM assignment.

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(ggplot2)
  library(stringr)
  library(forcats)
  library(purrr)
  library(broom)
  library(patchwork)
  library(scales)
})

use_ggrepel <- requireNamespace("ggrepel", quietly = TRUE)

m1_path <- Sys.getenv("M1_PATH", unset = "outputs/model_outputs/Model1.csv")
m2_path <- Sys.getenv("M2_PATH", unset = "outputs/model_outputs/Model2.csv")
m3_path <- Sys.getenv("M3_PATH", unset = "outputs/model_outputs/Model3.csv")
individual_path <- Sys.getenv("METWAS_INDIVIDUAL_PATH", unset = "data/private/analysis_ready_data.csv")
som_path <- Sys.getenv("SOM_ASSIGNMENT_PATH", unset = "outputs/som/som_assignments_k6.csv")
fig_dir <- Sys.getenv("FIGURE_OUT_DIR", unset = "outputs/figures")
table_dir <- Sys.getenv("TABLE_OUT_DIR", unset = "outputs/tables")

dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

# Read model outputs here so the remainder of the script is the audited final logic.
stopifnot(file.exists(m1_path), file.exists(m2_path), file.exists(m3_path))
m1_raw <- read_csv(m1_path, show_col_types = FALSE)
m2_raw <- read_csv(m2_path, show_col_types = FALSE)
m3_raw <- read_csv(m3_path, show_col_types = FALSE)

# ======================================================================
# 1. PRIMARY ANALYSIS DEFINITIONS
# ======================================================================

alpha <- 0.05

primary_symptoms <- c(
  "deprDepressedEver",
  "deprLossOfInterestEver",
  "deprHopelessness",
  "deprWorthlessness",
  "deprConfidence",
  "deprAttention",
  "deprThoughtsOfDeath",
  "deprTired",
  "deprSleepGain",
  "deprSleepLoss",
  "deprAppetite1",
  "deprGainWeight",
  "deprLostWeight",
  "deprSleep1"
)

symptom_labels <- c(
  deprDepressedEver      = "Depressed mood*",
  deprLossOfInterestEver = "Anhedonia*",
  deprHopelessness       = "Hopelessness",
  deprWorthlessness      = "Worthlessness",
  deprConfidence         = "Low self-confidence",
  deprAttention          = "Attention problems",
  deprThoughtsOfDeath    = "Thoughts of death",
  deprTired              = "Fatigue",
  deprSleepGain          = "Hypersomnia",
  deprSleepLoss          = "Insomnia",
  deprAppetite1          = "Appetite or weight change",
  deprGainWeight         = "Weight gain",
  deprLostWeight         = "Weight loss",
  deprSleep1             = "Sleep problems (any)"
)

model_cols <- c(
  "Model 1" = "#1F77B4",
  "Model 2" = "#FF7F0E",
  "Model 3" = "#9467BD"
)

som_cols <- c(
  "C1" = "#0072B2",
  "C2" = "#E69F00",
  "C3" = "#009E73",
  "C4" = "#D55E00",
  "C5" = "#CC79A7",
  "C6" = "#000000"
)


# ======================================================================
# 2. READ AND VALIDATE MODEL OUTPUTS
# ======================================================================

required_model_cols <- c(
  "variable", "feature", "p_value", "odds_ratio",
  "CI_low", "CI_high", "N", "N_cases", "N_controls"
)

check_required_cols <- function(df, model_name) {
  missing <- setdiff(required_model_cols, names(df))
  if (length(missing) > 0) {
    stop(model_name, " is missing columns: ", paste(missing, collapse = ", "))
  }
}
check_required_cols(m1_raw, "Model 1")
check_required_cols(m2_raw, "Model 2")
check_required_cols(m3_raw, "Model 3")

# MODEL 1 is the source of truth for the metabolite universe:
# the final analysis uses exactly 249 metabolite features; Model 1 is used
# as the source of truth so covariate coefficient rows can never enter counts.
metabolite_features <- unique(as.character(m1_raw$feature))

if (length(metabolite_features) != 249L) {
  stop(
    "Expected exactly 249 metabolite features in Model 1, found ",
    length(metabolite_features),
    ". Do not continue until this is resolved."
  )
}

bonf_alpha <- alpha / (length(primary_symptoms) * length(metabolite_features))
message("Global Bonferroni threshold = ", format(bonf_alpha, scientific = TRUE, digits = 8))

prep_primary_model <- function(df, model_label) {
  out <- df %>%
    filter(
      variable %in% primary_symptoms,
      feature %in% metabolite_features
    ) %>%
    mutate(
      variable = as.character(variable),
      feature = as.character(feature)
    )
  
  duplicates <- out %>%
    count(variable, feature) %>%
    filter(n != 1L)
  
  if (nrow(duplicates) > 0) {
    stop(model_label, " contains duplicate symptom–metabolite rows.")
  }
  
  expected_n <- length(primary_symptoms) * length(metabolite_features)
  if (nrow(out) != expected_n) {
    stop(
      model_label, ": expected ", expected_n,
      " primary symptom–metabolite rows, found ", nrow(out)
    )
  }
  
  per_symptom <- out %>% count(variable)
  if (any(per_symptom$n != length(metabolite_features))) {
    stop(model_label, ": at least one primary symptom does not have 249 metabolites.")
  }
  
  out %>%
    mutate(
      model = model_label,
      beta = log(odds_ratio),
      sig_bonf = !is.na(p_value) & p_value < bonf_alpha
    )
}

m1 <- prep_primary_model(m1_raw, "Model 1")
m2 <- prep_primary_model(m2_raw, "Model 2")
m3 <- prep_primary_model(m3_raw, "Model 3")


# ======================================================================
# 3. VERIFIED COUNTS AND MANUSCRIPT SUMMARY
# ======================================================================

count_sig <- function(df) {
  df %>%
    group_by(variable) %>%
    summarise(n_sig = sum(sig_bonf, na.rm = TRUE), .groups = "drop") %>%
    mutate(model = unique(df$model))
}

counts_long <- bind_rows(count_sig(m1), count_sig(m2), count_sig(m3)) %>%
  mutate(
    symptom = unname(symptom_labels[variable]),
    model = factor(model, levels = c("Model 1", "Model 2", "Model 3"))
  )

model_totals <- counts_long %>%
  group_by(model) %>%
  summarise(n_sig = sum(n_sig), .groups = "drop")

print(model_totals)

# These are the verified totals from the final attached model files.
expected_totals <- c("Model 1" = 1367L, "Model 2" = 660L, "Model 3" = 136L)
observed_totals <- setNames(model_totals$n_sig, as.character(model_totals$model))

if (!all(observed_totals[names(expected_totals)] == expected_totals)) {
  stop(
    "Primary result totals do not match the verified final-data values.\n",
    "Expected: M1=1367, M2=660, M3=136.\n",
    "Observed: ",
    paste(names(observed_totals), observed_totals, collapse = "; "),
    "\nCheck that the final model files are being used."
  )
}

write_csv(
  counts_long %>% select(model, variable, symptom, n_sig),
  file.path(table_dir, "Verified_significant_counts_by_symptom_models123.csv")
)


# ======================================================================
# 4. MODEL 2 -> MODEL 3 PAIRWISE COMPARISON
#    Critical distinction: retained != total significant in Model 3
# ======================================================================

cmp <- m2 %>%
  select(
    variable, feature,
    N_cases_M2 = N_cases, N_controls_M2 = N_controls,
    beta2 = beta, OR2 = odds_ratio,
    CI2_low = CI_low, CI2_high = CI_high,
    p2 = p_value, sig_m2 = sig_bonf
  ) %>%
  inner_join(
    m3 %>%
      select(
        variable, feature,
        N_cases_M3 = N_cases, N_controls_M3 = N_controls,
        beta3 = beta, OR3 = odds_ratio,
        CI3_low = CI_low, CI3_high = CI_high,
        p3 = p_value, sig_m3 = sig_bonf
      ),
    by = c("variable", "feature")
  ) %>%
  mutate(
    attenuation = ifelse(
      is.finite(beta2) & abs(beta2) > 1e-12,
      1 - abs(beta3) / abs(beta2),
      NA_real_
    ),
    same_direction = case_when(
      is.na(beta2) | is.na(beta3) ~ NA,
      beta2 == 0 | beta3 == 0 ~ NA,
      sign(beta2) == sign(beta3) ~ TRUE,
      TRUE ~ FALSE
    ),
    retained = sig_m2 & sig_m3,
    new_in_m3 = !sig_m2 & sig_m3,
    lost_after_bmi = sig_m2 & !sig_m3,
    symptom = unname(symptom_labels[variable])
  )

retention_summary <- cmp %>%
  group_by(variable, symptom) %>%
  summarise(
    n_sig_M2 = sum(sig_m2, na.rm = TRUE),
    n_sig_M3_total = sum(sig_m3, na.rm = TRUE),
    n_retained_M2_to_M3 = sum(retained, na.rm = TRUE),
    retained_pct = ifelse(
      n_sig_M2 > 0,
      100 * n_retained_M2_to_M3 / n_sig_M2,
      NA_real_
    ),
    n_new_in_M3 = sum(new_in_m3, na.rm = TRUE),
    n_lost_after_BMI = sum(lost_after_bmi, na.rm = TRUE),
    median_attenuation_M2sig = ifelse(
      n_sig_M2 > 0,
      median(attenuation[sig_m2], na.rm = TRUE),
      NA_real_
    ),
    mean_attenuation_M2sig = ifelse(
      n_sig_M2 > 0,
      mean(attenuation[sig_m2], na.rm = TRUE),
      NA_real_
    ),
    .groups = "drop"
  )

write_csv(
  retention_summary,
  file.path(table_dir, "Verified_Model2_to_Model3_retention_and_attenuation.csv")
)

print(retention_summary, n = Inf)


# ======================================================================
# 5. FIGURE 1A — SIGNIFICANT ASSOCIATION COUNTS
# ======================================================================

order_levels <- counts_long %>%
  filter(model == "Model 3") %>%
  arrange(desc(n_sig), symptom) %>%
  pull(symptom)

counts_plot <- counts_long %>%
  mutate(
    symptom = factor(symptom, levels = rev(order_levels)),
    model = factor(model, levels = c("Model 1", "Model 2", "Model 3"))
  )

p1A <- ggplot(counts_plot, aes(x = symptom, y = n_sig, fill = model)) +
  geom_col(
    position = position_dodge2(width = 0.8, preserve = "single"),
    width = 0.7,
    color = "white",
    linewidth = 0.25
  ) +
  geom_text(
    aes(label = n_sig),
    position = position_dodge2(width = 0.8, preserve = "single"),
    hjust = -0.12,
    size = 3.7
  ) +
  coord_flip() +
  scale_fill_manual(values = model_cols, drop = FALSE) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.16))) +
  labs(
    x = NULL,
    y = "Number of Bonferroni-significant metabolite associations",
    fill = NULL
  ) +
  theme_classic(base_size = 14) +
  theme(
    legend.position = "bottom",
    axis.text.y = element_text(size = 12),
    axis.text.x = element_text(size = 11),
    axis.title.x = element_text(size = 13)
  )

ggsave(
  file.path(fig_dir, "Figure1A_BonferroniCounts_models123.png"),
  p1A, width = 10.5, height = 6.2, dpi = 400
)
ggsave(
  file.path(fig_dir, "Figure1A_BonferroniCounts_models123.pdf"),
  p1A, width = 10.5, height = 6.2
)


# ======================================================================
# 6. FIGURE 1B/C — MODEL 2 AND MODEL 3 HEATMAPS
# ======================================================================

met_group <- function(f) {
  f_low <- tolower(f)
  
  if (grepl("glyca", f_low)) return("Inflammation")
  if (grepl("albumin|creatinine", f_low)) return("Proteins / fluid balance")
  if (grepl("^glucose$|^lactate$|^pyruvate$|^citrate$|^acetate$|acetoacetate|acetone|3hb", f_low))
    return("Small molecules / glycolysis")
  if (grepl("^val$|^leu$|^ile$|bcaa|^phe$|^tyr$|^trp$|^hist$|^gln$|^ala$|^gly$|^ser$|^pro$", f_low))
    return("Amino acids")
  if (grepl("omega|pufa|mufa|sfa|unsaturation|dha|epa|total_fa|^la$|fatty", f_low))
    return("Fatty acids")
  if (grepl("apoa|apob|apolipop", f_low))
    return("Apolipoproteins")
  if (grepl("vldl|idl|ldl|hdl", f_low))
    return("Lipoproteins")
  if (grepl("chol|_ce|_fc|total_c", f_low))
    return("Cholesterol")
  if (grepl("tg|triglycer", f_low))
    return("Triglycerides")
  if (grepl("pl|phosph", f_low))
    return("Phospholipids")
  "Other"
}

met_order_tbl <- tibble(feature = metabolite_features) %>%
  mutate(group = vapply(feature, met_group, character(1))) %>%
  arrange(group, feature) %>%
  mutate(feature_order = row_number())

met_order <- met_order_tbl$feature

heat_symptom_order <- c(
  "Depressed mood*", "Anhedonia*", "Hopelessness", "Worthlessness",
  "Low self-confidence", "Attention problems", "Thoughts of death",
  "Fatigue", "Hypersomnia", "Insomnia", "Sleep problems (any)",
  "Appetite or weight change", "Weight gain", "Weight loss"
)

heat_dat <- bind_rows(m2, m3) %>%
  mutate(
    symptom = unname(symptom_labels[variable]),
    symptom = factor(symptom, levels = heat_symptom_order),
    feature = factor(feature, levels = rev(met_order))
  )

# Common cap for comparability across Model 2 and Model 3.
heat_cap <- quantile(abs(heat_dat$beta), probs = 0.995, na.rm = TRUE)

make_heatmap <- function(df, model_name) {
  ggplot(df %>% filter(model == model_name),
         aes(x = symptom, y = feature, fill = beta)) +
    geom_tile() +
    geom_point(
      data = df %>% filter(model == model_name, sig_bonf),
      shape = 16, size = 0.35, color = "black"
    ) +
    scale_fill_gradient2(
      low = "#2166AC",
      mid = "white",
      high = "#B35806",
      midpoint = 0,
      limits = c(-heat_cap, heat_cap),
      oob = scales::squish,
      name = "log(OR)"
    ) +
    labs(x = NULL, y = NULL, title = model_name) +
    theme_minimal(base_size = 11) +
    theme(
      panel.grid = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      plot.title = element_text(face = "bold", size = 12)
    )
}

p1B <- make_heatmap(heat_dat, "Model 2")
p1C <- make_heatmap(heat_dat, "Model 3")

p1BC <- p1B | p1C

ggsave(
  file.path(fig_dir, "Figure1BC_heatmaps_Model2_Model3.png"),
  p1BC, width = 13.5, height = 8.0, dpi = 400
)
ggsave(
  file.path(fig_dir, "Figure1BC_heatmaps_Model2_Model3.pdf"),
  p1BC, width = 13.5, height = 8.0
)


# ======================================================================
# 7. FIGURE 2A — MODEL 2 vs MODEL 3 EFFECT ESTIMATES
# ======================================================================

cmp_sig <- cmp %>%
  filter(sig_m2, is.finite(OR2), is.finite(OR3)) %>%
  mutate(class = factor(vapply(feature, met_group, character(1))))

class_levels <- c(
  "Inflammation",
  "Proteins / fluid balance",
  "Small molecules / glycolysis",
  "Amino acids",
  "Fatty acids",
  "Apolipoproteins",
  "Lipoproteins",
  "Cholesterol",
  "Triglycerides",
  "Phospholipids",
  "Other"
)

cmp_sig$class <- factor(cmp_sig$class, levels = class_levels)

class_cols <- c(
  "Inflammation" = "#1F77B4",
  "Proteins / fluid balance" = "#6B6ECF",
  "Small molecules / glycolysis" = "#9467BD",
  "Amino acids" = "#17BECF",
  "Fatty acids" = "#FF7F0E",
  "Apolipoproteins" = "#7F7F7F",
  "Lipoproteins" = "#2CA02C",
  "Cholesterol" = "#E377C2",
  "Triglycerides" = "#BCBD22",
  "Phospholipids" = "#8C564B",
  "Other" = "#BDBDBD"
)

# Label the strongest absolute coefficient changes, but avoid repeatedly
# labeling the same metabolite code.
top_hits <- cmp_sig %>%
  mutate(abs_delta_beta = abs(beta3 - beta2)) %>%
  arrange(desc(abs_delta_beta)) %>%
  group_by(feature) %>%
  slice_head(n = 1) %>%
  ungroup() %>%
  slice_head(n = 10)

common_or_limits <- c(0.78, 1.29)

p2A <- ggplot(cmp_sig, aes(x = OR2, y = OR3, color = class)) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed",
              linewidth = 0.6, color = "grey40") +
  geom_point(alpha = 0.70, size = 2.2) +
  scale_x_log10(limits = common_or_limits) +
  scale_y_log10(limits = common_or_limits) +
  scale_color_manual(values = class_cols, drop = FALSE, name = "Metabolite group") +
  coord_fixed() +
  labs(
    x = "Odds ratio (Model 2)",
    y = "Odds ratio (Model 3)"
  ) +
  theme_classic(base_size = 13) +
  theme(
    legend.position = "right",
    legend.text = element_text(size = 9),
    legend.title = element_text(size = 10)
  )

if (use_ggrepel) {
  p2A <- p2A +
    ggrepel::geom_text_repel(
      data = top_hits,
      aes(label = feature),
      size = 3.0,
      color = "black",
      box.padding = 0.30,
      point.padding = 0.15,
      max.overlaps = Inf,
      min.segment.length = 0,
      show.legend = FALSE
    )
}

ggsave(
  file.path(fig_dir, "Figure2A_Model2_vs_Model3_ORs.png"),
  p2A, width = 8.2, height = 7.0, dpi = 400
)
ggsave(
  file.path(fig_dir, "Figure2A_Model2_vs_Model3_ORs.pdf"),
  p2A, width = 8.2, height = 7.0
)


# ======================================================================
# 8. FIGURE 2B — BMI ATTENUATION BY SYMPTOM
# ======================================================================

att_plot_dat <- cmp %>%
  filter(sig_m2, is.finite(attenuation)) %>%
  mutate(symptom = unname(symptom_labels[variable]))

att_order <- att_plot_dat %>%
  group_by(symptom) %>%
  summarise(median_att = median(attenuation, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(median_att)) %>%
  pull(symptom)

att_plot_dat <- att_plot_dat %>%
  mutate(symptom = factor(symptom, levels = rev(att_order)))

# Annotation is M2 significant -> SAME PAIRS retained in M3.
att_annotations <- retention_summary %>%
  mutate(
    symptom = factor(symptom, levels = rev(att_order)),
    label = paste0(n_sig_M2, " \u2192 ", n_retained_M2_to_M3)
  )

p2B <- ggplot(att_plot_dat, aes(x = attenuation, y = symptom)) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.45,
             color = "grey40") +
  geom_boxplot(
    width = 0.36,
    outlier.shape = NA,
    linewidth = 0.35,
    fill = "white",
    color = "grey25"
  ) +
  geom_point(
    position = position_jitter(height = 0.13, width = 0),
    alpha = 0.32,
    size = 1.4,
    color = "#1F77B4"
  ) +
  stat_summary(
    fun = median,
    geom = "point",
    size = 2.2,
    color = "black"
  ) +
  geom_text(
    data = att_annotations,
    aes(x = -0.47, y = symptom, label = label),
    inherit.aes = FALSE,
    size = 5,
    hjust = 0,
    fontface = "bold"
  ) +
  coord_cartesian(xlim = c(-0.5, 1.0), clip = "off") +
  labs(
    x = "BMI attenuation (Model 2 -> Model 3)",
    y = NULL
  ) +
  theme_classic(base_size = 13) +
  theme(
    axis.text.y = element_text(size = 16),
    axis.text.x = element_text(size = 16),
    axis.title      = element_text(size = 18),
    plot.margin = margin(10, 30, 10, 10)
  )

ggsave(
  file.path(fig_dir, "Figure2B_BMI_attenuation_by_symptom.png"),
  p2B, dpi = 400
)
ggsave(
  file.path(fig_dir, "Figure2B_BMI_attenuation_by_symptom.pdf"),
  p2B
)

fig2 <- p2A | p2B
ggsave(
  file.path(fig_dir, "Figure2_AB_combined.png"),
  fig2, width = 16.0, height = 7.0, dpi = 400
)
ggsave(
  file.path(fig_dir, "Figure2_AB_combined.pdf"),
  fig2, width = 16.0, height = 7.0
)


# ======================================================================
# 9. CORRECT SUPPLEMENTARY TABLES S2–S5
# ======================================================================

make_metwas_supp <- function(df, model_label) {
  df %>%
    transmute(
      model = model_label,
      variable,
      feature,
      N_cases,
      N_controls,
      OR = odds_ratio,
      OR_CI_low_95 = CI_low,
      OR_CI_high_95 = CI_high,
      p_value,
      bonf_alpha_249x14 = bonf_alpha,
      sig_bonf_249x14 = sig_bonf
    ) %>%
    arrange(variable, feature)
}

S2 <- make_metwas_supp(m1, "Model 1")
S3 <- make_metwas_supp(m2, "Model 2")
S4 <- make_metwas_supp(m3, "Model 3")

write_csv(S2, file.path(table_dir, "Supplementary_Table2_Model1_CORRECTED.csv"))
write_csv(S3, file.path(table_dir, "Supplementary_Table3_Model2_CORRECTED.csv"))
write_csv(S4, file.path(table_dir, "Supplementary_Table4_Model3_CORRECTED.csv"))

# By Methods, attenuation is defined for associations significant in Model 2.
S5 <- cmp %>%
  filter(sig_m2) %>%
  transmute(
    Symptom = variable,
    Metabolite = feature,
    Beta_Model2 = beta2,
    Beta_Model3 = beta3,
    OR_Model2 = OR2,
    OR_Model3 = OR3,
    p_Model2 = p2,
    p_Model3 = p3,
    Sig_Model2 = sig_m2,
    Sig_Model3 = sig_m3,
    Retained_in_Model3 = retained,
    Attenuation = attenuation,
    Direction_consistent = same_direction
  ) %>%
  arrange(Symptom, desc(Attenuation))

write_csv(S5, file.path(table_dir, "Supplementary_Table5_BMI_attenuation_CORRECTED.csv"))


# ======================================================================
# 10. BUILD THE CORRECT BMI-INTERACTION CANDIDATE SET
#     Main effect Bonferroni-significant in EITHER M2 or M3
# ======================================================================

bmi_candidate_pairs <- cmp %>%
  filter(sig_m2 | sig_m3) %>%
  transmute(
    symptom = variable,
    metabolite = feature,
    sig_M2 = sig_m2,
    sig_M3 = sig_m3,
    p_M2 = p2,
    p_M3 = p3,
    OR_M2 = OR2,
    OR_M3 = OR3
  ) %>%
  distinct(symptom, metabolite, .keep_all = TRUE) %>%
  arrange(symptom, metabolite)

if (nrow(bmi_candidate_pairs) != 691L) {
  stop(
    "Expected 691 unique BMI-interaction candidate pairs from M2/M3 union; found ",
    nrow(bmi_candidate_pairs),
    "."
  )
}

write_csv(
  bmi_candidate_pairs,
  file.path(table_dir, "BMI_interaction_candidate_pairs_691.csv")
)

message("Verified BMI-interaction candidate set: n = ", nrow(bmi_candidate_pairs))


# ======================================================================
# 11. OPTIONAL: RERUN ALL 691 BMI INTERACTIONS + FIGURE 3
#     Requires individual-level MetWAS file.
# ======================================================================

bmi_var <- "PersonPortrait.lastBmi"

covariates_model2 <- c(
  "age_at_sample",
  "Person.gender.code",
  "deltaTime",
  "PersonPortrait.lastSmokingStatus.name",
  "AlcoholConsumptionGroup",
  "PersonPortrait.educationGroup.code",
  "A10_before_blood",
  "C02_before_blood",
  "C10_before_blood"
)

coerce_binary01 <- function(x, name = "outcome") {
  if (is.logical(x)) return(as.integer(x))
  if (is.factor(x)) x <- as.character(x)
  y <- suppressWarnings(as.numeric(x))
  vals <- sort(unique(y[!is.na(y)]))
  if (!all(vals %in% c(0, 1))) {
    stop(name, " is not coded 0/1. Values: ", paste(vals, collapse = ", "))
  }
  y
}

mode_value <- function(x) {
  x_nonmiss <- x[!is.na(x)]
  if (length(x_nonmiss) == 0) return(NA)
  ux <- unique(x_nonmiss)
  ux[which.max(tabulate(match(x_nonmiss, ux)))]
}

reference_covariate_frame <- function(df, covars, n_rows = 1L) {
  out <- vector("list", length(covars))
  names(out) <- covars
  
  for (v in covars) {
    x <- df[[v]]
    
    if (is.numeric(x)) {
      val <- median(x, na.rm = TRUE)
      out[[v]] <- rep(val, n_rows)
    } else if (is.factor(x)) {
      val <- as.character(mode_value(x))
      out[[v]] <- factor(rep(val, n_rows), levels = levels(x))
    } else {
      val <- mode_value(x)
      out[[v]] <- rep(val, n_rows)
    }
  }
  
  as.data.frame(out, stringsAsFactors = FALSE)
}

prepare_interaction_data <- function(df, symptom, metabolite, include_profile = FALSE) {
  vars <- c(symptom, metabolite, bmi_var, covariates_model2)
  if (include_profile) vars <- c(vars, "final_cluster")
  
  missing <- setdiff(vars, names(df))
  if (length(missing) > 0) {
    stop(
      "Missing variables for ", symptom, " × ", metabolite, ": ",
      paste(missing, collapse = ", ")
    )
  }
  
  d <- df %>%
    select(all_of(vars)) %>%
    drop_na()
  
  d[[symptom]] <- coerce_binary01(d[[symptom]], symptom)
  
  met_mean <- mean(d[[metabolite]])
  met_sd <- sd(d[[metabolite]])
  bmi_mean <- mean(d[[bmi_var]])
  bmi_sd <- sd(d[[bmi_var]])
  
  if (!is.finite(met_sd) || met_sd <= 0) stop("Zero/invalid SD: ", metabolite)
  if (!is.finite(bmi_sd) || bmi_sd <= 0) stop("Zero/invalid BMI SD.")
  
  d <- d %>%
    mutate(
      met_z = (.data[[metabolite]] - met_mean) / met_sd,
      bmi_z = (.data[[bmi_var]] - bmi_mean) / bmi_sd
    )
  
  if (include_profile) {
    d <- d %>%
      mutate(
        final_cluster = factor(final_cluster, levels = paste0("C", 1:6))
      ) %>%
      filter(!is.na(final_cluster)) %>%
      droplevels()
  }
  
  attr(d, "met_mean") <- met_mean
  attr(d, "met_sd") <- met_sd
  attr(d, "bmi_mean") <- bmi_mean
  attr(d, "bmi_sd") <- bmi_sd
  d
}

run_one_bmi_interaction <- function(df, symptom, metabolite) {
  d <- prepare_interaction_data(df, symptom, metabolite, include_profile = FALSE)
  
  if (length(unique(d[[symptom]])) < 2) {
    return(tibble(
      Symptom = symptom, Metabolite = metabolite,
      beta_int = NA_real_, SE = NA_real_, p_value = NA_real_,
      OR_int = NA_real_, OR_CI_low = NA_real_, OR_CI_high = NA_real_,
      N_used = nrow(d), note = "Outcome has no variation"
    ))
  }
  
  f <- reformulate(
    c("met_z", "bmi_z", "met_z:bmi_z", covariates_model2),
    response = symptom
  )
  
  fit <- glm(f, data = d, family = binomial())
  
  tt <- broom::tidy(fit) %>% filter(term %in% c("met_z:bmi_z", "bmi_z:met_z"))
  if (nrow(tt) != 1L) {
    return(tibble(
      Symptom = symptom, Metabolite = metabolite,
      beta_int = NA_real_, SE = NA_real_, p_value = NA_real_,
      OR_int = NA_real_, OR_CI_low = NA_real_, OR_CI_high = NA_real_,
      N_used = nobs(fit), note = "Interaction term not uniquely found"
    ))
  }
  
  tibble(
    Symptom = symptom,
    Metabolite = metabolite,
    beta_int = tt$estimate,
    SE = tt$std.error,
    p_value = tt$p.value,
    OR_int = exp(tt$estimate),
    OR_CI_low = exp(tt$estimate - 1.96 * tt$std.error),
    OR_CI_high = exp(tt$estimate + 1.96 * tt$std.error),
    N_used = nobs(fit),
    note = NA_character_
  )
}

fit_bmi_panel <- function(df, symptom, metabolite,
                          met_grid = seq(-2.5, 2.5, length.out = 151),
                          bmi_levels = c(-1, 0, 1)) {
  d <- prepare_interaction_data(df, symptom, metabolite, include_profile = FALSE)
  
  f <- reformulate(
    c("met_z", "bmi_z", "met_z:bmi_z", covariates_model2),
    response = symptom
  )
  fit <- glm(f, data = d, family = binomial())
  
  grid <- tidyr::crossing(
    met_z = met_grid,
    bmi_z = bmi_levels
  )
  
  ref <- reference_covariate_frame(d, covariates_model2, nrow(grid))
  newdat <- bind_cols(grid, ref)
  
  pr <- predict(fit, newdata = newdat, type = "link", se.fit = TRUE)
  
  newdat %>%
    mutate(
      probability = plogis(pr$fit),
      low = plogis(pr$fit - 1.96 * pr$se.fit),
      high = plogis(pr$fit + 1.96 * pr$se.fit),
      BMI_level = factor(
        bmi_z,
        levels = c(-1, 0, 1),
        labels = c("Low BMI (-1 SD)", "Mean BMI", "High BMI (+1 SD)")
      ),
      Symptom = symptom,
      Metabolite = metabolite
    )
}

if (file.exists(individual_path)) {
  
  message("Individual-level data found. Running 691 BMI interaction models...")
  individual <- read.csv(individual_path, sep = ";", check.names = FALSE)
  
  # Verify all candidate symptom/metabolite columns are present.
  missing_indiv <- setdiff(
    unique(c(
      bmi_candidate_pairs$symptom,
      bmi_candidate_pairs$metabolite,
      bmi_var,
      covariates_model2
    )),
    names(individual)
  )
  if (length(missing_indiv) > 0) {
    stop(
      "Individual-level file is missing columns needed for the 691-pair scan:\n",
      paste(missing_indiv, collapse = ", ")
    )
  }
  
  S6 <- pmap_dfr(
    list(bmi_candidate_pairs$symptom, bmi_candidate_pairs$metabolite),
    ~ run_one_bmi_interaction(individual, ..1, ..2)
  ) %>%
    mutate(
      p_FDR = p.adjust(p_value, method = "BH"),
      Direction = case_when(
        beta_int > 0 ~ "positive",
        beta_int < 0 ~ "negative",
        TRUE ~ NA_character_
      )
    ) %>%
    arrange(p_FDR, p_value)
  
  if (nrow(S6) != 691L) {
    stop("S6 did not return 691 rows.")
  }
  if (any(!is.na(S6$note))) {
    warning("Some BMI interaction models were skipped; inspect S6$note.")
  }
  
  write_csv(S6, file.path(table_dir, "Supplementary_Table6_BMI_interactions_CORRECTED.csv"))
  
  bmi_sig <- S6 %>% filter(!is.na(p_FDR), p_FDR < 0.05)
  
  message("BMI interactions FDR < 0.05: ", nrow(bmi_sig))
  message(
    "Interaction beta range: ",
    signif(min(S6$beta_int, na.rm = TRUE), 4), " to ",
    signif(max(S6$beta_int, na.rm = TRUE), 4)
  )
  message(
    "Interaction OR range: ",
    signif(min(S6$OR_int, na.rm = TRUE), 4), " to ",
    signif(max(S6$OR_int, na.rm = TRUE), 4)
  )
  
  write_csv(
    bmi_sig,
    file.path(table_dir, "BMI_interactions_FDRsignificant.csv")
  )
  
  # Six manuscript Figure-3 examples.
  # All six are FDR-significant in the current corrected S6 results.
  fig3_pairs <- tribble(
    ~panel, ~Symptom,                  ~Metabolite,      ~title,
    "A",    "deprDepressedEver",       "Gln",            "Glutamine × BMI — Depressed mood",
    "B",    "deprLossOfInterestEver",  "Albumin",        "Albumin × BMI — Anhedonia",
    "C",    "deprLossOfInterestEver",  "Omega_3_pct",    "Omega-3 (%) × BMI — Anhedonia",
    "D",    "deprThoughtsOfDeath",     "Omega_3",        "Omega-3 × BMI — Thoughts of death",
    "E",    "deprTired",               "ApoB_by_ApoA1",  "ApoB/ApoA1 × BMI — Fatigue",
    "F",    "deprWorthlessness",       "GlycA",          "GlycA × BMI — Worthlessness"
  )
  
  fig3_check <- fig3_pairs %>%
    left_join(
      S6 %>% select(Symptom, Metabolite, p_FDR),
      by = c("Symptom", "Metabolite")
    )
  
  if (any(is.na(fig3_check$p_FDR)) || any(fig3_check$p_FDR >= 0.05)) {
    stop(
      "At least one predefined Figure 3 example is not FDR-significant in the rerun. ",
      "Inspect fig3_check before plotting."
    )
  }
  
  bmi_plot_cols <- c(
    "Low BMI (-1 SD)" = "#0072B2",
    "Mean BMI" = "#E69F00",
    "High BMI (+1 SD)" = "#CC79A7"
  )
  
  fig3_plots <- pmap(
    fig3_pairs,
    function(panel, Symptom, Metabolite, title) {
      pd <- fit_bmi_panel(individual, Symptom, Metabolite)
      
      ggplot(pd, aes(
        x = met_z, y = probability,
        color = BMI_level, fill = BMI_level
      )) +
        geom_ribbon(aes(ymin = low, ymax = high), alpha = 0.16, color = NA) +
        geom_line(linewidth = 1.0) +
        scale_color_manual(values = bmi_plot_cols) +
        scale_fill_manual(values = bmi_plot_cols) +
        labs(
          title = paste0(panel, "  ", title),
          x = paste0(Metabolite, " (z-score)"),
          y = "Predicted probability",
          color = "BMI",
          fill = "BMI"
        ) +
        theme_classic(base_size = 11) +
        theme(
          plot.title = element_text(face = "bold", size = 10.5),
          legend.position = "bottom"
        )
    }
  )
  
  fig3 <- wrap_plots(fig3_plots, ncol = 2, guides = "collect") &
    theme(legend.position = "bottom")
  
  ggsave(
    file.path(fig_dir, "Figure3_BMI_interactions_CORRECTED.png"),
    fig3, width = 11.0, height = 11.5, dpi = 400
  )
  ggsave(
    file.path(fig_dir, "Figure3_BMI_interactions_CORRECTED.pdf"),
    fig3, width = 11.0, height = 11.5
  )
  
} else {
  message(
    "Individual-level MetWAS file not found at: ", individual_path,
    "\nFigures 1–2 and S2–S5 are complete. Figure 3/S6 were skipped."
  )
  S6 <- NULL
  individual <- NULL
}


# ======================================================================
# 12. OPTIONAL: CORRECT SOM-PROFILE FOLLOW-UP + FIGURE 4 / S9
#
# Scientific estimand:
#   Among metabolite–symptom pairs with FDR-significant BMI interaction,
#   is there additional heterogeneity by SOM profile after retaining the
#   continuous metabolite × BMI interaction?
#
# Base:
#   symptom ~ met_z * bmi_z + SOM_profile + Model2 covariates
# Full:
#   Base + met_z:SOM_profile
#
# LRT tests only the 5-df met_z:SOM_profile contribution.
# Profile-specific ORs are evaluated at mean BMI (bmi_z = 0).
# ======================================================================

id_var <- "Person.skood"

harmonize_som_profile <- function(x) {
  x <- as.character(x)
  if (all(na.omit(x) %in% as.character(1:6))) {
    x <- paste0("C", x)
  } else if (all(grepl("^Cluster[1-6]$", na.omit(x)))) {
    x <- str_replace(x, "^Cluster", "C")
  }
  factor(x, levels = paste0("C", 1:6))
}

get_profile_slopes_at_mean_bmi <- function(fit, clusters = paste0("C", 1:6)) {
  cf <- coef(fit)
  V <- vcov(fit)
  base <- "met_z"
  
  map_dfr(clusters, function(cl) {
    if (cl == clusters[1]) {
      beta <- cf[[base]]
      se <- sqrt(V[base, base])
    } else {
      int1 <- paste0("met_z:final_cluster", cl)
      int2 <- paste0("final_cluster", cl, ":met_z")
      int <- if (int1 %in% names(cf)) int1 else int2
      
      if (!int %in% names(cf)) {
        return(tibble(
          SOM_Profile = cl, Beta_logOR_per_SD = NA_real_, SE = NA_real_,
          OR = NA_real_, OR_CI_low = NA_real_, OR_CI_high = NA_real_
        ))
      }
      
      beta <- cf[[base]] + cf[[int]]
      se <- sqrt(
        V[base, base] +
          V[int, int] +
          2 * V[base, int]
      )
    }
    
    tibble(
      SOM_Profile = cl,
      Beta_logOR_per_SD = beta,
      SE = se,
      OR = exp(beta),
      OR_CI_low = exp(beta - 1.96 * se),
      OR_CI_high = exp(beta + 1.96 * se)
    )
  })
}

run_one_som_interaction <- function(df, symptom, metabolite) {
  d <- prepare_interaction_data(df, symptom, metabolite, include_profile = TRUE)
  
  if (nlevels(d$final_cluster) != 6L) {
    return(tibble(
      Symptom = symptom, Metabolite = metabolite, SOM_Profile = NA_character_,
      Beta_logOR_per_SD = NA_real_, SE = NA_real_, OR = NA_real_,
      OR_CI_low = NA_real_, OR_CI_high = NA_real_,
      Interaction_LRT_p = NA_real_, N = nrow(d),
      note = "Not all 6 SOM profiles represented"
    ))
  }
  
  f0 <- reformulate(
    c(
      "met_z", "bmi_z", "met_z:bmi_z",
      "final_cluster",
      covariates_model2
    ),
    response = symptom
  )
  
  f1 <- reformulate(
    c(
      "met_z", "bmi_z", "met_z:bmi_z",
      "final_cluster", "met_z:final_cluster",
      covariates_model2
    ),
    response = symptom
  )
  
  fit0 <- glm(f0, data = d, family = binomial())
  fit1 <- glm(f1, data = d, family = binomial())
  
  if (nobs(fit0) != nobs(fit1)) {
    stop("Nested SOM models used different N for ", symptom, " × ", metabolite)
  }
  
  lrt <- anova(fit0, fit1, test = "LRT")
  p_lrt <- lrt$`Pr(>Chi)`[2]
  
  get_profile_slopes_at_mean_bmi(fit1, levels(d$final_cluster)) %>%
    mutate(
      Symptom = symptom,
      Metabolite = metabolite,
      Interaction_LRT_p = p_lrt,
      N = nobs(fit1),
      note = NA_character_
    )
}

fit_som_prediction_panel <- function(df, symptom, metabolite,
                                     met_grid = seq(-2.5, 2.5, length.out = 121)) {
  d <- prepare_interaction_data(df, symptom, metabolite, include_profile = TRUE)
  
  f <- reformulate(
    c(
      "met_z", "bmi_z", "met_z:bmi_z",
      "final_cluster", "met_z:final_cluster",
      covariates_model2
    ),
    response = symptom
  )
  fit <- glm(f, data = d, family = binomial())
  
  grid <- tidyr::crossing(
    met_z = met_grid,
    bmi_z = 0,
    final_cluster = factor(paste0("C", 1:6), levels = paste0("C", 1:6))
  )
  
  ref <- reference_covariate_frame(d, covariates_model2, nrow(grid))
  newdat <- bind_cols(grid, ref)
  
  pr <- predict(fit, newdata = newdat, type = "link", se.fit = TRUE)
  
  newdat %>%
    mutate(
      probability = plogis(pr$fit),
      low = plogis(pr$fit - 1.96 * pr$se.fit),
      high = plogis(pr$fit + 1.96 * pr$se.fit),
      Symptom = symptom,
      Metabolite = metabolite
    )
}

if (!is.null(individual) && !is.null(S6) && file.exists(som_path)) {
  
  som_assign <- read.csv(som_path, sep = ",", check.names = FALSE)
  
  if (!id_var %in% names(individual)) stop("Missing ID in individual data: ", id_var)
  if (!id_var %in% names(som_assign)) stop("Missing ID in SOM file: ", id_var)
  if (!"final_cluster" %in% names(som_assign)) stop("SOM file lacks final_cluster.")
  
  som_data <- individual %>%
    left_join(
      som_assign %>% select(all_of(id_var), final_cluster),
      by = id_var
    ) %>%
    mutate(final_cluster = harmonize_som_profile(final_cluster)) %>%
    filter(!is.na(final_cluster))
  
  # Follow up ONLY BMI interactions that survived global FDR.
  som_candidate_pairs <- S6 %>%
    filter(!is.na(p_FDR), p_FDR < 0.05) %>%
    transmute(symptom = Symptom, metabolite = Metabolite) %>%
    distinct()
  
  message(
    "SOM follow-up candidate pairs (BMI-interaction FDR < 0.05): ",
    nrow(som_candidate_pairs)
  )
  
  # Current corrected S6 is expected to give 211.
  if (nrow(som_candidate_pairs) != 211L) {
    warning(
      "Expected 211 FDR-significant BMI interactions from the current final analysis, found ",
      nrow(som_candidate_pairs),
      ". Continue only if this difference is intended."
    )
  }
  
  som_slopes <- pmap_dfr(
    list(som_candidate_pairs$symptom, som_candidate_pairs$metabolite),
    ~ run_one_som_interaction(som_data, ..1, ..2)
  )
  
  som_lrt <- som_slopes %>%
    distinct(Symptom, Metabolite, Interaction_LRT_p, N, note) %>%
    mutate(
      Interaction_LRT_FDR = p.adjust(Interaction_LRT_p, method = "BH")
    ) %>%
    arrange(Interaction_LRT_FDR, Interaction_LRT_p)
  
  S9 <- som_slopes %>%
    left_join(
      som_lrt %>%
        select(
          Symptom, Metabolite, Interaction_LRT_p, N,
          Interaction_LRT_FDR
        ),
      by = c("Symptom", "Metabolite", "Interaction_LRT_p", "N")
    ) %>%
    mutate(
      OR_95CI = paste0(
        sprintf("%.3f", OR),
        " (", sprintf("%.3f", OR_CI_low),
        "–", sprintf("%.3f", OR_CI_high), ")"
      )
    ) %>%
    arrange(Interaction_LRT_FDR, Symptom, Metabolite, SOM_Profile)
  
  write_csv(
    S9,
    file.path(table_dir, "Supplementary_Table9_SOM_interactions_CORRECTED.csv")
  )
  write_csv(
    som_lrt,
    file.path(table_dir, "SOM_interaction_LRT_summary_CORRECTED.csv")
  )
  
  som_sig_pairs <- som_lrt %>%
    filter(!is.na(Interaction_LRT_FDR), Interaction_LRT_FDR < 0.05)
  
  write_csv(
    som_sig_pairs,
    file.path(table_dir, "SOM_interactions_FDRsignificant_CORRECTED.csv")
  )
  
  message("Corrected SOM interactions FDR < 0.05: ", nrow(som_sig_pairs))
  
  # --------------------------------------------------------------------
  # Figure 4A — descriptive metabolic profile characteristics
  # Fixed markers, NOT data-driven top-marker selection.
  # --------------------------------------------------------------------
  
  profile_markers <- c(
    "XL_VLDL_TG",
    "ApoB_by_ApoA1",
    "L_HDL_CE",
    "PUFA_by_MUFA",
    "GlycA",
    "Total_BCAA"
  )
  profile_markers <- profile_markers[profile_markers %in% names(som_data)]
  
  if (length(profile_markers) < 4) {
    warning("Fewer than 4 fixed profile markers available for Figure 4A.")
  }
  
  profile_summary <- som_data %>%
    group_by(final_cluster) %>%
    summarise(
      Age = mean(age_at_sample, na.rm = TRUE),
      BMI = mean(.data[[bmi_var]], na.rm = TRUE),
      across(all_of(profile_markers), ~ mean(.x, na.rm = TRUE)),
      .groups = "drop"
    )
  
  z_across_profiles <- function(x) {
    s <- sd(x, na.rm = TRUE)
    if (!is.finite(s) || s == 0) return(rep(0, length(x)))
    (x - mean(x, na.rm = TRUE)) / s
  }
  
  profile_heat <- profile_summary %>%
    mutate(across(-final_cluster, z_across_profiles)) %>%
    pivot_longer(-final_cluster, names_to = "Marker", values_to = "Profile_z") %>%
    mutate(
      final_cluster = factor(final_cluster, levels = paste0("C", 1:6)),
      Marker = factor(Marker, levels = c("Age", "BMI", profile_markers))
    )
  
  p4A <- ggplot(profile_heat, aes(x = Marker, y = final_cluster, fill = Profile_z)) +
    geom_tile(color = "white", linewidth = 0.4) +
    scale_fill_gradient2(
      low = "#2166AC",
      mid = "white",
      high = "#B35806",
      midpoint = 0,
      name = "Profile\nz-score"
    ) +
    labs(
      x = NULL,
      y = "SOM metabolic profile"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      panel.grid = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.text.y = element_text(size = 11)
    )
  
  # --------------------------------------------------------------------
  # Figure 4B — strongest corrected profile interactions
  # All lines are solid. Significance refers to the GLOBAL interaction LRT,
  # not individual profile-specific slopes.
  # Predictions are shown at mean BMI (bmi_z = 0).
  # --------------------------------------------------------------------
  
  if (nrow(som_sig_pairs) > 0) {
    main_som_pairs <- som_sig_pairs %>%
      slice_head(n = min(6L, nrow(som_sig_pairs)))
    
    pred4 <- pmap_dfr(
      list(main_som_pairs$Symptom, main_som_pairs$Metabolite),
      ~ fit_som_prediction_panel(som_data, ..1, ..2)
    ) %>%
      left_join(
        main_som_pairs %>%
          select(Symptom, Metabolite, Interaction_LRT_FDR),
        by = c("Symptom", "Metabolite")
      ) %>%
      mutate(
        Symptom_label = unname(symptom_labels[Symptom]),
        pair = paste0(
          Symptom_label, " × ", Metabolite,
          "\nFDR=", signif(Interaction_LRT_FDR, 3)
        )
      )
    
    pair_order <- main_som_pairs %>%
      mutate(
        Symptom_label = unname(symptom_labels[Symptom]),
        pair = paste0(
          Symptom_label, " × ", Metabolite,
          "\nFDR=", signif(Interaction_LRT_FDR, 3)
        )
      ) %>%
      pull(pair)
    
    pred4 <- pred4 %>%
      mutate(
        pair = factor(pair, levels = pair_order),
        final_cluster = factor(final_cluster, levels = paste0("C", 1:6))
      )
    
    p4B <- ggplot(
      pred4,
      aes(
        x = met_z, y = probability,
        color = final_cluster, fill = final_cluster,
        group = final_cluster
      )
    ) +
      geom_ribbon(
        aes(ymin = low, ymax = high),
        alpha = 0.12, color = NA
      ) +
      geom_line(linewidth = 0.9) +
      facet_wrap(~ pair, ncol = 3, scales = "free_y") +
      scale_color_manual(values = som_cols, name = "Metabolic profile") +
      scale_fill_manual(values = som_cols, name = "Metabolic profile") +
      labs(
        x = "Metabolite (z-score)",
        y = "Predicted probability at mean BMI"
      ) +
      theme_classic(base_size = 11) +
      theme(
        legend.position = "bottom",
        strip.text = element_text(face = "bold", size = 9)
      )
    
    fig4 <- p4A / p4B + plot_layout(heights = c(0.75, 2.0))
    
    ggsave(
      file.path(fig_dir, "Figure4_SOM_profiles_and_interactions_CORRECTED.png"),
      fig4, width = 12.5, height = 10.5, dpi = 400
    )
    ggsave(
      file.path(fig_dir, "Figure4_SOM_profiles_and_interactions_CORRECTED.pdf"),
      fig4, width = 12.5, height = 10.5
    )
  } else {
    warning(
      "No SOM interaction survived the corrected conditional FDR analysis. ",
      "Figure 4B was not generated."
    )
  }
  
} else {
  message(
    "SOM follow-up skipped. It requires individual data, corrected S6, and SOM assignments."
  )
}


# ======================================================================
# 13. WRITE A MANUSCRIPT-READY NUMERICAL AUDIT
# ======================================================================

get_row <- function(var) {
  retention_summary %>% filter(variable == var) %>% slice(1)
}

wg <- get_row("deprGainWeight")
wl <- get_row("deprLostWeight")
ap <- get_row("deprAppetite1")
ins <- get_row("deprSleepLoss")
fat <- get_row("deprTired")
att <- get_row("deprAttention")
dm <- get_row("deprDepressedEver")
hop <- get_row("deprHopelessness")
tod <- get_row("deprThoughtsOfDeath")
anh <- get_row("deprLossOfInterestEver")
wor <- get_row("deprWorthlessness")
conf <- get_row("deprConfidence")

audit_lines <- c(
  "VERIFIED PRIMARY MetWAS NUMBERS",
  paste0("Bonferroni threshold: ", format(bonf_alpha, scientific = TRUE, digits = 8)),
  "Model 1 significant pairs: 1367",
  "Model 2 significant pairs: 660",
  "Model 3 significant pairs: 136",
  "BMI interaction candidate pairs (M2 OR M3 significant): 691",
  "",
  "MODEL 2 -> MODEL 3 RETENTION (same pair significant in both)",
  sprintf("Weight gain: %d/%d retained (%.1f%%); median attenuation %.3f",
          wg$n_retained_M2_to_M3, wg$n_sig_M2, wg$retained_pct, wg$median_attenuation_M2sig),
  sprintf("Weight loss: %d/%d retained (%.1f%%); median attenuation %.3f",
          wl$n_retained_M2_to_M3, wl$n_sig_M2, wl$retained_pct, wl$median_attenuation_M2sig),
  sprintf("Appetite/weight change: %d/%d retained (%.1f%%); median attenuation %.3f",
          ap$n_retained_M2_to_M3, ap$n_sig_M2, ap$retained_pct, ap$median_attenuation_M2sig),
  sprintf("Insomnia: %d/%d retained (%.1f%%); median attenuation %.3f",
          ins$n_retained_M2_to_M3, ins$n_sig_M2, ins$retained_pct, ins$median_attenuation_M2sig),
  sprintf("Fatigue: %d/%d retained (%.1f%%); median attenuation %.3f",
          fat$n_retained_M2_to_M3, fat$n_sig_M2, fat$retained_pct, fat$median_attenuation_M2sig),
  sprintf("Attention: %d/%d retained (%.1f%%); median attenuation %.3f",
          att$n_retained_M2_to_M3, att$n_sig_M2, att$retained_pct, att$median_attenuation_M2sig),
  sprintf("Depressed mood: %d/%d retained (%.1f%%); median attenuation %.3f",
          dm$n_retained_M2_to_M3, dm$n_sig_M2, dm$retained_pct, dm$median_attenuation_M2sig),
  sprintf("Hopelessness: %d/%d retained (%.1f%%); median attenuation %.3f",
          hop$n_retained_M2_to_M3, hop$n_sig_M2, hop$retained_pct, hop$median_attenuation_M2sig),
  sprintf("Thoughts of death: %d/%d retained (%.1f%%); median attenuation %.3f",
          tod$n_retained_M2_to_M3, tod$n_sig_M2, tod$retained_pct, tod$median_attenuation_M2sig),
  sprintf("Anhedonia: %d/%d retained (%.1f%%); median attenuation %.3f",
          anh$n_retained_M2_to_M3, anh$n_sig_M2, anh$retained_pct, anh$median_attenuation_M2sig),
  sprintf("Worthlessness: %d/%d retained; %d NEW in M3; median attenuation %.3f",
          wor$n_retained_M2_to_M3, wor$n_sig_M2, wor$n_new_in_M3, wor$median_attenuation_M2sig),
  sprintf("Low self-confidence: %d/%d retained; %d NEW in M3; median attenuation %.3f",
          conf$n_retained_M2_to_M3, conf$n_sig_M2, conf$n_new_in_M3, conf$median_attenuation_M2sig)
)

if (exists("S6") && !is.null(S6)) {
  audit_lines <- c(
    audit_lines,
    "",
    "BMI INTERACTION SCAN",
    paste0("Successful candidate rows: ", nrow(S6)),
    paste0("FDR < 0.05: ", sum(S6$p_FDR < 0.05, na.rm = TRUE)),
    paste0(
      "beta range: ",
      signif(min(S6$beta_int, na.rm = TRUE), 4),
      " to ",
      signif(max(S6$beta_int, na.rm = TRUE), 4)
    ),
    paste0(
      "interaction OR range: ",
      signif(min(S6$OR_int, na.rm = TRUE), 4),
      " to ",
      signif(max(S6$OR_int, na.rm = TRUE), 4)
    )
  )
}

if (exists("som_lrt")) {
  audit_lines <- c(
    audit_lines,
    "",
    "CORRECTED SOM FOLLOW-UP",
    paste0("Pairs tested: ", nrow(som_lrt)),
    paste0("FDR < 0.05: ", sum(som_lrt$Interaction_LRT_FDR < 0.05, na.rm = TRUE)),
    "NOTE: update all Figure 4/S9 manuscript p/FDR/OR claims from the corrected output."
  )
}

writeLines(
  audit_lines,
  con = file.path(table_dir, "MANUSCRIPT_NUMBERS_VERIFIED.txt")
)

cat(paste(audit_lines, collapse = "\n"), "\n")

message("\nDONE.")
message("Figures written to: ", fig_dir)
message("Tables + numerical audit written to: ", table_dir)
