#!/usr/bin/env Rscript

# ==============================================================================
# 03_som_training.R
#
# Train the metabolomic self-organizing map (SOM) used in the manuscript:
#
# "Body mass index modifies symptom-specific metabolomic associations with
# depressive symptoms in the Estonian Biobank"
#
# Final manuscript specification:
#   - 249 Nightingale NMR metabolite measures
#   - metabolite values are already QC'd and sex-standardized upstream
#   - metabolite-specific mean imputation for missing SOM inputs
#   - 21 x 21 hexagonal SOM
#   - bubble neighbourhood function
#   - online training
#   - rlen = 20,000
#   - alpha = c(0.05, 0.01)
#   - random seed = 2025
#
# This script trains the SOM only.
# Ward.D2 clustering and SOM quality-control analyses are performed in:
#   R/04_som_clustering_qc.R
#
# Individual-level Estonian Biobank data are not distributed with this
# repository.
# ==============================================================================


# ------------------------------------------------------------------------------
# 1. Packages
# ------------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(kohonen)
  library(readr)
  library(tibble)
})


# ------------------------------------------------------------------------------
# 2. Reproducible final SOM specification
# ------------------------------------------------------------------------------

# File locations may be changed using environment variables.
input_path <- Sys.getenv(
  "METWAS_DATA",
  unset = "data/private/analysis_ready_data.csv"
)

out_dir <- Sys.getenv(
  "SOM_OUT_DIR",
  unset = "outputs/som"
)

delim <- Sys.getenv(
  "METWAS_DELIM",
  unset = ";"
)

id_var <- Sys.getenv(
  "ID_VAR",
  unset = "Person.skood"
)

# IMPORTANT:
# These values are intentionally fixed because this script represents
# the SOM configuration reported in the manuscript.
grid_x <- 21L
grid_y <- 21L
rlen <- 20000L

training_mode <- "online"
neighbourhood_function <- "bubble"

alpha_start <- 0.05
alpha_end <- 0.01

seed <- 2025L

set.seed(seed)

dir.create(
  out_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# ------------------------------------------------------------------------------
# 3. Read analysis-ready data
# ------------------------------------------------------------------------------

if (!file.exists(input_path)) {
  stop(
    "Input data file was not found:\n",
    input_path,
    "\n\nSet METWAS_DATA to the analysis-ready Estonian Biobank dataset."
  )
}

message("Reading analysis-ready data: ", input_path)

dat <- readr::read_delim(
  input_path,
  delim = delim,
  show_col_types = FALSE,
  progress = FALSE
)

message(
  "Input dimensions: ",
  nrow(dat), " participants x ",
  ncol(dat), " variables"
)

if (nrow(dat) == 0L) {
  stop("Input dataset contains zero rows.")
}


# ------------------------------------------------------------------------------
# 4. Identify the 249 manuscript metabolites
# ------------------------------------------------------------------------------

# In the analysis-ready dataset, the manuscript Nightingale metabolite measures
# are stored contiguously from Total_C through S_HDL_TG_pct.
start_col <- match("Total_C", names(dat))
end_col <- match("S_HDL_TG_pct", names(dat))

if (
  is.na(start_col) ||
  is.na(end_col) ||
  start_col > end_col
) {
  stop(
    "Could not identify the metabolite block ",
    "'Total_C' ... 'S_HDL_TG_pct'."
  )
}

metabolites <- names(dat)[start_col:end_col]

if (length(metabolites) != 249L) {
  stop(
    "Expected exactly 249 metabolite measures, but found ",
    length(metabolites),
    ".\nCheck the analysis-ready input dataset before continuing."
  )
}

message("Identified 249 Nightingale metabolite measures.")


# ------------------------------------------------------------------------------
# 5. Construct SOM input matrix
# ------------------------------------------------------------------------------

# IMPORTANT:
# The analysis-ready input is expected to contain the metabolite measures
# already quality-controlled and sex-standardized as used in the manuscript.
#
# DO NOT call scale() here unless the manuscript analysis itself is changed.
# Performing another global scaling step would constitute a different
# preprocessing pipeline.

X <- as.data.frame(
  dat[, metabolites, drop = FALSE]
)

# Confirm that all metabolite inputs are numeric.
non_numeric <- metabolites[
  !vapply(X, is.numeric, logical(1))
]

if (length(non_numeric) > 0L) {
  stop(
    "The following metabolite columns are not numeric:\n",
    paste(non_numeric, collapse = ", ")
  )
}


# ------------------------------------------------------------------------------
# 6. Missing-data diagnostics and mean imputation
# ------------------------------------------------------------------------------

missing_before <- vapply(
  X,
  function(x) sum(is.na(x)),
  numeric(1)
)

missing_summary <- tibble::tibble(
  metabolite = metabolites,
  n_missing = as.integer(missing_before),
  pct_missing = 100 * n_missing / nrow(X)
)

readr::write_csv(
  missing_summary,
  file.path(out_dir, "som_missingness_before_imputation.csv")
)

if (any(missing_before == nrow(X))) {
  all_missing <- names(missing_before)[
    missing_before == nrow(X)
  ]

  stop(
    "At least one metabolite contains only missing values:\n",
    paste(all_missing, collapse = ", ")
  )
}

message(
  "Total missing metabolite values before SOM imputation: ",
  sum(missing_before)
)

# Metabolite-specific mean imputation.
#
# This is used only for SOM construction so that participants with occasional
# missing metabolite values can be retained. Mean imputation can reduce
# within-variable variance and is therefore reported as a limitation.
for (j in seq_along(X)) {

  x <- X[[j]]

  if (anyNA(x)) {

    x_mean <- mean(
      x,
      na.rm = TRUE
    )

    x[is.na(x)] <- x_mean

    X[[j]] <- x
  }
}

if (anyNA(X)) {
  stop(
    "Missing values remain after metabolite-specific mean imputation."
  )
}

X <- as.matrix(X)

storage.mode(X) <- "double"

rownames(X) <- NULL


# ------------------------------------------------------------------------------
# 7. Basic input checks
# ------------------------------------------------------------------------------

if (ncol(X) != 249L) {
  stop(
    "Internal error: SOM input matrix does not contain 249 metabolites."
  )
}

if (nrow(X) != nrow(dat)) {
  stop(
    "Internal error: SOM input row count differs from source dataset."
  )
}

if (any(!is.finite(X))) {
  stop(
    "SOM input matrix contains non-finite values after preprocessing."
  )
}

input_sd <- apply(
  X,
  2,
  stats::sd
)

if (any(input_sd == 0 | is.na(input_sd))) {

  bad <- metabolites[
    input_sd == 0 | is.na(input_sd)
  ]

  stop(
    "At least one metabolite has zero or undefined variance:\n",
    paste(bad, collapse = ", ")
  )
}

# Save descriptive input statistics as a reproducibility check.
input_summary <- tibble::tibble(
  metabolite = metabolites,
  mean = colMeans(X),
  sd = apply(X, 2, stats::sd),
  min = apply(X, 2, min),
  max = apply(X, 2, max)
)

readr::write_csv(
  input_summary,
  file.path(out_dir, "som_input_summary.csv")
)


# ------------------------------------------------------------------------------
# 8. Save participant row order
# ------------------------------------------------------------------------------

# SOM unit assignments are indexed according to input row order.
# Saving the identifier sequence provides an explicit check that downstream
# participant assignments use the same ordering.

if (id_var %in% names(dat)) {

  participant_order <- tibble::tibble(
    row_index = seq_len(nrow(dat)),
    participant_id = as.character(dat[[id_var]])
  )

  readr::write_csv(
    participant_order,
    file.path(out_dir, "som_participant_row_order.csv")
  )

} else {

  warning(
    "ID variable '",
    id_var,
    "' was not found. Participant row-order file was not written."
  )
}


# ------------------------------------------------------------------------------
# 9. Define final SOM grid
# ------------------------------------------------------------------------------

som_grid <- kohonen::somgrid(
  xdim = grid_x,
  ydim = grid_y,
  topo = "hexagonal",
  neighbourhood.fct = neighbourhood_function,
  toroidal = FALSE
)

message(
  "SOM grid: ",
  grid_x, " x ", grid_y,
  " hexagonal units (",
  grid_x * grid_y,
  " total units)"
)


# ------------------------------------------------------------------------------
# 10. Train final SOM
# ------------------------------------------------------------------------------

message("")
message("Training final manuscript SOM...")
message("  mode                 = ", training_mode)
message("  neighbourhood        = ", neighbourhood_function)
message("  grid                  = ", grid_x, " x ", grid_y)
message("  rlen                  = ", rlen)
message(
  "  alpha                 = ",
  alpha_start, " -> ", alpha_end
)
message("  seed                  = ", seed)
message("  participants          = ", nrow(X))
message("  metabolites           = ", ncol(X))
message("")

set.seed(seed)

som_model <- kohonen::som(
  X = X,
  grid = som_grid,
  rlen = rlen,
  alpha = c(alpha_start, alpha_end),
  keep.data = TRUE,
  mode = training_mode
)


# ------------------------------------------------------------------------------
# 11. Validate fitted SOM object
# ------------------------------------------------------------------------------

if (is.null(som_model$codes)) {
  stop("Fitted SOM object does not contain codebook vectors.")
}

if (is.null(som_model$unit.classif)) {
  stop("Fitted SOM object does not contain participant BMU assignments.")
}

if (length(som_model$unit.classif) != nrow(X)) {
  stop(
    "SOM assignment count does not match the number of participants."
  )
}

if (nrow(som_model$codes[[1]]) != grid_x * grid_y) {
  stop(
    "Unexpected number of SOM codebook vectors."
  )
}

message("SOM training completed successfully.")


# ------------------------------------------------------------------------------
# 12. Save trained model
# ------------------------------------------------------------------------------

model_filename <- sprintf(
  "som_%dx%d_rlen%d_%s.rds",
  grid_x,
  grid_y,
  rlen,
  training_mode
)

model_path <- file.path(
  out_dir,
  model_filename
)

saveRDS(
  som_model,
  model_path
)

message("Saved SOM model: ", model_path)


# ------------------------------------------------------------------------------
# 13. Save training metadata
# ------------------------------------------------------------------------------

training_metadata <- tibble::tibble(
  manuscript_metabolites = ncol(X),
  participants = nrow(X),

  grid_x = grid_x,
  grid_y = grid_y,
  n_units = grid_x * grid_y,

  topology = "hexagonal",
  neighbourhood_function = neighbourhood_function,
  toroidal = FALSE,

  training_mode = training_mode,
  rlen = rlen,

  alpha_start = alpha_start,
  alpha_end = alpha_end,

  random_seed = seed,

  mean_imputation = TRUE,
  additional_scaling_in_this_script = FALSE,

  kohonen_version = as.character(
    utils::packageVersion("kohonen")
  ),

  R_version = R.version.string
)

readr::write_csv(
  training_metadata,
  file.path(
    out_dir,
    "som_training_metadata.csv"
  )
)


# ------------------------------------------------------------------------------
# 14. Save training-change diagnostic
# ------------------------------------------------------------------------------

# kohonen stores training changes in the fitted object.
# The values are exported so the convergence-style diagnostic can be recreated
# without requiring the participant-level input data.

if (!is.null(som_model$changes)) {

  changes <- som_model$changes

  if (is.matrix(changes) || is.data.frame(changes)) {

    changes_df <- as.data.frame(changes)

    changes_df$training_step <- seq_len(
      nrow(changes_df)
    )

    changes_df <- changes_df[
      c(
        "training_step",
        setdiff(
          names(changes_df),
          "training_step"
        )
      )
    ]

  } else {

    changes_df <- tibble::tibble(
      training_step = seq_along(changes),
      mean_distance_change = as.numeric(changes)
    )
  }

  readr::write_csv(
    changes_df,
    file.path(
      out_dir,
      "som_training_changes.csv"
    )
  )
}


# ------------------------------------------------------------------------------
# 15. Save software environment
# ------------------------------------------------------------------------------

capture.output(
  sessionInfo(),
  file = file.path(
    out_dir,
    "som_sessionInfo.txt"
  )
)


# ------------------------------------------------------------------------------
# 16. Completion summary
# ------------------------------------------------------------------------------

message("")
message("============================================================")
message("Final SOM training completed")
message("============================================================")
message("Input participants: ", nrow(X))
message("Input metabolites:  ", ncol(X))
message("Grid:               ", grid_x, " x ", grid_y)
message("Training mode:      ", training_mode)
message("Neighbourhood:      ", neighbourhood_function)
message("rlen:               ", rlen)
message(
  "alpha:              ",
  alpha_start,
  " -> ",
  alpha_end
)
message("Random seed:         ", seed)
message("Model:               ", model_path)
message("")
message(
  "Next step: R/04_som_clustering_qc.R"
)
message("============================================================")
