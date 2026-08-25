#!/usr/bin/env Rscript

# Final SOM clustering and quality-control workflow.
# - Ward.D2 hierarchical clustering of SOM codebook vectors
# - fixed final K = 6
# - internal validation across candidate K values
# - quantization error and topographic error
# - participant-to-profile assignment for downstream analysis

suppressPackageStartupMessages({
  library(kohonen)
  library(clusterCrit)
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(readr)
  library(tibble)
})

som_path <- Sys.getenv("SOM_MODEL_PATH", unset = "outputs/som/som_21x21_rlen20000_online.rds")
data_path <- Sys.getenv("METWAS_DATA", unset = "data/private/analysis_ready_data.csv")
out_dir <- Sys.getenv("SOM_OUT_DIR", unset = "outputs/som")
id_var <- Sys.getenv("ID_VAR", unset = "Person.skood")
delim <- Sys.getenv("METWAS_DELIM", unset = ";")
final_k <- as.integer(Sys.getenv("SOM_FINAL_K", unset = "6"))
k_range <- 3:10

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
if (!file.exists(som_path)) stop("SOM model not found: ", som_path)
if (!file.exists(data_path)) stop("Analysis data not found: ", data_path)

som_model <- readRDS(som_path)
dat <- readr::read_delim(data_path, delim = delim, show_col_types = FALSE)
if (!id_var %in% names(dat)) stop("ID variable not found: ", id_var)
if (nrow(dat) != length(som_model$unit.classif)) {
  stop("Input rows do not match SOM training rows; assignment by row order is unsafe")
}

codes <- as.matrix(som_model$codes[[1]])
d_codes <- stats::dist(codes, method = "euclidean")
hc <- stats::hclust(d_codes, method = "ward.D2")

indices <- c(
  "silhouette", "dunn", "calinski_harabasz", "wemmert_gancarski",
  "davies_bouldin", "c_index"
)

safe_metrics <- function(k) {
  cl <- stats::cutree(hc, k = k)
  out <- tryCatch(
    clusterCrit::intCriteria(codes, as.integer(cl), indices),
    error = function(e) setNames(as.list(rep(NA_real_, length(indices))), indices)
  )
  tibble(k = k, !!!as.list(out))
}
validation <- purrr::map_dfr(k_range, safe_metrics)
readr::write_csv(validation, file.path(out_dir, "som_cluster_validation.csv"))

unit_cluster <- stats::cutree(hc, k = final_k)
sample_cluster <- unit_cluster[som_model$unit.classif]

assignments <- tibble(
  !!id_var := dat[[id_var]],
  som_unit = som_model$unit.classif,
  final_cluster = factor(paste0("C", sample_cluster), levels = paste0("C", 1:final_k))
)
readr::write_csv(assignments, file.path(out_dir, "som_assignments_k6.csv"))

compute_te <- function(som) {
  if (is.null(som$data)) return(NA_real_)
  x <- som$data[[1]]
  code <- som$codes[[1]]
  ud <- kohonen::unit.distances(som$grid)
  err <- logical(nrow(x))
  for (i in seq_len(nrow(x))) {
    xi <- as.numeric(x[i, ])
    d <- sqrt(rowSums((sweep(code, 2, xi, "-"))^2))
    bmus <- order(d)[1:2]
    err[i] <- ud[bmus[1], bmus[2]] > 1
  }
  mean(err)
}

qe <- if (!is.null(som_model$distances)) mean(som_model$distances, na.rm = TRUE) else NA_real_
te <- compute_te(som_model)
counts <- table(som_model$unit.classif)
n_units <- som_model$grid$xdim * som_model$grid$ydim
n_empty <- n_units - length(counts)

quality <- tibble(
  n = length(som_model$unit.classif),
  n_units = n_units,
  final_k = final_k,
  quantization_error = qe,
  topographic_error = te,
  n_empty_units = n_empty
)
readr::write_csv(quality, file.path(out_dir, "som_quality_metrics.csv"))
print(quality)
message("Saved K=6 assignments and SOM QC outputs to: ", out_dir)
