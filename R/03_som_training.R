#!/usr/bin/env Rscript

# Clean SOM implementation corresponding to the final manuscript design.


suppressPackageStartupMessages({
  library(kohonen)
  library(readr)
})

set.seed(as.integer(Sys.getenv("SOM_SEED", unset = "2025")))

input_path <- Sys.getenv("METWAS_DATA", unset = "data/private/analysis_ready_data.csv")
out_dir <- Sys.getenv("SOM_OUT_DIR", unset = "outputs/som")
delim <- Sys.getenv("METWAS_DELIM", unset = ";")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Current cleaned defaults reflect the manuscript's direct-metabolite 21×21 SOM.
grid_x <- as.integer(Sys.getenv("SOM_XDIM", unset = "21"))
grid_y <- as.integer(Sys.getenv("SOM_YDIM", unset = "21"))
rlen <- as.integer(Sys.getenv("SOM_RLEN", unset = "20000"))
mode <- Sys.getenv("SOM_MODE", unset = "online")
alpha_start <- as.numeric(Sys.getenv("SOM_ALPHA_START", unset = "0.05"))
alpha_end <- as.numeric(Sys.getenv("SOM_ALPHA_END", unset = "0.01"))

if (!file.exists(input_path)) stop("Input not found: ", input_path)
dat <- readr::read_delim(input_path, delim = delim, show_col_types = FALSE)

start_col <- match("Total_C", names(dat))
end_col <- match("S_HDL_TG_pct", names(dat))
if (is.na(start_col) || is.na(end_col) || start_col > end_col) {
  stop("Could not identify metabolite span Total_C ... S_HDL_TG_pct")
}
metabolites <- names(dat)[start_col:end_col]
if (length(metabolites) != 249L) stop("Expected 249 metabolites; found ", length(metabolites))

X <- as.data.frame(dat[, metabolites, drop = FALSE])
for (j in seq_along(X)) {
  x <- X[[j]]
  if (anyNA(x)) x[is.na(x)] <- mean(x, na.rm = TRUE)
  X[[j]] <- x
}
X <- as.matrix(X)



som_grid <- kohonen::somgrid(xdim = grid_x, ydim = grid_y, topo = "hexagonal")

args <- list(
  X = X,
  grid = som_grid,
  rlen = rlen,
  keep.data = TRUE,
  mode = mode
)
# alpha is used by online training; kohonen ignores it for batch/pbatch.
if (mode == "online") args$alpha <- c(alpha_start, alpha_end)

som_model <- do.call(kohonen::som, args)

model_path <- file.path(out_dir, sprintf("som_%dx%d_rlen%d_%s.rds", grid_x, grid_y, rlen, mode))
saveRDS(som_model, model_path)

meta <- data.frame(
  grid_x = grid_x,
  grid_y = grid_y,
  rlen = rlen,
  mode = batch,
  seed = as.integer(Sys.getenv("SOM_SEED", unset = "2025")),
  n = nrow(X),
  metabolites = ncol(X)
)
readr::write_csv(meta, file.path(out_dir, "som_training_metadata.csv"))
message("Saved SOM model: ", model_path)
