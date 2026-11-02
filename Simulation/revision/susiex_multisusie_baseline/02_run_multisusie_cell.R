#!/usr/bin/env Rscript
.infuse_repo <- Sys.getenv("INFUSE_ANALYSIS_ROOT", "/path/to/INFUSE_analysis")

args <- commandArgs(TRUE)
if (length(args) != 1L) stop("Usage: 02_run_multisusie_cell.R <task_id>", call. = FALSE)
this_file <- normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1L]),
                           mustWork = TRUE)
source(file.path(.infuse_repo, "Simulation/revision/susiex_multisusie_baseline/00_common.R"))
suppressPackageStartupMessages(library(reticulate))

cell <- decode_task(args[[1L]])
z <- read_z(cell, "d10")
out_file <- cell_output("multisusie", cell)
old <- safe_read_method_object(out_file, z, "MultiSuSiE")
if (!is.null(old)) {
  cat("Output exists; skipping:", out_file, "\n")
  quit(status = 0L)
}

R_list <- list(read_ld(cell, "EUR", z$RSID), read_ld(cell, "AFR", z$RSID))
freq_eur <- read_true_maf(cell$locus, "EUR", nrow(z))
freq_afr <- read_true_maf(cell$locus, "AFR", nrow(z))
maf_list <- list(freq_eur$maf, freq_afr$maf)
z_list <- list(as.numeric(z$zscore_1), as.numeric(z$zscore_2))

if (!file.exists(PYTHON_BIN)) stop("Configured Python executable does not exist: ", PYTHON_BIN, call. = FALSE)
use_python(PYTHON_BIN, required = TRUE)
ms <- import("MultiSuSiE", convert = FALSE)
np <- import("numpy", convert = FALSE)

R_list_py <- lapply(R_list, np$array)
z_list_py <- lapply(z_list, np$array)
maf_list_py <- lapply(maf_list, np$array)
start <- proc.time()[[3L]]
fit <- ms$multisusie_rss(
  z_list = z_list_py,
  R_list = R_list_py,
  rho = np$array(matrix(c(1, 0.8, 0.8, 1), 2L, byrow = TRUE)),
  population_sizes = np$array(as.integer(N_GWAS)),
  L = as.integer(10L),
  scaled_prior_variance = 0.2,
  low_memory_mode = FALSE,
  min_abs_corr = 0.5,
  single_population_mac_thresh = as.integer(20L),
  maf_list = maf_list_py,
  coverage = 0.95
)
runtime_minutes <- (proc.time()[[3L]] - start) / 60

pip <- as.numeric(py_to_r(fit$pip))
if (length(pip) != nrow(z) || any(!is.finite(pip)) || any(pip < -1e-10 | pip > 1 + 1e-10)) {
  stop("Invalid MultiSuSiE PIP output", call. = FALSE)
}
sets <- py_to_r(fit$sets)
component_sets <- sets[[1L]]
include <- as.logical(sets[[4L]])
if (length(component_sets) != length(include)) stop("Malformed MultiSuSiE CS output", call. = FALSE)
cs <- lapply(which(include), function(i) {
  idx0 <- as.integer(component_sets[[i]])
  if (!length(idx0)) return(integer())
  if (any(idx0 < 0L | idx0 >= nrow(z))) stop("MultiSuSiE CS index out of range", call. = FALSE)
  idx0 + 1L
})
cs <- unname(Filter(length, cs))

obj <- list(
  method = "MultiSuSiE", task_id = cell$task_id, locus = cell$locus,
  causal_index = cell$causal_index, causal_count = cell$causal_count,
  h2_index = cell$h2_index, h2 = cell$h2,
  rsid = z$RSID, truth = z$Signal, pip = pmin(1, pmax(0, pip)), cs = cs,
  runtime_minutes = runtime_minutes, L = 10L, coverage = 0.95,
  min_abs_corr = 0.5, rho = 0.8, scaled_prior_variance = 0.2,
  single_population_mac_thresh = 20L,
  input_frequency_semantics = "MAF = min(ALT frequency, 1 - ALT frequency)",
  eur_alt_frequency_above_half = freq_eur$n_above_half,
  afr_alt_frequency_above_half = freq_afr$n_above_half,
  package_path = normalizePath(py_to_r(ms$`__file__`), mustWork = TRUE)
)
if (!validate_method_object(obj, z, "MultiSuSiE")) stop("Internal MultiSuSiE output validation failed", call. = FALSE)
atomic_save_rds(obj, out_file)
cat(sprintf("wrote %s; n_snp=%d; n_cs=%d; minutes=%.4f\n",
            out_file, nrow(z), length(cs), runtime_minutes))

