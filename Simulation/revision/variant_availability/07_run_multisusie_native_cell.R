#!/usr/bin/env Rscript
.infuse_repo <- Sys.getenv("INFUSE_ANALYSIS_ROOT", "/path/to/INFUSE_analysis")

args <- commandArgs(TRUE)
if (length(args) != 4L) {
  stop("Usage: Rscript 07_run_multisusie_native_cell.R <locus> <base_index> <h2_index> <n_assv_each>")
}
locus <- as.integer(args[[1L]])
base_index <- as.integer(args[[2L]])
h2_index <- as.integer(args[[3L]])
n_assv_each <- as.integer(args[[4L]])
stopifnot(locus %in% 1:100, base_index %in% 1:3, h2_index %in% 1:2,
          n_assv_each %in% 0:2)

this_file <- normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE),
                                                    value = TRUE)[1L]),
                           mustWork = TRUE)
NATIVE_SRC_DIR <- dirname(this_file)
source(file.path(.infuse_repo, "Simulation/revision/variant_availability/native_sim_common.R"))
suppressMessages(library(reticulate))

cell_file <- native_cell_path(locus, base_index, h2_index, n_assv_each)
maf_file <- native_multisusie_maf_path(locus)
out_file <- native_multisusie_path(locus, base_index, h2_index, n_assv_each)
if (file.exists(out_file)) {
  cat("complete MultiSuSiE output exists; leaving untouched:", out_file, "\n")
  quit(status = 0L)
}
if (!file.exists(cell_file)) stop("Native simulation cell is missing: ", cell_file)
if (!file.exists(maf_file)) stop("MultiSuSiE MAF cache is missing: ", maf_file)
dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)

cell <- readRDS(cell_file)
lo <- load_native_locus(locus, load_genotypes = FALSE)
mc <- readRDS(maf_file)
stopifnot(identical(as.character(cell$truth$SNP), lo$snp),
          identical(unname(cell$availability), unname(lo$avail)),
          identical(mc$ancestry, ANCESTRY))

# Fit the MultiSuSiE intersection comparator using common variants only.
# Map its results to union coordinates, assigning PIP zero to excluded variants.
# The union analysis is implemented separately in multisusie_union/run_locus.R.
shared_idx <- which(rowSums(lo$avail) == length(ANCESTRY))
shared_ids <- lo$snp[shared_idx]
if (!identical(shared_ids, mc$snp)) stop("Shared-SNP/MAF-cache mismatch")
R_list <- lapply(ANCESTRY, function(anc) {
  R <- lo$R_obs[[anc]]
  hit <- match(shared_ids, rownames(R))
  if (anyNA(hit)) stop("Shared-SNP/LD mismatch for ", anc)
  unname(R[hit, hit, drop = FALSE])
})
z_list <- list(as.numeric(cell$zfile$zscore_1[shared_idx]),
               as.numeric(cell$zfile$zscore_2[shared_idx]))
maf_list <- lapply(ANCESTRY, function(anc) unname(mc$maf[[anc]][shared_ids]))

use_python(Sys.getenv("INFUSE_PYTHON_BIN", Sys.which("python")), required = TRUE)
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
  population_sizes = np$array(as.integer(c(NUM_EU, NUM_BB))),
  L = as.integer(10L),
  scaled_prior_variance = 0.2,
  low_memory_mode = FALSE,
  min_abs_corr = 0.5,
  single_population_mac_thresh = as.integer(20L),
  maf_list = maf_list_py,
  coverage = 0.95
)
runtime_minutes <- (proc.time()[[3L]] - start) / 60

pip_shared_universe <- as.numeric(py_to_r(fit$pip))
if (length(pip_shared_universe) != length(shared_idx) ||
    any(!is.finite(pip_shared_universe)) ||
    any(pip_shared_universe < -1e-12 | pip_shared_universe > 1 + 1e-12)) {
  stop("Invalid MultiSuSiE PIP output")
}
pip <- numeric(lo$p)
pip[shared_idx] <- pmin(1, pmax(0, pip_shared_universe))

sets <- py_to_r(fit$sets)
component_sets <- sets[[1L]]
include <- as.logical(sets[[4L]])
if (length(component_sets) != length(include)) stop("Malformed MultiSuSiE CS output")
cs <- lapply(which(include), function(i) {
  idx0 <- as.integer(component_sets[[i]])
  if (!length(idx0)) return(integer())
  if (any(idx0 < 0L | idx0 >= length(shared_idx))) stop("MultiSuSiE CS index out of range")
  shared_idx[idx0 + 1L]
})
cs <- Filter(length, cs)

obj <- list(
  metadata = list(locus = locus, base_count_index = base_index,
                  n_base = NATIVE_BASE_COUNTS[base_index], h2_index = h2_index,
                  h2_per_causal = NATIVE_H2[h2_index],
                  n_assv_each = n_assv_each, L = 10L,
                  variant_set = "intersection"),
  pip = pip, cs = cs, runtime_minutes = runtime_minutes,
  shared_idx = shared_idx, shared_snp = shared_ids,
  package_path = normalizePath(py_to_r(ms$`__file__`), mustWork = TRUE)
)
stopifnot(all(obj$pip[-shared_idx] == 0),
          all(unlist(obj$cs, use.names = FALSE) %in% shared_idx))
tmp <- paste0(out_file, ".partial_", Sys.getpid())
saveRDS(obj, tmp, compress = FALSE)
if (!file.rename(tmp, out_file)) stop("Atomic MultiSuSiE output rename failed: ", out_file)
cat(sprintf("wrote %s; p_shared=%d; n_cs=%d; minutes=%.4f\n",
            out_file, length(shared_idx), length(cs), runtime_minutes))
