#!/usr/bin/env Rscript
.infuse_repo <- Sys.getenv("INFUSE_ANALYSIS_ROOT", "/path/to/INFUSE_analysis")

args <- commandArgs(TRUE)
if (length(args) < 4L || length(args) > 5L) {
  stop("Usage: 08_run_infuse_d187_native_cell.R <locus> <base_index> <h2_index> <n_assv_each> [L]")
}
locus <- as.integer(args[1L]); base_index <- as.integer(args[2L])
h2_index <- as.integer(args[3L]); n_assv_each <- as.integer(args[4L])
L_effects <- if (length(args) == 5L) as.integer(args[5L]) else 10L

this_file <- normalizePath(sub("^--file=", "",
  grep("^--file=", commandArgs(FALSE), value = TRUE)[1L]), mustWork = TRUE)
NATIVE_SRC_DIR <- dirname(this_file)
source(file.path(.infuse_repo, "Simulation/revision/variant_availability/native_sim_common.R"))
suppressMessages(library(INFUSE))

stopifnot(locus %in% 1:100, base_index %in% 1:3, h2_index %in% 1:2,
          n_assv_each %in% NATIVE_ASSV_LEVELS, L_effects >= 1L)
out_file <- native_infuse_d187_path(locus, base_index, h2_index, n_assv_each)
if (file.exists(out_file)) {
  cat("complete d=187 output exists; leaving untouched:", out_file, "\n")
  quit(status = 0L)
}

cell_file <- native_cell_path(locus, base_index, h2_index, n_assv_each)
ann_file <- native_annotation187_path(locus)
if (!file.exists(cell_file)) stop("Missing native simulation cell: ", cell_file)
if (!file.exists(ann_file)) stop("Missing d=187 annotations: ", ann_file)
cell <- readRDS(cell_file)
lo <- load_native_locus(locus, load_genotypes = FALSE)
A187 <- readRDS(ann_file)
if (!is.matrix(A187) || !is.numeric(A187) || ncol(A187) != 187L ||
    !identical(rownames(A187), lo$snp) || any(!is.finite(A187))) {
  stop("Malformed or misaligned d=187 annotation matrix for locus ", locus)
}
if (!identical(as.character(cell$truth$SNP), lo$snp) ||
    !identical(unname(cell$availability), unname(lo$avail))) {
  stop("Saved native cell is not aligned to the current locus inputs")
}

R_list <- pad_ld_unavailable(lo$R_union, lo$avail)
ss <- make_sumstats(cell$zfile)
set.seed(20260903L + locus + 1009L * base_index + 7919L * h2_index +
           104729L * n_assv_each)
t0 <- proc.time()[["elapsed"]]
invisible(capture.output(
  fit <- INFUSE::infuse_core_avail(
    R_mat_list = R_list, summary_stat_list = ss, L = L_effects,
    variant_available = lo$avail, annot = A187, annot_method = "mlk",
    optim_method = "optim", est_annot_prior = "fixed",
    estimate_residual_variance = FALSE, max_iter = 100L, verbose = FALSE
  )
))
elapsed_minutes <- (proc.time()[["elapsed"]] - t0) / 60

cp <- config_pips(fit$pip_config, fit$availability$config_names)
res <- list(
  metadata = list(locus = locus, base_count_index = base_index,
                  n_base = NATIVE_BASE_COUNTS[base_index], h2_index = h2_index,
                  h2_per_causal = NATIVE_H2[h2_index], n_assv_each = n_assv_each,
                  L = L_effects, n_annotations = 187L,
                  annot_method = "mlk", est_annot_prior = "fixed",
                  variant_set = "union"),
  pip = as.numeric(fit$pip), pip_EUR = as.numeric(cp$EUR),
  pip_AFR = as.numeric(cp$AFR), pip_shared = as.numeric(cp$shared),
  cs = fit$cs$cs, timing_minutes = elapsed_minutes,
  annotation_names = colnames(A187)
)
if (length(res$pip) != lo$p || any(!is.finite(res$pip)) ||
    any(res$pip < -1e-12 | res$pip > 1 + 1e-12)) {
  stop(sprintf(paste0("Invalid d=187 PIP output: length=%d expected=%d; ",
                      "nonfinite=%d; range=[%.17g, %.17g]"),
               length(res$pip), lo$p, sum(!is.finite(res$pip)),
               suppressWarnings(min(res$pip, na.rm = TRUE)),
               suppressWarnings(max(res$pip, na.rm = TRUE))))
}
# Clip only machine-precision excursions at probability boundaries.
res$pip <- pmin(1, pmax(0, res$pip))
if (any(res$pip_shared[rowSums(lo$avail) < 2L] != 0) ||
    any(res$pip_EUR[!lo$avail[, "EUR"]] != 0) ||
    any(res$pip_AFR[!lo$avail[, "AFR"]] != 0)) {
  stop("d=187 result has posterior mass outside the availability support")
}

dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)
tmp <- paste0(out_file, ".partial_", Sys.getpid())
saveRDS(res, tmp, compress = FALSE)
if (!file.rename(tmp, out_file)) stop("Atomic output rename failed: ", out_file)
cat("wrote", out_file, "runtime_minutes=", elapsed_minutes, "\n")
