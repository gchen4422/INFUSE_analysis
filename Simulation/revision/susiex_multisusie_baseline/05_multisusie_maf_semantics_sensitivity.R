#!/usr/bin/env Rscript
.infuse_repo <- Sys.getenv("INFUSE_ANALYSIS_ROOT", "/path/to/INFUSE_analysis")

this_file <- normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1L]),
                           mustWork = TRUE)
source(file.path(.infuse_repo, "Simulation/revision/susiex_multisusie_baseline/00_common.R"))
suppressPackageStartupMessages(library(reticulate))

# The historical notebook passed ALT-frequency vectors through an argument named
# maf_list. At N=300,000 all variants are expected to pass the MAC filter under
# either representation. Verify that correcting the semantics does not alter the fit.
cell <- decode_task(1L)
z <- read_z(cell, "d10")
reference <- safe_read_method_object(cell_output("multisusie", cell), z, "MultiSuSiE")
if (is.null(reference)) stop("Validated task-1 MultiSuSiE result is required", call. = FALSE)

R_list <- list(read_ld(cell, "EUR", z$RSID), read_ld(cell, "AFR", z$RSID))
freq_eur <- read_true_maf(cell$locus, "EUR", nrow(z))
freq_afr <- read_true_maf(cell$locus, "AFR", nrow(z))
raw_frequency_list <- list(freq_eur$alt_frequency, freq_afr$alt_frequency)

use_python(PYTHON_BIN, required = TRUE)
ms <- import("MultiSuSiE", convert = FALSE)
np <- import("numpy", convert = FALSE)
fit <- ms$multisusie_rss(
  z_list = lapply(list(z$zscore_1, z$zscore_2), np$array),
  R_list = lapply(R_list, np$array),
  rho = np$array(matrix(c(1, 0.8, 0.8, 1), 2L, byrow = TRUE)),
  population_sizes = np$array(as.integer(N_GWAS)),
  L = as.integer(10L), scaled_prior_variance = 0.2,
  low_memory_mode = FALSE, min_abs_corr = 0.5,
  single_population_mac_thresh = as.integer(20L),
  maf_list = lapply(raw_frequency_list, np$array), coverage = 0.95
)
raw_pip <- as.numeric(py_to_r(fit$pip))
max_error <- max(abs(raw_pip - reference$pip))
out <- data.frame(
  task_id = 1L, n_snp = nrow(z),
  eur_alt_frequencies_above_half = freq_eur$n_above_half,
  afr_alt_frequencies_above_half = freq_afr$n_above_half,
  minimum_mac_using_true_maf = min(c(freq_eur$maf, freq_afr$maf) * 2 * N_GWAS[1L]),
  max_abs_pip_difference_raw_alt_frequency_vs_true_maf = max_error,
  pass = max_error < 1e-10
)
fwrite(out, file.path(RESULT_DIR, "multisusie_maf_semantics_sensitivity.tsv"), sep = "\t")
print(out)
if (!out$pass) stop("Raw ALT-frequency versus true-MAF sensitivity check changed the fit", call. = FALSE)

