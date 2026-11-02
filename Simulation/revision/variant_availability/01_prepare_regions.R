#!/usr/bin/env Rscript
.infuse_repo <- Sys.getenv("INFUSE_ANALYSIS_ROOT", "/path/to/INFUSE_analysis")

this_file <- normalizePath(sub("^--file=", "",
  grep("^--file=", commandArgs(FALSE), value = TRUE)[1]), mustWork = TRUE)
source(file.path(.infuse_repo, "Simulation/revision/variant_availability/native_assv_config.R"))
stop_if_missing_sources()

dir.create(dirname(REGION_FILE), recursive = TRUE, showWarnings = FALSE)
dir.create(PANEL_ROOT, recursive = TRUE, showWarnings = FALSE)
dir.create(WORK_ROOT, recursive = TRUE, showWarnings = FALSE)
dir.create(LOG_ROOT, recursive = TRUE, showWarnings = FALSE)

regions <- read.table(SOURCE_REGION_FILE, header = TRUE, stringsAsFactors = FALSE)
stopifnot(nrow(regions) == N_LOCI, all(c("CHR", "MinBP", "MaxBP") %in% names(regions)))
if (any(!is.finite(regions$CHR)) || any(!is.finite(regions$MinBP)) ||
    any(!is.finite(regions$MaxBP)) || any(regions$MinBP >= regions$MaxBP)) {
  stop("Invalid source region coordinates.")
}

out <- data.frame(
  locus = seq_len(nrow(regions)),
  chromosome = as.integer(regions$CHR),
  start_bp = as.integer(regions$MinBP),
  end_bp = as.integer(regions$MaxBP),
  source_region_file = SOURCE_REGION_FILE,
  source_EUR_bfile = unname(SOURCE_BFILE["EUR"]),
  source_AFR_bfile = unname(SOURCE_BFILE["AFR"]),
  locus_maf_filter = "none (preserve source panel)",
  locus_genotype_missingness_filter = "none (preserve source panel)",
  stringsAsFactors = FALSE
)
data.table::fwrite(out, REGION_FILE, sep = "\t")
cat("wrote", REGION_FILE, "with", nrow(out), "regions\n")
