#!/usr/bin/env Rscript
.infuse_repo <- Sys.getenv("INFUSE_ANALYSIS_ROOT", "/path/to/INFUSE_analysis")
.infuse_work <- Sys.getenv("INFUSE_REVISION", "/path/to/revision_output")

# Raw, signed Evo 2 external validation for BMI. This is a clean,
# self-contained replacement for the corresponding source plotting script.

module_root <- file.path(.infuse_work, "evo2_raw_score_external_validation")
source(file.path(.infuse_repo, "Real_data/revision/evo2/00_evo2_enrichment_common.R"))

bmi <- load_raw_evo2_result(
  file.path(module_root, "data", "res_BMI_snp_embedding_updated_raw_evo2.RData"),
  "BMI"
)

run_raw_evo2_validation(bmi, "bmi", module_root)
