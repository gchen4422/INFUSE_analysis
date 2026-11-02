#!/usr/bin/env Rscript
.infuse_repo <- Sys.getenv("INFUSE_ANALYSIS_ROOT", "/path/to/INFUSE_analysis")
.infuse_work <- Sys.getenv("INFUSE_REVISION", "/path/to/revision_output")

# Raw, signed Evo 2 external validation for SBP and DBP. This is a clean,
# self-contained replacement for the corresponding source plotting script.

module_root <- file.path(.infuse_work, "evo2_raw_score_external_validation")
source(file.path(.infuse_repo, "Real_data/revision/evo2/00_evo2_enrichment_common.R"))

sbp <- load_raw_evo2_result(
  file.path(module_root, "data", "res_sbp_snp_embedding_updated_raw_evo2.RData"),
  "SBP"
)
dbp <- load_raw_evo2_result(
  file.path(module_root, "data", "res_dbp_snp_embedding_updated_raw_evo2.RData"),
  "DBP"
)

# Pool SBP and DBP before ranking variants and calculating enrichment, following
# the combined BP analysis used in 9_plot_enrichment_either_bp.R. The final
# analysis treats BP as one outcome, so overlapping SBP and DBP regions are
# assigned the same genomic bootstrap cluster and are resampled together.
sbp[, Source_Trait := "SBP"]
dbp[, Source_Trait := "DBP"]
bp <- rbindlist(list(sbp, dbp), use.names = TRUE)
n_source_regions <- uniqueN(bp[, .(Source_Trait, Region)])
rm(sbp, dbp)
gc(verbose = FALSE)

clustered <- assign_overlapping_region_clusters(bp)
bp <- clustered$data
bp[, Trait := "BP"]

cluster_audit_file <- file.path(
  module_root, "results", "raw_evo2_bp_bootstrap_cluster_map.tsv"
)
dir.create(dirname(cluster_audit_file), recursive = TRUE, showWarnings = FALSE)
fwrite(clustered$cluster_audit, cluster_audit_file, sep = "\t")

message(
  "Pooled ", n_source_regions, " source trait-regions into ",
  uniqueN(bp$Region),
  " genomic bootstrap clusters; ",
  sum(clustered$cluster_audit$n_source_regions > 1L),
  " cluster(s) contain overlapping source regions."
)
message("Saved BP bootstrap cluster map: ", cluster_audit_file)
rm(clustered)

run_raw_evo2_validation(bp, "bp", module_root)
