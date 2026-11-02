#!/usr/bin/env Rscript
.infuse_work <- Sys.getenv("INFUSE_REVISION", "/path/to/revision_output")

# Add the signed, untransformed Evo 2 delta score to a SNP-level result table.
# The source analysis directory is read-only. All outputs are written to the
# self-contained revision module specified below.

suppressPackageStartupMessages(library(data.table))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) {
  stop("Usage: Rscript 01_add_raw_evo2_scores.R <SBP|DBP|BMI>", call. = FALSE)
}

trait <- toupper(args[[1L]])
if (!trait %in% c("SBP", "DBP", "BMI")) {
  stop("Trait must be SBP, DBP, or BMI.", call. = FALSE)
}

source_root <- Sys.getenv("INFUSE_AOU_SCORES", "/path/to/aou_scores")
embedding_root <- file.path(source_root, "didong_embedding")
module_root <- file.path(.infuse_work, "evo2_raw_score_external_validation")

input_file <- switch(
  trait,
  SBP = file.path(embedding_root, "summary_res", "res_sbp_snp_embedding_updated.RData"),
  DBP = file.path(embedding_root, "summary_res", "res_dbp_snp_embedding_updated.RData"),
  BMI = file.path(embedding_root, "summary_res", "res_BMI_snp_embedding_updated.RData")
)

output_file <- switch(
  trait,
  SBP = file.path(module_root, "data", "res_sbp_snp_embedding_updated_raw_evo2.RData"),
  DBP = file.path(module_root, "data", "res_dbp_snp_embedding_updated_raw_evo2.RData"),
  BMI = file.path(module_root, "data", "res_BMI_snp_embedding_updated_raw_evo2.RData")
)

score_dir <- file.path(source_root, paste0("results_", trait), "summary_data")
score_files <- list.files(
  score_dir,
  pattern = "^loci_[0-9]+_summ_scored[.]csv$",
  full.names = TRUE
)

if (!file.exists(input_file)) stop("Missing input: ", input_file, call. = FALSE)
if (!length(score_files)) stop("No Evo 2 score files found in: ", score_dir, call. = FALSE)

region_from_file <- function(path) {
  as.integer(sub("^loci_([0-9]+)_.*$", "\\1", basename(path)))
}

message("Reading raw Evo 2 scores for ", trait, " from ", length(score_files), " loci.")
score_table <- rbindlist(
  Map(
    function(path, region) {
      ans <- fread(path, select = c("SNP", "ALT", "REF", "evo2_delta_score"))
      ans[, Region := region]
      ans[, .(
        Region,
        SNP = as.character(SNP),
        source_ALT = as.character(ALT),
        source_REF = as.character(REF),
        evo2_delta_score
      )]
    },
    score_files,
    vapply(score_files, region_from_file, integer(1))
  ),
  use.names = TRUE
)

if (anyNA(score_table$evo2_delta_score)) {
  stop("The raw Evo 2 source contains missing scores.", call. = FALSE)
}
if (score_table[, anyDuplicated(paste(Region, SNP, sep = "\r"))] != 0L) {
  stop("The raw Evo 2 source contains duplicate Region/SNP keys.", call. = FALSE)
}

load_env <- new.env(parent = emptyenv())
loaded <- load(input_file, envir = load_env)
if (!identical(loaded, "res_all_new")) {
  stop("Expected the input RData to contain only res_all_new; found: ",
       paste(loaded, collapse = ", "), call. = FALSE)
}

res_all_new <- as.data.table(load_env$res_all_new)
rm(load_env)

required <- c("Region", "SNP", "ALT", "REF")
missing_required <- setdiff(required, names(res_all_new))
if (length(missing_required)) {
  stop("Input is missing key columns: ", paste(missing_required, collapse = ", "), call. = FALSE)
}
if ("evo2_delta_score" %in% names(res_all_new)) {
  stop("Input already contains evo2_delta_score; refusing to overwrite it.", call. = FALSE)
}
if (res_all_new[, anyDuplicated(paste(Region, SNP, sep = "\r"))] != 0L) {
  stop("Input contains duplicate Region/SNP keys.", call. = FALSE)
}

# This keyed lookup preserves the exact result-table row order. The assigned
# value is the source evo2_delta_score itself: no absolute value, logarithm,
# scaling, ranking, normalization, clipping, or other transformation is used.
raw_score <- score_table[
  res_all_new,
  on = .(Region, SNP),
  evo2_delta_score
]
source_alt <- score_table[res_all_new, on = .(Region, SNP), source_ALT]
source_ref <- score_table[res_all_new, on = .(Region, SNP), source_REF]

if (length(raw_score) != nrow(res_all_new)) {
  stop("Internal join error: joined score length differs from result rows.", call. = FALSE)
}
if (anyNA(raw_score)) {
  bad <- which(is.na(raw_score))
  example <- paste(
    paste(res_all_new$Region[head(bad, 5L)], res_all_new$SNP[head(bad, 5L)], sep = ":"),
    collapse = ", "
  )
  stop("Some result SNPs lack a raw Evo 2 score. Examples: ", example, call. = FALSE)
}
allele_match <- toupper(as.character(res_all_new$ALT)) == toupper(source_alt) &
  toupper(as.character(res_all_new$REF)) == toupper(source_ref)
if (anyNA(allele_match) || !all(allele_match)) {
  stop("ALT/REF mismatch between the result table and raw Evo 2 scoring input.",
       call. = FALSE)
}

res_all_new[, evo2_delta_score := raw_score]

# Exact equality check against a second keyed lookup guards against accidental
# modification of the signed source values.
check_score <- score_table[
  res_all_new,
  on = .(Region, SNP),
  evo2_delta_score
]
if (!identical(res_all_new$evo2_delta_score, check_score)) {
  stop("Raw Evo 2 equality check failed.", call. = FALSE)
}

dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
save(res_all_new, file = output_file, compress = "gzip")

input_keys <- unique(res_all_new[, .(Region, SNP)])
source_only <- score_table[!input_keys, on = .(Region, SNP), .N]
audit <- data.table(
  trait = trait,
  input_file = input_file,
  output_file = output_file,
  score_directory = score_dir,
  n_score_files = length(score_files),
  n_result_rows = nrow(res_all_new),
  n_result_regions = uniqueN(res_all_new$Region),
  n_result_snps = uniqueN(res_all_new$SNP),
  n_missing_scores = sum(is.na(res_all_new$evo2_delta_score)),
  n_alt_ref_mismatches = sum(!allele_match),
  n_source_rows_not_in_results = source_only,
  raw_score_min = min(res_all_new$evo2_delta_score),
  raw_score_max = max(res_all_new$evo2_delta_score),
  raw_score_mean = mean(res_all_new$evo2_delta_score),
  raw_score_median = median(res_all_new$evo2_delta_score),
  exact_source_equality = TRUE,
  transformation = "none_signed_raw_evo2_delta_score"
)

audit_file <- file.path(module_root, "results", paste0("raw_evo2_merge_audit_", tolower(trait), ".tsv"))
fwrite(audit, audit_file, sep = "\t")

message("Wrote: ", output_file)
message("Audit: ", audit_file)
