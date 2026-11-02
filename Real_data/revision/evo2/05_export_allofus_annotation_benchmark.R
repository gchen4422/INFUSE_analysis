#!/usr/bin/env Rscript
.infuse_work <- Sys.getenv("INFUSE_REVISION", "/path/to/revision_output")

# Export the minimal per-SNP information needed to benchmark annotations that
# are available only inside the All of Us Researcher Workbench.
#
# Deliberately excluded: Evo 2 scores, other annotations, and all other methods.
# Coordinates and both allele orientations are retained for VAT matching. The
# two PIP columns are copied without transformation.

suppressPackageStartupMessages(library(data.table))

module_root <- file.path(.infuse_work, "evo2_raw_score_external_validation")
output_dir <- file.path(module_root, "allofus_upload")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

input_files <- c(
  SBP = file.path(module_root, "data", "res_sbp_snp_embedding_updated_raw_evo2.RData"),
  DBP = file.path(module_root, "data", "res_dbp_snp_embedding_updated_raw_evo2.RData"),
  BMI = file.path(module_root, "data", "res_BMI_snp_embedding_updated_raw_evo2.RData")
)

source_columns <- c(
  MESuSiE_PIP_Either = "MESuSiE_PIP_Either",
  INFUSE_Classic_Annotations_PIP_Either = "MESuSiE_PIP_Either_all_mlk"
)

export_one <- function(trait, input_file) {
  if (!file.exists(input_file)) {
    stop("Missing input for ", trait, ": ", input_file, call. = FALSE)
  }

  e <- new.env(parent = emptyenv())
  loaded <- load(input_file, envir = e)
  if (!identical(loaded, "res_all_new")) {
    stop("Expected only res_all_new in ", input_file, "; found: ",
         paste(loaded, collapse = ", "), call. = FALSE)
  }
  d <- as.data.table(e$res_all_new)
  rm(e)

  required <- c("Region", "SNP", "CHR", "POS", "REF", "ALT", unname(source_columns))
  missing_columns <- setdiff(required, names(d))
  if (length(missing_columns)) {
    stop(trait, " input is missing: ", paste(missing_columns, collapse = ", "),
         call. = FALSE)
  }

  out <- d[, .(
    Trait = trait,
    Region = as.integer(Region),
    SNP = as.character(SNP),
    CHR = as.character(CHR),
    POS = as.integer(POS),
    REF = toupper(as.character(REF)),
    ALT = toupper(as.character(ALT)),
    vid_REF_ALT = paste(sub("^chr", "", as.character(CHR), ignore.case = TRUE),
                        as.integer(POS), toupper(as.character(REF)),
                        toupper(as.character(ALT)), sep = "-"),
    vid_ALT_REF = paste(sub("^chr", "", as.character(CHR), ignore.case = TRUE),
                        as.integer(POS), toupper(as.character(ALT)),
                        toupper(as.character(REF)), sep = "-"),
    MESuSiE_PIP_Either = as.numeric(get(source_columns[["MESuSiE_PIP_Either"]])),
    INFUSE_Classic_Annotations_PIP_Either = as.numeric(get(source_columns[["INFUSE_Classic_Annotations_PIP_Either"]]))
  )]

  if (anyNA(out$SNP) || any(!nzchar(out$SNP))) {
    stop(trait, " contains missing or empty SNP identifiers.", call. = FALSE)
  }
  if (anyDuplicated(out$SNP)) {
    stop(trait, " contains duplicate SNP identifiers; SNP alone is not a unique key.",
         call. = FALSE)
  }
  if (any(!grepl("^rs[0-9]+$", out$SNP))) {
    stop(trait, " contains SNP identifiers that are not rsIDs.", call. = FALSE)
  }
  if (anyNA(out$Region) || any(out$Region <= 0L)) {
    stop(trait, " contains invalid Region values.", call. = FALSE)
  }
  if (anyNA(out$CHR) || any(!nzchar(out$CHR)) || anyNA(out$POS) || any(out$POS <= 0L)) {
    stop(trait, " contains invalid chromosome or position values.", call. = FALSE)
  }
  if (out[, anyDuplicated(paste(CHR, POS, sep = ":"))] != 0L) {
    stop(trait, " contains duplicate CHR:POS coordinates.", call. = FALSE)
  }
  if (anyNA(out$REF) || anyNA(out$ALT) || any(!nzchar(out$REF)) ||
      any(!nzchar(out$ALT)) || any(out$REF == out$ALT)) {
    stop(trait, " contains invalid REF/ALT alleles.", call. = FALSE)
  }
  expected_vid_ref_alt <- paste(sub("^chr", "", out$CHR, ignore.case = TRUE),
                                out$POS, out$REF, out$ALT, sep = "-")
  expected_vid_alt_ref <- paste(sub("^chr", "", out$CHR, ignore.case = TRUE),
                                out$POS, out$ALT, out$REF, sep = "-")
  if (!identical(out$vid_REF_ALT, expected_vid_ref_alt) ||
      !identical(out$vid_ALT_REF, expected_vid_alt_ref) ||
      any(grepl("^chr", out$vid_REF_ALT, ignore.case = TRUE)) ||
      any(grepl("^chr", out$vid_ALT_REF, ignore.case = TRUE))) {
    stop(trait, " VAT vid columns are inconsistent.", call. = FALSE)
  }
  metric_columns <- names(source_columns)
  for (column in metric_columns) {
    value <- out[[column]]
    if (any(!is.finite(value))) {
      stop(trait, " has non-finite values in ", column, call. = FALSE)
    }
    if (any(value < 0 | value > 1)) {
      stop(trait, " has values outside [0,1] in ", column, call. = FALSE)
    }
  }

  setorder(out, Region, CHR, POS, SNP)
  output_file <- file.path(
    output_dir,
    paste0(tolower(trait), "_snp_coordinates_mesusie_infuse_classic_pip.tsv.gz")
  )
  fwrite(out, output_file, sep = "\t", quote = FALSE, compress = "gzip")

  data.table(
    trait = trait,
    input_file = input_file,
    output_file = output_file,
    n_rows = nrow(out),
    n_unique_snps = uniqueN(out$SNP),
    n_columns = ncol(out),
    mesusie_pip_min = min(out$MESuSiE_PIP_Either),
    mesusie_pip_max = max(out$MESuSiE_PIP_Either),
    infuse_classic_annotations_pip_min = min(out$INFUSE_Classic_Annotations_PIP_Either),
    infuse_classic_annotations_pip_max = max(out$INFUSE_Classic_Annotations_PIP_Either),
    snp_identifier = "dbSNP rsID",
    transformed = FALSE
  )
}

manifest <- rbindlist(Map(export_one, names(input_files), input_files))
fwrite(
  manifest,
  file.path(output_dir, "MANIFEST.tsv"),
  sep = "\t",
  quote = FALSE
)

message("Wrote upload-ready tables to: ", output_dir)
print(manifest)
