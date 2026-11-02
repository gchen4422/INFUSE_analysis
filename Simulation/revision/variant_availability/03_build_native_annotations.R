#!/usr/bin/env Rscript
.infuse_repo <- Sys.getenv("INFUSE_ANALYSIS_ROOT", "/path/to/INFUSE_analysis")

this_file <- normalizePath(sub("^--file=", "",
  grep("^--file=", commandArgs(FALSE), value = TRUE)[1]), mustWork = TRUE)
source(file.path(.infuse_repo, "Simulation/revision/variant_availability/native_assv_config.R"))
suppressMessages(library(data.table))

regions <- fread(REGION_FILE)
if (sum(file.exists(file.path(PANEL_ROOT, sprintf("locus_%03d", seq_len(N_LOCI)),
                              "COMPLETE"))) != N_LOCI) {
  stop("All native panels must be complete before annotation extraction.")
}

qc_out <- list()
for (chr in sort(unique(regions$chromosome))) {
  rr <- regions[chromosome == chr]
  source_file <- file.path(ANNOT_SOURCE_DIR,
                           sprintf("baselineLF2.2.UKB.%d.annot.gz", chr))
  if (!file.exists(source_file)) stop("Missing annotation source: ", source_file)

  range_expr <- paste(sprintf("($2 >= %d && $2 <= %d)", rr$start_bp, rr$end_bp),
                      collapse = " || ")
  cmd <- paste("gzip -dc", shQuote(source_file), "| awk -F '\\t'",
               shQuote(sprintf("NR == 1 || ($1 == %d && (%s))", chr, range_expr)))
  ann_chr <- fread(cmd = cmd)
  missing_cols <- setdiff(c("CHR", "BP", ANNOT_FIT_COLS), names(ann_chr))
  if (length(missing_cols)) stop("Annotation columns absent: ", paste(missing_cols, collapse = ", "))
  ann_chr <- ann_chr[!duplicated(BP)]

  for (locus_id in rr$locus) {
    locus_dir <- file.path(PANEL_ROOT, sprintf("locus_%03d", locus_id))
    vm <- fread(file.path(locus_dir, "variant_map.tsv"))
    hit <- match(vm$position, ann_chr$BP)
    ok <- !is.na(hit)
    # Follow the established INFUSE simulation: a variant without a row in the
    # baselineUKB2.2 annotation panel is excluded from the simulation universe.
    # It must not be encoded as an all-zero annotation vector because "not in the
    # annotation resource" is not evidence that every annotation is truly zero.
    kept_vm <- vm[ok]
    A <- as.matrix(ann_chr[hit[ok], ..ANNOT_FIT_COLS])
    dimnames(A) <- list(kept_vm$SNP, ANNOT_FIT_COLS)
    if (any(!is.finite(A))) stop("Non-finite annotation in locus ", locus_id)

    out <- data.table(SNP = kept_vm$SNP, A)
    fwrite(out, file.path(locus_dir, "annotation.tsv"), sep = "\t")
    saveRDS(A, file.path(locus_dir, "annotation.rds"), compress = FALSE)
    audit <- vm[, .(SNP, availability_pattern, annotation_matched = ok)]
    fwrite(audit, file.path(locus_dir, "annotation_match.tsv"), sep = "\t")
    if (any(!ok)) {
      fwrite(audit[annotation_matched == FALSE],
             file.path(locus_dir, "annotation_excluded.tsv"), sep = "\t")
    }
    count_pattern <- function(x, pattern) sum(x & vm$availability_pattern == pattern)
    qc_out[[length(qc_out) + 1L]] <- data.table(
      locus = locus_id, chromosome = chr, n_union = nrow(vm),
      n_annotation_matched = sum(ok), n_annotation_excluded = sum(!ok),
      kept_shared = count_pattern(ok, "EUR+AFR"),
      kept_EUR_only = count_pattern(ok, "EUR"),
      kept_AFR_only = count_pattern(ok, "AFR"),
      excluded_shared = count_pattern(!ok, "EUR+AFR"),
      excluded_EUR_only = count_pattern(!ok, "EUR"),
      excluded_AFR_only = count_pattern(!ok, "AFR")
    )
  }
  rm(ann_chr); invisible(gc())
}

qc <- rbindlist(qc_out)
setorder(qc, locus)
fwrite(qc, file.path(NATIVE_DIR, "metadata", "annotation_qc_summary.tsv"), sep = "\t")
cat("annotation-complete extraction:", nrow(qc), "loci; excluded variants:",
    sum(qc$n_annotation_excluded), "\n")
