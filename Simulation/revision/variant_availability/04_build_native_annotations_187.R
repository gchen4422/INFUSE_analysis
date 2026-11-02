#!/usr/bin/env Rscript
.infuse_repo <- Sys.getenv("INFUSE_ANALYSIS_ROOT", "/path/to/INFUSE_analysis")

args <- commandArgs(TRUE)
if (length(args) != 1L) stop("Usage: 04_build_native_annotations_187.R <locus>")
locus <- as.integer(args[[1L]])
if (is.na(locus) || !locus %in% 1:100) stop("locus must be in 1:100")

this_file <- normalizePath(sub("^--file=", "",
  grep("^--file=", commandArgs(FALSE), value = TRUE)[1L]), mustWork = TRUE)
NATIVE_SRC_DIR <- dirname(this_file)
source(file.path(.infuse_repo, "Simulation/revision/variant_availability/native_sim_common.R"))
suppressMessages(library(data.table))

lo <- load_native_locus(locus, load_genotypes = FALSE)
chr <- unique(lo$vm$chromosome)
if (length(chr) != 1L) stop("Expected exactly one chromosome for locus ", locus)

source_file <- file.path(ANNOT_SOURCE_DIR,
                         sprintf("baselineLF2.2.UKB.%d.annot.gz", chr))
if (!file.exists(source_file)) stop("Missing annotation source: ", source_file)
range_expr <- sprintf("NR == 1 || ($1 == %d && $2 >= %d && $2 <= %d)",
                      chr, min(lo$vm$position), max(lo$vm$position))
cmd <- paste("gzip -dc", shQuote(source_file), "| awk -F '\\t'", shQuote(range_expr))
ann <- fread(cmd = cmd)
if (ncol(ann) != 191L) {
  stop("Expected 191 source columns (4 metadata + 187 annotations); got ", ncol(ann))
}
# Preserve the exact rule used by 04_build_native_annotations.R. The source can
# contain more than one record at a BP; the established native d=10 universe was
# built by retaining the first source row at each position.
ann <- ann[!duplicated(BP)]

hit <- match(lo$vm$position, ann$BP)
if (anyNA(hit)) {
  stop("The annotation-complete native universe contains a position absent from the full source")
}
annotation_names <- names(ann)[5:191]
if (anyNA(annotation_names) || any(!nzchar(annotation_names)) ||
    anyDuplicated(annotation_names)) stop("Invalid full annotation column names")
A187 <- as.matrix(ann[hit, 5:191, with = FALSE])
storage.mode(A187) <- "double"
dimnames(A187) <- list(lo$snp, annotation_names)
if (any(!is.finite(A187))) stop("Non-finite full annotation values in locus ", locus)

# The ten-dimensional primary panel must be an exact column subset of the
# high-dimensional matrix. This guards both SNP-row and annotation-column order.
if (!all(colnames(lo$annot) %in% colnames(A187))) {
  stop("The primary ten annotations are not all present in the 187-column panel")
}
delta10 <- max(abs(A187[, colnames(lo$annot), drop = FALSE] - lo$annot))
if (!identical(delta10, 0)) stop("d=10/d=187 annotation mismatch: max delta ", delta10)

out_file <- native_annotation187_path(locus)
tmp <- paste0(out_file, ".partial_", Sys.getpid())
saveRDS(A187, tmp, compress = FALSE)
if (!file.rename(tmp, out_file)) stop("Atomic annotation rename failed: ", out_file)

qc_dir <- file.path(NATIVE_DIR, "metadata", "annotation_187_qc")
dir.create(qc_dir, recursive = TRUE, showWarnings = FALSE)
qc <- data.frame(locus = locus, chromosome = chr, n_variants = nrow(A187),
                 n_annotations = ncol(A187), max_delta_primary10 = delta10)
write.table(qc, file.path(qc_dir, sprintf("locus_%03d.tsv", locus)),
            sep = "\t", quote = FALSE, row.names = FALSE)
cat("wrote", out_file, nrow(A187), "x", ncol(A187), "\n")
