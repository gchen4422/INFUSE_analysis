#!/usr/bin/env Rscript
.infuse_repo <- Sys.getenv("INFUSE_ANALYSIS_ROOT", "/path/to/INFUSE_analysis")

args <- commandArgs(TRUE)
if (length(args) != 1L) {
  stop("Usage: Rscript 06_prepare_multisusie_maf.R <locus>")
}
locus <- as.integer(args[[1L]])
if (is.na(locus) || !locus %in% 1:100) stop("locus must be in 1:100")

this_file <- normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE),
                                                    value = TRUE)[1L]),
                           mustWork = TRUE)
NATIVE_SRC_DIR <- dirname(this_file)
source(file.path(.infuse_repo, "Simulation/revision/variant_availability/native_sim_common.R"))

out_file <- native_multisusie_maf_path(locus)
if (file.exists(out_file)) {
  cat("complete MAF cache exists; leaving untouched:", out_file, "\n")
  quit(status = 0L)
}
dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)

lo <- load_native_locus(locus, load_genotypes = FALSE)
shared_ids <- lo$snp[rowSums(lo$avail) == length(ANCESTRY)]
if (length(shared_ids) < 10L) stop("Too few shared variants at locus ", locus)

maf <- setNames(vector("list", length(ANCESTRY)), ANCESTRY)
for (anc in ANCESTRY) {
  prefix <- file.path(lo$dir, "plink", anc)
  pl <- snpStats::read.plink(paste0(prefix, ".bed"))
  hit <- match(shared_ids, as.character(pl$map$snp.name))
  if (anyNA(hit) || anyDuplicated(hit)) {
    stop("Shared-SNP/PLINK alignment failure for ", anc, " at locus ", locus)
  }
  cs <- snpStats::col.summary(pl$genotypes[, hit, drop = FALSE])
  maf[[anc]] <- setNames(as.numeric(cs$MAF), shared_ids)
  if (any(!is.finite(maf[[anc]]) | maf[[anc]] <= 0 | maf[[anc]] > 0.5)) {
    stop("Invalid shared-SNP MAF for ", anc, " at locus ", locus)
  }
  rm(pl, cs)
  invisible(gc())
}

obj <- list(locus = locus, snp = shared_ids, maf = maf,
            n_shared = length(shared_ids), ancestry = ANCESTRY)
tmp <- paste0(out_file, ".partial_", Sys.getpid())
saveRDS(obj, tmp, compress = FALSE)
if (!file.rename(tmp, out_file)) stop("Atomic MAF-cache rename failed: ", out_file)
cat("wrote", out_file, "with", length(shared_ids), "shared SNPs\n")
