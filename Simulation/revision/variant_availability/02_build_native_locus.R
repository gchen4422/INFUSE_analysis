#!/usr/bin/env Rscript
.infuse_repo <- Sys.getenv("INFUSE_ANALYSIS_ROOT", "/path/to/INFUSE_analysis")

args <- commandArgs(TRUE)
if (!length(args)) stop("Usage: Rscript 02_build_native_locus.R <locus>")
locus_id <- as.integer(args[1])
if (!is.finite(locus_id) || locus_id < 1L || locus_id > 100L) stop("locus must be 1..100")

this_file <- normalizePath(sub("^--file=", "",
  grep("^--file=", commandArgs(FALSE), value = TRUE)[1]), mustWork = TRUE)
source(file.path(.infuse_repo, "Simulation/revision/variant_availability/native_assv_config.R"))
stop_if_missing_sources()

suppressMessages({
  library(data.table)
  library(snpStats)
})

if (!file.exists(REGION_FILE)) stop("Run 01_prepare_regions.R first: missing ", REGION_FILE)
regions <- fread(REGION_FILE)
reg <- regions[regions$locus == locus_id, ]
stopifnot(nrow(reg) == 1L)

locus_dir <- file.path(PANEL_ROOT, sprintf("locus_%03d", locus_id))
complete_file <- file.path(locus_dir, "COMPLETE")
force <- identical(Sys.getenv("FORCE", "0"), "1")
if (file.exists(complete_file) && !force) {
  cat("locus", locus_id, "already complete; leaving outputs untouched\n")
  quit(status = 0L)
}
if (dir.exists(locus_dir) && !force) {
  stop("Partial locus output exists without COMPLETE: ", locus_dir,
       ". Inspect it, or rerun deliberately with FORCE=1.")
}
if (dir.exists(locus_dir) && force) unlink(locus_dir, recursive = TRUE)

task_id <- paste0("locus_", sprintf("%03d", locus_id), "_", Sys.getpid())
work <- file.path(WORK_ROOT, task_id)
dir.create(work, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(locus_dir, "plink"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(locus_dir, "ld"), recursive = TRUE, showWarnings = FALSE)

run_plink <- function(args, label) {
  log <- file.path(work, paste0(label, ".stdout.log"))
  err <- file.path(work, paste0(label, ".stderr.log"))
  st <- system2(PLINK2_BIN, args, stdout = log, stderr = err)
  if (!identical(st, 0L)) stop(label, " failed with exit ", st, "; see ", log, " and ", err)
}

read_bim <- function(prefix) {
  x <- fread(paste0(prefix, ".bim"), header = FALSE,
             col.names = c("CHR", "SNP", "CM", "BP", "A1", "A2"))
  x[, `:=`(SNP = as.character(SNP), A1 = toupper(as.character(A1)),
           A2 = toupper(as.character(A2)))]
  if (anyDuplicated(x$SNP)) stop("Duplicated SNP IDs in ", prefix)
  x
}

chr <- as.integer(reg$chromosome)
start_bp <- as.integer(reg$start_bp)
end_bp <- as.integer(reg$end_bp)

raw_prefix <- setNames(file.path(work, paste0("raw_", ANCESTRIES)), ANCESTRIES)
for (anc in ANCESTRIES) {
  run_plink(c(
    "--bfile", unname(SOURCE_BFILE[anc]),
    "--chr", as.character(chr), "--from-bp", as.character(start_bp),
    "--to-bp", as.character(end_bp),
    "--make-bed", "--threads", Sys.getenv("SLURM_CPUS_PER_TASK", "2"),
    "--out", raw_prefix[anc]
  ), paste0("extract_", anc))
}

raw_bim <- lapply(raw_prefix, read_bim)
if (any(vapply(raw_bim, nrow, integer(1)) < 10L)) {
  stop("Fewer than 10 QC-passing variants in an ancestry")
}

shared <- intersect(raw_bim$EUR$SNP, raw_bim$AFR$SNP)
ie <- match(shared, raw_bim$EUR$SNP); ia <- match(shared, raw_bim$AFR$SNP)
exact <- raw_bim$EUR$A1[ie] == raw_bim$AFR$A1[ia] &
         raw_bim$EUR$A2[ie] == raw_bim$AFR$A2[ia]
swapped <- raw_bim$EUR$A1[ie] == raw_bim$AFR$A2[ia] &
           raw_bim$EUR$A2[ie] == raw_bim$AFR$A1[ia]
compatible <- exact | swapped
incompatible_ids <- shared[!compatible]
shared_compatible <- shared[compatible]
shared_swapped <- shared[compatible & swapped]

keep <- list(
  EUR = setdiff(raw_bim$EUR$SNP, incompatible_ids),
  AFR = setdiff(raw_bim$AFR$SNP, incompatible_ids)
)
final_prefix <- setNames(file.path(locus_dir, "plink", ANCESTRIES), ANCESTRIES)
for (anc in ANCESTRIES) {
  keep_file <- file.path(work, paste0("keep_", anc, ".txt"))
  fwrite(data.table(SNP = keep[[anc]]), keep_file, col.names = FALSE)
  run_plink(c("--bfile", raw_prefix[anc], "--extract", keep_file,
              "--make-bed", "--threads", Sys.getenv("SLURM_CPUS_PER_TASK", "2"),
              "--out", final_prefix[anc]), paste0("finalize_", anc))
}

final_bim <- lapply(final_prefix, read_bim)
union_ids <- c(final_bim$EUR$SNP, setdiff(final_bim$AFR$SNP, final_bim$EUR$SNP))
idx_e <- match(union_ids, final_bim$EUR$SNP)
idx_a <- match(union_ids, final_bim$AFR$SNP)

variant_map <- data.table(
  union_index = seq_along(union_ids), SNP = union_ids,
  chromosome = chr,
  position = as.integer(sub("^[^:]+:", "", union_ids)),
  available_EUR = !is.na(idx_e), available_AFR = !is.na(idx_a),
  index_EUR = idx_e, index_AFR = idx_a,
  A1_EUR = ifelse(is.na(idx_e), NA_character_, final_bim$EUR$A1[idx_e]),
  A2_EUR = ifelse(is.na(idx_e), NA_character_, final_bim$EUR$A2[idx_e]),
  A1_AFR = ifelse(is.na(idx_a), NA_character_, final_bim$AFR$A1[idx_a]),
  A2_AFR = ifelse(is.na(idx_a), NA_character_, final_bim$AFR$A2[idx_a])
)
variant_map[, availability_pattern := fifelse(available_EUR & available_AFR, "EUR+AFR",
                                        fifelse(available_EUR, "EUR", "AFR"))]
variant_map[, flip_AFR_to_EUR := SNP %in% shared_swapped]

standardized_genotype <- function(prefix, anc) {
  pl <- read.plink(paste0(prefix, ".bed"))
  G <- as(pl$genotypes, "numeric")
  ids <- colnames(pl$genotypes)
  target <- final_bim[[anc]]$SNP
  hit <- match(target, ids)
  if (anyNA(hit)) stop("Genotype/BIM mismatch in ", anc)
  G <- G[, hit, drop = FALSE]
  cm <- colMeans(G, na.rm = TRUE)
  miss <- which(is.na(G), arr.ind = TRUE)
  if (nrow(miss)) G[miss] <- cm[miss[, 2]]
  G <- sweep(G, 2, colMeans(G), "-")
  s <- sqrt(colSums(G^2) / (nrow(G) - 1))
  if (any(!is.finite(s) | s <= 0)) stop("Zero/non-finite post-QC genotype SD in ", anc)
  G <- sweep(G, 2, s, "/")
  colnames(G) <- target
  if (anc == "AFR" && length(shared_swapped)) {
    j <- match(shared_swapped, colnames(G)); j <- j[!is.na(j)]
    G[, j] <- -G[, j, drop = FALSE]
  }
  G
}

ld_summary <- list()
for (anc in ANCESTRIES) {
  G <- standardized_genotype(final_prefix[anc], anc)
  R <- cov2cor(crossprod(G))
  R <- (R + t(R)) / 2
  diag(R) <- 1
  dimnames(R) <- list(colnames(G), colnames(G))
  if (any(!is.finite(R))) stop("Non-finite LD for ", anc)
  if (max(abs(R - t(R))) > 1e-12) stop("Nonsymmetric LD for ", anc)
  saveRDS(R, file.path(locus_dir, "ld", paste0(anc, ".rds")), compress = FALSE)
  ld_summary[[anc]] <- c(n_samples = nrow(G), n_variants = ncol(G),
                         max_asymmetry = max(abs(R - t(R))))
  rm(G, R); invisible(gc())
}

fwrite(variant_map, file.path(locus_dir, "variant_map.tsv"), sep = "\t", na = "NA")
if (length(incompatible_ids)) {
  fwrite(data.table(SNP = incompatible_ids),
         file.path(locus_dir, "incompatible_shared_alleles.tsv"), sep = "\t")
}

qc <- data.table(
  locus = locus_id, chromosome = chr, start_bp = start_bp, end_bp = end_bp,
  locus_maf_filter = "none", locus_genotype_missingness_filter = "none",
  raw_EUR = nrow(raw_bim$EUR), raw_AFR = nrow(raw_bim$AFR),
  final_EUR = nrow(final_bim$EUR), final_AFR = nrow(final_bim$AFR),
  shared = sum(variant_map$availability_pattern == "EUR+AFR"),
  EUR_only = sum(variant_map$availability_pattern == "EUR"),
  AFR_only = sum(variant_map$availability_pattern == "AFR"),
  union = nrow(variant_map), incompatible_shared_alleles = length(incompatible_ids),
  shared_allele_swaps = length(shared_swapped),
  n_EUR_samples = as.integer(ld_summary$EUR["n_samples"]),
  n_AFR_samples = as.integer(ld_summary$AFR["n_samples"])
)
fwrite(qc, file.path(locus_dir, "qc.tsv"), sep = "\t")

if (!file.create(complete_file)) {
  stop("Cannot create completion marker: ", complete_file, call. = FALSE)
}

unlink(work, recursive = TRUE)
cat("completed locus", locus_id, "union", nrow(variant_map),
    "shared", qc$shared, "EUR-only", qc$EUR_only, "AFR-only", qc$AFR_only, "\n")
