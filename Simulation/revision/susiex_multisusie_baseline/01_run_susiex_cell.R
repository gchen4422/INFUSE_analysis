#!/usr/bin/env Rscript
.infuse_repo <- Sys.getenv("INFUSE_ANALYSIS_ROOT", "/path/to/INFUSE_analysis")

args <- commandArgs(TRUE)
if (length(args) != 1L) stop("Usage: 01_run_susiex_cell.R <task_id>", call. = FALSE)
this_file <- normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1L]),
                           mustWork = TRUE)
source(file.path(.infuse_repo, "Simulation/revision/susiex_multisusie_baseline/00_common.R"))

cell <- decode_task(args[[1L]])
z <- read_z(cell, "d10")
out_file <- cell_output("susiex", cell)
old <- safe_read_method_object(out_file, z, "SuSiEx")
if (!is.null(old)) {
  cat("Output exists; skipping:", out_file, "\n")
  quit(status = 0L)
}

b_eur <- read_bim(cell$locus, "EUR")
b_afr <- read_bim(cell$locus, "AFR")
if (!identical(z$RSID, b_eur$SNP) || !identical(z$RSID, b_afr$SNP)) {
  stop("BIM and summary-statistic SNP order differ", call. = FALSE)
}
if (!identical(b_eur[, c("SNP", "BP", "A1", "A2")],
               b_afr[, c("SNP", "BP", "A1", "A2")])) {
  stop("EUR/AFR reference alleles or IDs differ", call. = FALSE)
}

raw_dir <- file.path(RESULT_DIR, "susiex", "raw", sprintf("task_%03d", cell$task_id))
dir.create(raw_dir, recursive = TRUE, showWarnings = FALSE)
eur_sst <- file.path(raw_dir, "eur.tsv")
afr_sst <- file.path(raw_dir, "afr.tsv")

make_sst <- function(zscore, sample_size, bim) {
  data.frame(CHR = z$CHR, POS = z$POS, RSID = z$RSID,
             A1 = bim$A1, A2 = bim$A2,
             Beta = zscore / sqrt(sample_size), Se = 1 / sqrt(sample_size),
             Z = zscore, N = sample_size,
             P = 2 * pnorm(abs(zscore), lower.tail = FALSE),
             check.names = FALSE)
}
fwrite(make_sst(z$zscore_1, z$N_1, b_eur), eur_sst, sep = "\t")
fwrite(make_sst(z$zscore_2, z$N_2, b_afr), afr_sst, sep = "\t")

locus_info <- fread(LOCUS_INFO_FILE, data.table = FALSE)
if (nrow(locus_info) < cell$locus) stop("Locus-info table is incomplete", call. = FALSE)
chr <- locus_info$CHR[cell$locus]
bp_start <- locus_info$MinBP[cell$locus]
bp_end <- locus_info$MaxBP[cell$locus]
out_name <- "fit"
threads <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = "4"))
threads <- max(1L, threads)

cmd_args <- c(
  sprintf("--sst_file=%s,%s", eur_sst, afr_sst),
  sprintf("--ref_file=%s,%s", reference_prefix(cell$locus, "EUR"),
          reference_prefix(cell$locus, "AFR")),
  sprintf("--ld_file=%s,%s", susiex_ld_prefix(cell$locus, "EUR"),
          susiex_ld_prefix(cell$locus, "AFR")),
  sprintf("--n_gwas=%d,%d", N_GWAS[1L], N_GWAS[2L]),
  sprintf("--out_dir=%s", raw_dir), sprintf("--out_name=%s", out_name),
  sprintf("--chr=%s", chr), sprintf("--bp=%s,%s", bp_start, bp_end),
  "--chr_col=1,1", "--bp_col=2,2", "--snp_col=3,3",
  "--a1_col=4,4", "--a2_col=5,5", "--eff_col=6,6", "--se_col=7,7",
  "--pval_col=10,10", sprintf("--plink=%s", PLINK_BIN),
  "--keep-ambig=True", "--mult-step=True", "--maf=0.001", "--level=0.95",
  "--min_purity=0.5", "--pval_thresh=1e-05", "--tol=0.0001",
  "--n_sig=10", "--max_iter=200", sprintf("--threads=%d", threads)
)

start <- proc.time()[[3L]]
status <- system2(SUSIEX_BIN, cmd_args,
                  stdout = file.path(raw_dir, "stdout.log"),
                  stderr = file.path(raw_dir, "stderr.log"))
runtime_minutes <- (proc.time()[[3L]] - start) / 60
if (!identical(as.integer(status), 0L)) stop("SuSiEx exited with status ", status, call. = FALSE)

snp_file <- file.path(raw_dir, paste0(out_name, ".snp"))
cs_file <- file.path(raw_dir, paste0(out_name, ".cs"))
summary_file <- file.path(raw_dir, paste0(out_name, ".summary"))
if (!all(file.exists(c(snp_file, cs_file, summary_file)))) stop("SuSiEx output contract incomplete", call. = FALSE)

snp_header <- readLines(snp_file, n = 1L, warn = FALSE)
pip_columns <- grep("^PIP(\\(CS[0-9]+\\)|\\.CS[0-9]+\\.)$",
                    names(fread(snp_file, nrows = 0L, data.table = FALSE)), value = TRUE)
cs_text <- trimws(readLines(cs_file, warn = FALSE))
summary_text <- trimws(readLines(summary_file, warn = FALSE))
explicit_no_cs <- !length(pip_columns) && any(cs_text == "NULL") && any(summary_text == "NULL")

if (explicit_no_cs) {
  pip <- numeric(nrow(z))
  cs <- list()
} else {
  snp <- fread(snp_file, data.table = FALSE, check.names = FALSE)
  if (!"SNP" %in% names(snp) || anyDuplicated(snp$SNP) || !setequal(snp$SNP, z$RSID)) {
    stop("SuSiEx SNP output cannot be aligned to the simulated variants", call. = FALSE)
  }
  pip_columns <- grep("^PIP(\\(CS[0-9]+\\)|\\.CS[0-9]+\\.)$", names(snp), value = TRUE)
  if (!length(pip_columns)) stop("SuSiEx returned no PIP columns without explicit no-CS contract", call. = FALSE)
  M <- as.matrix(snp[match(z$RSID, snp$SNP), pip_columns, drop = FALSE])
  storage.mode(M) <- "double"
  if (any(!is.finite(M)) || any(M < -1e-10 | M > 1 + 1e-10)) stop("Invalid SuSiEx component PIPs", call. = FALSE)
  pip <- 1 - apply(1 - M, 1L, prod)

  cs_table <- fread(cs_file, data.table = FALSE)
  if (!all(c("CS_ID", "SNP") %in% names(cs_table))) stop("Malformed SuSiEx CS file", call. = FALSE)
  cs <- lapply(split(cs_table$SNP, cs_table$CS_ID), function(ids) {
    idx <- match(ids, z$RSID)
    if (anyNA(idx)) stop("SuSiEx CS contains an unknown SNP", call. = FALSE)
    as.integer(idx)
  })
  cs <- unname(Filter(length, cs))
}

obj <- list(
  method = "SuSiEx", task_id = cell$task_id, locus = cell$locus,
  causal_index = cell$causal_index, causal_count = cell$causal_count,
  h2_index = cell$h2_index, h2 = cell$h2,
  rsid = z$RSID, truth = z$Signal, pip = pmin(1, pmax(0, pip)), cs = cs,
  runtime_minutes = runtime_minutes, L = 10L, coverage = 0.95,
  min_abs_corr = 0.5, explicit_no_cs = explicit_no_cs,
  raw_dir = raw_dir
)
if (!validate_method_object(obj, z, "SuSiEx")) stop("Internal SuSiEx output validation failed", call. = FALSE)
atomic_save_rds(obj, out_file)
cat(sprintf("wrote %s; n_snp=%d; n_cs=%d; minutes=%.4f\n",
            out_file, nrow(z), length(cs), runtime_minutes))

