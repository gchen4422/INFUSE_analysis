.infuse_repo <- Sys.getenv("INFUSE_ANALYSIS_ROOT", "/path/to/INFUSE_analysis")
.infuse_project <- Sys.getenv("INFUSE_PROJECT", "/path/to/baseline_simulation")
.infuse_work <- Sys.getenv("INFUSE_REVISION", "/path/to/revision_output")
options(stringsAsFactors = FALSE)

suppressPackageStartupMessages(library(data.table))

PROJECT_ROOT <- .infuse_project
REVISION_ROOT <- .infuse_work
source(file.path(.infuse_repo, "Simulation/revision/common/baseline_sources.R"))
MODULE_DIR <- file.path(REVISION_ROOT, "baseline_susiex_multisusie_benchmark")
RESULT_DIR <- file.path(MODULE_DIR, "results")
FIGURE_DIR <- file.path(REVISION_ROOT, "figures", "baseline_susiex_multisusie_benchmark")
LOG_DIR <- file.path(MODULE_DIR, "logs")
WORK_DIR <- file.path(MODULE_DIR, ".work")

# External paths are read-only inputs. This module never writes outside REVISION_ROOT.
BASELINE_LD_ROOT <- Sys.getenv("INFUSE_SHARED_LD", file.path(dirname(.infuse_project), "shared_50"))
SIMULATION_ROOT <- Sys.getenv("INFUSE_SIMULATION", dirname(.infuse_project))
EUR_REF_ROOT <- file.path(SIMULATION_ROOT, "risk_loci_ld_eur")
AFR_REF_ROOT <- file.path(SIMULATION_ROOT, "risk_loci_ld_afr")
SUSIEX_LD_EUR_ROOT <- file.path(BASELINE_LD_ROOT, "susiex_ld_eur")
SUSIEX_LD_AFR_ROOT <- file.path(BASELINE_LD_ROOT, "susiex_ld_afr")
LOCUS_INFO_FILE <- file.path(SIMULATION_ROOT, "simulation_locus_info.txt")
SUSIEX_BIN <- Sys.getenv("INFUSE_SUSIEX_BIN", unname(Sys.which("SuSiEx")))
PLINK_BIN <- Sys.getenv("INFUSE_PLINK_BIN", unname(Sys.which("plink")))
PYTHON_BIN <- Sys.getenv("INFUSE_PYTHON_BIN", unname(Sys.which("python3")))

for (d in c(RESULT_DIR, FIGURE_DIR, LOG_DIR, WORK_DIR,
            file.path(RESULT_DIR, "susiex", "cells"),
            file.path(RESULT_DIR, "susiex", "raw"),
            file.path(RESULT_DIR, "multisusie", "cells"))) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

N_LOCI <- 100L
N_CAUSAL_INDEX <- 3L
N_H2 <- 2L
N_CELLS <- N_LOCI * N_CAUSAL_INDEX * N_H2
CAUSAL_COUNTS <- c(2L, 3L, 5L)
H2_VALUES <- c(1e-4, 2e-4)
N_GWAS <- c(300000L, 300000L)

decode_task <- function(task_id) {
  task_id <- as.integer(task_id)
  if (length(task_id) != 1L || is.na(task_id) || task_id < 1L || task_id > N_CELLS) {
    stop("task_id must be an integer in 1:", N_CELLS, call. = FALSE)
  }
  x <- task_id - 1L
  locus <- x %% N_LOCI + 1L
  x <- x %/% N_LOCI
  causal_index <- x %% N_CAUSAL_INDEX + 1L
  x <- x %/% N_CAUSAL_INDEX
  h2_index <- x %% N_H2 + 1L
  list(task_id = task_id, locus = locus, causal_index = causal_index,
       causal_count = CAUSAL_COUNTS[causal_index], h2_index = h2_index,
       h2 = H2_VALUES[h2_index])
}

cell_stem <- function(cell) {
  sprintf("CAUSAL_%d_LOCI_%d_h2_%d",
          cell$causal_index, cell$locus, cell$h2_index)
}

z_path <- function(cell, dimension = c("d10", "d187")) {
  dimension <- match.arg(dimension)
  spec <- figure2_baseline_spec(dimension)
  file.path(spec$data_root, sprintf("causal_num_%d", cell$causal_index),
            spec$data_subdir, cell_stem(cell))
}

ld_path <- function(cell, ancestry = c("EUR", "AFR")) {
  ancestry <- match.arg(ancestry)
  suffix <- if (ancestry == "EUR") ".LD1" else ".LD2"
  file.path(BASELINE_LD_ROOT, sprintf("causal_num_%d", cell$causal_index),
            "summary_data", paste0(cell_stem(cell), suffix))
}

reference_prefix <- function(locus, ancestry = c("EUR", "AFR")) {
  ancestry <- match.arg(ancestry)
  file.path(if (ancestry == "EUR") EUR_REF_ROOT else AFR_REF_ROOT,
            sprintf("loci_%d", as.integer(locus)))
}

susiex_ld_prefix <- function(locus, ancestry = c("EUR", "AFR")) {
  ancestry <- match.arg(ancestry)
  root <- if (ancestry == "EUR") SUSIEX_LD_EUR_ROOT else SUSIEX_LD_AFR_ROOT
  file.path(root, sprintf("loci_%d.ld", as.integer(locus)))
}

frequency_path <- function(locus, ancestry = c("EUR", "AFR")) {
  ancestry <- match.arg(ancestry)
  root <- if (ancestry == "EUR") EUR_REF_ROOT else AFR_REF_ROOT
  file.path(root, sprintf("maf_loci_%d_maf.txt", as.integer(locus)))
}

cell_output <- function(method = c("susiex", "multisusie"), cell) {
  method <- match.arg(method)
  file.path(RESULT_DIR, method, "cells", sprintf("task_%03d.rds", cell$task_id))
}

read_z <- function(cell, dimension = "d10") {
  f <- z_path(cell, dimension)
  if (!file.exists(f)) stop("Missing summary-statistic input: ", f, call. = FALSE)
  z <- fread(f, data.table = FALSE)
  required <- c("CHR", "POS", "RSID", "zscore_1", "zscore_2", "Signal", "N_1", "N_2")
  if (!identical(names(z), required)) stop("Unexpected columns in ", f, call. = FALSE)
  if (!nrow(z) || anyDuplicated(z$RSID) || anyNA(z)) stop("Invalid summary statistics: ", f, call. = FALSE)
  if (!all(z$Signal %in% 0:3) || sum(z$Signal != 0L) != cell$causal_count) {
    stop("Truth count/configuration mismatch in ", f, call. = FALSE)
  }
  if (!all(z$N_1 == N_GWAS[1L]) || !all(z$N_2 == N_GWAS[2L])) {
    stop("Unexpected GWAS sample size in ", f, call. = FALSE)
  }
  z
}

read_bim <- function(locus, ancestry = c("EUR", "AFR")) {
  ancestry <- match.arg(ancestry)
  f <- paste0(reference_prefix(locus, ancestry), ".bim")
  if (!file.exists(f)) stop("Missing BIM: ", f, call. = FALSE)
  b <- fread(f, header = FALSE, data.table = FALSE)
  if (ncol(b) != 6L || anyDuplicated(b[[2L]])) stop("Malformed BIM: ", f, call. = FALSE)
  names(b) <- c("CHR", "SNP", "CM", "BP", "A1", "A2")
  b
}

read_ld <- function(cell, ancestry = c("EUR", "AFR"), expected_ids = NULL) {
  ancestry <- match.arg(ancestry)
  f <- ld_path(cell, ancestry)
  if (!file.exists(f)) stop("Missing LD matrix: ", f, call. = FALSE)
  R <- as.matrix(fread(f, header = FALSE, data.table = FALSE))
  storage.mode(R) <- "double"
  if (nrow(R) != ncol(R)) stop("LD matrix is not square: ", f, call. = FALSE)
  if (!is.null(expected_ids) && nrow(R) != length(expected_ids)) {
    stop("LD/SNP dimension mismatch: ", f, call. = FALSE)
  }
  if (any(!is.finite(R))) stop("LD matrix contains non-finite values: ", f, call. = FALSE)
  if (max(abs(diag(R) - 1)) > 1e-6) stop("LD diagonal is not one: ", f, call. = FALSE)
  if (max(abs(R - t(R))) > 1e-6) stop("LD matrix is not symmetric: ", f, call. = FALSE)
  if (!is.null(expected_ids)) dimnames(R) <- list(expected_ids, expected_ids)
  R
}

read_true_maf <- function(locus, ancestry = c("EUR", "AFR"), expected_n) {
  ancestry <- match.arg(ancestry)
  f <- frequency_path(locus, ancestry)
  if (!file.exists(f)) stop("Missing allele-frequency input: ", f, call. = FALSE)
  alt_frequency <- scan(f, quiet = TRUE)
  if (length(alt_frequency) != expected_n || any(!is.finite(alt_frequency)) ||
      any(alt_frequency <= 0 | alt_frequency >= 1)) {
    stop("Invalid allele-frequency input: ", f, call. = FALSE)
  }
  # The historical files are named MAF but contain ALT frequencies in some panels.
  # MultiSuSiE requests minor allele frequency, so convert explicitly.
  maf <- pmin(alt_frequency, 1 - alt_frequency)
  if (any(maf <= 0 | maf > 0.5)) stop("MAF conversion failed: ", f, call. = FALSE)
  list(maf = maf, alt_frequency = alt_frequency,
       n_above_half = sum(alt_frequency > 0.5))
}

atomic_save_rds <- function(object, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  tmp <- paste0(path, ".partial_", Sys.getpid())
  saveRDS(object, tmp, compress = FALSE)
  if (!file.rename(tmp, path)) stop("Atomic rename failed for ", path, call. = FALSE)
  invisible(path)
}

validate_method_object <- function(x, z, method) {
  if (!is.list(x) || !identical(x$method, method)) return(FALSE)
  if (length(x$pip) != nrow(z) || any(!is.finite(x$pip)) ||
      any(x$pip < -1e-10 | x$pip > 1 + 1e-10)) return(FALSE)
  if (!is.list(x$cs)) return(FALSE)
  idx <- unlist(x$cs, use.names = FALSE)
  if (length(idx) && any(idx < 1L | idx > nrow(z))) return(FALSE)
  identical(as.character(x$rsid), as.character(z$RSID))
}

safe_read_method_object <- function(path, z, method) {
  if (!file.exists(path)) return(NULL)
  x <- tryCatch(readRDS(path), error = function(e) NULL)
  if (is.null(x) || !validate_method_object(x, z, method)) return(NULL)
  x
}
