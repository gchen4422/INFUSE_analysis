.infuse_repo <- Sys.getenv("INFUSE_ANALYSIS_ROOT", "/path/to/INFUSE_analysis")
.infuse_project <- Sys.getenv("INFUSE_PROJECT", "/path/to/baseline_simulation")
.infuse_work <- Sys.getenv("INFUSE_REVISION", "/path/to/revision_output")
# Fine-mapping analysis for the ancestry-specific causal-prior sensitivity study.
#
# This is a focused counterpart of 4_shared_50_mesusie_inf.R and
# More_annotes/4_shared_50_mesusie_inf.R. It runs:
#   1. MESuSiE without annotations (multi-ancestry control);
#   2. ancestry-specific SuSiE (compatibility controls); and
#   3. the current GitHub INFUSE model. The package directly optimizes the
#      annotation coefficients for 1--14 inputs and uses sparse PCA/BIC for
#      15 or more inputs.
#
# PAINTOR and secondary INFUSE fitting variants are intentionally excluded: the
# revision question is the robustness of the primary INFUSE fit when the true
# causal prior differs across ancestries.

suppressPackageStartupMessages({
  library(data.table)
  library(INFUSE)
  library(MESuSiE)
  library(susieR)
})

args <- commandArgs(trailingOnly = TRUE)
if (!length(args) %in% 5:6) {
  stop(
    "Usage: Rscript 02_shared_50_mesusie_inf_ancestry_specific_prior.R ",
    "<scenario> <base_d> <LD_BLOCK> <h2_num> <num_causal> [result_tag]",
    call. = FALSE
  )
}

scenario <- args[[1L]]
base_d <- as.integer(args[[2L]])
LD_BLOCK <- as.integer(args[[3L]])
h2_num <- as.integer(args[[4L]])
num_causal <- as.integer(args[[5L]])
result_tag <- if (length(args) == 6L) args[[6L]] else ""

valid_scenarios <- c(
  "baseline",
  "different_annotation",
  "different_weight",
  "combined"
)
if (!scenario %in% valid_scenarios) {
  stop("Unknown scenario: ", scenario, call. = FALSE)
}
if (!base_d %in% c(10L, 187L)) {
  stop("base_d must be 10 or 187.", call. = FALSE)
}
if (!LD_BLOCK %in% seq_len(100L)) {
  stop("LD_BLOCK must be between 1 and 100.", call. = FALSE)
}
if (!h2_num %in% 1:2) {
  stop("h2_num must be 1 or 2.", call. = FALSE)
}
if (!num_causal %in% 1:3) {
  stop("num_causal must be 1, 2, or 3.", call. = FALSE)
}
if (nzchar(result_tag) && !grepl("^[A-Za-z0-9_]+$", result_tag)) {
  stop("result_tag may contain only letters, numbers, and underscores.", call. = FALSE)
}

# INFUSE runtime/reproducibility requirement. These variables must be set
# before R starts, so this check detects an incorrectly launched analysis.
thread_variables <- c(
  "OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS",
  "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS"
)
thread_settings <- Sys.getenv(thread_variables, unset = NA_character_)
if (anyNA(thread_settings) || any(thread_settings != "1")) {
  stop(
    "All INFUSE numerical thread variables must equal 1: ",
    paste(paste0(thread_variables, "=", thread_settings), collapse = ", "),
    call. = FALSE
  )
}

project_dir <- .infuse_project
module_dir <- file.path(.infuse_work, "ancestry_specific_eaf_prior_simulation")
script_dir <- file.path(module_dir, "scripts")
github_infuse_dir <- Sys.getenv("INFUSE_PACKAGE_SOURCE", file.path(dirname(.infuse_repo), "INFUSE"))
input_d <- base_d + 2L

condition_dir <- file.path(
  module_dir,
  paste0("scenario_", scenario),
  paste0("base_d_", base_d, "_input_d_", input_d)
)
data_dir <- file.path(condition_dir, paste0("causal_num_", num_causal), "summary_data")
annot_dir <- file.path(condition_dir, "annotations")
result_dir_name <- if (nzchar(result_tag)) paste0("result_", result_tag) else "result"
result_dir <- file.path(
  condition_dir, paste0("causal_num_", num_causal), result_dir_name
)
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)

stem <- paste0("CAUSAL_", num_causal, "_LOCI_", LD_BLOCK, "_h2_", h2_num)
summary_path <- file.path(data_dir, stem)
annotation_path <- file.path(annot_dir, paste0("Region_", LD_BLOCK, ".annot"))
result_path <- file.path(result_dir, paste0("Res_", stem, ".RData"))
timing_path <- file.path(result_dir, paste0("time_", stem, ".tsv"))

overwrite <- identical(Sys.getenv("INFUSE_OVERWRITE", unset = "0"), "1")
if (file.exists(result_path) && !overwrite) {
  if (!file.exists(timing_path)) {
    saved_result <- new.env(parent = emptyenv())
    load(result_path, envir = saved_result)
    if (!exists("timing", envir = saved_result, inherits = FALSE)) {
      stop("Existing result has no timing object: ", result_path, call. = FALSE)
    }
    fwrite(get("timing", envir = saved_result, inherits = FALSE), timing_path, sep = "\t")
    message("Recovered missing timing file: ", timing_path)
  }
  message("Result already exists; skipping: ", result_path)
  quit(save = "no", status = 0L)
}

required_paths <- c(summary_path, annotation_path)
if (any(!file.exists(required_paths))) {
  stop(
    "Missing required input(s): ",
    paste(required_paths[!file.exists(required_paths)], collapse = ", "),
    call. = FALSE
  )
}

zfile <- fread(summary_path)
required_columns <- c("RSID", "zscore_1", "zscore_2", "Signal", "N_1", "N_2")
if (!all(required_columns %in% names(zfile))) {
  stop("Summary-statistic file is missing required columns.", call. = FALSE)
}
if (any(!is.finite(zfile$zscore_1)) || any(!is.finite(zfile$zscore_2)) ||
    any(!is.finite(zfile$N_1)) || any(!is.finite(zfile$N_2))) {
  stop("Summary statistics contain non-finite values.", call. = FALSE)
}

annot_file_subset <- fread(annotation_path)
if (nrow(annot_file_subset) != nrow(zfile)) {
  stop("Annotation and summary-statistic row counts differ.", call. = FALSE)
}
if (ncol(annot_file_subset) != input_d) {
  stop(
    "Expected ", input_d, " annotation columns but found ",
    ncol(annot_file_subset), ".",
    call. = FALSE
  )
}
annot_matrix <- as.matrix(annot_file_subset)
storage.mode(annot_matrix) <- "double"
if (any(!is.finite(annot_matrix))) {
  stop("Annotation matrix contains non-finite values.", call. = FALSE)
}

summary_stat_1 <- data.frame(
  SNP = zfile$RSID,
  Beta = zfile$zscore_1 / sqrt(zfile$N_1),
  Se = 1 / sqrt(zfile$N_1),
  Z = zfile$zscore_1,
  N = zfile$N_1
)
summary_stat_2 <- data.frame(
  SNP = zfile$RSID,
  Beta = zfile$zscore_2 / sqrt(zfile$N_2),
  Se = 1 / sqrt(zfile$N_2),
  Z = zfile$zscore_2,
  N = zfile$N_2
)
summary_stat_sd_list <- list(EUR = summary_stat_1, AFR = summary_stat_2)

ld_dir_eur <- file.path(Sys.getenv("INFUSE_SIMULATION", dirname(.infuse_project)), "risk_loci_ld_eur")
ld_dir_afr <- file.path(Sys.getenv("INFUSE_SIMULATION", dirname(.infuse_project)), "risk_loci_ld_afr")
EU_cov <- as.matrix(fread(file.path(ld_dir_eur, paste0("loci_", LD_BLOCK, ".ld"))))
BB_cov <- as.matrix(fread(file.path(ld_dir_afr, paste0("loci_", LD_BLOCK, ".ld"))))
storage.mode(EU_cov) <- storage.mode(BB_cov) <- "double"

n_snp <- nrow(zfile)
if (!identical(dim(EU_cov), c(n_snp, n_snp)) ||
    !identical(dim(BB_cov), c(n_snp, n_snp))) {
  stop("LD and summary-statistic dimensions differ.", call. = FALSE)
}
dimnames(EU_cov) <- list(zfile$RSID, zfile$RSID)
dimnames(BB_cov) <- list(zfile$RSID, zfile$RSID)
R_mat_list <- list(EUR = EU_cov, AFR = BB_cov)

bayes_fac <- 3
ancestry_weight <- c(
  bayes_fac / (bayes_fac * 2 + 1),
  bayes_fac / (bayes_fac * 2 + 1),
  1 / (bayes_fac * 2 + 1)
)

start_time <- proc.time()[["elapsed"]]
MESuSiE_res <- MESuSiE::meSuSie_core(
  R_mat_list,
  summary_stat_sd_list,
  L = 10
)
MESuSiE_time <- proc.time()[["elapsed"]] - start_time

start_time <- proc.time()[["elapsed"]]
susie_EU <- susieR::susie_rss(zfile$zscore_1, EU_cov, check_prior = FALSE)
susie_BB <- susieR::susie_rss(zfile$zscore_2, BB_cov, check_prior = FALSE)
susie_time <- proc.time()[["elapsed"]] - start_time

# Verify that the installed package used for production is the same R-level
# implementation as the current local GitHub checkout. This prevents a stale
# installation from silently changing the low-dimensional annotation path.
github_core_path <- file.path(github_infuse_dir, "R", "infuse_core.R")
if (!file.exists(github_core_path)) {
  stop("Current GitHub INFUSE source is missing: ", github_core_path, call. = FALSE)
}
github_env <- new.env(parent = globalenv())
sys.source(github_core_path, envir = github_env)
installed_body <- paste(deparse(body(INFUSE::infuse_core)), collapse = "\n")
github_body <- paste(deparse(body(github_env$infuse_core)), collapse = "\n")
if (!identical(installed_body, github_body)) {
  stop("Installed INFUSE::infuse_core() differs from the current GitHub checkout.",
       call. = FALSE)
}
github_commit <- trimws(system2(
  "git", c("-C", github_infuse_dir, "rev-parse", "HEAD"), stdout = TRUE
))
if (length(github_commit) != 1L || !grepl("^[0-9a-f]{40}$", github_commit)) {
  stop("Could not resolve the current GitHub INFUSE commit.", call. = FALSE)
}

start_time <- proc.time()[["elapsed"]]
INFUSE_res <- INFUSE::infuse_core(
  R_mat_list,
  summary_stat_sd_list,
  L = 10,
  residual_variance = NULL,
  prior_weights = NULL,
  ancestry_weight = ancestry_weight,
  optim_method = "optim",
  estimate_residual_variance = FALSE,
  max_iter = 100,
  annot = annot_matrix,
  annot_method = "mlk"
)
INFUSE_time <- proc.time()[["elapsed"]] - start_time

if (length(INFUSE_res$pip) != n_snp || any(!is.finite(INFUSE_res$pip))) {
  stop("INFUSE returned invalid marginal PIPs.", call. = FALSE)
}

# Preserve the primary historical object name for straightforward adaptation of
# the existing result-summarizing code.
MESuSiE_res_mlk_annot_subset <- INFUSE_res
MAFPC.time <- INFUSE_time
MESuSiE.time <- MESuSiE_time
susie.time <- susie_time

timing <- data.frame(
  scenario = scenario,
  base_d = base_d,
  input_d = input_d,
  LD_BLOCK = LD_BLOCK,
  h2_num = h2_num,
  num_causal = num_causal,
  MESuSiE_seconds = MESuSiE_time,
  SuSiE_seconds = susie_time,
  INFUSE_seconds = INFUSE_time
)
analysis_metadata <- list(
  scenario = scenario,
  base_d = base_d,
  input_d = input_d,
  LD_BLOCK = LD_BLOCK,
  h2_num = h2_num,
  num_causal = num_causal,
  result_tag = result_tag,
  annotation_names = colnames(annot_matrix),
  infuse_package_version = as.character(utils::packageVersion("INFUSE")),
  infuse_package_path = find.package("INFUSE"),
  infuse_github_checkout = github_infuse_dir,
  infuse_github_commit = github_commit,
  installed_core_matches_github = TRUE,
  numerical_threads = setNames(as.list(thread_settings), thread_variables),
  infuse_annotation_path = if (input_d < 15L) {
    "direct_annotation_optimization"
  } else {
    "sparse_pca_bic"
  },
  paintor_run = FALSE
)

temporary_result <- tempfile(
  pattern = paste0(".", stem, "_"),
  tmpdir = result_dir,
  fileext = ".RData"
)
save(
  susie_EU,
  susie_BB,
  MESuSiE_res,
  INFUSE_res,
  MESuSiE_res_mlk_annot_subset,
  MESuSiE.time,
  susie.time,
  MAFPC.time,
  timing,
  analysis_metadata,
  file = temporary_result
)
if (!file.rename(temporary_result, result_path)) {
  stop("Failed to atomically move the result into place.", call. = FALSE)
}
fwrite(timing, timing_path, sep = "\t")

message("Completed: ", result_path)
