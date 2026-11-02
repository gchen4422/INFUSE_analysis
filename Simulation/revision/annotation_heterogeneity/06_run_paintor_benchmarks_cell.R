#!/usr/bin/env Rscript
.infuse_repo <- Sys.getenv("INFUSE_ANALYSIS_ROOT", "/path/to/INFUSE_analysis")
.infuse_project <- Sys.getenv("INFUSE_PROJECT", "/path/to/baseline_simulation")
.infuse_work <- Sys.getenv("INFUSE_REVISION", "/path/to/revision_output")

# Run the three PAINTOR benchmarks for one ancestry-prior sensitivity cell:
#   Paintor       : no functional annotations;
#   Paintor+      : the complete INFUSE 12-column informative panel;
#   Paintor+H     : the complete INFUSE-H 189-column panel.
#
# The command and file formats follow 4_shared_50_mesusie_all_annotes_paintor_only.R.
# Unlike that legacy script, this runner uses a deterministic MCMC seed and
# validates every result before declaring the cell complete.

suppressPackageStartupMessages(library(data.table))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4L) {
  stop(
    "Usage: Rscript 06_run_paintor_benchmarks_cell.R ",
    "<scenario> <locus> <h2_index> <causal_index>",
    call. = FALSE
  )
}

scenario <- args[[1L]]
locus <- as.integer(args[[2L]])
h2_index <- as.integer(args[[3L]])
causal_index <- as.integer(args[[4L]])

valid_scenarios <- c(
  "baseline", "different_annotation", "different_weight", "combined"
)
if (!scenario %in% valid_scenarios) {
  stop("Unknown scenario: ", scenario, call. = FALSE)
}
if (is.na(locus) || !locus %in% 1:100) {
  stop("locus must be an integer from 1 to 100.", call. = FALSE)
}
if (is.na(h2_index) || !h2_index %in% 1:2) {
  stop("h2_index must be 1 or 2.", call. = FALSE)
}
if (is.na(causal_index) || !causal_index %in% 1:3) {
  stop("causal_index must be 1, 2, or 3.", call. = FALSE)
}

project_dir <- .infuse_project
module_dir <- file.path(.infuse_work, "ancestry_specific_eaf_prior_simulation")
paintor_executable <- Sys.getenv("INFUSE_PAINTOR_BIN", unname(Sys.which("PAINTOR")))
source(file.path(.infuse_repo, "Simulation/revision/annotation_heterogeneity/paintor_output_contract.R"))
ld_eur_path <- file.path(
  dirname(project_dir), "risk_loci_ld_eur", paste0("loci_", locus, ".ld")
)
ld_afr_path <- file.path(
  dirname(project_dir), "risk_loci_ld_afr", paste0("loci_", locus, ".ld")
)

if (!file.exists(paintor_executable)) {
  stop("PAINTOR executable is missing: ", paintor_executable, call. = FALSE)
}

condition_dir <- function(base_d) {
  file.path(
    module_dir,
    paste0("scenario_", scenario),
    paste0("base_d_", base_d, "_input_d_", base_d + 2L)
  )
}

stem <- paste0(
  "CAUSAL_", causal_index,
  "_LOCI_", locus,
  "_h2_", h2_index
)
summary_d10 <- file.path(
  condition_dir(10L), paste0("causal_num_", causal_index),
  "summary_data", stem
)
summary_d187 <- file.path(
  condition_dir(187L), paste0("causal_num_", causal_index),
  "summary_data", stem
)
annotation_d10 <- file.path(
  condition_dir(10L), "annotations", paste0("Region_", locus, ".annot")
)
annotation_d187 <- file.path(
  condition_dir(187L), "annotations", paste0("Region_", locus, ".annot")
)

required_inputs <- c(
  summary_d10, summary_d187, annotation_d10, annotation_d187,
  ld_eur_path, ld_afr_path
)
if (any(!file.exists(required_inputs))) {
  stop(
    "Missing input(s): ",
    paste(required_inputs[!file.exists(required_inputs)], collapse = ", "),
    call. = FALSE
  )
}

# d=10 and d=187 use the same simulated data. Enforce this rather than merely
# assuming it, because all three PAINTOR fits must benchmark the same cell.
summary_hashes <- unname(tools::md5sum(c(summary_d10, summary_d187)))
if (!identical(summary_hashes[[1L]], summary_hashes[[2L]])) {
  stop("d=10 and d=187 summary statistics differ for ", stem, call. = FALSE)
}

zfile <- fread(summary_d10)
required_z_columns <- c("RSID", "zscore_1", "zscore_2", "Signal")
if (!all(required_z_columns %in% names(zfile))) {
  stop("Summary file lacks required columns: ", summary_d10, call. = FALSE)
}
if (anyDuplicated(zfile$RSID)) {
  stop("Summary file contains duplicated RSIDs: ", summary_d10, call. = FALSE)
}
n_snp <- nrow(zfile)

annotations10 <- fread(annotation_d10)
annotations187 <- fread(annotation_d187)
if (nrow(annotations10) != n_snp || ncol(annotations10) != 12L) {
  stop("Expected the complete nSNP x 12 Paintor+/INFUSE annotation matrix.", call. = FALSE)
}
if (nrow(annotations187) != n_snp || ncol(annotations187) != 189L) {
  stop("Expected the complete nSNP x 189 Paintor+H/INFUSE-H annotation matrix.", call. = FALSE)
}
expected_eaf <- c("EAF_EUR", "EAF_AFR")
if (!identical(tail(names(annotations10), 2L), expected_eaf) ||
    !identical(tail(names(annotations187), 2L), expected_eaf)) {
  stop("Both annotation matrices must end with EAF_EUR and EAF_AFR.", call. = FALSE)
}
if (any(!vapply(annotations10, is.numeric, logical(1))) ||
    any(!vapply(annotations187, is.numeric, logical(1))) ||
    anyNA(annotations10) || anyNA(annotations187)) {
  stop("PAINTOR annotation matrices must be complete and numeric.", call. = FALSE)
}

# The source LD files have SNP identifiers in their header. fread removes that
# header, matching the legacy runner's conversion to a numeric PAINTOR matrix.
ld_eur <- as.matrix(fread(ld_eur_path))
ld_afr <- as.matrix(fread(ld_afr_path))
storage.mode(ld_eur) <- "double"
storage.mode(ld_afr) <- "double"
if (!identical(dim(ld_eur), c(n_snp, n_snp)) ||
    !identical(dim(ld_afr), c(n_snp, n_snp))) {
  stop("LD and summary-statistic dimensions differ for locus ", locus, call. = FALSE)
}

benchmark_root <- file.path(module_dir, "paintor_benchmarks")
work_parent <- file.path(benchmark_root, ".work")
dir.create(work_parent, recursive = TRUE, showWarnings = FALSE)
work_dir <- file.path(
  work_parent,
  paste0(
    "scenario_", scenario, "_", stem, "_",
    Sys.getenv("SLURM_JOB_ID", unset = Sys.getpid()), "_",
    Sys.getenv("SLURM_ARRAY_TASK_ID", unset = "local")
  )
)
if (dir.exists(work_dir)) unlink(work_dir, recursive = TRUE, force = TRUE)
dir.create(work_dir, recursive = TRUE, showWarnings = FALSE)
on.exit(unlink(work_dir, recursive = TRUE, force = TRUE), add = TRUE)

file.copy(summary_d10, file.path(work_dir, stem), overwrite = TRUE)
fwrite(
  as.data.table(ld_eur), file.path(work_dir, paste0(stem, ".LD1")),
  sep = " ", col.names = FALSE
)
fwrite(
  as.data.table(ld_afr), file.path(work_dir, paste0(stem, ".LD2")),
  sep = " ", col.names = FALSE
)
input_file <- file.path(work_dir, paste0(stem, ".input"))
writeLines(stem, input_file)

# PAINTOR requires the locus annotation file to exist even when no annotation
# names are supplied. For the unannotated fit it reads this file but fits only
# its intercept; this is also how the reference runner is structured.
fwrite(
  annotations10, file.path(work_dir, paste0(stem, ".annotations")),
  sep = " ", col.names = TRUE
)

validate_result <- function(path) {
  if (!file.exists(path) || file.info(path)$size <= 0) return(FALSE)
  ans <- tryCatch(fread(path), error = function(e) NULL)
  if (is.null(ans) || nrow(ans) != n_snp ||
      !all(c("RSID", "Posterior_Prob") %in% names(ans))) return(FALSE)
  if (!identical(as.character(ans$RSID), as.character(zfile$RSID))) return(FALSE)
  pip <- ans$Posterior_Prob
  paintor_valid_pip(pip) && !paintor_log_diagnostics(sub("\\.mcmc\\.paintor$",".paintor.log",path))$optimizer_error_logged
}

run_one <- function(method, annotation_table = NULL, base_d = 0L) {
  method_label <- c(
    paintor = "Paintor",
    paintor_plus_10 = "Paintor+",
    paintor_plus_187 = "Paintor+H"
  )[[method]]
  expected_input_d <- if (is.null(annotation_table)) 0L else ncol(annotation_table)
  method_dir <- file.path(
    benchmark_root, paste0("scenario_", scenario), method,
    paste0("causal_num_", causal_index), "result"
  )
  dir.create(method_dir, recursive = TRUE, showWarnings = FALSE)
  result_path <- file.path(method_dir, paste0(stem, ".mcmc.paintor"))
  metadata_path <- file.path(method_dir, paste0(stem, ".metadata.tsv"))
  log_path <- file.path(method_dir, paste0(stem, ".paintor.log"))

  if (validate_result(result_path) && file.exists(metadata_path)) {
    existing_metadata <- tryCatch(fread(metadata_path), error = function(e) NULL)
    metadata_matches <- !is.null(existing_metadata) && nrow(existing_metadata) == 1L &&
      identical(as.character(existing_metadata$method), method) &&
      identical(as.character(existing_metadata$method_label), method_label) &&
      identical(as.integer(existing_metadata$base_annotation_dimension), base_d) &&
      identical(as.integer(existing_metadata$fitted_annotation_dimension), expected_input_d)
    if (metadata_matches) {
      if(paintor_pip_status(fread(result_path)$Posterior_Prob)=="all_zero_no_selection")
        message("Existing all-zero output retained as no selection; convergence not certified: ",result_path)
      return(invisible(result_path))
    }
  }

  annotation_args <- character()
  input_d <- 0L
  if (!is.null(annotation_table)) {
    input_d <- ncol(annotation_table)
    fwrite(
      annotation_table, file.path(work_dir, paste0(stem, ".annotations")),
      sep = " ", col.names = TRUE
    )
    annotation_args <- c(
      "-annotations", paste(names(annotation_table), collapse = ",")
    )
  }

  scenario_index <- match(scenario, valid_scenarios)
  method_index <- match(method, c("paintor", "paintor_plus_10", "paintor_plus_187"))
  seed <- as.integer(
    10000000L + scenario_index * 1000000L + causal_index * 100000L +
      h2_index * 10000L + locus * 10L + method_index
  )

  command_args <- c(
    "-input", input_file,
    "-Zhead", "zscore_1,zscore_2",
    "-LDname", "LD1,LD2",
    "-in", work_dir,
    "-out", method_dir,
    "-mcmc",
    annotation_args,
    "-set_seed", seed,
    "-Gname", paste0(stem, ".enrich"),
    "-Lname", paste0(stem, "_bayes_factor"),
    "-RESname", "mcmc.paintor"
  )

  started <- Sys.time()
  status <- system2(
    paintor_executable, args = command_args,
    stdout = log_path, stderr = log_path
  )
  elapsed_seconds <- as.numeric(difftime(Sys.time(), started, units = "secs"))
  if (!identical(status, 0L) || !validate_result(result_path)) {
    stop(
      "PAINTOR failed validation for ", scenario, "/", method, "/", stem,
      ". See ", log_path, call. = FALSE
    )
  }

  metadata <- data.table(
    scenario = scenario,
    method = method,
    method_label = method_label,
    base_annotation_dimension = base_d,
    fitted_annotation_dimension = input_d,
    includes_EAF_EUR = input_d > 0L,
    includes_EAF_AFR = input_d > 0L,
    locus = locus,
    h2_index = h2_index,
    causal_index = causal_index,
    n_snp = n_snp,
    seed = seed,
    summary_md5 = summary_hashes[[1L]],
    elapsed_seconds = elapsed_seconds,
    paintor_executable = paintor_executable,
    created_at = format(Sys.time(), tz = "America/Indianapolis", usetz = TRUE)
  )
  pip_saved <- fread(result_path)$Posterior_Prob
  metadata[, `:=`(output_status=paintor_pip_status(pip_saved),total_pip=sum(pip_saved),
                   n_positive_pip=sum(pip_saved>0),exit_status=status)]
  diagnostic <- paintor_log_diagnostics(log_path)
  for(nm in names(diagnostic)) metadata[,(nm):=diagnostic[[nm]]]
  fwrite(metadata, metadata_path, sep = "\t")
  if(metadata$output_status=="all_zero_no_selection")
    message("All-zero output: no selection; convergence not certified: ",result_path)
  invisible(result_path)
}

run_one("paintor", base_d = 0L)
run_one("paintor_plus_10", annotations10, base_d = 10L)
run_one("paintor_plus_187", annotations187, base_d = 187L)

message(
  "Completed all PAINTOR benchmarks for scenario=", scenario,
  ", locus=", locus, ", h2_index=", h2_index,
  ", causal_index=", causal_index
)
