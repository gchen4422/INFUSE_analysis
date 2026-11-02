#!/usr/bin/env Rscript
.infuse_repo <- Sys.getenv("INFUSE_ANALYSIS_ROOT", "/path/to/INFUSE_analysis")
################################################################################
# One task evaluates one LD-mismatch scenario at one locus over all six
# (num_causal, h2) cells.  Z scores and causal truth are read unchanged from the
# existing w2_startover simulation.
#
# Usage:
#   Rscript 01_run_sparse.R <locus:1-100> <scenario>
################################################################################

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) {
  stop("Usage: Rscript 01_run_sparse.R <locus:1-100> <scenario>", call. = FALSE)
}
locus <- as.integer(args[[1]])
scenario_id <- args[[2]]
if (is.na(locus) || locus < 1L || locus > 100L) stop("locus must be 1-100", call. = FALSE)

this_file <- normalizePath(sub("^--file=", "",
  grep("^--file=", commandArgs(FALSE), value = TRUE)[1]), mustWork = FALSE)
script_dir <- dirname(this_file)
source(file.path(.infuse_repo, "Simulation/revision/ld_mismatch/ld_mismatch_common.R"))

suppressMessages({
  library(data.table)
  library(dplyr)
  library(susieR)
  library(MESuSiE)
  library(Matrix)
  library(mixsqp)
  library(ggplot2)
})
source(KRIGING_SOURCE)
source(KRIGING_SEPARATE_SOURCE)

grid <- scenario_grid()
scenario_row <- grid[grid$scenario == scenario_id, , drop = FALSE]
if (nrow(scenario_row) != 1L) {
  stop("Unknown scenario: ", scenario_id, ". Valid: ",
       paste(grid$scenario, collapse = ", "), call. = FALSE)
}

data.table::setDTthreads(1L)

dir.create(file.path(RESULT_DIR, "per_task"), recursive = TRUE, showWarnings = FALSE)
out_file <- file.path(RESULT_DIR, "per_task",
                      sprintf("ld_mismatch_%s_locus_%03d.RData", scenario_id, locus))

if (file.exists(out_file)) stop("Refusing to overwrite existing task: ", out_file)
metrics_all <- list()
matrix_all <- list()
scores_all <- list()
s_all <- list()
targets_all <- list()
errors_all <- list()
timing_all <- list()
anchor_detection_all <- list()

cat(sprintf("LD mismatch worker: locus=%d scenario=%s\n", locus, scenario_id))

for (causal_num in 1:3) {
  for (h2_num in 1:2) {
    cell_start <- Sys.time()
    zfile <- load_zfile(causal_num, locus, h2_num)
    R_true <- load_true_ld(causal_num, locus, h2_num, zfile$RSID)
    validate_correlation_matrix(R_true$EUR, "true EUR LD")
    validate_correlation_matrix(R_true$AFR, "true AFR LD")

    seed_base <- 1729L + 100000L * locus + 1000L * causal_num + 100L * h2_num
    scen <- apply_ld_scenario(scenario_row, R_true, zfile, seed_base)
    validate_correlation_matrix(scen$R$EUR, paste0(scenario_id, " EUR LD"))
    validate_correlation_matrix(scen$R$AFR, paste0(scenario_id, " AFR LD"))

    tag <- data.frame(
      scenario = scenario_id,
      scenario_class = scenario_row$scenario_class,
      delta = scenario_row$delta,
      lambda = scenario_row$lambda,
      anchor_scale = scenario_row$anchor_scale,
      mismatch_ancestry = scenario_row$mismatch_ancestry,
      locus = locus,
      causal_num = causal_num,
      h2_num = h2_num,
      n_snp = nrow(zfile),
      truth_scope = scen$truth_scope,
      target_rule = scen$target_rule,
      stringsAsFactors = FALSE
    )

    matrix_all[[length(matrix_all) + 1L]] <- cbind(
      tag,
      rbind(ld_difference_summary(R_true$EUR, scen$R$EUR, "EUR"),
            ld_difference_summary(R_true$AFR, scen$R$AFR, "AFR"))
    )

    target_rows <- function(ancestry, idx) {
      if (!length(idx)) {
        out <- tag[0, , drop = FALSE]
        out$ancestry <- character(0)
        out$SNP_index <- integer(0)
        out$RSID <- character(0)
        out$Signal <- integer(0)
        out$anchor_z <- numeric(0)
        out$abs_anchor_z <- numeric(0)
        out$max_abs_r <- numeric(0)
        return(out)
      }
      z_ancestry <- if (ancestry == "EUR") zfile$zscore_1 else zfile$zscore_2
      R_ancestry <- if (ancestry == "EUR") R_true$EUR else R_true$AFR
      abs_R <- abs(R_ancestry)
      diag(abs_R) <- 0
      max_abs_r <- apply(abs_R, 1L, max)
      cbind(tag, ancestry = ancestry, SNP_index = idx,
            RSID = zfile$RSID[idx], Signal = zfile$Signal[idx],
            anchor_z = z_ancestry[idx], abs_anchor_z = abs(z_ancestry[idx]),
            max_abs_r = max_abs_r[idx])
    }
    targets_all[[length(targets_all) + 1L]] <- rbind(
      target_rows("EUR", scen$truth$EUR),
      target_rows("AFR", scen$truth$AFR)
    )

    diag_res <- tryCatch(
      run_kriging_diagnostics(
        zfile$zscore_1, zfile$zscore_2,
        scen$R$EUR, scen$R$AFR,
        stats::median(zfile$N_1), stats::median(zfile$N_2)
      ),
      error = function(e) e
    )

    if (inherits(diag_res, "error")) {
      errors_all[[length(errors_all) + 1L]] <- data.frame(
        tag, error = conditionMessage(diag_res), stringsAsFactors = FALSE,
        check.names = FALSE
      )
      warning("Diagnostic failed: ", conditionMessage(diag_res))
    } else {
      metrics_all[[length(metrics_all) + 1L]] <- cbind(
        tag, build_metric_rows(diag_res, scen$truth, nrow(zfile), scen$truth_scope)
      )
      scores_all[[length(scores_all) + 1L]] <- cbind(tag, diag_res$score_summary)
      s_all[[length(s_all) + 1L]] <- cbind(tag, diag_res$s_summary)
      anchor_detection <- build_anchor_detection_rows(
        diag_res, scen$truth, zfile$zscore_1, zfile$zscore_2
      )
      if (nrow(anchor_detection)) {
        anchor_detection$RSID <- zfile$RSID[anchor_detection$SNP_index]
        anchor_detection$Signal <- zfile$Signal[anchor_detection$SNP_index]
        anchor_detection_all[[length(anchor_detection_all) + 1L]] <- cbind(
          tag, anchor_detection
        )
      }
    }

    timing_all[[length(timing_all) + 1L]] <- cbind(
      tag, minutes = as.numeric(difftime(Sys.time(), cell_start, units = "mins"))
    )
    cat(sprintf("  causal=%d h2=%d completed in %.2f min\n", causal_num, h2_num,
                as.numeric(difftime(Sys.time(), cell_start, units = "mins"))))
  }
}

bind_or_empty <- function(x) if (length(x)) data.table::rbindlist(x, fill = TRUE) else data.frame()
metrics <- bind_or_empty(metrics_all)
matrix_summary <- bind_or_empty(matrix_all)
score_summary <- bind_or_empty(scores_all)
s_summary <- bind_or_empty(s_all)
targets <- bind_or_empty(targets_all)
errors <- bind_or_empty(errors_all)
timing <- bind_or_empty(timing_all)
anchor_detection <- bind_or_empty(anchor_detection_all)
if (nrow(errors)) stop("Diagnostic failures; refusing to publish incomplete task")
scenario_metadata <- scenario_row
scenario_metadata$ld_connected_min_abs_r <- 0.1

save(metrics, matrix_summary, score_summary, s_summary, targets, errors, timing,
     anchor_detection, scenario_metadata, file = out_file)
cat("Wrote ", out_file, "\n", sep = "")
