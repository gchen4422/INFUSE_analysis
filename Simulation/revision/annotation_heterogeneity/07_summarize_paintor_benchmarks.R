#!/usr/bin/env Rscript
.infuse_repo <- Sys.getenv("INFUSE_ANALYSIS_ROOT", "/path/to/INFUSE_analysis")
.infuse_project <- Sys.getenv("INFUSE_PROJECT", "/path/to/baseline_simulation")
.infuse_work <- Sys.getenv("INFUSE_REVISION", "/path/to/revision_output")

# Summarize Paintor, Paintor+, and Paintor+H across the four ancestry-prior
# sensitivity scenarios. Metric definitions follow the root
# 03_result_summarizing.R workflow. In particular, PAINTOR credible sets are
# constructed after normalizing locus PIPs and retaining variants while the
# ordered cumulative PIP is strictly below 0.95.

suppressPackageStartupMessages(library(data.table))

project_dir <- .infuse_project
revision_dir <- .infuse_work
module_dir <- file.path(
  revision_dir, "ancestry_specific_eaf_prior_simulation"
)
benchmark_root <- file.path(module_dir, "paintor_benchmarks")
source(file.path(.infuse_repo, "Simulation/revision/annotation_heterogeneity/paintor_output_contract.R"))
res_dir <- file.path(
  revision_dir, "res", "ancestry_specific_prior_sensitivity",
  "across_scenarios", "paintor_benchmarks"
)
dir.create(res_dir, recursive = TRUE, showWarnings = FALSE)

scenarios <- c(
  "baseline", "different_annotation", "different_weight", "combined"
)
method_specs <- data.table(
  method_dir = c("paintor", "paintor_plus_10", "paintor_plus_187"),
  Method = c("Paintor", "Paintor+", "Paintor+H"),
  base_d = c(0L, 10L, 187L),
  input_d = c(0L, 12L, 189L)
)
causal_counts <- c(2L, 3L, 5L)
h2_values <- c(1e-4, 2e-4)
cutoffs <- seq(1, 0.01, by = -0.01)
target_fdrs <- c(0.05, 0.10, 0.50)

clamp_pip <- function(x, label) {
  if (any(!is.finite(x))) stop(label, " contains non-finite PIPs.")
  if (any(x < -1e-8 | x > 1 + 1e-8)) {
    stop(label, " contains PIPs outside [0, 1].")
  }
  pmin(pmax(as.numeric(x), 0), 1)
}

# legacy_paintor_set is shared with the runner/checker contract. Zero mass is
# explicitly an empty selection, not an estimated 95% credible set.

# Exact strict PIP > cutoff rule used in the current ancestry-prior summary.
compute_roc <- function(pip, truth) {
  pip <- clamp_pip(pip, "PAINTOR ROC input")
  truth <- as.logical(truth)
  n_truth <- sum(truth)
  order_index <- order(pip)
  sorted_pip <- pip[order_index]
  cumulative_truth <- cumsum(truth[order_index])
  n_at_or_below <- findInterval(cutoffs, sorted_pip)
  truth_at_or_below <- ifelse(
    n_at_or_below > 0L,
    cumulative_truth[pmax(n_at_or_below, 1L)],
    0L
  )
  n_selected <- length(pip) - n_at_or_below
  true_positive <- n_truth - truth_at_or_below
  false_positive <- n_selected - true_positive
  data.table(
    Cutoff = cutoffs,
    Power = if (n_truth) true_positive / n_truth else 0,
    FDR = fifelse(n_selected > 0L, false_positive / n_selected, 0)
  )
}

add_metadata <- function(x, scenario, method_row, locus, h2_index,
                         causal_index, ancestry = NA_character_) {
  x[, `:=`(
    scenario = scenario,
    base_d = method_row$base_d,
    input_d = method_row$input_d,
    Method = method_row$Method,
    Ancestry = ancestry,
    locus = locus,
    h2_index = h2_index,
    h2 = h2_values[h2_index],
    causal_index = causal_index,
    causal_count = causal_counts[causal_index]
  )]
}

set_list <- vector("list", 4L * 600L * 3L)
either_list <- vector("list", 4L * 600L * 3L)
shared_list <- vector("list", 4L * 600L * 3L)
ancestry_list <- vector("list", 4L * 600L * 3L * 2L)
timing_list <- vector("list", 4L * 600L * 3L)
set_counter <- either_counter <- shared_counter <- ancestry_counter <- 0L
timing_counter <- 0L

for (scenario in scenarios) {
  condition_dir <- file.path(
    module_dir, paste0("scenario_", scenario), "base_d_10_input_d_12"
  )
  for (causal_index in seq_len(3L)) {
    data_dir <- file.path(
      condition_dir, paste0("causal_num_", causal_index), "summary_data"
    )
    for (h2_index in seq_len(2L)) {
      for (locus in seq_len(100L)) {
        stem <- paste0(
          "CAUSAL_", causal_index,
          "_LOCI_", locus,
          "_h2_", h2_index
        )
        summary_path <- file.path(data_dir, stem)
        if (!file.exists(summary_path)) stop("Missing summary: ", summary_path)
        zfile <- fread(summary_path, select = c("RSID", "Signal"))
        if (anyDuplicated(zfile$RSID)) stop("Duplicated RSIDs in ", summary_path)
        truth_any <- zfile$Signal != 0L
        truth_shared <- zfile$Signal == 3L
        truth_eur <- zfile$Signal == 1L
        truth_afr <- zfile$Signal == 2L

        for (method_i in seq_len(nrow(method_specs))) {
          method_row <- method_specs[method_i]
          result_dir <- file.path(
            benchmark_root, paste0("scenario_", scenario),
            method_row$method_dir, paste0("causal_num_", causal_index),
            "result"
          )
          result_path <- file.path(result_dir, paste0(stem, ".mcmc.paintor"))
          metadata_path <- file.path(result_dir, paste0(stem, ".metadata.tsv"))
          if (!file.exists(result_path) || !file.exists(metadata_path)) {
            stop("Missing PAINTOR result or metadata for: ", result_path)
          }
          result <- fread(result_path)
          if (!all(c("RSID", "Posterior_Prob") %in% names(result)) ||
              nrow(result) != nrow(zfile) || anyDuplicated(result$RSID) ||
              !setequal(result$RSID, zfile$RSID)) {
            stop("Malformed or misaligned PAINTOR result: ", result_path)
          }
          result <- result[match(zfile$RSID, result$RSID)]
          if (!identical(as.character(result$RSID), as.character(zfile$RSID))) {
            stop("Failed to align PAINTOR result: ", result_path)
          }
          pip <- clamp_pip(result$Posterior_Prob, result_path)
          diagnostic <- paintor_log_diagnostics(file.path(result_dir,paste0(stem,".paintor.log")))
          if(diagnostic$optimizer_error_logged) stop("Explicit PAINTOR optimizer error: ",result_path)
          current_output_status <- paintor_pip_status(pip)

          set_counter <- set_counter + 1L
          set_list[[set_counter]] <- legacy_paintor_set(pip, truth_any)[
            , `:=`(
              scenario = scenario,
              base_d = method_row$base_d,
              input_d = method_row$input_d,
              Method = method_row$Method,
              locus = locus,
              h2_index = h2_index,
              h2 = h2_values[h2_index],
              causal_index = causal_index,
              causal_count = causal_counts[causal_index],
              output_status = current_output_status
            )
          ]

          either_counter <- either_counter + 1L
          either_list[[either_counter]] <- add_metadata(
            compute_roc(pip, truth_any), scenario, method_row,
            locus, h2_index, causal_index
          )
          if (any(truth_shared)) {
            shared_counter <- shared_counter + 1L
            shared_list[[shared_counter]] <- add_metadata(
              compute_roc(pip, truth_shared), scenario, method_row,
              locus, h2_index, causal_index
            )
          }
          if (any(truth_eur)) {
            ancestry_counter <- ancestry_counter + 1L
            ancestry_list[[ancestry_counter]] <- add_metadata(
              compute_roc(pip, truth_eur), scenario, method_row,
              locus, h2_index, causal_index, "EUR"
            )
          }
          if (any(truth_afr)) {
            ancestry_counter <- ancestry_counter + 1L
            ancestry_list[[ancestry_counter]] <- add_metadata(
              compute_roc(pip, truth_afr), scenario, method_row,
              locus, h2_index, causal_index, "AFR"
            )
          }

          metadata <- fread(metadata_path)
          expected_metadata <- list(
            scenario = scenario,
            method = method_row$method_dir,
            base_annotation_dimension = method_row$base_d,
            fitted_annotation_dimension = method_row$input_d,
            locus = locus,
            h2_index = h2_index,
            causal_index = causal_index
          )
          for (field in names(expected_metadata)) {
            if (nrow(metadata) != 1L ||
                !identical(as.character(metadata[[field]]),
                           as.character(expected_metadata[[field]]))) {
              stop("Metadata mismatch for ", field, " in ", metadata_path)
            }
          }
          timing_counter <- timing_counter + 1L
          metadata[, `:=`(output_status=current_output_status,total_pip=sum(pip),
                           n_positive_pip=sum(pip>0))]
          for(nm in names(diagnostic)) metadata[,(nm):=diagnostic[[nm]]]
          timing_list[[timing_counter]] <- metadata
        }
      }
    }
  }
}

all_Set_data_dataframe <- rbindlist(set_list[seq_len(set_counter)])
either_locus_ROC <- rbindlist(either_list[seq_len(either_counter)])
shared_locus_ROC <- rbindlist(shared_list[seq_len(shared_counter)])
ancestry_locus_ROC <- rbindlist(ancestry_list[seq_len(ancestry_counter)])
timing_all <- rbindlist(timing_list[seq_len(timing_counter)], fill = TRUE)
output_status_summary <- timing_all[,.(N=.N),by=.(scenario,method,output_status)]
fwrite(timing_all,file.path(res_dir,"paintor_output_diagnostics.tsv"),sep="\t")
fwrite(output_status_summary,file.path(res_dir,"paintor_output_status_summary.tsv"),sep="\t")

roc_group <- c(
  "scenario", "base_d", "input_d", "Method", "h2_index", "h2",
  "causal_index", "causal_count", "Cutoff"
)
either_all_ROC_data_dataframe <- either_locus_ROC[
  , .(Power = mean(Power), FDR = mean(FDR), N_loci = uniqueN(locus)),
  by = roc_group
]
shared_all_ROC_data_dataframe <- shared_locus_ROC[
  , .(Power = mean(Power), FDR = mean(FDR), N_loci = uniqueN(locus)),
  by = roc_group
]
ancestry_all_ROC_data_dataframe <- ancestry_locus_ROC[
  , .(Power = mean(Power), FDR = mean(FDR), N_loci = uniqueN(locus)),
  by = c(roc_group, "Ancestry")
]

set_power_summary <- all_Set_data_dataframe[
  , .(
    Power = mean(Power),
    FDR = mean(FDR),
    Median_size = median(Size),
    Mean_size = mean(Size),
    N_loci = .N
  ),
  by = .(
    scenario, base_d, input_d, Method,
    h2_index, h2, causal_index, causal_count
  )
]
set_power_summary[, Power_name := Power]

fdr_power_summary <- function(roc_data, ancestry = FALSE) {
  group_columns <- c(
    "scenario", "base_d", "input_d", "Method",
    "h2_index", "h2", "causal_index", "causal_count"
  )
  if (ancestry) group_columns <- c(group_columns, "Ancestry")
  groups <- unique(roc_data[, ..group_columns])
  output <- vector("list", nrow(groups) * length(target_fdrs))
  counter <- 0L
  for (i in seq_len(nrow(groups))) {
    key <- groups[i]
    curve <- roc_data[key, on = group_columns]
    for (target in target_fdrs) {
      distance <- abs(curve$FDR - target)
      candidates <- curve[distance == min(distance)]
      selected <- candidates[which.max(Power)]
      counter <- counter + 1L
      output[[counter]] <- cbind(
        key,
        data.table(
          Target_FDR = target,
          Power = selected$Power,
          Achieved_FDR = selected$FDR,
          Cutoff = selected$Cutoff
        )
      )
    }
  }
  rbindlist(output[seq_len(counter)])
}

FDR_Power_either <- fdr_power_summary(either_all_ROC_data_dataframe)
FDR_Power_shared <- fdr_power_summary(shared_all_ROC_data_dataframe)
FDR_Power_ancestry <- fdr_power_summary(
  ancestry_all_ROC_data_dataframe, ancestry = TRUE
)


summary_path <- file.path(
  res_dir, "paintor_benchmarks_across_scenarios.RData"
)
save(
  all_Set_data_dataframe,
  set_power_summary,
  either_locus_ROC,
  shared_locus_ROC,
  ancestry_locus_ROC,
  either_all_ROC_data_dataframe,
  shared_all_ROC_data_dataframe,
  ancestry_all_ROC_data_dataframe,
  FDR_Power_either,
  FDR_Power_shared,
  FDR_Power_ancestry,
  timing_all,
  output_status_summary,
  file = summary_path,
  compress = "xz"
)

fwrite(all_Set_data_dataframe, file.path(res_dir, "credible_set_locus.tsv"), sep = "\t")
fwrite(set_power_summary, file.path(res_dir, "credible_set_summary.tsv"), sep = "\t")
fwrite(either_all_ROC_data_dataframe, file.path(res_dir, "pip_either_roc.tsv"), sep = "\t")
fwrite(shared_all_ROC_data_dataframe, file.path(res_dir, "pip_shared_roc.tsv"), sep = "\t")
fwrite(ancestry_all_ROC_data_dataframe, file.path(res_dir, "pip_ancestry_roc.tsv"), sep = "\t")
fwrite(FDR_Power_either, file.path(res_dir, "pip_either_fdr_power.tsv"), sep = "\t")
fwrite(FDR_Power_shared, file.path(res_dir, "pip_shared_fdr_power.tsv"), sep = "\t")
fwrite(FDR_Power_ancestry, file.path(res_dir, "pip_ancestry_fdr_power.tsv"), sep = "\t")
fwrite(timing_all, file.path(res_dir, "paintor_run_metadata.tsv"), sep = "\t")

message("Saved PAINTOR benchmark summary: ", summary_path)
