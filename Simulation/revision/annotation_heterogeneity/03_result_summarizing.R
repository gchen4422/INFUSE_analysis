#!/usr/bin/env Rscript
.infuse_repo <- Sys.getenv("INFUSE_ANALYSIS_ROOT", "/path/to/INFUSE_analysis")
.infuse_project <- Sys.getenv("INFUSE_PROJECT", "/path/to/baseline_simulation")
.infuse_work <- Sys.getenv("INFUSE_REVISION", "/path/to/revision_output")

# Result summary for one ancestry-specific causal-prior sensitivity setting.
# This is a direct adaptation of the original 03_result_summarizing.R workflow.
# The necessary changes are current revision paths and result-object names,
# INFUSE in place of MAFPC, EUR/AFR labels, and omission of PAINTOR. The four
# output figures retain the original plotting functions and layout:
#   1. 95% credible-set size and power;
#   2. either-causal precision-recall and power at target FDR;
#   3. shared-causal precision-recall and power at target FDR;
#   4. ancestry-specific precision-recall and power at target FDR.
#
# PAINTOR is omitted because it was not run in this focused sensitivity study.

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(forcats)
  library(ggplot2)
  library(patchwork)
  library(ggpattern)
})

args <- commandArgs(trailingOnly = TRUE)
if (!length(args) %in% 2:3) {
  stop(
    "Usage: Rscript 03_result_summarizing.R ",
    "<scenario> <base_d> [result_tag]",
    call. = FALSE
  )
}

scenario <- args[[1L]]
base_d <- as.integer(args[[2L]])
result_tag <- if (length(args) == 3L) args[[3L]] else ""
valid_scenarios <- c(
  "baseline", "different_annotation", "different_weight", "combined"
)
if (!scenario %in% valid_scenarios) {
  stop("Unknown scenario: ", scenario, call. = FALSE)
}
if (!base_d %in% c(10L, 187L)) {
  stop("base_d must be 10 or 187.", call. = FALSE)
}
if (nzchar(result_tag) && !grepl("^[A-Za-z0-9_]+$", result_tag)) {
  stop("result_tag may contain only letters, numbers, and underscores.", call. = FALSE)
}

project_dir <- .infuse_project
source(file.path(.infuse_repo, "Simulation/revision/annotation_heterogeneity/utility.R"))
revision_dir <- .infuse_work
analysis_name <- "ancestry_specific_prior_sensitivity"
module_dir <- file.path(revision_dir, "ancestry_specific_eaf_prior_simulation")
input_d <- base_d + 2L
condition_dir <- file.path(
  module_dir,
  paste0("scenario_", scenario),
  paste0("base_d_", base_d, "_input_d_", input_d)
)
res_dir <- file.path(
  revision_dir, "res", analysis_name,
  paste0("scenario_", scenario), paste0("base_d_", base_d)
)
figure_dir <- file.path(
  revision_dir, "figures", analysis_name,
  paste0("scenario_", scenario), paste0("base_d_", base_d)
)
dir.create(res_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

scenario_labels <- c(
  baseline = "Baseline: shared annotation and weight",
  different_annotation = "Ancestry-specific EAF; shared weight",
  different_weight = "Shared annotation; ancestry-specific weight",
  combined = "Ancestry-specific EAF and weight"
)
scenario_label <- unname(scenario_labels[[scenario]])
causal_counts <- c(2L, 3L, 5L)
h2_values <- c(1e-4, 2e-4)
cutoffs <- seq(1, 0.01, by = -0.01)
target_fdrs <- c(0.05, 0.10, 0.50)
all_loci <- seq_len(100L)
result_dir_name <- if (nzchar(result_tag)) paste0("result_", result_tag) else "result"

clamp_pip <- function(x, label) {
  if (any(!is.finite(x))) {
    stop(label, " contains non-finite PIPs.", call. = FALSE)
  }
  if (any(x < -1e-8 | x > 1 + 1e-8)) {
    stop(label, " contains PIPs outside [0, 1].", call. = FALSE)
  }
  pmin(pmax(as.numeric(x), 0), 1)
}

get_config_pip <- function(fit, label, n_snp) {
  out <- as.matrix(fit$pip_config)
  if (!identical(dim(out), c(n_snp, 3L))) {
    stop(label, " must contain an nSNP x 3 configuration-PIP matrix.", call. = FALSE)
  }
  out
}

get_cs_indices <- function(x, n_snp, label) {
  selected <- unique(as.integer(unlist(x, use.names = FALSE)))
  selected <- selected[is.finite(selected)]
  if (length(selected) && any(selected < 1L | selected > n_snp)) {
    stop(label, " credible set contains an invalid SNP index.", call. = FALSE)
  }
  selected
}

credible_set_metrics <- function(selected, truth) {
  n_selected <- length(selected)
  n_causal <- sum(truth)
  n_true_selected <- if (n_selected) sum(truth[selected]) else 0L
  data.table(
    Power = if (n_causal) n_true_selected / n_causal else 0,
    FDR = if (n_selected) (n_selected - n_true_selected) / n_selected else 0,
    True_selected = n_true_selected,
    Size = n_selected
  )
}

# Exact implementation of the legacy strict threshold rule PIP > cutoff,
# computed from a sorted PIP vector rather than a variant-by-cutoff matrix.
compute_roc <- function(pip, truth) {
  pip <- clamp_pip(pip, "ROC input")
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

# Legacy SuSiE proxy for an ancestry-specific configuration: the target-
# ancestry PIP must exceed the cutoff while the other-ancestry PIP is below it.
compute_roc_susie_specific <- function(target_pip, other_pip, truth) {
  target_pip <- clamp_pip(target_pip, "SuSiE target-ancestry PIP")
  other_pip <- clamp_pip(other_pip, "SuSiE other-ancestry PIP")
  truth <- as.logical(truth)
  n_truth <- sum(truth)
  rows <- lapply(cutoffs, function(cutoff) {
    selected <- target_pip > cutoff & other_pip < cutoff
    true_positive <- sum(selected & truth)
    n_selected <- sum(selected)
    data.table(
      Cutoff = cutoff,
      Power = if (n_truth) true_positive / n_truth else 0,
      FDR = if (n_selected) (n_selected - true_positive) / n_selected else 0
    )
  })
  rbindlist(rows)
}

add_roc_metadata <- function(x, method, locus, h2_index, causal_index,
                             ancestry = NA_character_) {
  x[
    , `:=`(
      scenario = scenario,
      base_d = base_d,
      input_d = input_d,
      Method = method,
      Ancestry = ancestry,
      locus = locus,
      h2_index = h2_index,
      h2 = h2_values[h2_index],
      causal_index = causal_index,
      causal_count = causal_counts[causal_index]
    )
  ]
}

set_list <- vector("list", 600L * 3L)
either_list <- vector("list", 600L * 3L)
shared_list <- vector("list", 600L * 3L)
ancestry_list <- vector("list", 600L * 6L)
timing_list <- vector("list", 600L)
data_all_list <- vector("list", 600L)
compare_list <- vector("list", 600L)
set_counter <- either_counter <- shared_counter <- ancestry_counter <- 0L
result_counter <- 0L
infuse_metadata_reference <- NULL

for (causal_index in seq_len(3L)) {
  data_dir <- file.path(
    condition_dir, paste0("causal_num_", causal_index), "summary_data"
  )
  result_dir <- file.path(
    condition_dir, paste0("causal_num_", causal_index), result_dir_name
  )

  for (h2_index in seq_len(2L)) {
    for (locus in all_loci) {
      stem <- paste0(
        "CAUSAL_", causal_index,
        "_LOCI_", locus,
        "_h2_", h2_index
      )
      summary_path <- file.path(data_dir, stem)
      result_path <- file.path(result_dir, paste0("Res_", stem, ".RData"))
      required_paths <- c(summary_path, result_path)
      if (any(!file.exists(required_paths))) {
        stop(
          "Incomplete setting; missing: ",
          paste(required_paths[!file.exists(required_paths)], collapse = ", "),
          call. = FALSE
        )
      }

      zfile <- fread(summary_path, select = c("RSID", "Signal"))
      n_snp <- nrow(zfile)
      truth_any <- zfile$Signal != 0L
      truth_shared <- zfile$Signal == 3L
      truth_eur <- zfile$Signal == 1L
      truth_afr <- zfile$Signal == 2L

      result_env <- new.env(parent = emptyenv())
      load(result_path, envir = result_env)
      required_objects <- c(
        "MESuSiE_res", "INFUSE_res", "susie_EU", "susie_BB", "timing",
        "analysis_metadata"
      )
      absent <- required_objects[
        !vapply(required_objects, exists, logical(1), envir = result_env,
                inherits = FALSE)
      ]
      if (length(absent)) {
        stop(
          "Result is missing objects: ", paste(absent, collapse = ", "),
          " in ", result_path, call. = FALSE
        )
      }

      mesusie <- result_env$MESuSiE_res
      infuse <- result_env$INFUSE_res
      susie_eur <- result_env$susie_EU
      susie_afr <- result_env$susie_BB
      fit_metadata <- result_env$analysis_metadata
      if (!identical(fit_metadata$result_tag, result_tag)) {
        stop("Result-tag mismatch in ", result_path, call. = FALSE)
      }
      if (identical(result_tag, "current_github")) {
        required_threads <- c(
          OPENBLAS_NUM_THREADS = "1", OMP_NUM_THREADS = "1",
          MKL_NUM_THREADS = "1", VECLIB_MAXIMUM_THREADS = "1",
          NUMEXPR_NUM_THREADS = "1"
        )
        observed_threads <- unlist(fit_metadata$numerical_threads)
        if (!isTRUE(fit_metadata$installed_core_matches_github) ||
            !identical(fit_metadata$infuse_annotation_path,
                       "direct_annotation_optimization") ||
            !identical(observed_threads[names(required_threads)], required_threads)) {
          stop("Current-GitHub provenance check failed in ", result_path,
               call. = FALSE)
        }
      }
      metadata_signature <- fit_metadata[c(
        "result_tag", "infuse_package_version", "infuse_package_path",
        "infuse_github_checkout", "infuse_github_commit",
        "installed_core_matches_github", "infuse_annotation_path"
      )]
      if (is.null(infuse_metadata_reference)) {
        infuse_metadata_reference <- metadata_signature
      } else if (!identical(metadata_signature, infuse_metadata_reference)) {
        stop("INFUSE implementation metadata changed across results: ", result_path,
             call. = FALSE)
      }
      if (length(mesusie$pip) != n_snp || length(infuse$pip) != n_snp ||
          length(susie_eur$pip) != n_snp || length(susie_afr$pip) != n_snp) {
        stop("PIP and summary-statistic dimensions differ: ", result_path, call. = FALSE)
      }
      mesusie_config <- get_config_pip(mesusie, "MESuSiE", n_snp)
      infuse_config <- get_config_pip(infuse, "INFUSE", n_snp)

      cs_specs <- list(
        MESuSiE = get_cs_indices(mesusie$cs$cs, n_snp, "MESuSiE"),
        INFUSE = get_cs_indices(infuse$cs$cs, n_snp, "INFUSE"),
        SuSiE = get_cs_indices(
          c(susie_eur$sets$cs, susie_afr$sets$cs), n_snp, "SuSiE"
        )
      )
      for (method in names(cs_specs)) {
        set_counter <- set_counter + 1L
        set_list[[set_counter]] <- credible_set_metrics(
          cs_specs[[method]], truth_any
        )[
          , `:=`(
            scenario = scenario,
            base_d = base_d,
            input_d = input_d,
            Method = method,
            locus = locus,
            h2_index = h2_index,
            h2 = h2_values[h2_index],
            causal_index = causal_index,
            causal_count = causal_counts[causal_index]
          )
        ]
      }

      mesusie_set <- credible_set_metrics(cs_specs$MESuSiE, truth_any)
      infuse_set <- credible_set_metrics(cs_specs$INFUSE, truth_any)
      comparison_category <- if (
        infuse_set$Power >= mesusie_set$Power && infuse_set$FDR <= mesusie_set$FDR
      ) {
        1L
      } else if (
        infuse_set$Power >= mesusie_set$Power && infuse_set$FDR > mesusie_set$FDR
      ) {
        2L
      } else if (
        infuse_set$Power < mesusie_set$Power && infuse_set$FDR <= mesusie_set$FDR
      ) {
        3L
      } else {
        4L
      }

      either_specs <- list(
        MESuSiE = mesusie$pip,
        INFUSE = infuse$pip,
        SuSiE = pmax(susie_eur$pip, susie_afr$pip)
      )
      for (method in names(either_specs)) {
        either_counter <- either_counter + 1L
        either_list[[either_counter]] <- add_roc_metadata(
          compute_roc(either_specs[[method]], truth_any),
          method, locus, h2_index, causal_index
        )
      }

      if (any(truth_shared)) {
        shared_specs <- list(
          MESuSiE = mesusie_config[, 3L],
          INFUSE = infuse_config[, 3L],
          SuSiE = pmin(susie_eur$pip, susie_afr$pip)
        )
        for (method in names(shared_specs)) {
          shared_counter <- shared_counter + 1L
          shared_list[[shared_counter]] <- add_roc_metadata(
            compute_roc(shared_specs[[method]], truth_shared),
            method, locus, h2_index, causal_index
          )
        }
      }

      if (any(truth_eur)) {
        eur_specs <- list(
          MESuSiE = compute_roc(mesusie_config[, 1L], truth_eur),
          INFUSE = compute_roc(infuse_config[, 1L], truth_eur),
          SuSiE = compute_roc_susie_specific(
            susie_eur$pip, susie_afr$pip, truth_eur
          )
        )
        for (method in names(eur_specs)) {
          ancestry_counter <- ancestry_counter + 1L
          ancestry_list[[ancestry_counter]] <- add_roc_metadata(
            eur_specs[[method]], method, locus, h2_index, causal_index, "EUR"
          )
        }
      }

      if (any(truth_afr)) {
        afr_specs <- list(
          MESuSiE = compute_roc(mesusie_config[, 2L], truth_afr),
          INFUSE = compute_roc(infuse_config[, 2L], truth_afr),
          SuSiE = compute_roc_susie_specific(
            susie_afr$pip, susie_eur$pip, truth_afr
          )
        )
        for (method in names(afr_specs)) {
          ancestry_counter <- ancestry_counter + 1L
          ancestry_list[[ancestry_counter]] <- add_roc_metadata(
            afr_specs[[method]], method, locus, h2_index, causal_index, "AFR"
          )
        }
      }

      result_counter <- result_counter + 1L
      timing_list[[result_counter]] <- as.data.table(result_env$timing)
      compare_list[[result_counter]] <- data.table(
        scenario = scenario,
        base_d = base_d,
        input_d = input_d,
        causal_index = causal_index,
        causal_count = causal_counts[causal_index],
        h2_index = h2_index,
        h2 = h2_values[h2_index],
        locus = locus,
        Category = comparison_category
      )
      data_all_list[[result_counter]] <- data.table(
        SNP = zfile$RSID,
        Signal = zfile$Signal,
        scenario = scenario,
        base_d = base_d,
        input_d = input_d,
        h2_index = h2_index,
        h2 = h2_values[h2_index],
        causal_index = causal_index,
        causal_count = causal_counts[causal_index],
        locus = locus,
        MESuSiE_Either = clamp_pip(mesusie$pip, "MESuSiE either PIP"),
        MESuSiE_Shared = clamp_pip(mesusie_config[, 3L], "MESuSiE shared PIP"),
        MESuSiE_EUR = clamp_pip(mesusie_config[, 1L], "MESuSiE EUR-only PIP"),
        MESuSiE_AFR = clamp_pip(mesusie_config[, 2L], "MESuSiE AFR-only PIP"),
        INFUSE_Either = clamp_pip(infuse$pip, "INFUSE either PIP"),
        INFUSE_Shared = clamp_pip(infuse_config[, 3L], "INFUSE shared PIP"),
        INFUSE_EUR = clamp_pip(infuse_config[, 1L], "INFUSE EUR-only PIP"),
        INFUSE_AFR = clamp_pip(infuse_config[, 2L], "INFUSE AFR-only PIP"),
        SuSiE_Either = clamp_pip(
          pmax(susie_eur$pip, susie_afr$pip), "SuSiE either PIP"
        ),
        SuSiE_Shared = clamp_pip(
          pmin(susie_eur$pip, susie_afr$pip), "SuSiE shared PIP"
        ),
        SuSiE_EUR = clamp_pip(susie_eur$pip, "SuSiE EUR PIP"),
        SuSiE_AFR = clamp_pip(susie_afr$pip, "SuSiE AFR PIP")
      )
      rm(result_env, mesusie, infuse, susie_eur, susie_afr)
      if (result_counter %% 25L == 0L) {
        message(
          "Processed ", result_counter, "/600 results for scenario=", scenario,
          ", base_d=", base_d
        )
        invisible(gc())
      }
    }
  }
}

all_Set_data_dataframe <- rbindlist(set_list[seq_len(set_counter)])
either_locus_ROC <- rbindlist(either_list[seq_len(either_counter)])
shared_locus_ROC <- rbindlist(shared_list[seq_len(shared_counter)])
ancestry_locus_ROC <- rbindlist(ancestry_list[seq_len(ancestry_counter)])
time_all <- rbindlist(timing_list[seq_len(result_counter)], fill = TRUE)
data_all <- rbindlist(data_all_list[seq_len(result_counter)])
compare_table <- rbindlist(compare_list[seq_len(result_counter)])
summary_table <- compare_table[, .(Count = .N), by = Category][order(Category)]
all_Set_data_dataframe[, causal_num := causal_index]
either_locus_ROC[, causal_num := causal_index]
shared_locus_ROC[, causal_num := causal_index]
ancestry_locus_ROC[, causal_num := causal_index]
data_all[, causal_num := causal_index]
rm(
  set_list, either_list, shared_list, ancestry_list, timing_list,
  data_all_list, compare_list
)
invisible(gc())

roc_group <- c(
  "scenario", "base_d", "input_d", "Method", "h2_index", "h2",
  "causal_index", "causal_count", "Cutoff"
)
either_all_ROC_data_dataframe <- either_locus_ROC[
  , .(
    Power = mean(Power),
    FDR = mean(FDR),
    N_loci = uniqueN(locus)
  ),
  by = roc_group
]
shared_all_ROC_data_dataframe <- shared_locus_ROC[
  , .(
    Power = mean(Power),
    FDR = mean(FDR),
    N_loci = uniqueN(locus)
  ),
  by = roc_group
]
ancestry_all_ROC_data_dataframe <- ancestry_locus_ROC[
  , .(
    Power = mean(Power),
    FDR = mean(FDR),
    N_loci = uniqueN(locus)
  ),
  by = c(roc_group, "Ancestry")
]
either_all_ROC_data_dataframe[, causal_num := causal_index]
shared_all_ROC_data_dataframe[, causal_num := causal_index]
ancestry_all_ROC_data_dataframe[, causal_num := causal_index]

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

# Match the established FDR-power summary rule: choose the point on the mean
# locus-level curve closest to each target FDR, breaking ties by higher power.
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
FDR_Power_either[, causal_num := causal_index]
FDR_Power_shared[, causal_num := causal_index]
FDR_Power_ancestry[, causal_num := causal_index]

method_levels <- c("SuSiE", "MESuSiE", "INFUSE")
legacy_h2_labels <- c(
  "~h^2 == 10^-4",
  "~h^2 == 2%*%10^-4"
)
legacy_causal_labels <- c(
  "Num~Causal == 2",
  "Num~Causal == 3",
  "Num~Causal == 5"
)

# Plot through the same functions used by the original number-5 workflow.
# Only the method subset (SuSiE, MESuSiE, INFUSE) and ancestry names differ.
as_legacy_plot_data <- function(x) {
  x <- as.data.frame(copy(x))
  x$Method <- factor(x$Method, levels = method_levels)
  x$h2 <- factor(
    x$h2_index, levels = 1:2, labels = legacy_h2_labels
  )
  x$causal_num <- factor(
    x$causal_index, levels = 1:3, labels = legacy_causal_labels
  )
  x
}

set_plot_data <- as_legacy_plot_data(all_Set_data_dataframe)
power_plot_data <- as_legacy_plot_data(set_power_summary)
set_plot_data$Size <- log2(set_plot_data$Size + 1)
upper_limit <- max(set_plot_data$Size, na.rm = TRUE)

p_set_size <- Set_Size_fun(set_plot_data, upper_limit) +
  ylab("log2(Set Size + 1)")
p_set_power <- Set_Power_fun(power_plot_data)
p_cs <- (p_set_size / p_set_power) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(size = 7, face = "bold"))

make_legacy_fdr_data <- function(x) {
  out <- as_legacy_plot_data(x)
  out$FDR <- factor(
    out$Target_FDR,
    levels = target_fdrs,
    labels = c("0.05", "0.1", "0.5")
  )
  out$Power <- round(out$Power, 2)
  out
}

either_roc_plot_data <- as_legacy_plot_data(either_all_ROC_data_dataframe)
either_fdr_plot_data <- make_legacy_fdr_data(FDR_Power_either)
p_either_roc <- ROC_shared_fun(either_roc_plot_data)
p_either_fdr <- FDR_Power_shared_fun(either_fdr_plot_data) +
  ylim(0, min(ceiling(max(either_fdr_plot_data$Power) * 10) / 10 + 0.1, 1))
p_either <- (p_either_roc / p_either_fdr) +
  plot_annotation(tag_levels = "a") +
  plot_layout(guides = "collect", heights = c(2, 1)) &
  theme(
    legend.position = "bottom",
    plot.tag = element_text(size = 7, face = "bold")
  )

shared_roc_plot_data <- as_legacy_plot_data(shared_all_ROC_data_dataframe)
shared_fdr_plot_data <- make_legacy_fdr_data(FDR_Power_shared)
p_shared_roc <- ROC_shared_fun(shared_roc_plot_data)
p_shared_fdr <- FDR_Power_shared_fun(shared_fdr_plot_data) +
  ylim(0, min(ceiling(max(shared_fdr_plot_data$Power) * 10) / 10 + 0.1, 1))
p_shared <- (p_shared_roc / p_shared_fdr) +
  plot_annotation(tag_levels = "a") +
  plot_layout(guides = "collect", heights = c(2, 1)) &
  theme(
    legend.position = "bottom",
    plot.tag = element_text(size = 7, face = "bold")
  )

ancestry_roc_plot_data <- as_legacy_plot_data(ancestry_all_ROC_data_dataframe)
ancestry_roc_plot_data$Ancestry <- factor(
  ancestry_roc_plot_data$Ancestry, levels = c("EUR", "AFR")
)
p_ancestry_roc <- ROC_ancestry_fun(ancestry_roc_plot_data)

ancestry_fdr_plot_data <- make_legacy_fdr_data(FDR_Power_ancestry)
ancestry_fdr_plot_data$Ancestry <- factor(
  ancestry_fdr_plot_data$Ancestry, levels = c("EUR", "AFR")
)
p_ancestry_fdr <- FDR_Power_ancestry_fun(ancestry_fdr_plot_data) +
  ylim(0, min(ceiling(max(ancestry_fdr_plot_data$Power) * 10) / 10 + 0.1, 1))

p_ancestry <- (p_ancestry_roc / p_ancestry_fdr) +
  plot_annotation(tag_levels = "a") +
  plot_layout(heights = c(2, 1)) &
  theme(plot.tag = element_text(size = 7, face = "bold"))


res_out <- paste0(
  "shared_50_sample_mlk_ancestry_specific_prior_", scenario,
  "_base_d_", base_d,
  if (nzchar(result_tag)) paste0("_", result_tag) else "_spca",
  ".RData"
)
summary_path <- file.path(res_dir, res_out)
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
  time_all,
  data_all,
  compare_table,
  summary_table,
  result_tag,
  p_cs,
  p_either,
  p_shared,
  p_ancestry,
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

save_figure <- function(plot, stem, width, height) {
  pdf_path <- file.path(figure_dir, paste0(stem, ".pdf"))
  png_path <- file.path(figure_dir, paste0(stem, ".png"))
  ggsave(pdf_path, plot, width = width, height = height, units = "in")
  ggsave(png_path, plot, width = width, height = height, units = "in", dpi = 500)
}

save_figure(p_cs, "Set_size_power", 270 / 25.4, 180 / 25.4)
save_figure(p_either, "PIP_ROC_FDR_Power_Either", 270 / 25.4, 270 / 25.4)
save_figure(p_shared, "PIP_ROC_FDR_Power_Shared", 270 / 25.4, 270 / 25.4)
save_figure(p_ancestry, "PIP_ROC_FDR_Power_Ancestry", 270 / 25.4, 270 / 25.4)

message("Saved manuscript-style summary: ", summary_path)
message("Saved four figures under: ", figure_dir)
