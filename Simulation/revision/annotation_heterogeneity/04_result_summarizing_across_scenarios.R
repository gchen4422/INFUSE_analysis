#!/usr/bin/env Rscript
.infuse_repo <- Sys.getenv("INFUSE_ANALYSIS_ROOT", "/path/to/INFUSE_analysis")
.infuse_project <- Sys.getenv("INFUSE_PROJECT", "/path/to/baseline_simulation")
.infuse_work <- Sys.getenv("INFUSE_REVISION", "/path/to/revision_output")

# Cross-scenario performance summary for the ancestry-specific causal-prior
# sensitivity analysis. Run d=10 and d=187 separately. This script reads the
# four per-setting number-5 summaries, retains their metric definitions, and
# uses the same colors and visual theme as the original simulation figures.

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(ggpattern)
})

args <- commandArgs(trailingOnly = TRUE)
if (!length(args) %in% 1:2) {
  stop("Usage: Rscript 04_result_summarizing_across_scenarios.R <base_d> [result_tag]")
}

base_d <- as.integer(args[[1L]])
result_tag <- if (length(args) == 2L) args[[2L]] else ""
if (!base_d %in% c(10L, 187L)) stop("base_d must be 10 or 187.")
if (nzchar(result_tag) && !grepl("^[A-Za-z0-9_]+$", result_tag)) {
  stop("result_tag may contain only letters, numbers, and underscores.")
}

simulation_dir <- .infuse_project
revision_dir <- .infuse_work
analysis_dir <- file.path(revision_dir, "ancestry_specific_eaf_prior_simulation")
source(file.path(.infuse_repo, "Simulation/revision/annotation_heterogeneity/utility.R"))

analysis_name <- "ancestry_specific_prior_sensitivity"
summary_root <- file.path(revision_dir, "res", analysis_name)
plot_dir <- file.path(
  revision_dir, "figures", analysis_name, "across_scenarios",
  paste0("base_d_", base_d)
)
res_dir <- file.path(
  revision_dir, "res", analysis_name, "across_scenarios",
  paste0("base_d_", base_d)
)
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(res_dir, recursive = TRUE, showWarnings = FALSE)

scenarios <- c(
  "baseline", "different_annotation", "different_weight", "combined"
)
scenario_labels <- c(
  baseline = "Baseline",
  different_annotation = "Different\nannotation",
  different_weight = "Different\nweight",
  combined = "Both"
)
method_levels <- c("SuSiE", "MESuSiE", "INFUSE")
method_colors <- c(
  "SuSiE" = "#fbe79e",
  "MESuSiE" = "#84c3b7",
  "INFUSE" = "#71b8ed"
)
legacy_h2_labels <- c(
  "~h^2 == 10^-4",
  "~h^2 == 2%*%10^-4"
)
legacy_causal_labels <- c(
  "Num~Causal == 2",
  "Num~Causal == 3",
  "Num~Causal == 5"
)
report_fdr <- 0.05

set_list <- list()
set_summary_list <- list()
either_list <- list()
shared_list <- list()
ancestry_list <- list()
either_roc_list <- list()

for (scenario in scenarios) {
  summary_file <- file.path(
    summary_root,
    paste0("scenario_", scenario),
    paste0("base_d_", base_d),
    paste0(
      "shared_50_sample_mlk_ancestry_specific_prior_", scenario,
      "_base_d_", base_d,
      if (nzchar(result_tag)) paste0("_", result_tag) else "_spca",
      ".RData"
    )
  )
  if (!file.exists(summary_file)) stop("Missing number-5 summary: ", summary_file)

  result_env <- new.env(parent = emptyenv())
  load(summary_file, envir = result_env)
  required <- c(
    "all_Set_data_dataframe", "set_power_summary",
    "either_all_ROC_data_dataframe",
    "FDR_Power_either", "FDR_Power_shared", "FDR_Power_ancestry",
    "result_tag"
  )
  missing_objects <- setdiff(required, ls(result_env))
  if (length(missing_objects)) {
    stop("Missing objects in ", summary_file, ": ", paste(missing_objects, collapse = ", "))
  }
  if (!identical(result_env$result_tag, result_tag)) {
    stop("Unexpected result tag in ", summary_file)
  }
  cells <- unique(as.data.table(result_env$all_Set_data_dataframe)[,
    .(locus, h2_index, causal_index)])
  if (nrow(cells) != 600L) {
    stop("Expected 600 results in ", summary_file)
  }

  set_list[[scenario]] <- copy(as.data.table(result_env$all_Set_data_dataframe))
  set_summary_list[[scenario]] <- copy(as.data.table(result_env$set_power_summary))
  either_list[[scenario]] <- copy(as.data.table(result_env$FDR_Power_either))
  shared_list[[scenario]] <- copy(as.data.table(result_env$FDR_Power_shared))
  ancestry_list[[scenario]] <- copy(as.data.table(result_env$FDR_Power_ancestry))
  either_roc_list[[scenario]] <- copy(
    as.data.table(result_env$either_all_ROC_data_dataframe)
  )
  rm(result_env)
  invisible(gc())
}

all_Set_data_dataframe <- rbindlist(set_list, idcol = "source_scenario")
set_power_summary <- rbindlist(set_summary_list, idcol = "source_scenario")
FDR_Power_either <- rbindlist(either_list, idcol = "source_scenario")
FDR_Power_shared <- rbindlist(shared_list, idcol = "source_scenario")
FDR_Power_ancestry <- rbindlist(ancestry_list, idcol = "source_scenario")
either_all_ROC_data_dataframe <- rbindlist(
  either_roc_list, idcol = "source_scenario"
)
rm(
  set_list, set_summary_list, either_list, shared_list, ancestry_list,
  either_roc_list
)

add_plot_factors <- function(x) {
  x <- as.data.frame(copy(x))
  x$Scenario <- factor(
    x$source_scenario, levels = scenarios,
    labels = unname(scenario_labels[scenarios])
  )
  x$Method <- factor(x$Method, levels = method_levels)
  x$h2 <- factor(x$h2_index, levels = 1:2, labels = legacy_h2_labels)
  x$causal_num <- factor(
    x$causal_index, levels = 1:3, labels = legacy_causal_labels
  )
  x
}

cross_scenario_theme <- function() {
  list(
    theme_bw(),
    custom_theme(),
    theme(
      axis.text.x = element_text(size = 8, face = "bold"),
      legend.position = "bottom"
    ),
    annotate(
      "rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = -Inf,
      colour = "black", size = 0.2
    )
  )
}

################################################################
#
#        95% credible-set size and power across scenarios
#
################################################################
set_plot_data <- add_plot_factors(all_Set_data_dataframe)
set_plot_data$log2_Size <- log2(set_plot_data$Size + 1)
power_plot_data <- add_plot_factors(set_power_summary)

p_size_box <- ggplot(
  set_plot_data,
  aes(x = Scenario, y = log2_Size, fill = Method)
) +
  geom_boxplot(
    position = position_dodge2(width = 0.85, preserve = "single"),
    width = 0.72, outlier.size = 0.1, fatten = 0.5, color = "darkgray"
  ) +
  scale_fill_manual(values = method_colors, guide = "none") +
  ylab("log2(Set Size + 1)") +
  xlab("") +
  facet_grid(h2 ~ causal_num, labeller = label_parsed) +
  cross_scenario_theme()

p_power_bar <- ggplot(
  power_plot_data,
  aes(x = Scenario, y = Power_name, fill = Method)
) +
  geom_bar(
    stat = "identity", position = position_dodge(width = 0.84),
    width = 0.72
  ) +
  geom_text(
    aes(label = round(Power_name, 2), group = Method),
    position = position_dodge(width = 0.84), vjust = -0.4,
    size = 5 * 4 / 14
  ) +
  scale_fill_manual(values = method_colors) +
  ylab("Power") +
  xlab("") +
  scale_y_continuous(
    limits = c(0, 1.15), breaks = c(0, 0.25, 0.5, 0.75, 1),
    expand = expansion(mult = c(0, 0))
  ) +
  facet_grid(h2 ~ causal_num, labeller = label_parsed) +
  cross_scenario_theme()

size_power <- (p_size_box / p_power_bar) +
  plot_annotation(tag_levels = "a") +
  plot_layout(guides = "collect") &
  theme(
    plot.tag = element_text(size = 7, face = "bold"),
    legend.position = "bottom"
  )

################################################################
#
#        Power at FDR = 0.05 across scenarios
#
################################################################
make_signal_power_plot <- function(x) {
  plot_data <- add_plot_factors(x)
  plot_data <- plot_data[abs(plot_data$Target_FDR - report_fdr) < 1e-10, ]
  plot_data$Power <- round(plot_data$Power, 2)
  ggplot(plot_data, aes(x = Scenario, y = Power, fill = Method)) +
    geom_bar(
      stat = "identity", position = position_dodge(width = 0.84),
      width = 0.72
    ) +
    geom_text(
      aes(label = Power, group = Method),
    position = position_dodge(width = 0.84), vjust = -0.2,
      size = 5 * 4 / 14
    ) +
    scale_fill_manual(values = method_colors) +
    ylab("Power") +
    xlab("FDR = 0.05") +
    scale_y_continuous(
      limits = c(0, 1.15), breaks = c(0, 0.25, 0.5, 0.75, 1),
      expand = expansion(mult = c(0, 0))
    ) +
    facet_grid(h2 ~ causal_num, labeller = label_parsed) +
    cross_scenario_theme()
}

p_power_either <- make_signal_power_plot(FDR_Power_either)
p_power_shared <- make_signal_power_plot(FDR_Power_shared)

ancestry_plot_data <- add_plot_factors(FDR_Power_ancestry)
ancestry_plot_data <- ancestry_plot_data[
  abs(ancestry_plot_data$Target_FDR - report_fdr) < 1e-10,
]
ancestry_plot_data$Ancestry <- factor(
  ancestry_plot_data$Ancestry, levels = c("EUR", "AFR")
)
p_power_ancestry <- ggplot(
  ancestry_plot_data,
  aes(
    x = Scenario, y = Power, fill = Method, pattern = Ancestry,
    group = interaction(Method, Ancestry)
  )
) +
  geom_col_pattern(
    position = position_dodge2(width = 0.9, preserve = "single"),
    width = 0.78, pattern_density = 0.01, pattern_spacing = 0.05,
    pattern_key_scale_factor = 0.3, pattern_fill = "#DEDBD2",
    pattern_colour = "ivory3"
  ) +
  scale_fill_manual(name = "Method", values = method_colors) +
  scale_pattern_manual(values = c("EUR" = "stripe", "AFR" = "circle")) +
  ylab("Power") +
  xlab("FDR = 0.05") +
  scale_y_continuous(
  limits = c(0, 1.15), breaks = c(0, 0.25, 0.5, 0.75, 1),
    expand = expansion(mult = c(0, 0))
  ) +
  facet_grid(h2 ~ causal_num, labeller = label_parsed) +
  cross_scenario_theme() +
  guides(pattern = guide_legend(override.aes = list(fill = "white", order = 1)))

################################################################
#
#        Save result and figures
#
################################################################
performance_at_fdr_005 <- rbindlist(
  list(
    copy(FDR_Power_either)[, Signal := "Either"],
    copy(FDR_Power_shared)[, Signal := "Shared"],
    copy(FDR_Power_ancestry)[, Signal := "Ancestry-specific"]
  ),
  fill = TRUE
)[abs(Target_FDR - report_fdr) < 1e-10]

################################################################
#
#        Averaged summary across six h2 x causal-count settings
#
################################################################
overall_set_plot_data <- add_plot_factors(all_Set_data_dataframe)
overall_set_plot_data$log2_Size <- log2(overall_set_plot_data$Size + 1)

overall_set_power <- set_power_summary[
  , .(Power = mean(Power_name)),
  by = .(source_scenario, Method)
]
overall_set_power_plot <- as.data.frame(copy(overall_set_power))
overall_set_power_plot$Scenario <- factor(
  overall_set_power_plot$source_scenario,
  levels = scenarios,
  labels = unname(scenario_labels[scenarios])
)
overall_set_power_plot$Method <- factor(
  overall_set_power_plot$Method, levels = method_levels
)

p_size_box_overall <- ggplot(
  overall_set_plot_data,
  aes(x = Scenario, y = log2_Size, fill = Method)
) +
  geom_boxplot(
    position = position_dodge2(width = 0.85, preserve = "single"),
    width = 0.72, outlier.size = 0.1, fatten = 0.5, color = "darkgray"
  ) +
  scale_fill_manual(values = method_colors, guide = "none") +
  ylab("log2(Set Size + 1)") +
  xlab("") +
  cross_scenario_theme()

p_power_bar_overall <- ggplot(
  overall_set_power_plot,
  aes(x = Scenario, y = Power, fill = Method)
) +
  geom_bar(
    stat = "identity", position = position_dodge(width = 0.84),
    width = 0.72
  ) +
  geom_text(
    aes(label = round(Power, 2), group = Method),
    position = position_dodge(width = 0.84), vjust = -0.2,
    size = 5 * 5 / 14
  ) +
  scale_fill_manual(values = method_colors, guide = "none") +
  scale_y_continuous(
    limits = c(0, 1.15), breaks = c(0, 0.25, 0.5, 0.75, 1),
    expand = expansion(mult = c(0, 0))
  ) +
  ylab("Power") +
  xlab("") +
  cross_scenario_theme()

overall_either_roc <- either_all_ROC_data_dataframe[
  , .(Power = mean(Power), FDR = mean(FDR)),
  by = .(source_scenario, Method, Cutoff)
]
overall_either_roc$Scenario <- factor(
  overall_either_roc$source_scenario,
  levels = scenarios,
  labels = unname(scenario_labels[scenarios])
)
overall_either_roc$Method <- factor(
  overall_either_roc$Method, levels = method_levels
)

p_roc_either_overall <- ggplot(
  overall_either_roc, aes(x = Power, y = 1 - FDR, color = Method)
) +
  geom_line(size = 0.5) +
  geom_point(
    data = overall_either_roc[abs(Cutoff - 0.5) < 1e-10],
    aes(shape = "PIP > 0.5"), size = 2.5,
    show.legend = c(color = FALSE)
  ) +
  scale_color_manual(values = method_colors) +
  scale_shape_manual(name = " ", values = c("PIP > 0.5" = 21)) +
  scale_x_continuous(
    limits = c(0, 1), breaks = c(0, 0.25, 0.5, 0.75, 1)
  ) +
  scale_y_continuous(
    limits = c(0, 1), breaks = c(0, 0.25, 0.5, 0.75, 1)
  ) +
  facet_wrap(~Scenario, ncol = 2) +
  xlab("Recall") +
  ylab("Precision") +
  cross_scenario_theme() +
  guides(color = guide_legend(order = 1), shape = guide_legend(order = 2))

overall_fdr_power <- rbindlist(
  list(
    FDR_Power_either[
      abs(Target_FDR - report_fdr) < 1e-10,
      .(Power = mean(Power)), by = .(source_scenario, Method)
    ][, Signal := "Either"],
    FDR_Power_shared[
      abs(Target_FDR - report_fdr) < 1e-10,
      .(Power = mean(Power)), by = .(source_scenario, Method)
    ][, Signal := "Shared"],
    FDR_Power_ancestry[
      abs(Target_FDR - report_fdr) < 1e-10,
      .(Power = mean(Power)), by = .(source_scenario, Method)
    ][, Signal := "Ancestry-specific"]
  )
)
overall_fdr_power[, Power := round(Power, 2)]
overall_fdr_power[, Scenario := factor(
  source_scenario, levels = scenarios,
  labels = unname(scenario_labels[scenarios])
)]
overall_fdr_power[, Method := factor(Method, levels = method_levels)]
overall_fdr_power[, Signal := factor(
  Signal, levels = c("Either", "Shared", "Ancestry-specific")
)]

p_fdr_power_overall <- ggplot(
  overall_fdr_power,
  aes(x = Signal, y = Power, fill = Method)
) +
  geom_bar(
    stat = "identity", position = position_dodge(width = 0.84),
    width = 0.72
  ) +
  geom_text(
    aes(label = Power, group = Method),
    position = position_dodge(width = 0.84), vjust = -0.2,
    size = 5 * 5 / 14
  ) +
  scale_fill_manual(values = method_colors) +
  scale_y_continuous(
    limits = c(0, 1.05), breaks = c(0, 0.25, 0.5, 0.75, 1),
    expand = expansion(mult = c(0, 0))
  ) +
  facet_wrap(~Scenario, ncol = 2) +
  ylab("Power") +
  xlab("FDR = 0.05") +
  cross_scenario_theme()

size_power_overall <- p_size_box_overall / p_power_bar_overall
roc_fdr_power_overall <- p_roc_either_overall / p_fdr_power_overall
averaged_summary <- (size_power_overall | roc_fdr_power_overall) +
  plot_annotation(tag_levels = "a") +
  plot_layout(guides = "collect") &
  theme(
    plot.tag = element_text(size = 7, face = "bold"),
    legend.position = "bottom"
  )

res_out <- paste0(
  "ancestry_specific_prior_across_scenarios_base_d_", base_d,
  if (nzchar(result_tag)) paste0("_", result_tag) else "",
  ".RData"
)
save(
  all_Set_data_dataframe,
  set_power_summary,
  FDR_Power_either,
  FDR_Power_shared,
  FDR_Power_ancestry,
  either_all_ROC_data_dataframe,
  performance_at_fdr_005,
  overall_set_power,
  overall_either_roc,
  overall_fdr_power,
  size_power,
  p_power_either,
  p_power_shared,
  p_power_ancestry,
  averaged_summary,
  result_tag,
  file = file.path(res_dir, res_out),
  compress = "xz"
)

fwrite(set_power_summary, file.path(res_dir, "credible_set_summary.tsv"), sep = "\t")
fwrite(
  performance_at_fdr_005,
  file.path(res_dir, "performance_at_fdr_0.05.tsv"), sep = "\t"
)
fwrite(overall_set_power, file.path(res_dir, "averaged_cs_power.tsv"), sep = "\t")
fwrite(overall_either_roc, file.path(res_dir, "averaged_pip_either_roc.tsv"), sep = "\t")
fwrite(overall_fdr_power, file.path(res_dir, "averaged_pip_power_fdr_0.05.tsv"), sep = "\t")

save_figure <- function(plot, stem, width, height) {
  ggsave(
    file.path(plot_dir, paste0(stem, ".pdf")), plot,
    height = height, width = width, units = "mm", dpi = 500
  )
  ggsave(
    file.path(plot_dir, paste0(stem, ".png")), plot,
    height = height, width = width, units = "mm", dpi = 500
  )
}

save_figure(size_power, "Sensitivity_Set_size_power", 320, 190)
save_figure(p_power_either, "Sensitivity_Power_Either_FDR_0.05", 320, 190)
save_figure(p_power_shared, "Sensitivity_Power_Shared_FDR_0.05", 320, 190)
save_figure(
  p_power_ancestry, "Sensitivity_Power_Ancestry_FDR_0.05", 340, 190
)
save_figure(averaged_summary, "Sensitivity_Averaged_Summary", 360, 240)

message("Saved cross-scenario d=", base_d, " summary: ", file.path(res_dir, res_out))
message("Saved cross-scenario figures under: ", plot_dir)
