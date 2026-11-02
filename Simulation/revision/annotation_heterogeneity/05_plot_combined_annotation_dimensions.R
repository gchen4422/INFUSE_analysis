#!/usr/bin/env Rscript
.infuse_repo <- Sys.getenv("INFUSE_ANALYSIS_ROOT", "/path/to/INFUSE_analysis")
.infuse_project <- Sys.getenv("INFUSE_PROJECT", "/path/to/baseline_simulation")
.infuse_work <- Sys.getenv("INFUSE_REVISION", "/path/to/revision_output")

# Combine the corrected d = 10 cross-scenario results, the d = 187 INFUSE
# result, and the three PAINTOR benchmarks. The d = 10 analysis remains
# labelled "INFUSE", the d = 187 analysis is labelled "INFUSE-H", and the
# annotation-informed Paintor fits are labelled "+" or "+H" for their 10- or
# 187-column baseline panels. Both fits additionally include EAF_EUR and EAF_AFR,
# giving actual fitted dimensions of 12 and 189.

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(ggpattern)
})

project_dir <- .infuse_project
revision_dir <- .infuse_work
analysis_dir <- file.path(
  revision_dir, "ancestry_specific_eaf_prior_simulation"
)
source(file.path(.infuse_repo, "Simulation/revision/annotation_heterogeneity/utility.R"))

summary_root <- file.path(
  revision_dir, "res", "ancestry_specific_prior_sensitivity",
  "across_scenarios"
)
figure_dir <- file.path(
  revision_dir, "figures", "ancestry_specific_prior_sensitivity",
  "across_scenarios", "combined_d10_d187"
)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

summary_file <- function(base_d) {
  # The corrected low-dimensional analysis follows the current INFUSE code
  # path and directly optimizes the 12 annotation coefficients.  Its summary
  # is intentionally kept separate from the obsolete sparse-PCA/BIC rerun.
  suffix <- if (identical(as.integer(base_d), 10L)) {
    "_current_github"
  } else {
    ""
  }

  file.path(
    summary_root, paste0("base_d_", base_d),
    paste0(
      "ancestry_specific_prior_across_scenarios_base_d_",
      base_d, suffix, ".RData"
    )
  )
}

load_summary <- function(base_d) {
  path <- summary_file(base_d)
  if (!file.exists(path)) stop("Missing cross-scenario summary: ", path)
  message("Loading cross-scenario summary: ", path)
  env <- new.env(parent = emptyenv())
  load(path, envir = env)
  required <- c(
    "all_Set_data_dataframe", "set_power_summary",
    "FDR_Power_either", "FDR_Power_shared", "FDR_Power_ancestry",
    "either_all_ROC_data_dataframe"
  )
  missing_objects <- setdiff(required, ls(env))
  if (length(missing_objects)) {
    stop(
      "Missing objects in ", path, ": ",
      paste(missing_objects, collapse = ", ")
    )
  }
  observed_base_d <- unique(as.integer(env$set_power_summary$base_d))
  if (!identical(observed_base_d, as.integer(base_d))) {
    stop(
      "Unexpected base_d value in ", path, ": ",
      paste(observed_base_d, collapse = ", ")
    )
  }
  if (identical(as.integer(base_d), 10L)) {
    if (!exists("result_tag", envir = env, inherits = FALSE) ||
        !identical(env$result_tag, "current_github")) {
      stop(
        "The d = 10 combined figure requires the corrected current_github ",
        "summary; unexpected result tag in ", path
      )
    }
  }
  env
}

low <- load_summary(10L)
high <- load_summary(187L)

load_paintor_summary <- function() {
  path <- file.path(
    summary_root, "paintor_benchmarks",
    "paintor_benchmarks_across_scenarios.RData"
  )
  if (!file.exists(path)) stop("Missing PAINTOR benchmark summary: ", path)
  message("Loading PAINTOR benchmark summary: ", path)
  env <- new.env(parent = emptyenv())
  load(path, envir = env)
  required <- c(
    "all_Set_data_dataframe", "set_power_summary",
    "FDR_Power_either", "FDR_Power_shared", "FDR_Power_ancestry",
    "either_all_ROC_data_dataframe"
  )
  missing_objects <- setdiff(required, ls(env))
  if (length(missing_objects)) {
    stop(
      "Missing objects in ", path, ": ",
      paste(missing_objects, collapse = ", ")
    )
  }
  expected_methods <- c("Paintor", "Paintor+", "Paintor+H")
  observed_methods <- sort(unique(env$set_power_summary$Method))
  if (!identical(observed_methods, sort(expected_methods))) {
    stop("Unexpected PAINTOR methods in ", path)
  }
  cells <- as.data.table(env$all_Set_data_dataframe)[,
    .(scenario, Method, locus, h2_index, causal_index)]
  if (nrow(cells) != 7200L || anyDuplicated(cells)) {
    stop("Expected 7,200 PAINTOR method-by-cell results in ", path)
  }
  env
}

paintor <- load_paintor_summary()

scenarios <- c(
  "baseline", "different_annotation", "different_weight", "combined"
)
scenario_labels <- c(
  baseline = "Baseline",
  different_annotation = "Different\nannotation",
  different_weight = "Different\nweight",
  combined = "Both"
)
method_levels <- c(
  "SuSiE", "Paintor", "Paintor+", "Paintor+H",
  "MESuSiE", "INFUSE", "INFUSE-H"
)
method_colors <- c(
  "SuSiE" = "#fbe79e",
  "Paintor" = "#f2b56e",
  "Paintor+" = "#f57c6e",
  "Paintor+H" = "#c94f6d",
  "MESuSiE" = "#84c3b7",
  "INFUSE" = "#71b8ed",
  "INFUSE-H" = "#4a90e2"
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

# Keep every method from d = 10, and add only high-dimensional INFUSE.
# SuSiE and MESuSiE are intentionally not duplicated because they do not use
# the annotation panel and their results are unchanged between d = 10 and 187.
combine_with_infuse_h_and_paintor <- function(
    low_object, high_object, paintor_object) {
  low_dt <- copy(as.data.table(low_object))
  high_dt <- copy(as.data.table(high_object))[Method == "INFUSE"]
  high_dt[, Method := "INFUSE-H"]
  paintor_dt <- copy(as.data.table(paintor_object))
  if (!"source_scenario" %in% names(paintor_dt)) {
    if (!"scenario" %in% names(paintor_dt)) {
      stop("PAINTOR object lacks scenario information.")
    }
    paintor_dt[, source_scenario := scenario]
  }
  rbindlist(
    list(low_dt, paintor_dt, high_dt), use.names = TRUE, fill = TRUE
  )
}

all_set <- combine_with_infuse_h_and_paintor(
  low$all_Set_data_dataframe, high$all_Set_data_dataframe,
  paintor$all_Set_data_dataframe
)
set_power <- combine_with_infuse_h_and_paintor(
  low$set_power_summary, high$set_power_summary,
  paintor$set_power_summary
)
fdr_either <- combine_with_infuse_h_and_paintor(
  low$FDR_Power_either, high$FDR_Power_either,
  paintor$FDR_Power_either
)
fdr_shared <- combine_with_infuse_h_and_paintor(
  low$FDR_Power_shared, high$FDR_Power_shared,
  paintor$FDR_Power_shared
)
fdr_ancestry <- combine_with_infuse_h_and_paintor(
  low$FDR_Power_ancestry, high$FDR_Power_ancestry,
  paintor$FDR_Power_ancestry
)
either_roc <- combine_with_infuse_h_and_paintor(
  low$either_all_ROC_data_dataframe,
  high$either_all_ROC_data_dataframe,
  paintor$either_all_ROC_data_dataframe
)

add_plot_factors <- function(x) {
  x <- as.data.frame(copy(as.data.table(x)))
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
# 95% credible-set size and power: existing panel layout
################################################################

set_plot_data <- add_plot_factors(all_set)
set_plot_data$log2_Size <- log2(set_plot_data$Size + 1)
power_plot_data <- add_plot_factors(set_power)

p_size_box <- ggplot(
  set_plot_data,
  aes(x = Scenario, y = log2_Size, fill = Method)
) +
  geom_boxplot(
    position = position_dodge2(width = 0.85, preserve = "single"),
    width = 0.72, outlier.size = 0.1, fatten = 0.5,
    color = "darkgray"
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
  scale_fill_manual(
    values = method_colors,
    guide = guide_legend(
      title = "Method", order = 1, nrow = 1, byrow = TRUE
    )
  ) +
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
    legend.position = "bottom",
    legend.box = "horizontal"
  )

################################################################
# PIP power at FDR = 0.05: existing panel layouts
################################################################

make_signal_power_plot <- function(x) {
  plot_data <- add_plot_factors(x)
  plot_data <- plot_data[
    abs(plot_data$Target_FDR - report_fdr) < 1e-10,
  ]
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
    scale_fill_manual(
      values = method_colors,
      guide = guide_legend(
        title = "Method", order = 1, nrow = 1, byrow = TRUE
      )
    ) +
    ylab("Power") +
    xlab("FDR = 0.05") +
    scale_y_continuous(
      limits = c(0, 1.15), breaks = c(0, 0.25, 0.5, 0.75, 1),
      expand = expansion(mult = c(0, 0))
    ) +
    facet_grid(h2 ~ causal_num, labeller = label_parsed) +
    cross_scenario_theme()
}

p_power_either <- make_signal_power_plot(fdr_either)
p_power_shared <- make_signal_power_plot(fdr_shared)

ancestry_plot_data <- add_plot_factors(fdr_ancestry)
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
  scale_fill_manual(
    name = "Method", values = method_colors,
    guide = guide_legend(order = 2, nrow = 1, byrow = TRUE)
  ) +
  scale_pattern_manual(
    name = "Ancestry",
    values = c("EUR" = "stripe", "AFR" = "circle"),
    breaks = c("EUR", "AFR"),
    labels = c("White British", "Black British")
  ) +
  ylab("Power") +
  xlab("FDR = 0.05") +
  scale_y_continuous(
    limits = c(0, 1.15), breaks = c(0, 0.25, 0.5, 0.75, 1),
    expand = expansion(mult = c(0, 0))
  ) +
  facet_grid(h2 ~ causal_num, labeller = label_parsed) +
  cross_scenario_theme() +
  guides(
    pattern = guide_legend(
      order = 1, nrow = 1, byrow = TRUE,
      override.aes = list(fill = "white")
    )
  ) +
  theme(
    legend.box = "horizontal",
    legend.box.just = "center"
  )

################################################################
# Averaged summary: existing four-panel layout
################################################################

overall_set_plot_data <- add_plot_factors(all_set)
overall_set_plot_data$log2_Size <- log2(overall_set_plot_data$Size + 1)

overall_set_power <- set_power[
  , .(Power = mean(Power_name)),
  by = .(source_scenario, Method)
]
overall_set_power_plot <- as.data.frame(copy(overall_set_power))
overall_set_power_plot$Scenario <- factor(
  overall_set_power_plot$source_scenario,
  levels = scenarios, labels = unname(scenario_labels[scenarios])
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
    width = 0.72, outlier.size = 0.1, fatten = 0.5,
    color = "darkgray"
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

overall_either_roc <- either_roc[
  , .(Power = mean(Power), FDR = mean(FDR)),
  by = .(source_scenario, Method, Cutoff)
]
overall_either_roc$Scenario <- factor(
  overall_either_roc$source_scenario,
  levels = scenarios, labels = unname(scenario_labels[scenarios])
)
overall_either_roc$Method <- factor(
  overall_either_roc$Method, levels = method_levels
)

p_roc_either_overall <- ggplot(
  overall_either_roc,
  aes(x = Power, y = 1 - FDR, color = Method)
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
  guides(
    color = "none",
    shape = guide_legend(title = NULL, order = 1,
      override.aes = list(colour = "black", fill = "white", shape = 21,
                          size = 2.5, linetype = 0))
  )

overall_fdr_power <- rbindlist(
  list(
    fdr_either[
      abs(Target_FDR - report_fdr) < 1e-10,
      .(Power = mean(Power)), by = .(source_scenario, Method)
    ][, Signal := "Either"],
    fdr_shared[
      abs(Target_FDR - report_fdr) < 1e-10,
      .(Power = mean(Power)), by = .(source_scenario, Method)
    ][, Signal := "Shared"],
    fdr_ancestry[
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
p_fdr_power_overall <- p_fdr_power_overall + guides(
  fill = guide_legend(title = "Method", order = 2, nrow = 1, byrow = TRUE)
)
roc_fdr_power_overall <- p_roc_either_overall / p_fdr_power_overall
averaged_summary <- (size_power_overall | roc_fdr_power_overall) +
  plot_annotation(tag_levels = "a") +
  plot_layout(guides = "collect") &
  theme(
    plot.tag = element_text(size = 7, face = "bold"),
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.key.width = grid::unit(4, "mm"),
    legend.key.height = grid::unit(4, "mm")
  )

save_figure <- function(plot, stem, width, height, dpi = 500) {
  ggsave(
    file.path(figure_dir, paste0(stem, ".pdf")), plot,
    height = height, width = width, units = "mm", dpi = dpi,
    device = cairo_pdf
  )
  ggsave(
    file.path(figure_dir, paste0(stem, ".png")), plot,
    height = height, width = width, units = "mm", dpi = dpi
  )
}

save_figure(
  size_power, "Sensitivity_Set_size_power_INFUSE_H", 440, 210
)
save_figure(
  p_power_either, "Sensitivity_Power_Either_FDR_0.05_INFUSE_H", 440, 200
)
save_figure(
  p_power_shared, "Sensitivity_Power_Shared_FDR_0.05_INFUSE_H", 440, 200
)
save_figure(
  p_power_ancestry,
  "Sensitivity_Power_Ancestry_FDR_0.05_INFUSE_H", 460, 210
)
save_figure(
  averaged_summary, "Sensitivity_Averaged_Summary_INFUSE_H",
  170, 145, dpi = 600
)

# Save the numerical inputs used by the combined benchmark figures.
combined_res_dir <- file.path(
  summary_root, "combined_d10_d187_paintor"
)
dir.create(combined_res_dir, recursive = TRUE, showWarnings = FALSE)


save(
  all_set, set_power, fdr_either, fdr_shared, fdr_ancestry,
  either_roc, overall_set_power, overall_either_roc, overall_fdr_power,
  file = file.path(
    combined_res_dir, "combined_d10_d187_paintor_figure_data.RData"
  ),
  compress = "xz"
)
fwrite(
  set_power,
  file.path(combined_res_dir, "credible_set_summary.tsv"), sep = "\t"
)
fwrite(
  rbindlist(
    list(
      copy(fdr_either)[, Signal := "Either"],
      copy(fdr_shared)[, Signal := "Shared"],
      copy(fdr_ancestry)[, Signal := "Ancestry-specific"]
    ),
    fill = TRUE
  )[abs(Target_FDR - report_fdr) < 1e-10],
  file.path(combined_res_dir, "performance_at_fdr_0.05.tsv"), sep = "\t"
)

message(
  "Saved combined d=10, d=187, and PAINTOR benchmark figures under: ",
  figure_dir
)
