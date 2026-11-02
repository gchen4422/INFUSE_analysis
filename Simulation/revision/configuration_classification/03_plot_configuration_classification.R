#!/usr/bin/env Rscript
.infuse_repo <- Sys.getenv("INFUSE_ANALYSIS_ROOT", "/path/to/INFUSE_analysis")

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
})

file_arg <- commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))][1L]
script_dir <- dirname(normalizePath(sub("^--file=", "", file_arg)))
source(file.path(.infuse_repo, "Simulation/revision/configuration_classification/00_common.R"))

analysis_file <- file.path(result_out_dir, "configuration_classification_analysis.rds")
if (!file.exists(analysis_file)) stop("Analysis output is missing.", call. = FALSE)
x <- readRDS(analysis_file)

blue <- "#0072B2"
orange <- "#D55E00"
green <- "#009E73"
purple <- "#CC79A7"
dim_cols <- c("d = 10" = blue, "d = 187" = orange)
class_cols <- c("EUR-specific" = blue, "AFR-specific" = orange, "Shared" = green)

theme_pub <- theme_bw(base_size = 10) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(linewidth = 0.25, colour = "grey90"),
    strip.background = element_rect(fill = "grey95", colour = "grey40"),
    strip.text = element_text(face = "bold"),
    legend.position = "bottom",
    plot.title = element_text(face = "bold", size = 11),
    plot.subtitle = element_text(size = 9),
    axis.title = element_text(face = "bold")
  )

# Match the manuscript-ready supplementary LD-mismatch figures.  This theme is
# used only for the reviewer configuration-classification figure; the remaining
# diagnostic figures below retain their existing presentation.
theme_supp <- function(base_size = 10) {
  theme_classic(base_size = base_size) +
    theme(
      axis.text = element_text(color = "#4D4D4D", size = 8),
      axis.text.x = element_text(margin = margin(t = 3)),
      axis.title = element_text(face = "bold", size = 10),
      strip.background = element_blank(),
      strip.text = element_text(face = "bold", size = 9),
      panel.spacing = unit(0.35, "lines"),
      legend.position = "bottom",
      legend.title = element_blank(),
      legend.text = element_text(size = 9),
      legend.key.width = unit(1.2, "lines"),
      plot.title = element_text(
        face = "bold", size = 12, hjust = 0, margin = margin(b = 2)
      ),
      plot.margin = margin(4, 4, 4, 4)
    )
}

inf_curve <- x$binary_curves[method == "INFUSE"]
auc_labels <- x$metric_summary[method == "INFUSE" & metric == "shared_auroc",
  .(dimension, label = sprintf("%s: AUROC %.3f", dimension, estimate))]
p_roc <- ggplot(inf_curve, aes(false_positive_rate, sensitivity, colour = dimension)) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "grey65") +
  geom_path(aes(linetype = dimension), linewidth = 0.9) +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE) +
  scale_colour_manual(values = dim_cols) +
  scale_linetype_manual(values = c("d = 10" = "solid", "d = 187" = "longdash")) +
  annotate("text", x = 0.98, y = c(0.12, 0.04), hjust = 1,
           label = auc_labels$label, colour = unname(dim_cols[auc_labels$dimension]), size = 3) +
  labs(
    title = "a",
    x = "False-positive rate", y = "Shared-effect sensitivity",
    colour = NULL, linetype = NULL
  ) + theme_supp()

cal <- x$calibration[method == "INFUSE" & n > 0]
p_cal <- ggplot(cal, aes(mean_predicted, observed_shared, colour = dimension)) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "grey65") +
  geom_line(aes(linetype = dimension), linewidth = 0.7) +
  geom_point(aes(size = n, shape = dimension), alpha = 0.9) +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE) +
  scale_colour_manual(values = dim_cols) +
  scale_linetype_manual(values = c("d = 10" = "solid", "d = 187" = "longdash")) +
  scale_shape_manual(values = c("d = 10" = 16, "d = 187" = 17)) +
  scale_size_continuous(range = c(2, 5)) +
  labs(
    title = "b",
    x = "Mean normalized shared support", y = "Observed fraction shared",
    colour = NULL, linetype = NULL, shape = NULL, size = "Causal SNPs"
  ) + theme_supp() + guides(size = "none", colour = "none",
                            linetype = "none", shape = "none")

make_confusion <- function(dd, panel) {
  d <- x$confusion[method == "INFUSE" & dimension == dd]
  d[, truth_class := factor(truth_class, levels = rev(config_levels))]
  d[, pred_class := factor(pred_class, levels = config_levels)]
  ggplot(d, aes(pred_class, truth_class, fill = proportion)) +
    geom_tile(colour = "white", linewidth = 0.8) +
    geom_text(aes(label = sprintf("%.1f%%\n(n=%d)", 100 * proportion, N)), size = 3) +
    scale_fill_gradient(low = "white", high = purple, limits = c(0, 1), labels = scales::percent) +
    scale_x_discrete(labels = c(
      "EUR-specific" = "WB-\nspecific",
      "AFR-specific" = "BB-\nspecific",
      "Shared" = "Shared"
    )) +
    scale_y_discrete(labels = c(
      "EUR-specific" = "WB-\nspecific",
      "AFR-specific" = "BB-\nspecific",
      "Shared" = "Shared"
    )) +
    labs(
      title = paste0(panel, "   ", dd),
      x = NULL, y = NULL, fill = "Row fraction"
    ) + theme_supp() + guides(fill = "none") +
    theme(
      axis.line = element_blank(),
      axis.ticks = element_blank()
    )
}
p_cm10 <- make_confusion("d = 10", "c")
p_cm187 <- make_confusion("d = 187", "d")

primary <- wrap_plots(
  p_roc,
  p_cm10,
  p_cal,
  p_cm187,
  ncol = 2,
  guides = "collect"
) & theme(legend.position = "bottom")

ggsave(file.path(figure_out_dir, "REVIEWER_CONFIGURATION_CLASSIFICATION_PRIMARY.pdf"),
       primary, width = 10.2, height = 8.0, units = "in", device = cairo_pdf)
ggsave(file.path(figure_out_dir, "REVIEWER_CONFIGURATION_CLASSIFICATION_PRIMARY.png"),
       primary, width = 10.2, height = 8.0, units = "in", dpi = 600,
       device = "png", bg = "white")

if (identical(Sys.getenv("CONFIG_PRIMARY_ONLY"), "1")) {
  message("Primary configuration-classification figure written to ", figure_out_dir)
  quit(save = "no", status = 0L)
}

# Class-specific precision and recall with locus-bootstrap confidence intervals.
metric_map <- data.table(
  metric = c(paste0("recall_", c("EUR", "AFR", "Shared")),
             paste0("precision_", c("EUR", "AFR", "Shared"))),
  measure = rep(c("Recall", "Precision"), each = 3),
  configuration = rep(config_levels, 2)
)
class_perf <- merge(x$metric_summary[method == "INFUSE"], metric_map, by = "metric")
p_class <- ggplot(class_perf,
                  aes(configuration, estimate, colour = dimension, group = dimension)) +
  geom_hline(yintercept = 1 / 3, linetype = 3, colour = "grey70") +
  geom_errorbar(aes(ymin = conf_low, ymax = conf_high),
                width = 0.12, linewidth = 0.6, position = position_dodge(width = 0.35)) +
  geom_point(size = 2.6, position = position_dodge(width = 0.35)) +
  facet_wrap(~measure) +
  scale_colour_manual(values = dim_cols) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  labs(
    title = "Configuration-specific classification performance",
    subtitle = "Points are pooled estimates; bars are 95% locus-cluster bootstrap intervals",
    x = NULL, y = "Proportion", colour = NULL
  ) + theme_pub + theme(axis.text.x = element_text(angle = 20, hjust = 1))
ggsave(file.path(figure_out_dir, "Configuration_class_precision_recall.pdf"),
       p_class, width = 8, height = 4.6, units = "in")
ggsave(file.path(figure_out_dir, "Configuration_class_precision_recall.png"),
       p_class, width = 8, height = 4.6, units = "in", dpi = 350)

# Confidence/coverage: confidence is relative configuration support, not causal PIP.
cc <- x$confidence_coverage[method == "INFUSE"]
p_cc <- ggplot(cc, aes(coverage, accuracy, colour = dimension)) +
  geom_path(linewidth = 0.9) + geom_point(size = 1.6) +
  scale_colour_manual(values = dim_cols) +
  scale_x_continuous(labels = scales::percent, limits = c(0, 1)) +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
  labs(
    title = "Accuracy-coverage tradeoff for configuration calls",
    subtitle = "Coverage decreases as the required maximum normalized configuration support increases",
    x = "Causal SNPs receiving a call", y = "Three-class accuracy", colour = NULL
  ) + theme_pub
ggsave(file.path(figure_out_dir, "Configuration_accuracy_coverage.pdf"),
       p_cc, width = 6.5, height = 4.7, units = "in")
ggsave(file.path(figure_out_dir, "Configuration_accuracy_coverage.png"),
       p_cc, width = 6.5, height = 4.7, units = "in", dpi = 350)

# Supporting end-to-end result: raw configuration PIP must localize and classify.
e2e <- x$end_to_end[method == "INFUSE" & cutoff <= 1]
e2e[, configuration := factor(configuration, levels = config_levels)]
p_e2e <- ggplot(e2e, aes(recall, precision, colour = dimension)) +
  geom_path(linewidth = 0.8) +
  facet_wrap(~configuration) +
  scale_colour_manual(values = dim_cols) +
  scale_x_continuous(limits = c(0, 1)) + scale_y_continuous(limits = c(0, 1)) +
  labs(
    title = "Supporting end-to-end localization and configuration",
    subtitle = "Raw configuration PIPs evaluated over all variants; counts pooled across simulation cells",
    x = "Recall", y = "Precision", colour = NULL
  ) + theme_pub
ggsave(file.path(figure_out_dir, "SUPPLEMENT_End_to_end_configuration_precision_recall.pdf"),
       p_e2e, width = 9.5, height = 4.3, units = "in")
ggsave(file.path(figure_out_dir, "SUPPLEMENT_End_to_end_configuration_precision_recall.png"),
       p_e2e, width = 9.5, height = 4.3, units = "in", dpi = 350)

# Context comparison with the unannotated model.
keep_metrics <- c("accuracy_3class", "balanced_accuracy_3class", "macro_f1",
                  "shared_auroc", "shared_average_precision")
context <- x$metric_summary[metric %in% keep_metrics]
context[, metric_label := factor(metric,
  levels = keep_metrics,
  labels = c("3-class accuracy", "3-class balanced accuracy", "Macro-F1",
             "Shared-vs-specific AUROC", "Shared-vs-specific average precision"))]
p_context <- ggplot(context, aes(metric_label, estimate, colour = method, shape = dimension)) +
  geom_errorbar(aes(ymin = conf_low, ymax = conf_high), width = 0.1,
                position = position_dodge(width = 0.45)) +
  geom_point(size = 2.6, position = position_dodge(width = 0.45)) +
  scale_colour_manual(values = c("INFUSE" = green, "MESuSiE" = "#666666")) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(
    title = "Context benchmark: INFUSE and unannotated MESuSiE",
    subtitle = "95% locus-cluster bootstrap intervals",
    x = NULL, y = "Metric value", colour = NULL, shape = NULL
  ) + theme_pub + theme(axis.text.x = element_text(angle = 25, hjust = 1))
ggsave(file.path(figure_out_dir, "SUPPLEMENT_INFUSE_vs_MESuSiE_configuration_metrics.pdf"),
       p_context, width = 9, height = 5.2, units = "in")
ggsave(file.path(figure_out_dir, "SUPPLEMENT_INFUSE_vs_MESuSiE_configuration_metrics.png"),
       p_context, width = 9, height = 5.2, units = "in", dpi = 350)

message("Figures written to ", figure_out_dir)
