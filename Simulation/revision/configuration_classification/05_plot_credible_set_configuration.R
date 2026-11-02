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

analysis_file <- file.path(
  result_out_dir,
  "credible_set_configuration_classification_analysis.rds"
)
if (!file.exists(analysis_file)) stop("Credible-set analysis output is missing.", call. = FALSE)
x <- readRDS(analysis_file)

blue <- "#0072B2"
orange <- "#D55E00"
purple <- "#CC79A7"
dim_cols <- c("d = 10" = blue, "d = 187" = orange)

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
auc_labels <- x$metric_summary[
  method == "INFUSE" & metric == "shared_auroc",
  .(dimension, label = sprintf("%s: AUROC %.3f", dimension, estimate))
]
setorder(auc_labels, dimension)

p_roc <- ggplot(
  inf_curve,
  aes(false_positive_rate, sensitivity, colour = dimension)
) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "grey65") +
  geom_path(aes(linetype = dimension), linewidth = 0.9) +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE) +
  scale_colour_manual(values = dim_cols) +
  scale_linetype_manual(values = c("d = 10" = "solid", "d = 187" = "longdash")) +
  annotate(
    "text", x = 0.98, y = c(0.12, 0.04), hjust = 1,
    label = auc_labels$label,
    colour = unname(dim_cols[auc_labels$dimension]),
    size = 3
  ) +
  labs(
    title = "a",
    x = "False-positive rate",
    y = "Shared-effect sensitivity",
    colour = NULL,
    linetype = NULL
  ) +
  theme_supp()

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
    x = "Mean credible-set shared probability",
    y = "Observed fraction shared",
    colour = NULL,
    linetype = NULL,
    shape = NULL,
    size = "Credible sets"
  ) +
  theme_supp() +
  guides(size = "none", colour = "none", linetype = "none", shape = "none")

make_confusion <- function(dd, panel) {
  d <- x$confusion[method == "INFUSE" & dimension == dd]
  d[, truth_class := factor(truth_class, levels = rev(config_levels))]
  d[, pred_class := factor(pred_class, levels = config_levels)]

  ggplot(d, aes(pred_class, truth_class, fill = proportion)) +
    geom_tile(colour = "white", linewidth = 0.8) +
    geom_text(aes(label = sprintf("%.1f%%\n(n=%d)", 100 * proportion, N)), size = 3) +
    scale_fill_gradient(
      low = "white", high = purple, limits = c(0, 1), labels = scales::percent
    ) +
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
      x = NULL,
      y = NULL,
      fill = "Row fraction"
    ) +
    theme_supp() +
    guides(fill = "none") +
    theme(axis.line = element_blank(), axis.ticks = element_blank())
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

output_stem <- file.path(
  figure_out_dir,
  "REVIEWER_CREDIBLE_SET_CONFIGURATION_CLASSIFICATION_PRIMARY"
)

ggsave(
  paste0(output_stem, ".pdf"),
  primary,
  width = 10.2,
  height = 8.0,
  units = "in",
  device = cairo_pdf
)
ggsave(
  paste0(output_stem, ".png"),
  primary,
  width = 10.2,
  height = 8.0,
  units = "in",
  dpi = 600,
  device = "png",
  bg = "white"
)

message("Credible-set configuration-classification figure written to ", figure_out_dir)
