#!/usr/bin/env Rscript
.infuse_work <- Sys.getenv("INFUSE_REVISION", "/path/to/revision_output")

# Create the primary supplementary figure for the raw Evo 2 external-
# validation analysis. This script reads only the external-validation summary
# tables and deliberately excludes the signed_evo2_finemapping analysis.

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
})

this_file <- normalizePath(sub(
  "^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]
), mustWork = FALSE)
script_dir <- dirname(this_file)
module_root <- file.path(.infuse_work, "evo2_raw_score_external_validation")
result_dir <- file.path(module_root, "results")
figure_dir <- file.path(
  dirname(module_root), "figures", "evo2_raw_score_external_validation"
)

threshold_files <- file.path(
  result_dir,
  c(
    "raw_evo2_enrichment_by_pip_threshold_bp.tsv",
    "raw_evo2_enrichment_by_pip_threshold_bmi.tsv"
  )
)
top_n_files <- file.path(
  result_dir,
  c(
    "raw_evo2_enrichment_by_top_n_bp.tsv",
    "raw_evo2_enrichment_by_top_n_bmi.tsv"
  )
)

stopifnot(all(file.exists(threshold_files)), all(file.exists(top_n_files)))

threshold_results <- rbindlist(lapply(threshold_files, fread), fill = TRUE)
top_n_results <- rbindlist(lapply(top_n_files, fread), fill = TRUE)

method_names <- c(
  MESuSiE = "MESuSiE",
  Classic_Annotations = "INFUSE"
)
method_colors <- c(
  MESuSiE = "#84c3b7",
  INFUSE = "#71b8ed"
)
method_shapes <- c(MESuSiE = 16, INFUSE = 17)
method_linetypes <- c(MESuSiE = "solid", INFUSE = "22")
trait_order <- c("BMI", "BP")

prepare_results <- function(x) {
  x <- copy(x)
  x[, `:=`(
    Trait = factor(Trait, levels = trait_order),
    Method = factor(
      unname(method_names[method]),
      levels = c("MESuSiE", "INFUSE")
    )
  )]
  x
}

threshold_results <- prepare_results(threshold_results)
top_n_results <- prepare_results(top_n_results)
top_n_results[, cutoff_factor := factor(
  cutoff,
  levels = c(10, 20, 50, 100, 200),
  labels = c("10", "20", "50", "100", "200")
)]

theme_supplement <- function(base_size = 11) {
  theme_classic(base_size = base_size) +
    theme(
      axis.text = element_text(color = "#4D4D4D", size = 9),
      axis.text.x = element_text(margin = margin(t = 3)),
      axis.title = element_text(face = "bold", size = 11),
      strip.background = element_blank(),
      strip.text = element_text(face = "bold", size = 11),
      panel.spacing = unit(0.8, "lines"),
      legend.position = "bottom",
      legend.title = element_blank(),
      legend.text = element_text(size = 10),
      legend.key.width = unit(1.4, "lines"),
      plot.title = element_text(
        face = "bold", size = 14, hjust = 0,
        margin = margin(b = 3)
      ),
      plot.margin = margin(8, 8, 8, 8)
    )
}

make_shift_plot <- function(data, x_var, panel_label, x_label,
                            x_is_threshold = FALSE) {
  x_symbol <- as.name(x_var)
  p <- ggplot(
    data,
    aes(
      x = !!x_symbol,
      y = deleterious_shift,
      color = Method,
      shape = Method,
      linetype = Method,
      group = Method
    )
  ) +
    geom_hline(
      yintercept = 0, color = "#777777", linetype = "dashed",
      linewidth = 0.45
    ) +
    geom_line(linewidth = 0.65) +
    geom_point(size = 2.0) +
    facet_wrap(~Trait, nrow = 1, scales = "free_y") +
    scale_color_manual(values = method_colors, drop = FALSE) +
    scale_shape_manual(values = method_shapes, drop = FALSE) +
    scale_linetype_manual(values = method_linetypes, drop = FALSE) +
    scale_y_continuous(labels = scales::label_scientific(digits = 2)) +
    labs(
      title = panel_label,
      x = x_label,
      y = "Deleterious Evo 2 shift\n(background mean - prioritized mean)"
    ) +
    theme_supplement()

  if (x_is_threshold) {
    p <- p + scale_x_continuous(
      breaks = seq(0.5, 0.9, 0.1),
      limits = c(0.49, 0.91)
    )
  }
  p
}

p_bmi_threshold <- make_shift_plot(
  threshold_results[Trait == "BMI"],
  "cutoff", "a", "PIP threshold",
  x_is_threshold = TRUE
)
p_bp_threshold <- make_shift_plot(
  threshold_results[Trait == "BP"],
  "cutoff", "b", "PIP threshold",
  x_is_threshold = TRUE
) + theme(axis.title.y = element_blank())
p_bmi_top_n <- make_shift_plot(
  top_n_results[Trait == "BMI"],
  "cutoff_factor", "c", "Number of top-ranked SNPs"
)
p_bp_top_n <- make_shift_plot(
  top_n_results[Trait == "BP"],
  "cutoff_factor", "d", "Number of top-ranked SNPs"
) + theme(axis.title.y = element_blank())

primary_figure <- (
  p_bmi_threshold | p_bp_threshold
) / (
  p_bmi_top_n | p_bp_top_n
) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
output_stem <- file.path(
  figure_dir, "Evo2_external_validation_primary_supp_figure"
)

ggsave(
  paste0(output_stem, ".pdf"), primary_figure,
  width = 10.2, height = 8.0, units = "in", device = cairo_pdf
)
ggsave(
  paste0(output_stem, ".png"), primary_figure,
  width = 10.2, height = 8.0, units = "in", dpi = 600,
  device = "png", bg = "white"
)

message("Created: ", output_stem, ".pdf")
message("Created: ", output_stem, ".png")
