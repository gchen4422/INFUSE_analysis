#!/usr/bin/env Rscript
.infuse_repo <- Sys.getenv("INFUSE_ANALYSIS_ROOT", "/path/to/INFUSE_analysis")

# Extend the established combined low/high-dimensional Figure 2 with SuSiEx and
# MultiSuSiE. Panels a-e are constructed from the current canonical d=10 and
# d=187 summary results. This is important because the saved low-dimensional
# plot cache predates the corrected w2_startover summary and is stale. Panel f
# retains the established enrichment-MSE data and appends the two new methods.

this_file <- normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1L]),
                           mustWork = TRUE)
source(file.path(.infuse_repo, "Simulation/revision/susiex_multisusie_baseline/00_common.R"))
suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
})

required_files <- c(
  file.path(RESULT_DIR, "combined_credible_set_cell_metrics.tsv.gz"),
  file.path(RESULT_DIR, "combined_mean_threshold_curves.tsv.gz"),
  file.path(PROJECT_ROOT, "plot", "MSE_low_dim_w2.Rdata"),
  file.path(PROJECT_ROOT, "More_annotes", "plot", "MSE_high_dim_w2.Rdata"),
  file.path(RESULT_DIR, "new_method_enrichment_MSE_PIP_gt_0.5.tsv")
)
missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files)) {
  stop("Missing required figure input(s): ", paste(missing_files, collapse = ", "),
       call. = FALSE)
}

load_objects <- function(path) {
  e <- new.env(parent = emptyenv())
  load(path, envir = e)
  e
}
combined_sets <- fread(required_files[[1L]])
combined_curves <- fread(required_files[[2L]])
low_mse <- load_objects(required_files[[3L]])
high_mse <- load_objects(required_files[[4L]])
if (!"p_mse_box" %in% ls(low_mse) || !"p_mse_box" %in% ls(high_mse)) {
  stop("A reference enrichment-MSE RData file lacks p_mse_box", call. = FALSE)
}

relabel_reference_method <- function(method, high_dimension = FALSE) {
  method <- as.character(method)
  out <- method
  if (high_dimension) {
    out[method == "Paintor+annot"] <- "Paintor+H"
    out[method == "MAFPC"] <- "INFUSE-H"
  } else {
    out[method == "Paintor+annot"] <- "Paintor+"
    out[method == "MAFPC"] <- "INFUSE"
  }
  out
}

method_order <- c("SuSiE", "SuSiEx", "MultiSuSiE", "Paintor", "Paintor+",
                  "Paintor+H", "MESuSiE", "INFUSE", "INFUSE-H")
method_colors <- c(
  "SuSiE" = "#FBE79E", "SuSiEx" = "#E69F00", "MultiSuSiE" = "#CC79A7",
  "Paintor" = "#F2B56E", "Paintor+" = "#F57C6E", "Paintor+H" = "#D63447",
  "MESuSiE" = "#84C3B7", "INFUSE" = "#71B8ED", "INFUSE-H" = "#4A90E2"
)
high_keep <- c("Paintor+annot", "MAFPC")

new_mse_detail <- fread(required_files[[5L]])

# Match the established combined Figure 2: all d=10 methods, plus only the two
# annotation-aware d=187 methods. Annotation-free SuSiEx and MultiSuSiE are
# shown once because their input data and results are identical across d.
d10_keep <- c("SuSiE", "SuSiEx", "MultiSuSiE", "PAINTOR", "PAINTOR+",
              "MESuSiE", "INFUSE")
d187_keep <- c("PAINTOR+H", "INFUSE-H")
figure_sets <- rbindlist(list(
  combined_sets[dimension == "d = 10" & method %in% d10_keep],
  combined_sets[dimension == "d = 187" & method %in% d187_keep]
), use.names = TRUE, fill = TRUE)
figure_curves <- rbindlist(list(
  combined_curves[dimension == "d = 10" & method %in% d10_keep],
  combined_curves[dimension == "d = 187" & method %in% d187_keep]
), use.names = TRUE, fill = TRUE)
# The canonical combined tables use the internal label PAINTOR. Use the
# manuscript-facing capitalization without changing any upstream result.
figure_sets[method == "PAINTOR", method := "Paintor"]
figure_curves[method == "PAINTOR", method := "Paintor"]
figure_sets[method == "PAINTOR+", method := "Paintor+"]
figure_curves[method == "PAINTOR+", method := "Paintor+"]
figure_sets[method == "PAINTOR+H", method := "Paintor+H"]
figure_curves[method == "PAINTOR+H", method := "Paintor+H"]
if (any(figure_sets[, .N, by = method]$N != N_CELLS) ||
    !setequal(figure_sets$method, method_order)) {
  stop("Canonical credible-set data do not contain 600 rows for every method",
       call. = FALSE)
}

# a. Credible-set size: use the current canonical cell-level summaries.
size_data <- figure_sets[, .(
  Power = set_power, FDR = set_fdr, Size = set_size, Method = method,
  h2 = h2_index, causal_num = causal_index,
  source = fifelse(method %in% c("SuSiEx", "MultiSuSiE"), "new run",
                   fifelse(dimension == "d = 10", "canonical d=10", "canonical d=187"))
)]

# b. Credible-set power.
power_data <- size_data[, .(Power_name = mean(Power), source = first(source)), by = Method]

# c/d. Match utility.R::FDR_Power for the new methods.
canonical_fdr_bars <- function(target) {
  selected <- figure_curves[, {
    delta <- abs(FDR - target)
    candidate <- .SD[delta == min(delta)]
    candidate[which.max(Power)][1L]
  }, by = .(method, causal_index, h2_index, endpoint, ancestry)]
  selected[, Power := round(Power, 2)]
  out <- selected[, .(Power = round(mean(Power), 2)), by = .(Method = method, endpoint)]
  out[, Cat := fifelse(endpoint == "Either", "Either",
                       fifelse(endpoint == "Shared", "Shared", "Ancestry"))]
  out[, source := fifelse(Method %in% c("SuSiEx", "MultiSuSiE"), "new run",
                          fifelse(Method %in% c("Paintor+H", "INFUSE-H"),
                                  "canonical d=187", "canonical d=10"))]
  out[, FDR := target]
  out[, .(FDR, Method, Power, Cat, source)]
}
fdr05_data <- canonical_fdr_bars(0.05)
fdr10_data <- canonical_fdr_bars(0.10)

# e. Overall effect-in-either-ancestry precision-recall curve.
pr_data <- figure_curves[endpoint == "Either", .(
  Power = mean(Power), FDR = mean(FDR)
), by = .(Method = method, Cutoff)]
pr_data[, source := fifelse(Method %in% c("SuSiEx", "MultiSuSiE"), "new run",
                            fifelse(Method %in% c("Paintor+H", "INFUSE-H"),
                                    "canonical d=187", "canonical d=10"))]

# f. Existing enrichment-MSE values are copied; new values use the compatible
# definition from 05b. The complement-background sensitivity remains separate.
mse_low <- as.data.table(low_mse$p_mse_box$data)
mse_low[, `:=`(method = relabel_reference_method(method), source = "reference d=10")]
mse_high <- as.data.table(high_mse$p_mse_box$data)[as.character(method) %in% high_keep]
mse_high[, `:=`(method = relabel_reference_method(method, TRUE), source = "reference d=187")]
mse_new <- new_mse_detail[, .(method, causal, h2, annotation,
                              mse = mse_figure2, source = "new run")]
mse_data <- rbindlist(list(mse_low, mse_high, mse_new), use.names = TRUE, fill = TRUE)

figure_method_tables <- list(size_data, power_data, fdr05_data, fdr10_data, pr_data)
if (any(vapply(figure_method_tables, function(d) {
  any(!as.character(d$Method) %in% method_order)
}, logical(1)))) {
  stop("Unexpected figure method", call. = FALSE)
}
size_data[, Method := factor(as.character(Method), levels = method_order)]
power_data[, Method := factor(as.character(Method), levels = method_order)]
fdr05_data[, Method := factor(as.character(Method), levels = method_order)]
fdr10_data[, Method := factor(as.character(Method), levels = method_order)]
pr_data[, Method := factor(as.character(Method), levels = method_order)]
mse_data[, method := factor(as.character(method), levels = method_order)]
mse_data[, annotation := factor(as.character(annotation),
                                levels = c("Non-synonymous", "Conserved", "DHS",
                                           "H3K27ac", "H3K4me3"))]
fdr05_data[, Cat := factor(as.character(Cat), levels = c("Either", "Shared", "Ancestry"))]
fdr10_data[, Cat := factor(as.character(Cat), levels = c("Either", "Shared", "Ancestry"))]

fwrite(size_data, file.path(RESULT_DIR, "figure2_extended_credible_set_data.tsv.gz"), sep = "\t")
fwrite(power_data, file.path(RESULT_DIR, "figure2_extended_credible_set_power.tsv"), sep = "\t")
fwrite(fdr05_data, file.path(RESULT_DIR, "figure2_extended_power_FDR_0.05.tsv"), sep = "\t")
fwrite(fdr10_data, file.path(RESULT_DIR, "figure2_extended_power_FDR_0.10.tsv"), sep = "\t")
fwrite(pr_data, file.path(RESULT_DIR, "figure2_extended_precision_recall.tsv.gz"), sep = "\t")
fwrite(mse_data, file.path(RESULT_DIR, "figure2_extended_enrichment_MSE.tsv"), sep = "\t")

theme_figure2 <- theme_bw(base_size = 9) +
  theme(axis.text.x = element_text(size = 7, angle = 34, hjust = 1),
        axis.text.y = element_text(size = 7.5),
        axis.title = element_text(size = 10, face = "bold"),
        strip.text = element_text(size = 9, face = "bold"),
        strip.background = element_blank(), panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(colour = "black", linewidth = 0.35),
        plot.title = element_text(size = 12, hjust = 0.5),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 9, face = "bold"),
        legend.position = "bottom")

p_size <- ggplot(size_data, aes(Method, log2(as.numeric(Size) + 1), fill = Method)) +
  geom_boxplot(outlier.size = 0.12, linewidth = 0.32, colour = "darkgray") +
  scale_fill_manual(values = method_colors, drop = FALSE) +
  labs(x = NULL, y = expression(log[2]("95% credible-set size + 1")),
       title = "Credible-set size") +
  theme_figure2 + theme(legend.position = "none")

p_power <- ggplot(power_data, aes(Method, Power_name, fill = Method)) +
  geom_col(width = 0.78) +
  geom_text(aes(label = sprintf("%.2f", Power_name)), vjust = -0.45, size = 2.15) +
  scale_fill_manual(values = method_colors, drop = FALSE) +
  scale_y_continuous(limits = c(0, 1), expand = expansion(mult = c(0, 0.04))) +
  labs(x = NULL, y = "Power", title = "Credible-set power") +
  theme_figure2 + theme(legend.position = "none")

plot_fdr_panel <- function(d, target) {
  ggplot(d, aes(Cat, Power, fill = Method)) +
    geom_col(position = position_dodge(width = 0.9), width = 0.82) +
    geom_text(aes(label = sprintf("%.2f", Power), group = Method),
              position = position_dodge(width = 0.9), vjust = -0.4, size = 1.75) +
    scale_fill_manual(values = method_colors, drop = FALSE) +
    scale_y_continuous(limits = c(0, 0.65), expand = expansion(mult = c(0, 0.03))) +
    labs(x = NULL, y = "Power", title = sprintf("Power at FDR = %.2g", target)) +
    theme_figure2 +
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5), legend.position = "none")
}
p_fdr05 <- plot_fdr_panel(fdr05_data, 0.05)
p_fdr10 <- plot_fdr_panel(fdr10_data, 0.10)

p_pr <- ggplot(pr_data, aes(Power, 1 - FDR, colour = Method)) +
  geom_line(linewidth = 0.62, na.rm = TRUE) +
  geom_point(data = pr_data[abs(Cutoff - 0.5) < 1e-10],
             aes(shape = "PIP > 0.5"), size = 2.0, stroke = 0.45) +
  scale_colour_manual(values = method_colors, drop = FALSE) +
  scale_shape_manual(name = NULL, values = c("PIP > 0.5" = 21)) +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
  labs(x = "Recall", y = "Precision", title = "Effect in either ancestry") +
  guides(colour = guide_legend(order = 1, nrow = 1), shape = guide_legend(order = 2)) +
  theme_figure2 +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5),
        legend.direction = "horizontal", legend.box = "horizontal")

p_mse <- ggplot(mse_data, aes(method, log10(as.numeric(mse) + 1), fill = method)) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.35) +
  geom_boxplot(width = 0.70, outlier.size = 0.45, linewidth = 0.32) +
  facet_wrap(~ annotation, scales = "free_y", nrow = 1) +
  scale_fill_manual(values = method_colors, drop = FALSE) +
  labs(x = NULL, y = expression(log[10]("enrichment MSE + 1")),
       title = "Functional-enrichment recovery at PIP > 0.5") +
  theme_figure2 + theme(legend.position = "none")

# Only p_pr contributes guides. guide_area() receives that one method legend;
# the other panels keep legend.position = "none" and cannot generate duplicate
# fill legends when patchwork assembles the figure.
figure2_extended <-
  ((p_size | p_power | p_fdr05 | p_fdr10 | p_pr) /
     p_mse /
     guide_area()) +
  plot_annotation(tag_levels = "a") +
  plot_layout(heights = c(1, 1, 0.12), guides = "collect")

output_stem <- file.path(FIGURE_DIR, "Baseline_Figure2_with_SuSiEx_MultiSuSiE")
ggsave(paste0(output_stem, ".pdf"), figure2_extended,
       width = 22, height = 9, units = "in", device = cairo_pdf)
ggsave(paste0(output_stem, ".png"), figure2_extended,
       width = 22, height = 9, units = "in", dpi = 600, bg = "white")
saveRDS(list(plot = figure2_extended, size_data = size_data, power_data = power_data,
             fdr05_data = fdr05_data, fdr10_data = fdr10_data,
             pr_data = pr_data, mse_data = mse_data),
        file.path(RESULT_DIR, "Baseline_Figure2_with_SuSiEx_MultiSuSiE.rds"),
        compress = "xz")

cat("Extended Figure 2 saved to:\n", paste0(output_stem, ".pdf"), "\n",
    paste0(output_stem, ".png"), "\n", sep = "")
