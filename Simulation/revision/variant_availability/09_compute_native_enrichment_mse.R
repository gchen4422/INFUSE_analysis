#!/usr/bin/env Rscript
# Figure labels use the manuscript name; stored result keys remain compatible.
availability_method_labels <- function(x) {
  sub("^MESuSiE(?=([ ~]|$))", "INFUSE\u2212", as.character(x), perl = TRUE)
}

.infuse_work <- Sys.getenv("INFUSE_REVISION", "/path/to/revision_output")

# Functional-enrichment recovery for the native AS-SV simulation.
#
# Usage:
#   Rscript 09_compute_native_enrichment_mse.R <assv_each:0|1|2>
#
# Summarizes enrichment estimates from saved fine-mapping results.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) {
  stop("Usage: Rscript 09_compute_native_enrichment_mse.R <assv_each:0|1|2>",
       call. = FALSE)
}
assv_each <- suppressWarnings(as.integer(args[[1L]]))
if (is.na(assv_each) || !assv_each %in% 0:2) {
  stop("assv_each must be 0, 1, or 2.", call. = FALSE)
}

suppressPackageStartupMessages(library(ggplot2))

module_dir <- file.path(.infuse_work, "INFUSE_variant_availability")
revision_dir <- dirname(module_dir)
native_dir <- file.path(module_dir, "simulation_native_assv")
scenario <- paste0("native_ASSV_each_", assv_each)
plot_dir <- file.path(revision_dir, "figures", "INFUSE_variant_availability", scenario)
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

base_counts <- c(2L, 3L, 5L)
h2_values <- c(`1` = 1e-4, `2` = 2e-4)
pip_cutoff <- 0.5

annotation_columns <- c(
  "H3K27ac_Hnisz_common",
  "DHS_Trynka_common",
  "H3K4me3_Trynka_common",
  "Conserved_LindbladToh_common",
  "non_synonymous_common"
)
annotation_labels <- c(
  H3K27ac_Hnisz_common = "H3K27ac",
  DHS_Trynka_common = "DHS",
  H3K4me3_Trynka_common = "H3K4me3",
  Conserved_LindbladToh_common = "Conserved",
  non_synonymous_common = "Non-synonymous"
)
annotation_order <- c("Non-synonymous", "Conserved", "DHS", "H3K27ac", "H3K4me3")

method_order <- c(
  "SuSiE", "SuSiEx intersect", "SuSiEx union", "MultiSuSiE",
  "MESuSiE intersect", "MESuSiE union",
  "INFUSE intersect", "INFUSE union", "INFUSE-H"
)
method_colors <- c(
  "SuSiE" = "#FBE79E",
  "SuSiEx intersect" = "#F2B56E",
  "SuSiEx union" = "#F57C6E",
  "MultiSuSiE" = "#CC79A7",
  "MESuSiE intersect" = "#A8DDD3",
  "MESuSiE union" = "#3F8C7F",
  "INFUSE intersect" = "#A5C8E1",
  "INFUSE union" = "#2B6CB0",
  "INFUSE-H" = "#6A3D9A"
)

main_file <- function(base_count, h2_index, locus) {
  file.path(native_dir, "results_native_sim", "cells",
            paste0("base_", base_count), paste0("h2_", h2_index),
            paste0("assv_each_", assv_each), sprintf("locus_%03d.rds", locus))
}
multisusie_file <- function(base_count, h2_index, locus) {
  file.path(native_dir, "results_native_sim", "multisusie_intersect",
            paste0("base_", base_count), paste0("h2_", h2_index),
            paste0("assv_each_", assv_each), sprintf("locus_%03d.rds", locus))
}
infuse_h_file <- function(base_count, h2_index, locus) {
  file.path(native_dir, "results_native_sim", "infuse_union_d187",
            paste0("base_", base_count), paste0("h2_", h2_index),
            paste0("assv_each_", assv_each), sprintf("locus_%03d.rds", locus))
}
annotation_file <- function(locus) {
  file.path(native_dir, "panels_source_preserved", sprintf("locus_%03d", locus),
            "annotation.rds")
}

extract_pips <- function(cell, multisusie, infuse_h) {
  list(
    "SuSiE" = cell$methods$susie$pip,
    "SuSiEx intersect" = cell$methods$susiex_intersect$pip,
    "SuSiEx union" = cell$methods$susiex_union$pip,
    "MultiSuSiE" = multisusie$pip,
    "MESuSiE intersect" = cell$methods$mesusie_intersect$pip,
    "MESuSiE union" = cell$methods$mesusie_union$pip,
    "INFUSE intersect" = cell$methods$infuse_intersect$pip,
    "INFUSE union" = cell$methods$infuse_union$pip,
    "INFUSE-H" = infuse_h$pip
  )
}

safe_ratio <- function(numerator_frequency, denominator_frequency) {
  if (!is.finite(numerator_frequency) || !is.finite(denominator_frequency) ||
      denominator_frequency <= 0) return(NA_real_)
  numerator_frequency / denominator_frequency
}

all_rows <- list()
row_i <- 0L

for (base_count in base_counts) {
  for (h2_index in 1:2) {
    total_n <- 0L
    total_ann <- setNames(numeric(length(annotation_columns)), annotation_columns)
    causal_n <- 0L
    causal_ann <- setNames(numeric(length(annotation_columns)), annotation_columns)
    selected_n <- setNames(integer(length(method_order)), method_order)
    selected_ann <- setNames(
      replicate(length(method_order),
                setNames(numeric(length(annotation_columns)), annotation_columns),
                simplify = FALSE),
      method_order
    )

    for (locus in 1:100) {
      paths <- c(main_file(base_count, h2_index, locus),
                 multisusie_file(base_count, h2_index, locus),
                 infuse_h_file(base_count, h2_index, locus),
                 annotation_file(locus))
      missing_paths <- paths[!file.exists(paths)]
      if (length(missing_paths)) {
        stop("Missing required input(s): ", paste(missing_paths, collapse = ", "),
             call. = FALSE)
      }

      cell <- readRDS(paths[[1L]])
      multisusie <- readRDS(paths[[2L]])
      infuse_h <- readRDS(paths[[3L]])
      annot <- readRDS(paths[[4L]])

      snp <- as.character(cell$truth$SNP)
      if (is.null(rownames(annot)) || !identical(rownames(annot), snp)) {
        stop("Annotation/SNP row-order mismatch at locus ", locus, ".", call. = FALSE)
      }
      if (!all(annotation_columns %in% colnames(annot))) {
        stop("Missing signal annotation column(s) at locus ", locus, ".", call. = FALSE)
      }
      A <- annot[, annotation_columns, drop = FALSE]
      if (!is.numeric(A) || anyNA(A) || any(!is.finite(A)) ||
          any(!A %in% c(0, 1))) {
        stop("Signal annotations must be finite binary values at locus ", locus, ".",
             call. = FALSE)
      }

      truth <- as.integer(cell$truth$Signal != 0L)
      if (length(truth) != nrow(A)) {
        stop("Truth/annotation length mismatch at locus ", locus, ".", call. = FALSE)
      }
      pips <- extract_pips(cell, multisusie, infuse_h)
      if (!identical(names(pips), method_order)) {
        stop("Internal method-order mismatch.", call. = FALSE)
      }
      for (method in method_order) {
        z <- as.numeric(pips[[method]])
        pip_tolerance <- 1e-12
        if (length(z) != nrow(A) || anyNA(z) || any(!is.finite(z)) ||
            any(z < -pip_tolerance | z > 1 + pip_tolerance)) {
          stop("Invalid ", method, " PIPs at locus ", locus, ".", call. = FALSE)
        }
        # Posterior probabilities can exceed an endpoint by machine roundoff
        # (observed maximum: 1 + 2.22e-16). Clamp only after the strict
        # tolerance check; this cannot change membership at PIP > 0.5.
        z <- pmin(1, pmax(0, z))
        selected <- z > pip_cutoff
        selected_n[[method]] <- selected_n[[method]] + sum(selected)
        if (any(selected)) {
          selected_ann[[method]] <- selected_ann[[method]] + colSums(A[selected, , drop = FALSE])
        }
      }

      total_n <- total_n + nrow(A)
      total_ann <- total_ann + colSums(A)
      causal_n <- causal_n + sum(truth)
      causal_ann <- causal_ann + colSums(A[truth == 1L, , drop = FALSE])
    }

    noncausal_n <- total_n - causal_n
    noncausal_ann <- total_ann - causal_ann
    if (causal_n <= 0L || noncausal_n <= 0L) {
      stop("Degenerate truth counts for base=", base_count, ", h2=", h2_index, ".",
           call. = FALSE)
    }

    for (method in method_order) {
      n_sel <- selected_n[[method]]
      n_unsel <- total_n - n_sel
      for (ann in annotation_columns) {
        causal_freq <- causal_ann[[ann]] / causal_n
        noncausal_freq <- noncausal_ann[[ann]] / noncausal_n
        true_enrichment <- safe_ratio(causal_freq, noncausal_freq)

        if (n_sel > 0L && n_unsel > 0L) {
          selected_freq <- selected_ann[[method]][[ann]] / n_sel
          unselected_freq <- (total_ann[[ann]] - selected_ann[[method]][[ann]]) / n_unsel
          estimated_enrichment <- safe_ratio(selected_freq, unselected_freq)
        } else {
          selected_freq <- NA_real_
          unselected_freq <- NA_real_
          estimated_enrichment <- NA_real_
        }

        row_i <- row_i + 1L
        all_rows[[row_i]] <- data.frame(
          assv_each = assv_each,
          base_causal_count = base_count,
          h2_index = h2_index,
          h2_per_causal = unname(h2_values[as.character(h2_index)]),
          method = method,
          annotation = unname(annotation_labels[[ann]]),
          annotation_column = ann,
          pip_cutoff = pip_cutoff,
          n_variants = total_n,
          n_causal = causal_n,
          n_noncausal = noncausal_n,
          n_selected = n_sel,
          n_unselected = n_unsel,
          causal_annotation_frequency = causal_freq,
          noncausal_annotation_frequency = noncausal_freq,
          selected_annotation_frequency = selected_freq,
          unselected_annotation_frequency = unselected_freq,
          true_enrichment = true_enrichment,
          estimated_enrichment = estimated_enrichment,
          squared_error = (estimated_enrichment - true_enrichment)^2,
          stringsAsFactors = FALSE
        )
      }
    }
  }
}

detail <- do.call(rbind, all_rows)
expected_rows <- length(base_counts) * 2L * length(method_order) *
  length(annotation_columns)
if (nrow(detail) != expected_rows) {
  stop("Unexpected number of result rows: ", nrow(detail),
       "; expected ", expected_rows, ".", call. = FALSE)
}
if (anyNA(detail$true_enrichment) || any(!is.finite(detail$true_enrichment))) {
  stop("Ground-truth enrichment contains undefined values.", call. = FALSE)
}

summary_split <- split(detail, list(detail$method, detail$annotation), drop = TRUE)
summary_rows <- lapply(summary_split, function(d) {
  ok <- is.finite(d$squared_error)
  data.frame(
    method = d$method[[1L]], annotation = d$annotation[[1L]],
    mean_mse = if (any(ok)) mean(d$squared_error[ok]) else NA_real_,
    median_mse = if (any(ok)) median(d$squared_error[ok]) else NA_real_,
    sd_mse = if (sum(ok) > 1L) sd(d$squared_error[ok]) else NA_real_,
    n_parameter_cells = nrow(d), n_defined = sum(ok),
    stringsAsFactors = FALSE
  )
})
summary_table <- do.call(rbind, summary_rows)
summary_table$method <- factor(summary_table$method, levels = method_order)
summary_table$annotation <- factor(summary_table$annotation, levels = annotation_order)
summary_table <- summary_table[order(summary_table$annotation, summary_table$method), ]
summary_table$method <- as.character(summary_table$method)
summary_table$annotation <- as.character(summary_table$annotation)

overall_split <- split(detail, detail$method)
overall_table <- do.call(rbind, lapply(overall_split, function(d) {
  ok <- is.finite(d$squared_error)
  data.frame(
    method = d$method[[1L]],
    mean_mse = if (any(ok)) mean(d$squared_error[ok]) else NA_real_,
    median_mse = if (any(ok)) median(d$squared_error[ok]) else NA_real_,
    n_annotation_parameter_cells = nrow(d), n_defined = sum(ok),
    stringsAsFactors = FALSE
  )
}))
overall_table$method <- factor(overall_table$method, levels = method_order)
overall_table <- overall_table[order(overall_table$method), ]
overall_table$method <- as.character(overall_table$method)

detail_path <- file.path(plot_dir, "functional_enrichment_estimates_PIP_gt_0.5.tsv")
summary_path <- file.path(plot_dir, "functional_enrichment_MSE_summary_PIP_gt_0.5.tsv")
overall_path <- file.path(plot_dir, "functional_enrichment_MSE_overall_PIP_gt_0.5.tsv")
write.table(detail, detail_path, sep = "\t", quote = FALSE, row.names = FALSE)
write.table(summary_table, summary_path, sep = "\t", quote = FALSE, row.names = FALSE)
write.table(overall_table, overall_path, sep = "\t", quote = FALSE, row.names = FALSE)

plot_data <- detail[is.finite(detail$squared_error), , drop = FALSE]
if (!nrow(plot_data)) stop("No finite MSE values are available for plotting.", call. = FALSE)
plot_data$method <- factor(plot_data$method, levels = method_order)
plot_data$annotation <- factor(plot_data$annotation, levels = annotation_order)

p <- ggplot(plot_data, aes(x = method, y = squared_error, fill = method)) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.35,
             colour = "black") +
  geom_boxplot(width = 0.64, linewidth = 0.42, outlier.size = 0.75,
               outlier.stroke = 0.25, colour = "#333333") +
  facet_wrap(~ annotation, scales = "free_y", nrow = 1L) +
  scale_x_discrete(labels = availability_method_labels) +
  scale_fill_manual(values = method_colors, name = "Method", drop = FALSE, labels = availability_method_labels) +
  scale_y_continuous(expand = expansion(mult = c(0.03, 0.10))) +
  labs(
    x = NULL,
    y = "MSE of Enrichment",
    title = "MSE Distribution by Annotation",
    subtitle = paste0("Native AS-SV each = ", assv_each,
                      "; putative causal signals: overall PIP > 0.5")
  ) +
  guides(fill = guide_legend(nrow = 1L, byrow = TRUE,
                             override.aes = list(colour = NA))) +
  theme_bw(base_size = 8, base_family = "sans") +
  theme(
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.title = element_text(size = 8.5, face = "bold"),
    legend.text = element_text(size = 7.2),
    legend.key.height = grid::unit(3.2, "mm"),
    legend.key.width = grid::unit(4.2, "mm"),
    legend.spacing.x = grid::unit(0.6, "mm"),
    axis.text.x = element_text(angle = 25, hjust = 1, vjust = 1, size = 6.3),
    axis.text.y = element_text(size = 7),
    axis.title.y = element_text(size = 8.5, face = "bold"),
    strip.background = element_blank(),
    strip.text = element_text(size = 9, face = "bold"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(colour = "black", linewidth = 0.35),
    plot.title = element_text(size = 11, hjust = 0.5),
    plot.subtitle = element_text(size = 7.5, hjust = 0.5),
    panel.spacing.x = grid::unit(2.0, "mm"),
    plot.margin = margin(5.5, 7.5, 5.5, 7.5)
  )

pdf_path <- file.path(plot_dir, "functional_enrichment_MSE_PIP_gt_0.5.pdf")
png_path <- file.path(plot_dir, "functional_enrichment_MSE_PIP_gt_0.5.png")
ggsave(pdf_path, p, width = 370, height = 125, units = "mm", device = cairo_pdf)
ggsave(png_path, p, width = 370, height = 125, units = "mm", dpi = 400,
       bg = "white", limitsize = FALSE)

saveRDS(list(detail = detail, summary = summary_table,
             overall = overall_table, plot = p),
        file.path(plot_dir, "functional_enrichment_MSE_PIP_gt_0.5.rds"))

cat("Completed ", scenario, "\n", sep = "")
cat("  ", pdf_path, "\n", sep = "")
cat("  ", png_path, "\n", sep = "")
cat("  selected-count range: ", min(detail$n_selected), " to ",
    max(detail$n_selected), "\n", sep = "")
cat("  undefined MSE rows: ", sum(!is.finite(detail$squared_error)), " / ",
    nrow(detail), "\n", sep = "")
