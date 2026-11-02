#!/usr/bin/env Rscript
.infuse_repo <- Sys.getenv("INFUSE_ANALYSIS_ROOT", "/path/to/INFUSE_analysis")

# Summarize the matched control and four sparse-LD scenarios from the rerun in
# which anchor sampling is independent of z scores and causal status.

this_file <- normalizePath(sub(
  "^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]
), mustWork = FALSE)
script_dir <- dirname(this_file)
source(file.path(.infuse_repo, "Simulation/revision/ld_mismatch/ld_mismatch_common.R"))

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
})

required_scenarios <- c(
  "matched",
  "sparse_attenuation_025",
  "sparse_attenuation_050",
  "sparse_negative_050",
  "sparse_negative_100"
)
expected_files <- unlist(lapply(required_scenarios, function(scenario) {
  sprintf("ld_mismatch_%s_locus_%03d.RData", scenario, seq_len(100L))
}), use.names = FALSE)
per_task_dir <- file.path(RESULT_DIR, "per_task")
files <- file.path(per_task_dir, expected_files)
dir.create(RESULT_DIR, recursive=TRUE, showWarnings=FALSE)
dir.create(FIGURE_DIR, recursive=TRUE, showWarnings=FALSE)
data.table::setDTthreads(1L)
missing_files <- files[!file.exists(files)]
if (length(missing_files)) {
  stop("Sparse no-Z-filter rerun is incomplete: ", length(missing_files),
       " of ", length(files), " files are missing.", call. = FALSE)
}

metrics_list <- targets_list <- detection_list <- s_list <- error_list <- timing_list <- list()
for (file in files) {
  e <- new.env(parent = emptyenv())
  load(file, envir = e)
  metrics_list[[length(metrics_list) + 1L]] <- e$metrics
  targets_list[[length(targets_list) + 1L]] <- e$targets
  s_list[[length(s_list) + 1L]] <- e$s_summary
  timing_list[[length(timing_list) + 1L]] <- e$timing
  if (exists("anchor_detection", e, inherits = FALSE) && nrow(e$anchor_detection)) {
    detection_list[[length(detection_list) + 1L]] <- e$anchor_detection
  }
  if (exists("errors", e, inherits = FALSE) && nrow(e$errors)) {
    error_list[[length(error_list) + 1L]] <- e$errors
  }
}
bind <- function(x) if (length(x)) rbindlist(x, fill = TRUE) else data.table()
metrics <- bind(metrics_list)
targets <- bind(targets_list)
anchor_detection <- bind(detection_list)
s_summary <- bind(s_list)
errors <- bind(error_list)
timing <- bind(timing_list)

if (nrow(errors)) {
  fwrite(errors, file.path(RESULT_DIR, "failed_cells.csv"))
  stop("At least one diagnostic cell failed.", call. = FALSE)
}
if (nrow(metrics) != length(files) * 6L * 6L) {
  stop("Unexpected number of metric rows: ", nrow(metrics), call. = FALSE)
}
if (nrow(anchor_detection) == 0L) {
  stop("No per-anchor detection records were found.", call. = FALSE)
}

# Preserve the pooled, variant-weighted point estimates while allowing all
# SNP-level outcomes from the same locus to move together in the uncertainty
# calculation. The same locus resamples are used across groups so that paired
# comparisons remain possible.
bootstrap_reps <- 2000L
bootstrap_seed <- 20260904L
bootstrap_loci <- sort(unique(metrics$locus))
if (length(bootstrap_loci) != 100L) {
  stop("Expected 100 independent loci for the cluster bootstrap.", call. = FALSE)
}
set.seed(bootstrap_seed)
bootstrap_index <- replicate(
  bootstrap_reps,
  sample.int(length(bootstrap_loci), length(bootstrap_loci), replace = TRUE)
)

locus_ratio_ci <- function(numerator, denominator, locus, scale = 1) {
  if (length(numerator) != length(denominator) ||
      length(numerator) != length(locus)) {
    stop("Cluster-bootstrap inputs have incompatible lengths.", call. = FALSE)
  }
  if (any(!is.finite(numerator)) || any(!is.finite(denominator)) ||
      any(denominator < 0)) {
    stop("Cluster-bootstrap counts must be finite and nonnegative.", call. = FALSE)
  }

  by_locus <- data.table(
    locus = locus,
    numerator = numerator,
    denominator = denominator
  )[, .(
    numerator = sum(numerator),
    denominator = sum(denominator)
  ), by = locus]

  locus_position <- match(by_locus$locus, bootstrap_loci)
  if (anyNA(locus_position)) {
    stop("Unexpected locus in cluster-bootstrap input.", call. = FALSE)
  }
  numerator_by_locus <- denominator_by_locus <- numeric(length(bootstrap_loci))
  numerator_by_locus[locus_position] <- by_locus$numerator
  denominator_by_locus[locus_position] <- by_locus$denominator

  total_denominator <- sum(denominator_by_locus)
  estimate <- if (total_denominator > 0) {
    scale * sum(numerator_by_locus) / total_denominator
  } else {
    NA_real_
  }

  boot_numerator <- colSums(matrix(
    numerator_by_locus[bootstrap_index], nrow = length(bootstrap_loci)
  ))
  boot_denominator <- colSums(matrix(
    denominator_by_locus[bootstrap_index], nrow = length(bootstrap_loci)
  ))
  boot_ratio <- rep(NA_real_, bootstrap_reps)
  valid <- boot_denominator > 0
  boot_ratio[valid] <- scale * boot_numerator[valid] / boot_denominator[valid]

  if (!any(valid)) {
    return(c(estimate = estimate, low = NA_real_, high = NA_real_))
  }
  interval <- quantile(
    boot_ratio[valid], c(0.025, 0.975), names = FALSE, na.rm = TRUE
  )
  c(estimate = estimate, low = interval[[1L]], high = interval[[2L]])
}

primary_tests <- c("separate_union", "joint", "joint_separate_s")
test_labels <- c(
  separate_union = "Separate tests (union)",
  joint = "Joint test (shared s)",
  joint_separate_s = "Joint test (separate s)"
)
test_colors <- c(
  "Separate tests (union)" = "#E41A1C",
  "Joint test (shared s)" = "#377EB8",
  "Joint test (separate s)" = "#4DAF4A"
)
scenario_labels <- c(
  sparse_negative_100 = "Full\nreversal",
  sparse_negative_050 = "Partial\nreversal",
  sparse_attenuation_050 = "Moderate\nattenuation",
  sparse_attenuation_025 = "Mild\nattenuation"
)
scenario_order <- names(scenario_labels)

sparse_metrics <- metrics |>
  filter(scenario_class == "sparse", test %in% primary_tests) |>
  group_by(scenario, anchor_scale, test, causal_num, h2_num) |>
  group_modify(~{
    sensitivity_ci <- locus_ratio_ci(.x$TP, .x$n_truth, .x$locus)
    fdr_ci <- locus_ratio_ci(.x$FP, .x$TP + .x$FP, .x$locus)
    fp_ci <- locus_ratio_ci(.x$FP, .x$FP + .x$TN, .x$locus, scale = 1e5)
    tibble(
      TP = sum(.x$TP),
      FP = sum(.x$FP),
      TN = sum(.x$TN),
      n_truth = sum(.x$n_truth),
      Sensitivity = sensitivity_ci[["estimate"]],
      sensitivity_low = sensitivity_ci[["low"]],
      sensitivity_high = sensitivity_ci[["high"]],
      FDR = fdr_ci[["estimate"]],
      fdr_low = fdr_ci[["low"]],
      fdr_high = fdr_ci[["high"]],
      # FDR is not informative when a test has zero true detections: one or two
      # false calls make FDR equal to one despite a negligible genome-wide error
      # burden. Report false calls per 100,000 non-anchor variants in the primary
      # figure instead.
      FP_per_100k = fp_ci[["estimate"]],
      fp_per_100k_low = fp_ci[["low"]],
      fp_per_100k_high = fp_ci[["high"]]
    )
  }) |>
  ungroup() |>
  mutate(
    Test = factor(unname(test_labels[test]), levels = unname(test_labels)),
    Perturbation = factor(unname(scenario_labels[scenario]),
                          levels = unname(scenario_labels[scenario_order]))
  )

matched_metrics <- metrics |>
  filter(scenario == "matched", test %in% primary_tests) |>
  group_by(test, causal_num, h2_num) |>
  group_modify(~{
    flag_ci <- locus_ratio_ci(.x$n_flagged, .x$n_snp, .x$locus)
    tibble(
      n_flagged = sum(.x$n_flagged),
      n_snp = sum(.x$n_snp),
      flag_rate = flag_ci[["estimate"]],
      ci_low = flag_ci[["low"]],
      ci_high = flag_ci[["high"]]
    )
  }) |>
  ungroup() |>
  mutate(Test = factor(unname(test_labels[test]), levels = unname(test_labels)))

z_stratum_levels <- c("|Z|<1", "1<=|Z|<2", "2<=|Z|<4", "|Z|>=4")
z_stratified <- anchor_detection |>
  filter(scenario_class == "sparse", test %in% primary_tests) |>
  group_by(scenario, test, z_stratum) |>
  group_modify(~{
    by_locus <- .x |>
      group_by(locus) |>
      summarise(n_anchor = n(), n_detected = sum(detected), .groups = "drop")
    sensitivity_ci <- locus_ratio_ci(
      by_locus$n_detected, by_locus$n_anchor, by_locus$locus
    )
    tibble(
      n_anchor = sum(by_locus$n_anchor),
      n_detected = sum(by_locus$n_detected),
      sensitivity = sensitivity_ci[["estimate"]],
      ci_low = sensitivity_ci[["low"]],
      ci_high = sensitivity_ci[["high"]]
    )
  }) |>
  ungroup() |>
  mutate(
    Test = factor(unname(test_labels[test]), levels = unname(test_labels)),
    z_stratum = factor(z_stratum, levels = z_stratum_levels),
    Perturbation = factor(unname(scenario_labels[scenario]),
                          levels = unname(scenario_labels[scenario_order]))
  )

# The same anchors recur across the four perturbation geometries. Deduplicate
# them before characterizing the sampling distribution.
unique_targets <- targets |>
  filter(scenario_class == "sparse") |>
  distinct(locus, causal_num, h2_num, ancestry, SNP_index, .keep_all = TRUE) |>
  mutate(z_stratum = factor(as.character(anchor_signal_bin(abs_anchor_z)),
                            levels = z_stratum_levels))
anchor_distribution <- unique_targets |>
  count(ancestry, z_stratum, name = "n_anchor") |>
  group_by(ancestry) |>
  mutate(proportion = n_anchor / sum(n_anchor)) |>
  ungroup()
anchor_sampling_summary <- unique_targets |>
  group_by(ancestry) |>
  summarise(
    n_anchor = n(),
    median_abs_z = median(abs_anchor_z),
    fraction_abs_z_gt_2 = mean(abs_anchor_z > 2),
    fraction_abs_z_gt_4 = mean(abs_anchor_z >= 4),
    median_max_abs_r = median(max_abs_r),
    fraction_causal = mean(Signal != 0),
    .groups = "drop"
  )

cell_sensitivity <- anchor_detection |>
  filter(scenario_class == "sparse", test %in% primary_tests) |>
  group_by(scenario, locus, causal_num, h2_num, test) |>
  summarise(sensitivity = mean(detected), .groups = "drop") |>
  pivot_wider(names_from = test, values_from = sensitivity) |>
  mutate(
    separate_minus_joint_shared = separate_union - joint,
    joint_separate_s_minus_shared_s = joint_separate_s - joint
  )
paired_difference_summary <- cell_sensitivity |>
  group_by(scenario, causal_num, h2_num) |>
  summarise(
    mean_separate_minus_joint_shared = mean(separate_minus_joint_shared),
    mean_joint_separate_s_minus_shared_s = mean(joint_separate_s_minus_shared_s),
    sd_separate_minus_joint_shared = sd(separate_minus_joint_shared),
    sd_joint_separate_s_minus_shared_s = sd(joint_separate_s_minus_shared_s),
    n_loci = n_distinct(locus),
    .groups = "drop"
  )

s_parameter_summary <- s_summary |>
  filter(scenario %in% required_scenarios) |>
  group_by(scenario, causal_num, h2_num) |>
  summarise(
    mean_s_shared = mean(s_shared),
    mean_s_EUR = mean(s_EUR),
    mean_s_AFR = mean(s_AFR),
    .groups = "drop"
  )

dir.create(RESULT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FIGURE_DIR, recursive = TRUE, showWarnings = FALSE)
fwrite(as.data.table(sparse_metrics), file.path(RESULT_DIR, "sparse_metrics_no_z_filter.csv"))
fwrite(as.data.table(matched_metrics), file.path(RESULT_DIR, "matched_metrics_no_z_filter.csv"))
fwrite(as.data.table(z_stratified), file.path(RESULT_DIR, "sensitivity_by_post_sampling_z_stratum.csv"))
fwrite(as.data.table(anchor_sampling_summary), file.path(RESULT_DIR, "anchor_sampling_summary.csv"))
fwrite(as.data.table(anchor_distribution), file.path(RESULT_DIR, "anchor_z_stratum_distribution.csv"))
fwrite(as.data.table(paired_difference_summary), file.path(RESULT_DIR, "paired_test_sensitivity_differences.csv"))
fwrite(as.data.table(s_parameter_summary), file.path(RESULT_DIR, "s_parameter_summary.csv"))

theme_revision <- function() {
  theme_bw(base_size = 9) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_blank(),
      strip.background = element_blank(),
      strip.text = element_text(face = "bold"),
      axis.title = element_text(face = "bold"),
      legend.position = "bottom",
      legend.title = element_blank()
    )
}

p_sensitivity <- ggplot(
  sparse_metrics,
  aes(Perturbation, Sensitivity, color = Test, group = Test)
) +
  geom_errorbar(aes(ymin = sensitivity_low, ymax = sensitivity_high),
                width = 0.12, position = position_dodge(width = 0.25)) +
  geom_line(position = position_dodge(width = 0.25), linewidth = 0.55) +
  geom_point(position = position_dodge(width = 0.25), size = 1.5) +
  facet_grid(h2_num ~ causal_num, labeller = label_both) +
  scale_color_manual(values = test_colors) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(x = "Selected-anchor LD perturbation", y = "Sensitivity") +
  theme_revision() +
  theme(axis.text.x = element_text(size = 6.5))

p_false_positive_burden <- ggplot(
  sparse_metrics,
  aes(Perturbation, FP_per_100k, color = Test, group = Test)
) +
  geom_errorbar(aes(ymin = fp_per_100k_low, ymax = fp_per_100k_high),
                width = 0.12, position = position_dodge(width = 0.25)) +
  geom_line(position = position_dodge(width = 0.25), linewidth = 0.55) +
  geom_point(position = position_dodge(width = 0.25), size = 1.5) +
  facet_grid(h2_num ~ causal_num, labeller = label_both) +
  scale_color_manual(values = test_colors) +
  scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.08))) +
  labs(
    x = "Selected-anchor LD perturbation",
    y = "False-positive calls per 100,000 non-anchor variants"
  ) +
  theme_revision() +
  theme(axis.text.x = element_text(size = 6.5))

p_z_strata <- ggplot(
  z_stratified,
  aes(z_stratum, sensitivity, color = Test, group = Test)
) +
  geom_errorbar(aes(ymin = ci_low, ymax = ci_high), width = 0.12,
                position = position_dodge(width = 0.2)) +
  geom_line(position = position_dodge(width = 0.2), linewidth = 0.55) +
  geom_point(position = position_dodge(width = 0.2), size = 1.5) +
  facet_wrap(~Perturbation, ncol = 2) +
  scale_color_manual(values = test_colors) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(
    x = "Anchor signal stratum (defined after random sampling)",
    y = "Sensitivity"
  ) +
  theme_revision() +
  theme(axis.text.x = element_text(angle = 25, hjust = 1))

p_anchor_distribution <- ggplot(
  anchor_distribution,
  aes(z_stratum, proportion, fill = ancestry)
) +
  geom_col(position = "dodge", width = 0.72) +
  scale_fill_manual(values = c(EUR = "#377EB8", AFR = "#E41A1C")) +
  labs(
    x = "Anchor |Z| stratum",
    y = "Proportion of randomly selected anchors"
  ) +
  theme_revision()

p_primary <- (p_sensitivity / p_false_positive_burden) +
  plot_layout(guides = "collect") +
  plot_annotation(
    title = "Primary result: random LD-connected anchors (no Z-score selection)",
    subtitle = paste(
      "Reversal multipliers: d = -1 and -0.5;",
      "moderate and mild sign-preserving attenuation:",
      "d = sqrt(0.50) and sqrt(0.75), respectively."
    ),
    tag_levels = "a"
  ) &
  theme(legend.position = "bottom")
p_interpretation <- (p_anchor_distribution / p_z_strata) +
  plot_layout(heights = c(1, 2), guides = "collect") +
  plot_annotation(tag_levels = "a") & theme(legend.position = "bottom")

primary_figure_dir <- FIGURE_DIR
dir.create(primary_figure_dir, recursive = TRUE, showWarnings = FALSE)
ggsave(file.path(primary_figure_dir,
                 "PRIMARY_Sparse_LD_joint_vs_separate_no_Z_filter.pdf"),
       p_primary, width = 10, height = 8, device = cairo_pdf)
ggsave(file.path(primary_figure_dir,
                 "PRIMARY_Sparse_LD_joint_vs_separate_no_Z_filter.png"),
       p_primary, width = 10, height = 8, dpi = 400, bg = "white")
ggsave(file.path(FIGURE_DIR, "Sparse_LD_post_sampling_Z_strata.pdf"),
       p_interpretation, width = 9, height = 8, device = cairo_pdf)
ggsave(file.path(FIGURE_DIR, "Sparse_LD_post_sampling_Z_strata.png"),
       p_interpretation, width = 9, height = 8, dpi = 400, bg = "white")

save(
  metrics, targets, anchor_detection, s_summary, timing, sparse_metrics,
  matched_metrics, z_stratified, unique_targets, anchor_distribution,
  anchor_sampling_summary, cell_sensitivity, paired_difference_summary,
  s_parameter_summary,
  file = file.path(RESULT_DIR, "LD_mismatch_sparse_no_z_filter_summary.RData")
)

cat("Completed sparse no-Z-filter summary in ", RESULT_DIR, "\n", sep = "")
