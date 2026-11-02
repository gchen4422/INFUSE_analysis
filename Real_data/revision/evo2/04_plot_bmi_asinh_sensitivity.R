#!/usr/bin/env Rscript
.infuse_repo <- Sys.getenv("INFUSE_ANALYSIS_ROOT", "/path/to/INFUSE_analysis")
.infuse_work <- Sys.getenv("INFUSE_REVISION", "/path/to/revision_output")

# BMI-only sensitivity analysis for the signed Evo 2 delta score.
#
# The primary analysis remains the untransformed raw-score analysis produced by
# 03_plot_enrichment_either_bmi_evo2.R. This script applies the odd, monotone
# transformation
#
#   asinh(raw Evo 2 delta score / c),
#
# where c is the median absolute non-zero score across the full BMI SNP
# universe.  The transformation preserves every score's sign and ordering,
# remains approximately linear near zero, and reduces the influence of the
# long tails.  The scale c is fixed before selecting SNPs and is therefore the
# same for MESuSiE and Classic Annotations.

module_root <- file.path(.infuse_work, "evo2_raw_score_external_validation")
source(file.path(.infuse_repo, "Real_data/revision/evo2/00_evo2_enrichment_common.R"))

bmi <- load_raw_evo2_result(
  file.path(module_root, "data", "res_BMI_snp_embedding_updated_raw_evo2.RData"),
  "BMI"
)

raw_nonzero <- bmi$evo2_delta_score[
  is.finite(bmi$evo2_delta_score) & bmi$evo2_delta_score != 0
]
asinh_scale <- median(abs(raw_nonzero))
if (!is.finite(asinh_scale) || asinh_scale <= 0) {
  stop("Could not determine a positive asinh scale.", call. = FALSE)
}
bmi[, evo2_asinh_score := asinh(evo2_delta_score / asinh_scale)]

cluster_bootstrap_score <- function(score, region, selected,
                                    n_boot = 1000L, seed = 1L) {
  z <- data.table(score = score, region = region, selected = selected)
  by_region <- z[, .(
    selected_sum = sum(score[selected]),
    selected_n = sum(selected)
  ), by = region]

  if (sum(by_region$selected_n) == 0L) {
    return(c(mean = NA_real_, ci_low = NA_real_, ci_high = NA_real_))
  }

  estimate <- sum(by_region$selected_sum) / sum(by_region$selected_n)
  set.seed(seed)
  n_region <- nrow(by_region)
  draw <- matrix(
    sample.int(n_region, n_region * n_boot, replace = TRUE),
    nrow = n_region,
    ncol = n_boot
  )
  boot_n <- colSums(matrix(by_region$selected_n[draw], nrow = n_region))
  boot_mean <- colSums(matrix(by_region$selected_sum[draw], nrow = n_region)) /
    boot_n

  c(
    mean = estimate,
    ci_low = unname(quantile(boot_mean, 0.025, na.rm = TRUE)),
    ci_high = unname(quantile(boot_mean, 0.975, na.rm = TRUE))
  )
}

paired_cluster_bootstrap <- function(score, region, selected_mesusie,
                                     selected_classic, n_boot = 1000L,
                                     seed = 1L) {
  z <- data.table(
    score = score,
    region = region,
    selected_mesusie = selected_mesusie,
    selected_classic = selected_classic
  )
  by_region <- z[, .(
    mesusie_sum = sum(score[selected_mesusie]),
    mesusie_n = sum(selected_mesusie),
    classic_sum = sum(score[selected_classic]),
    classic_n = sum(selected_classic)
  ), by = region]

  mesusie_mean <- sum(by_region$mesusie_sum) / sum(by_region$mesusie_n)
  classic_mean <- sum(by_region$classic_sum) / sum(by_region$classic_n)
  estimate <- mesusie_mean - classic_mean

  set.seed(seed)
  n_region <- nrow(by_region)
  draw <- matrix(
    sample.int(n_region, n_region * n_boot, replace = TRUE),
    nrow = n_region,
    ncol = n_boot
  )
  boot_mesusie <-
    colSums(matrix(by_region$mesusie_sum[draw], nrow = n_region)) /
    colSums(matrix(by_region$mesusie_n[draw], nrow = n_region))
  boot_classic <-
    colSums(matrix(by_region$classic_sum[draw], nrow = n_region)) /
    colSums(matrix(by_region$classic_n[draw], nrow = n_region))
  boot_difference <- boot_mesusie - boot_classic

  c(
    mesusie_mean = mesusie_mean,
    classic_mean = classic_mean,
    classic_advantage = estimate,
    ci_low = unname(quantile(boot_difference, 0.025, na.rm = TRUE)),
    ci_high = unname(quantile(boot_difference, 0.975, na.rm = TRUE))
  )
}

top_n_selection <- function(pip, snp, n_take) {
  ord <- order(-pip, snp, na.last = NA)
  selected <- rep(FALSE, length(pip))
  if (n_take > 0L) {
    selected[ord[seq_len(min(n_take, length(ord)))]] <- TRUE
  }
  selected
}

summarize_rule <- function(d, cutoffs, rule, n_boot = 1000L,
                           seed = 20260828L) {
  method_rows <- vector("list", length(cutoffs) * length(method_columns))
  contrast_rows <- vector("list", length(cutoffs))
  row_id <- 1L

  for (i in seq_along(cutoffs)) {
    if (rule == "matched_top_N_from_Classic_PIP_threshold") {
      classic_pip <- d[[method_columns[["Classic_Annotations"]]]]
      n_take <- sum(is.finite(classic_pip) & classic_pip > cutoffs[i])
    } else {
      n_take <- as.integer(cutoffs[i])
    }
    if (n_take == 0L) next

    selected_by_method <- lapply(method_columns, function(column) {
      top_n_selection(d[[column]], d$SNP, n_take)
    })

    for (m in seq_along(method_columns)) {
      method <- names(method_columns)[m]
      selected <- selected_by_method[[method]]
      stat <- cluster_bootstrap_score(
        d$evo2_asinh_score,
        d$Region,
        selected,
        n_boot = n_boot,
        seed = seed + 100L * i + m
      )
      method_rows[[row_id]] <- data.table(
        Trait = "BMI",
        method = method,
        selection_rule = rule,
        cutoff = cutoffs[i],
        n_selected = sum(selected),
        n_regions_selected = uniqueN(d$Region[selected]),
        asinh_scale_raw_evo2 = asinh_scale,
        asinh_evo2_mean = stat[["mean"]],
        asinh_evo2_ci_low = stat[["ci_low"]],
        asinh_evo2_ci_high = stat[["ci_high"]]
      )
      row_id <- row_id + 1L
    }

    contrast <- paired_cluster_bootstrap(
      d$evo2_asinh_score,
      d$Region,
      selected_by_method[["MESuSiE"]],
      selected_by_method[["Classic_Annotations"]],
      n_boot = n_boot,
      seed = seed + 10000L + i
    )
    contrast_rows[[i]] <- data.table(
      Trait = "BMI",
      selection_rule = rule,
      cutoff = cutoffs[i],
      n_selected_per_method = n_take,
      asinh_scale_raw_evo2 = asinh_scale,
      mesusie_asinh_mean = contrast[["mesusie_mean"]],
      classic_asinh_mean = contrast[["classic_mean"]],
      classic_deleterious_advantage = contrast[["classic_advantage"]],
      contrast_ci_low = contrast[["ci_low"]],
      contrast_ci_high = contrast[["ci_high"]]
    )
  }

  list(
    methods = rbindlist(method_rows, use.names = TRUE, fill = TRUE),
    contrast = rbindlist(contrast_rows, use.names = TRUE, fill = TRUE)
  )
}

make_asinh_mean_plot <- function(summary, x_label, title) {
  ggplot(summary, aes(
    x = cutoff,
    y = asinh_evo2_mean,
    color = method,
    fill = method,
    group = method
  )) +
    geom_hline(yintercept = 0, color = "grey55", linewidth = 0.4,
               linetype = "dashed") +
    geom_ribbon(
      aes(ymin = asinh_evo2_ci_low, ymax = asinh_evo2_ci_high),
      color = NA,
      alpha = 0.13
    ) +
    geom_line(linewidth = 0.8) +
    geom_point(size = 1.9) +
    scale_color_manual(
      values = method_colors,
      breaks = names(method_columns),
      labels = method_legend_labels,
      drop = FALSE
    ) +
    scale_fill_manual(values = method_colors, drop = FALSE, guide = "none") +
    labs(
      title = title,
      subtitle = paste0(
        "Signed asinh transformation; c = median |raw score| = ",
        format(asinh_scale, scientific = TRUE, digits = 3),
        ".\nLower values indicate greater predicted disruption"
      ),
      x = x_label,
      y = expression(paste("Mean asinh(raw Evo 2 score / ", italic(c), ")"))
    ) +
    publication_theme()
}

make_contrast_plot <- function(contrast, x_label, title) {
  ggplot(contrast, aes(x = cutoff, y = classic_deleterious_advantage)) +
    geom_hline(yintercept = 0, color = "grey45", linewidth = 0.45,
               linetype = "dashed") +
    geom_ribbon(
      aes(ymin = contrast_ci_low, ymax = contrast_ci_high),
      fill = method_colors[["Classic_Annotations"]],
      color = NA,
      alpha = 0.15
    ) +
    geom_line(color = method_colors[["Classic_Annotations"]], linewidth = 0.85) +
    geom_point(color = method_colors[["Classic_Annotations"]], size = 2) +
    labs(
      title = title,
      subtitle = paste(
        "Positive values mean Classic Annotations selects lower, more deleterious-direction Evo 2 scores.",
        "Bands are locus-bootstrap 95% CIs.",
        sep = "\n"
      ),
      x = x_label,
      y = "MESuSiE mean - Classic Annotations mean"
    ) +
    publication_theme() +
    theme(legend.position = "none")
}

thresholds <- summarize_rule(
  bmi,
  seq(0.50, 0.90, 0.05),
  "matched_top_N_from_Classic_PIP_threshold"
)
top_n <- summarize_rule(bmi, c(10L, 20L, 50L, 100L, 200L), "top_N")

result_dir <- file.path(module_root, "results")
figure_dir <- file.path(dirname(module_root), "figures",
                        "evo2_raw_score_external_validation")
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)

fwrite(
  thresholds$methods,
  file.path(result_dir, "asinh_evo2_by_pip_threshold_bmi.tsv"),
  sep = "\t"
)
fwrite(
  thresholds$contrast,
  file.path(result_dir, "asinh_evo2_method_contrast_by_pip_threshold_bmi.tsv"),
  sep = "\t"
)
fwrite(
  top_n$methods,
  file.path(result_dir, "asinh_evo2_by_top_n_bmi.tsv"),
  sep = "\t"
)
fwrite(
  top_n$contrast,
  file.path(result_dir, "asinh_evo2_method_contrast_by_top_n_bmi.tsv"),
  sep = "\t"
)

save_plot_pair(
  make_asinh_mean_plot(
    thresholds$methods,
    "PIP cutoff",
    "Sign-preserving Evo 2 sensitivity analysis (BMI)"
  ),
  "asinh_evo2_mean_by_pip_threshold_bmi",
  figure_dir,
  7.8,
  4.2
)
save_plot_pair(
  make_contrast_plot(
    thresholds$contrast,
    "PIP cutoff",
    "Classic-Annotations Evo 2 contrast (BMI)"
  ),
  "asinh_evo2_classic_advantage_by_pip_threshold_bmi",
  figure_dir,
  7.8,
  4.2
)
save_plot_pair(
  make_asinh_mean_plot(
    top_n$methods,
    "Number of top-ranked SNPs",
    "Sign-preserving Evo 2 top-N sensitivity analysis (BMI)"
  ),
  "asinh_evo2_mean_by_top_n_bmi",
  figure_dir,
  7.8,
  4.2
)
save_plot_pair(
  make_contrast_plot(
    top_n$contrast,
    "Number of top-ranked SNPs",
    "Classic-Annotations top-N Evo 2 contrast (BMI)"
  ),
  "asinh_evo2_classic_advantage_by_top_n_bmi",
  figure_dir,
  7.8,
  4.2
)

message("BMI asinh sensitivity analysis complete. asinh scale c = ", asinh_scale)
