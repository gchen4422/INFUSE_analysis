#!/usr/bin/env Rscript
.infuse_repo <- Sys.getenv("INFUSE_ANALYSIS_ROOT", "/path/to/INFUSE_analysis")

suppressPackageStartupMessages(library(data.table))

file_arg <- commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))][1L]
script_dir <- dirname(normalizePath(sub("^--file=", "", file_arg)))
source(file.path(.infuse_repo, "Simulation/revision/configuration_classification/00_common.R"))

inputs <- file.path(result_out_dir, c(
  "configuration_predictions_d10.rds",
  "configuration_predictions_d187.rds"
))
if (!all(file.exists(inputs))) stop("Both extraction outputs are required.", call. = FALSE)
x <- lapply(inputs, readRDS)
causal <- rbindlist(lapply(x, `[[`, "causal"), use.names = TRUE)
roc_counts <- rbindlist(lapply(x, `[[`, "roc_counts"), use.names = TRUE)
cells <- rbindlist(lapply(x, `[[`, "cells"), use.names = TRUE)

# The two dimensions must use exactly the same causal truth in every cell.
truth_key <- c("causal_num", "h2", "locus", "snp_index", "SNP", "Signal", "truth_class")
truth_d10 <- unique(causal[dimension == "d = 10" & method == "INFUSE", ..truth_key])
truth_d187 <- unique(causal[dimension == "d = 187" & method == "INFUSE", ..truth_key])
setorder(truth_d10, causal_num, h2, locus, snp_index)
setorder(truth_d187, causal_num, h2, locus, snp_index)
if (!identical(truth_d10, truth_d187)) {
  stop("d=10 and d=187 do not have identical causal truth rows.", call. = FALSE)
}

metric_long <- function(d) {
  v <- classification_metrics(d)
  data.table(metric = names(v), estimate = as.numeric(v))
}

point_metrics <- causal[, metric_long(.SD), by = .(dimension, method)]

# Locus-cluster bootstrap. Each resampled locus retains every h2 and causal-number
# setting; the same sampled locus sequence is used across dimensions and methods.
set.seed(20260828)
B <- 1000L
loci <- sort(unique(causal$locus))
bootstrap_rows <- vector("list", B * 4L)
ii <- 0L
for (b in seq_len(B)) {
  sampled <- sample(loci, length(loci), replace = TRUE)
  boot <- rbindlist(lapply(sampled, function(ll) causal[locus == ll]), use.names = TRUE)
  for (dd in c("d = 10", "d = 187")) {
    for (mm in c("INFUSE", "MESuSiE")) {
      ii <- ii + 1L
      vals <- classification_metrics(boot[dimension == dd & method == mm])
      bootstrap_rows[[ii]] <- data.table(
        bootstrap = b, dimension = dd, method = mm,
        metric = names(vals), estimate = as.numeric(vals)
      )
    }
  }
  if (b %% 100L == 0L) message("Bootstrap ", b, "/", B)
}
bootstrap <- rbindlist(bootstrap_rows, use.names = TRUE)
metric_ci <- bootstrap[, .(
  conf_low = quantile(estimate, 0.025, na.rm = TRUE),
  conf_high = quantile(estimate, 0.975, na.rm = TRUE)
), by = .(dimension, method, metric)]
metric_summary <- merge(point_metrics, metric_ci, by = c("dimension", "method", "metric"), all.x = TRUE)

# Three-class confusion matrices, normalized within the true class.
confusion <- causal[, .N, by = .(dimension, method, truth_class, pred_class)]
confusion[, proportion := N / sum(N), by = .(dimension, method, truth_class)]

# Fixed-bin calibration of P(shared) among true causal variants.
calibration <- causal[, {
  breaks <- seq(0, 1, by = 0.1)
  bin <- cut(q_shared, breaks = breaks, include.lowest = TRUE, right = TRUE)
  d <- data.table(bin = bin, q_shared = q_shared, observed = as.integer(truth_class == "Shared"))
  d[, .(
    n = .N, mean_predicted = mean(q_shared), observed_shared = mean(observed)
  ), by = bin]
}, by = .(dimension, method)]

# Binary ROC and PR curves conditional on causality.
binary_thresholds <- seq(1, 0, by = -0.01)
binary_curves <- causal[, {
  y <- truth_class == "Shared"
  rbindlist(lapply(binary_thresholds, function(tt) {
    called <- q_shared >= tt
    tp <- sum(called & y); fp <- sum(called & !y)
    fn <- sum(!called & y); tn <- sum(!called & !y)
    data.table(
      threshold = tt,
      sensitivity = if ((tp + fn) > 0) tp / (tp + fn) else NA_real_,
      false_positive_rate = if ((fp + tn) > 0) fp / (fp + tn) else NA_real_,
      precision = if ((tp + fp) > 0) tp / (tp + fp) else 1,
      recall = if ((tp + fn) > 0) tp / (tp + fn) else NA_real_
    )
  }))
}, by = .(dimension, method)]

# Classification accuracy as increasingly confident causal-SNP calls are kept.
confidence_thresholds <- unique(c(1 / 3, seq(0.4, 0.95, by = 0.025)))
confidence_coverage <- causal[, rbindlist(lapply(confidence_thresholds, function(tt) {
  keep <- q_max >= tt
  data.table(
    threshold = tt,
    coverage = mean(keep),
    accuracy = if (any(keep)) mean(pred_class[keep] == truth_class[keep]) else NA_real_,
    n_kept = sum(keep)
  )
})), by = .(dimension, method)]

# Supporting end-to-end curves: pool discovery counts over all simulation cells.
end_to_end <- roc_counts[, .(
  TP = sum(TP), FP = sum(FP), n_true = sum(n_true), n_variant = sum(n_variant)
), by = .(dimension, method, configuration, cutoff)]
end_to_end[, `:=`(
  recall = TP / n_true,
  precision = fifelse(TP + FP > 0, TP / (TP + FP), 1),
  FDR = fifelse(TP + FP > 0, FP / (TP + FP), 0)
)]
power_at_fdr05 <- end_to_end[FDR <= 0.05,
  .SD[which.max(recall)], by = .(dimension, method, configuration)]

# Counts and stratified summaries make hidden class imbalance visible.
class_counts <- causal[, .N, by = .(dimension, method, truth_class, causal_num, h2)]
localization_strata <- causal[, .(
  n = .N,
  accuracy_3class = mean(pred_class == truth_class),
  shared_binary_accuracy = mean((q_shared >= 0.5) == (truth_class == "Shared"))
), by = .(
  dimension, method,
  pip_stratum = cut(pip_either, c(-Inf, 0.1, 0.5, Inf),
                    labels = c("PIP < 0.1", "0.1 <= PIP < 0.5", "PIP >= 0.5"))
)]

support_diagnostics <- causal[, .(
  n_causal = .N,
  zero_support = sum(support_sum <= 0),
  fraction_zero_support = mean(support_sum <= 0),
  minimum_support_sum = min(support_sum),
  fraction_marginal_pip_below_001 = mean(pip_either < 0.01),
  fraction_marginal_pip_below_01 = mean(pip_either < 0.1),
  median_abs_support_sum_minus_pip = median(abs(support_sum - pip_either)),
  max_abs_support_sum_minus_pip = max(abs(support_sum - pip_either))
), by = .(dimension, method)]

fwrite(metric_summary, file.path(result_out_dir, "classification_metric_summary.tsv"), sep = "\t")
fwrite(confusion, file.path(result_out_dir, "three_class_confusion.tsv"), sep = "\t")
fwrite(calibration, file.path(result_out_dir, "shared_probability_calibration.tsv"), sep = "\t")
fwrite(confidence_coverage, file.path(result_out_dir, "confidence_coverage.tsv"), sep = "\t")
fwrite(end_to_end, file.path(result_out_dir, "end_to_end_configuration_curves.tsv"), sep = "\t")
fwrite(power_at_fdr05, file.path(result_out_dir, "end_to_end_power_at_fdr05.tsv"), sep = "\t")
fwrite(class_counts, file.path(result_out_dir, "class_counts.tsv"), sep = "\t")
fwrite(localization_strata, file.path(result_out_dir, "classification_by_localization_pip.tsv"), sep = "\t")
fwrite(support_diagnostics, file.path(result_out_dir, "normalized_support_diagnostics.tsv"), sep = "\t")

saveRDS(list(
  causal = causal, cells = cells, metric_summary = metric_summary,
  bootstrap = bootstrap, confusion = confusion, calibration = calibration,
  binary_curves = binary_curves, confidence_coverage = confidence_coverage,
  end_to_end = end_to_end, power_at_fdr05 = power_at_fdr05,
  class_counts = class_counts, localization_strata = localization_strata
  , support_diagnostics = support_diagnostics
), file.path(result_out_dir, "configuration_classification_analysis.rds"), compress = "xz")

message("Analysis complete. Causal rows: ", nrow(causal), "; bootstrap replicates: ", B)
