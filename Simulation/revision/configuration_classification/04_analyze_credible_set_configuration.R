#!/usr/bin/env Rscript
.infuse_repo <- Sys.getenv("INFUSE_ANALYSIS_ROOT", "/path/to/INFUSE_analysis")

suppressPackageStartupMessages(library(data.table))

file_arg <- commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))][1L]
script_dir <- dirname(normalizePath(sub("^--file=", "", file_arg)))
source(file.path(.infuse_repo, "Simulation/revision/configuration_classification/00_common.R"))

# This analysis treats each retained 95% credible set as the unit of
# classification.  A credible set inherits the configuration probabilities of
# its originating single-effect component.  Because alpha_l[j, c] is the joint
# posterior probability that component l selects SNP j and configuration c,
# summing alpha_l over SNPs gives a proper three-category component posterior.

object_map <- c(MESuSiE = "MESuSiE_res", INFUSE = "MESuSiE_res_mlk_annot_subset")
cs_rows <- list()
row_counter <- 0L
cell_counter <- 0L

for (d_label in c("d10", "d187")) {
  spec <- dimension_spec(d_label)

  for (causal_num in 1:3) {
    for (h2 in 1:2) {
      for (locus in 1:100) {
        cell_counter <- cell_counter + 1L
        stem <- sprintf("CAUSAL_%d_LOCI_%d_h2_%d", causal_num, locus, h2)
        result_file <- file.path(
          spec$result_root, paste0("causal_num_", causal_num), spec$result_subdir,
          paste0("Res_", stem, ".RData")
        )
        data_file <- file.path(
          spec$data_root, paste0("causal_num_", causal_num), spec$data_subdir, stem
        )

        if (!file.exists(result_file)) stop("Missing result: ", result_file, call. = FALSE)
        if (!file.exists(data_file)) stop("Missing truth data: ", data_file, call. = FALSE)

        z <- fread(data_file, data.table = FALSE)
        if (!all(c("RSID", "Signal") %in% names(z))) {
          stop("Truth columns missing in ", data_file, call. = FALSE)
        }

        e <- new.env(parent = emptyenv())
        load(result_file, envir = e)
        if (!all(unname(object_map) %in% ls(e))) {
          stop("Required objects missing in ", result_file, call. = FALSE)
        }

        for (method in names(object_map)) {
          fit <- e[[object_map[[method]]]]
          cs_obj <- fit$cs

          if (is.null(cs_obj) || is.null(cs_obj$cs) || length(cs_obj$cs) == 0L) next
          if (is.null(cs_obj$cs_index) || length(cs_obj$cs_index) != length(cs_obj$cs)) {
            stop("Credible-set/component index mismatch in ", result_file,
                 " object ", method, call. = FALSE)
          }

          for (k in seq_along(cs_obj$cs)) {
            component <- as.integer(cs_obj$cs_index[[k]])
            cs_idx <- as.integer(cs_obj$cs[[k]])

            if (component < 1L || component > length(fit$alpha)) {
              stop("Invalid credible-set component index in ", result_file,
                   " object ", method, call. = FALSE)
            }
            if (!length(cs_idx) || any(cs_idx < 1L | cs_idx > nrow(z))) {
              stop("Invalid credible-set SNP indices in ", result_file,
                   " object ", method, call. = FALSE)
            }

            alpha_l <- as.matrix(fit$alpha[[component]])
            if (!identical(dim(alpha_l), c(nrow(z), 3L))) {
              stop("Unexpected alpha dimensions in ", result_file,
                   " object ", method, " component ", component, call. = FALSE)
            }

            q <- colSums(alpha_l)
            if (any(!is.finite(q)) || any(q < -1e-10) || abs(sum(q) - 1) > 1e-6) {
              stop("Invalid component configuration posterior in ", result_file,
                   " object ", method, " component ", component, call. = FALSE)
            }
            q <- pmax(q, 0)
            q <- q / sum(q)

            causal_in_cs <- cs_idx[z$Signal[cs_idx] != 0]
            causal_signals <- z$Signal[causal_in_cs]
            unique_signals <- sort(unique(causal_signals))
            truth_status <- if (!length(causal_in_cs)) {
              "No causal SNP"
            } else if (length(causal_in_cs) == 1L) {
              "Exactly one causal SNP"
            } else if (length(unique_signals) == 1L) {
              "Multiple causal SNPs, same configuration"
            } else {
              "Multiple causal SNPs, mixed configurations"
            }

            truth_class <- if (length(causal_in_cs) == 1L) {
              as.character(signal_to_class(causal_signals))
            } else {
              NA_character_
            }

            purity_row <- if (!is.null(cs_obj$purity) && nrow(cs_obj$purity) >= k) {
              as.numeric(cs_obj$purity[k, c("min.abs.corr", "mean.abs.corr", "median.abs.corr")])
            } else {
              rep(NA_real_, 3L)
            }

            row_counter <- row_counter + 1L
            cs_rows[[row_counter]] <- data.table(
              dimension = spec$dimension,
              method = method,
              causal_num = causal_num,
              h2 = h2,
              locus = locus,
              component = component,
              cs_name = if (!is.null(names(cs_obj$cs))) names(cs_obj$cs)[k] else paste0("L", component),
              cs_size = length(cs_idx),
              claimed_coverage = if (!is.null(cs_obj$coverage)) as.numeric(cs_obj$coverage[k]) else NA_real_,
              min_abs_corr = purity_row[1L],
              mean_abs_corr = purity_row[2L],
              median_abs_corr = purity_row[3L],
              n_causal_in_cs = length(causal_in_cs),
              n_truth_configurations = length(unique_signals),
              truth_status = truth_status,
              causal_snp_index = if (length(causal_in_cs) == 1L) causal_in_cs else NA_integer_,
              causal_snp = if (length(causal_in_cs) == 1L) as.character(z$RSID[causal_in_cs]) else NA_character_,
              truth_class = truth_class,
              q_eur = q[1L],
              q_afr = q[2L],
              q_shared = q[3L],
              q_max = max(q),
              pred_class = config_levels[which.max(q)]
            )
          }
        }

        rm(e, z)
        if (cell_counter %% 100L == 0L) {
          message("Processed ", cell_counter, "/1200 dimension-specific cells")
          invisible(gc(FALSE))
        }
      }
    }
  }
}

cs_all <- rbindlist(cs_rows, use.names = TRUE)
cs_primary <- cs_all[truth_status == "Exactly one causal SNP"]

# One row is one retained credible set (single-effect component).  A simulation
# evaluation may therefore contribute several rows.  This explicit key check
# prevents accidental region-level collapsing or repeated extraction of a CS.
cs_key <- c("dimension", "method", "causal_num", "h2", "locus", "component")
if (anyDuplicated(cs_all[, ..cs_key])) {
  stop("Credible-set extraction produced duplicate component rows.", call. = FALSE)
}

cell_grid <- CJ(
  dimension = c("d = 10", "d = 187"),
  method = names(object_map),
  causal_num = 1:3,
  h2 = 1:2,
  locus = 1:100,
  unique = TRUE
)
cs_count_per_evaluation <- cs_all[, .(
  n_credible_sets = .N,
  n_unique_components = uniqueN(component)
), by = .(dimension, method, causal_num, h2, locus)]
cs_count_per_evaluation <- merge(
  cell_grid,
  cs_count_per_evaluation,
  by = c("dimension", "method", "causal_num", "h2", "locus"),
  all.x = TRUE
)
cs_count_per_evaluation[is.na(n_credible_sets), `:=`(
  n_credible_sets = 0L,
  n_unique_components = 0L
)]
if (any(cs_count_per_evaluation$n_credible_sets !=
        cs_count_per_evaluation$n_unique_components)) {
  stop("A simulation evaluation contains repeated credible-set components.", call. = FALSE)
}

# Write the extraction table before downstream summaries so a completed raw
# pass remains available for diagnosis if a later summary step fails.
fwrite(
  cs_all,
  file.path(result_out_dir, "credible_set_configuration_predictions.tsv.gz"),
  sep = "\t"
)
fwrite(
  cs_count_per_evaluation,
  file.path(result_out_dir, "credible_set_count_per_evaluation.tsv"),
  sep = "\t"
)

if (!nrow(cs_primary)) stop("No evaluable credible sets were found.", call. = FALSE)
if (any(!cs_primary$truth_class %in% config_levels)) {
  stop("Invalid primary credible-set truth labels.", call. = FALSE)
}

cs_metrics <- function(d) {
  truth <- factor(d$truth_class, levels = config_levels)
  pred <- factor(d$pred_class, levels = config_levels)
  cm <- table(truth, pred)
  recall <- diag(cm) / rowSums(cm)
  precision <- diag(cm) / colSums(cm)
  f1 <- 2 * precision * recall / (precision + recall)
  y_shared <- as.integer(truth == "Shared")
  pred_shared <- d$q_shared >= 0.5

  c(
    n_credible_sets = nrow(d),
    accuracy_3class = mean(truth == pred),
    balanced_accuracy_3class = mean(recall, na.rm = TRUE),
    macro_f1 = mean(f1, na.rm = TRUE),
    shared_auroc = roc_auc(y_shared, d$q_shared),
    shared_average_precision = average_precision(y_shared, d$q_shared),
    shared_brier = mean((d$q_shared - y_shared)^2),
    binary_balanced_accuracy = mean(c(
      mean(pred_shared[y_shared == 1L]),
      mean(!pred_shared[y_shared == 0L])
    )),
    shared_sensitivity = mean(pred_shared[y_shared == 1L]),
    ancestry_specific_specificity = mean(!pred_shared[y_shared == 0L]),
    recall_EUR = unname(recall["EUR-specific"]),
    recall_AFR = unname(recall["AFR-specific"]),
    recall_Shared = unname(recall["Shared"]),
    precision_EUR = unname(precision["EUR-specific"]),
    precision_AFR = unname(precision["AFR-specific"]),
    precision_Shared = unname(precision["Shared"])
  )
}

metric_long <- function(d) {
  v <- cs_metrics(d)
  data.table(metric = names(v), estimate = as.numeric(v))
}

point_metrics <- cs_primary[, metric_long(.SD), by = .(dimension, method)]

# Locus-cluster bootstrap.  Each sampled locus retains every causal-number and
# heritability cell, and the same sampled loci are used across dimensions and
# methods.
set.seed(20260903)
B <- 1000L
loci <- sort(unique(cs_primary$locus))
bootstrap_rows <- vector("list", B * 4L)
ii <- 0L

for (b in seq_len(B)) {
  sampled <- sample(loci, length(loci), replace = TRUE)
  boot <- rbindlist(lapply(sampled, function(ll) cs_primary[locus == ll]), use.names = TRUE)

  for (dd in c("d = 10", "d = 187")) {
    for (mm in c("INFUSE", "MESuSiE")) {
      ii <- ii + 1L
      vals <- cs_metrics(boot[dimension == dd & method == mm])
      bootstrap_rows[[ii]] <- data.table(
        bootstrap = b,
        dimension = dd,
        method = mm,
        metric = names(vals),
        estimate = as.numeric(vals)
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
metric_summary <- merge(
  point_metrics,
  metric_ci,
  by = c("dimension", "method", "metric"),
  all.x = TRUE
)

confusion_counts <- cs_primary[, .N, by = .(dimension, method, truth_class, pred_class)]
confusion_grid <- CJ(
  dimension = sort(unique(cs_primary$dimension)),
  method = sort(unique(cs_primary$method)),
  truth_class = config_levels,
  pred_class = config_levels,
  unique = TRUE
)
confusion <- merge(
  confusion_grid,
  confusion_counts,
  by = c("dimension", "method", "truth_class", "pred_class"),
  all.x = TRUE
)
confusion[is.na(N), N := 0L]
confusion[, proportion := N / sum(N), by = .(dimension, method, truth_class)]

calibration <- cs_primary[, {
  breaks <- seq(0, 1, by = 0.1)
  bin <- cut(q_shared, breaks = breaks, include.lowest = TRUE, right = TRUE)
  d <- data.table(
    bin = bin,
    q_shared = q_shared,
    observed = as.integer(truth_class == "Shared")
  )
  d[, .(
    n = .N,
    mean_predicted = mean(q_shared),
    observed_shared = mean(observed)
  ), by = bin]
}, by = .(dimension, method)]

binary_thresholds <- seq(1, 0, by = -0.01)
binary_curves <- cs_primary[, {
  y <- truth_class == "Shared"
  rbindlist(lapply(binary_thresholds, function(tt) {
    called <- q_shared >= tt
    tp <- sum(called & y)
    fp <- sum(called & !y)
    fn <- sum(!called & y)
    tn <- sum(!called & !y)
    data.table(
      threshold = tt,
      sensitivity = if ((tp + fn) > 0) tp / (tp + fn) else NA_real_,
      false_positive_rate = if ((fp + tn) > 0) fp / (fp + tn) else NA_real_,
      precision = if ((tp + fp) > 0) tp / (tp + fp) else 1,
      recall = if ((tp + fn) > 0) tp / (tp + fn) else NA_real_
    )
  }))
}, by = .(dimension, method)]

audit <- cs_all[, .(
  n_credible_sets = .N,
  median_cs_size = as.numeric(median(cs_size)),
  n_no_causal = sum(truth_status == "No causal SNP"),
  n_exactly_one_causal = sum(truth_status == "Exactly one causal SNP"),
  n_multiple_same_configuration = sum(
    truth_status == "Multiple causal SNPs, same configuration"
  ),
  n_multiple_mixed_configurations = sum(
    truth_status == "Multiple causal SNPs, mixed configurations"
  ),
  fraction_primary_evaluable = mean(truth_status == "Exactly one causal SNP")
), by = .(dimension, method)]

posterior_diagnostics <- cs_all[, .(
  n_credible_sets = .N,
  min_posterior_sum = min(q_eur + q_afr + q_shared),
  max_posterior_sum = max(q_eur + q_afr + q_shared),
  min_probability = min(q_eur, q_afr, q_shared),
  max_probability = max(q_eur, q_afr, q_shared)
), by = .(dimension, method)]

fwrite(audit, file.path(result_out_dir, "credible_set_configuration_audit.tsv"), sep = "\t")
fwrite(metric_summary, file.path(result_out_dir, "credible_set_classification_metric_summary.tsv"), sep = "\t")
fwrite(confusion, file.path(result_out_dir, "credible_set_three_class_confusion.tsv"), sep = "\t")
fwrite(calibration, file.path(result_out_dir, "credible_set_shared_probability_calibration.tsv"), sep = "\t")
fwrite(posterior_diagnostics, file.path(result_out_dir, "credible_set_posterior_diagnostics.tsv"), sep = "\t")

saveRDS(list(
  cs_all = cs_all,
  cs_primary = cs_primary,
  cs_count_per_evaluation = cs_count_per_evaluation,
  audit = audit,
  metric_summary = metric_summary,
  bootstrap = bootstrap,
  confusion = confusion,
  calibration = calibration,
  binary_curves = binary_curves,
  posterior_diagnostics = posterior_diagnostics,
  primary_truth_rule = "Retained 95% credible sets containing exactly one simulated causal SNP"
), file.path(result_out_dir, "credible_set_configuration_classification_analysis.rds"), compress = "xz")

message(
  "Credible-set configuration analysis complete. Total retained sets: ",
  nrow(cs_all), "; primary evaluable sets: ", nrow(cs_primary)
)
