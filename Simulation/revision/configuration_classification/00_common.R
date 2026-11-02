.infuse_repo <- Sys.getenv("INFUSE_ANALYSIS_ROOT", "/path/to/INFUSE_analysis")
.infuse_project <- Sys.getenv("INFUSE_PROJECT", "/path/to/baseline_simulation")
.infuse_work <- Sys.getenv("INFUSE_REVISION", "/path/to/revision_output")
options(stringsAsFactors = FALSE)

project_root <- .infuse_project
module_dir <- file.path(.infuse_work, "configuration_classification_w2")
result_out_dir <- file.path(module_dir, "results")
figure_out_dir <- file.path(.infuse_work, "figures", "configuration_classification_w2")
dir.create(result_out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_out_dir, recursive = TRUE, showWarnings = FALSE)

config_levels <- c("EUR-specific", "AFR-specific", "Shared")
config_short <- c("EUR", "AFR", "Shared")

source(file.path(.infuse_repo, "Simulation/revision/common/baseline_sources.R"))
dimension_spec <- figure2_baseline_spec

signal_to_class <- function(x) {
  out <- rep(NA_character_, length(x))
  out[x == 1] <- config_levels[1]
  out[x == 2] <- config_levels[2]
  out[x == 3] <- config_levels[3]
  factor(out, levels = config_levels)
}

safe_divide <- function(num, den, zero = NA_real_) {
  out <- rep(zero, length(num))
  ok <- is.finite(num) & is.finite(den) & den != 0
  out[ok] <- num[ok] / den[ok]
  out
}

roc_auc <- function(y, score) {
  y <- as.integer(y)
  ok <- is.finite(score) & y %in% c(0L, 1L)
  y <- y[ok]
  score <- score[ok]
  n1 <- sum(y == 1L)
  n0 <- sum(y == 0L)
  if (!n1 || !n0) return(NA_real_)
  ranks <- rank(score, ties.method = "average")
  (sum(ranks[y == 1L]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

average_precision <- function(y, score) {
  y <- as.integer(y)
  ok <- is.finite(score) & y %in% c(0L, 1L)
  y <- y[ok]
  score <- score[ok]
  n1 <- sum(y == 1L)
  if (!n1) return(NA_real_)
  ord <- order(score, decreasing = TRUE)
  y <- y[ord]
  mean(cumsum(y)[y == 1L] / which(y == 1L))
}

classification_metrics <- function(d) {
  truth <- factor(d$truth_class, levels = config_levels)
  pred <- factor(d$pred_class, levels = config_levels)
  cm <- table(truth, pred)
  recall <- diag(cm) / rowSums(cm)
  precision <- diag(cm) / colSums(cm)
  f1 <- 2 * precision * recall / (precision + recall)
  y_shared <- as.integer(truth == "Shared")
  pred_shared <- d$q_shared >= 0.5

  c(
    n_causal = nrow(d),
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
