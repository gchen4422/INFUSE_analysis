#!/usr/bin/env Rscript
.infuse_repo <- Sys.getenv("INFUSE_ANALYSIS_ROOT", "/path/to/INFUSE_analysis")

suppressPackageStartupMessages(library(data.table))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L || !args[[1L]] %in% c("d10", "d187")) {
  stop("Usage: Rscript 01_extract_configuration_predictions.R <d10|d187>", call. = FALSE)
}
d_label <- args[[1L]]
file_arg <- commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))][1L]
script_dir <- dirname(normalizePath(sub("^--file=", "", file_arg)))
source(file.path(.infuse_repo, "Simulation/revision/configuration_classification/00_common.R"))

spec <- dimension_spec(d_label)
thresholds <- c(1.000001, seq(1, 0.01, by = -0.01))
causal_rows <- list()
roc_rows <- list()
cell_rows <- list()
counter <- 0L

for (causal_num in 1:3) {
  for (h2 in 1:2) {
    for (locus in 1:100) {
      counter <- counter + 1L
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
      required <- c("RSID", "Signal")
      if (!all(required %in% names(z))) stop("Truth columns missing in ", data_file, call. = FALSE)

      e <- new.env(parent = emptyenv())
      load(result_file, envir = e)
      object_map <- c(MESuSiE = "MESuSiE_res", INFUSE = "MESuSiE_res_mlk_annot_subset")
      if (!all(unname(object_map) %in% ls(e))) {
        stop("Required objects missing in ", result_file, call. = FALSE)
      }

      truth_class <- signal_to_class(z$Signal)
      causal_idx <- which(z$Signal != 0)
      cell_rows[[length(cell_rows) + 1L]] <- data.frame(
        dimension = spec$dimension, causal_num = causal_num, h2 = h2,
        locus = locus, n_snp = nrow(z), n_causal = length(causal_idx),
        result_file = result_file, stringsAsFactors = FALSE
      )

      for (method in names(object_map)) {
        fit <- e[[object_map[[method]]]]
        pip <- as.numeric(fit$pip)
        pc <- as.matrix(fit$pip_config)
        if (length(pip) != nrow(z) || !identical(dim(pc), c(nrow(z), 3L))) {
          stop("PIP dimensions do not match truth in ", result_file, call. = FALSE)
        }
        if (any(!is.finite(pip)) || any(!is.finite(pc)) || any(pc < -1e-12)) {
          stop("Invalid PIPs in ", result_file, " object ", method, call. = FALSE)
        }
        support_sum <- rowSums(pc)
        q <- matrix(1 / 3, nrow = nrow(pc), ncol = 3L)
        ok <- is.finite(support_sum) & support_sum > 0
        q[ok, ] <- pc[ok, , drop = FALSE] / support_sum[ok]
        colnames(q) <- c("q_eur", "q_afr", "q_shared")
        pred <- config_levels[max.col(q, ties.method = "first")]

        causal_rows[[length(causal_rows) + 1L]] <- data.frame(
          dimension = spec$dimension, method = method,
          causal_num = causal_num, h2 = h2, locus = locus,
          snp_index = causal_idx, SNP = as.character(z$RSID[causal_idx]),
          Signal = z$Signal[causal_idx],
          truth_class = as.character(truth_class[causal_idx]),
          pip_either = pip[causal_idx],
          pip_eur = pc[causal_idx, 1], pip_afr = pc[causal_idx, 2],
          pip_shared = pc[causal_idx, 3],
          support_sum = support_sum[causal_idx],
          q_eur = q[causal_idx, 1], q_afr = q[causal_idx, 2],
          q_shared = q[causal_idx, 3],
          q_max = apply(q[causal_idx, , drop = FALSE], 1, max),
          pred_class = pred[causal_idx],
          stringsAsFactors = FALSE
        )

        for (cc in seq_along(config_levels)) {
          is_true <- z$Signal == cc
          score <- pc[, cc]
          called <- outer(score, thresholds, `>`)
          tp <- colSums(called & is_true)
          fp <- colSums(called & !is_true)
          roc_rows[[length(roc_rows) + 1L]] <- data.frame(
            dimension = spec$dimension, method = method,
            causal_num = causal_num, h2 = h2, locus = locus,
            configuration = config_levels[cc], cutoff = thresholds,
            TP = tp, FP = fp, n_true = sum(is_true), n_variant = nrow(z),
            stringsAsFactors = FALSE
          )
        }
      }
      rm(e, z)
      if (counter %% 25L == 0L) {
        message(spec$dimension, ": processed ", counter, "/600 cells")
        invisible(gc(FALSE))
      }
    }
  }
}

causal <- rbindlist(causal_rows, use.names = TRUE)
roc_counts <- rbindlist(roc_rows, use.names = TRUE)
cells <- rbindlist(cell_rows, use.names = TRUE)

if (nrow(cells) != 600L || uniqueN(cells[, .(causal_num, h2, locus)]) != 600L) {
  stop("Extraction did not produce the complete 600-cell grid.", call. = FALSE)
}
expected_causal_per_method <- 2L * 100L * sum(c(2L, 3L, 5L))
if (nrow(causal) != expected_causal_per_method * 2L) {
  stop("Unexpected causal-row count: ", nrow(causal), call. = FALSE)
}

saveRDS(
  list(
    dimension = spec$dimension,
    source_spec = spec,
    config_order = config_levels,
    causal = causal,
    roc_counts = roc_counts,
    cells = cells,
    extracted_at = Sys.time()
  ),
  file.path(result_out_dir, paste0("configuration_predictions_", d_label, ".rds")),
  compress = "xz"
)
message("Completed ", spec$dimension, " extraction.")
