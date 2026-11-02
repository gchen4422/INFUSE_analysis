#!/usr/bin/env Rscript
.infuse_repo <- Sys.getenv("INFUSE_ANALYSIS_ROOT", "/path/to/INFUSE_analysis")

# Add SuSiEx and MultiSuSiE to the functional-enrichment MSE analysis used in
# the established baseline Figure 2. Existing methods are never recomputed by
# this script; 04_plot.R reads their values directly from the saved Figure 2
# ggplot objects.

this_file <- normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1L]),
                           mustWork = TRUE)
source(file.path(.infuse_repo, "Simulation/revision/susiex_multisusie_baseline/00_common.R"))

pip_cutoff <- 0.5
annotation_root <- file.path(SIMULATION_ROOT, "risk_loci_annot")
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

pip_file <- file.path(RESULT_DIR, "new_method_variant_pips.rds")
if (!file.exists(pip_file)) stop("Run 03_summarize.R first", call. = FALSE)
pips <- as.data.table(readRDS(pip_file))

required_pip_columns <- c("method", "locus", "causal_index", "causal_count",
                          "h2_index", "RSID", "Signal", "PIP")
if (!all(required_pip_columns %in% names(pips))) {
  stop("new_method_variant_pips.rds is missing required columns", call. = FALSE)
}
if (!setequal(unique(pips$method), c("SuSiEx", "MultiSuSiE"))) {
  stop("Unexpected methods in new_method_variant_pips.rds", call. = FALSE)
}

annotation_parts <- lapply(seq_len(N_LOCI), function(locus) {
  f <- file.path(annotation_root, sprintf("loci_%d.annot", locus))
  if (!file.exists(f)) stop("Missing annotation file: ", f, call. = FALSE)
  a <- fread(f, select = c("CHR", "BP", annotation_columns))
  a[, `:=`(locus = locus, RSID = paste0(CHR, ":", BP))]
  a <- a[, c("locus", "RSID", annotation_columns), with = FALSE]
  if (anyDuplicated(a$RSID)) stop("Duplicated annotation RSIDs at locus ", locus, call. = FALSE)
  for (nm in annotation_columns) {
    if (anyNA(a[[nm]]) || any(!a[[nm]] %in% c(0, 1))) {
      stop("Signal annotation is not finite binary: ", nm, ", locus ", locus,
           call. = FALSE)
    }
  }
  a
})
annotations <- rbindlist(annotation_parts)

setkey(annotations, locus, RSID)
setkey(pips, locus, RSID)
joined <- annotations[pips]
if (nrow(joined) != nrow(pips) || anyNA(joined[[annotation_columns[[1L]]]])) {
  stop("PIP-to-annotation alignment was not complete and one-to-one", call. = FALSE)
}

safe_ratio <- function(x, y) {
  if (!is.finite(x) || !is.finite(y) || y <= 0) return(NA_real_)
  x / y
}

result_rows <- list()
row_i <- 0L
for (method_value in c("SuSiEx", "MultiSuSiE")) {
  for (causal_index_value in seq_len(N_CAUSAL_INDEX)) {
    for (h2_index_value in seq_len(N_H2)) {
      d <- joined[method == method_value & causal_index == causal_index_value &
                    h2_index == h2_index_value]
      if (!nrow(d)) stop("Missing method/condition rows", call. = FALSE)
      observed_locus_counts <- d[, .N, by = locus][order(locus)]
      expected_locus_counts <- annotations[, .N, by = locus][order(locus)]
      if (nrow(d) != nrow(annotations) || uniqueN(d$locus) != N_LOCI ||
          !identical(observed_locus_counts$locus, expected_locus_counts$locus) ||
          !identical(observed_locus_counts$N, expected_locus_counts$N)) {
        stop("Incomplete method/condition variant set for ", method_value, call. = FALSE)
      }

      causal <- d$Signal != 0L
      selected <- d$PIP > pip_cutoff
      n_total <- nrow(d)
      n_causal <- sum(causal)
      n_selected <- sum(selected)
      if (!n_causal || !n_selected || n_selected >= n_total) {
        stop("Degenerate truth/selection counts for ", method_value,
             ", causal index ", causal_index_value, ", h2 index ", h2_index_value,
             call. = FALSE)
      }

      for (annotation_column in annotation_columns) {
        values <- d[[annotation_column]]
        total_annotation_count <- sum(values)
        causal_annotation_count <- sum(values[causal])
        selected_annotation_count <- sum(values[selected])

        true_enrichment_exact <- safe_ratio(
          causal_annotation_count / n_causal,
          (total_annotation_count - causal_annotation_count) / (n_total - n_causal)
        )
        # The established Figure 2 MSE joined against True_enrichment_long,
        # whose enrichment values were rounded to two decimal places.
        true_enrichment_figure2 <- round(true_enrichment_exact, 2)

        selected_frequency <- selected_annotation_count / n_selected

        # Exact legacy Figure 2 definition. The denominator follows
        # 9_1_enrichmen_results_pip_cutoff_either_pip_method.R so the two new
        # methods are directly comparable to the saved existing-method values.
        figure2_background_frequency <- total_annotation_count / (n_total - n_selected)
        estimated_enrichment_figure2 <- safe_ratio(
          selected_frequency, figure2_background_frequency
        )

        # Enrichment definition used by 09_compute_native_enrichment_mse.R:
        # the background is the actual complement of the selected set.
        complement_background_frequency <-
          (total_annotation_count - selected_annotation_count) / (n_total - n_selected)
        estimated_enrichment_complement <- safe_ratio(
          selected_frequency, complement_background_frequency
        )

        row_i <- row_i + 1L
        result_rows[[row_i]] <- data.frame(
          method = method_value,
          causal = causal_index_value,
          causal_count = CAUSAL_COUNTS[causal_index_value],
          h2 = h2_index_value,
          h2_per_causal = H2_VALUES[h2_index_value],
          annotation = unname(annotation_labels[[annotation_column]]),
          annotation_column = annotation_column,
          pip_cutoff = pip_cutoff,
          n_variants = n_total,
          n_causal = n_causal,
          n_selected = n_selected,
          total_annotation_count = total_annotation_count,
          selected_annotation_count = selected_annotation_count,
          true_enrichment_exact = true_enrichment_exact,
          true_enrichment_figure2 = true_enrichment_figure2,
          estimated_enrichment_figure2 = estimated_enrichment_figure2,
          mse_figure2 = (estimated_enrichment_figure2 - true_enrichment_figure2)^2,
          estimated_enrichment_complement = estimated_enrichment_complement,
          mse_complement = (estimated_enrichment_complement - true_enrichment_exact)^2,
          stringsAsFactors = FALSE
        )
      }
    }
  }
}

detail <- rbindlist(result_rows)
expected_rows <- 2L * N_CAUSAL_INDEX * N_H2 * length(annotation_columns)
if (nrow(detail) != expected_rows ||
    any(!is.finite(detail$mse_figure2)) || any(!is.finite(detail$mse_complement))) {
  stop("Unexpected or non-finite enrichment MSE output", call. = FALSE)
}

fwrite(detail, file.path(RESULT_DIR, "new_method_enrichment_MSE_PIP_gt_0.5.tsv"), sep = "\t")
fwrite(detail[, .(
  mean_mse = mean(mse_figure2), median_mse = median(mse_figure2),
  n_parameter_cells = .N
), by = .(method, annotation)],
file.path(RESULT_DIR, "new_method_enrichment_MSE_figure2_summary.tsv"), sep = "\t")
fwrite(detail[, .(
  mean_mse = mean(mse_complement), median_mse = median(mse_complement),
  n_parameter_cells = .N
), by = .(method, annotation)],
file.path(RESULT_DIR, "new_method_enrichment_MSE_complement_sensitivity_summary.tsv"),
sep = "\t")

cat("New-method enrichment MSE complete: ", nrow(detail), " rows.\n", sep = "")
