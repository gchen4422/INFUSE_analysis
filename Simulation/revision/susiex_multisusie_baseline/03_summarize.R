#!/usr/bin/env Rscript
.infuse_repo <- Sys.getenv("INFUSE_ANALYSIS_ROOT", "/path/to/INFUSE_analysis")

this_file <- normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1L]),
                           mustWork = TRUE)
source(file.path(.infuse_repo, "Simulation/revision/susiex_multisusie_baseline/00_common.R"))

cutoffs <- seq(1, 0.01, by = -0.01)
method_keys <- c(susiex = "SuSiEx", multisusie = "MultiSuSiE")

compute_curve <- function(pip, truth) {
  do.call(rbind, lapply(cutoffs, function(cutoff) {
    selected <- pip > cutoff
    tp <- sum(selected & truth)
    fp <- sum(selected) - tp
    power <- if (sum(truth)) tp / sum(truth) else 0
    fdr <- if (tp + fp) fp / (tp + fp) else 0
    c(Cutoff = cutoff, Power = power, FDR = fdr)
  }))
}

cell_metrics <- list()
curve_parts <- list()
variant_parts <- list()

for (method_key in names(method_keys)) {
  method <- unname(method_keys[[method_key]])
  for (task_id in seq_len(N_CELLS)) {
    cell <- decode_task(task_id)
    z <- read_z(cell, "d10")
    f <- cell_output(method_key, cell)
    obj <- safe_read_method_object(f, z, method)
    if (is.null(obj)) stop("Missing or invalid ", method, " output for task ", task_id, call. = FALSE)

    selected <- sort(unique(as.integer(unlist(obj$cs, use.names = FALSE))))
    truth_any <- z$Signal != 0L
    identified <- if (length(selected)) sum(truth_any[selected]) else 0L
    set_size <- length(selected)
    cell_metrics[[length(cell_metrics) + 1L]] <- data.frame(
      method = method, task_id = task_id, locus = cell$locus,
      causal_index = cell$causal_index, causal_count = cell$causal_count,
      h2_index = cell$h2_index, h2 = cell$h2,
      set_power = identified / sum(truth_any),
      set_fdr = if (set_size) (set_size - identified) / set_size else 0,
      set_size = set_size, n_cs = length(obj$cs),
      runtime_minutes = obj$runtime_minutes
    )

    endpoint_truth <- list(
      Either = truth_any,
      Shared = z$Signal == 3L,
      `Ancestry-specific|White British` = z$Signal == 1L,
      `Ancestry-specific|Black British` = z$Signal == 2L
    )
    for (endpoint_name in names(endpoint_truth)) {
      endpoint_split <- strsplit(endpoint_name, "\\|", fixed = FALSE)[[1L]]
      endpoint <- endpoint_split[1L]
      ancestry <- if (length(endpoint_split) == 2L) endpoint_split[2L] else NA_character_
      truth <- endpoint_truth[[endpoint_name]]
      # Match the established baseline summaries: shared curves condition on a
      # locus containing a shared effect; ancestry-specific curves retain every locus.
      if (endpoint == "Shared" && !any(truth)) next
      cv <- as.data.frame(compute_curve(obj$pip, truth))
      cv$method <- method
      cv$task_id <- task_id
      cv$locus <- cell$locus
      cv$causal_index <- cell$causal_index
      cv$causal_count <- cell$causal_count
      cv$h2_index <- cell$h2_index
      cv$h2 <- cell$h2
      cv$endpoint <- endpoint
      cv$ancestry <- ancestry
      curve_parts[[length(curve_parts) + 1L]] <- cv
    }

    variant_parts[[length(variant_parts) + 1L]] <- data.frame(
      method = method, task_id = task_id, locus = cell$locus,
      causal_index = cell$causal_index, causal_count = cell$causal_count,
      h2_index = cell$h2_index, h2 = cell$h2,
      RSID = z$RSID, Signal = z$Signal, PIP = obj$pip
    )
  }
}

new_cell_metrics <- rbindlist(cell_metrics)
new_cell_curves <- rbindlist(curve_parts)
new_variant_pips <- rbindlist(variant_parts)
fwrite(new_cell_metrics, file.path(RESULT_DIR, "new_method_cell_metrics.tsv"), sep = "\t")
fwrite(new_cell_curves, file.path(RESULT_DIR, "new_method_cell_threshold_curves.tsv.gz"), sep = "\t")
saveRDS(new_variant_pips, file.path(RESULT_DIR, "new_method_variant_pips.rds"), compress = "xz")

new_mean_curves <- new_cell_curves[, .(
  Power = mean(Power), FDR = mean(FDR), n_loci = .N
), by = .(method, causal_index, causal_count, h2_index, h2,
          endpoint, ancestry, Cutoff)]

load_existing <- function(dimension) {
  setting <- if (dimension == "d = 10") "d10" else "d187"
  spec <- figure2_baseline_spec(setting)
  f <- spec$legacy_cache
  e <- new.env(parent = emptyenv())
  load(f, envir = e)
  e
}

decode_h2 <- function(x) ifelse(grepl("2%\\*%10", as.character(x)), 2L, 1L)
decode_causal <- function(x) {
  ans <- as.integer(sub(".*==[[:space:]]*([235]).*", "\\1", as.character(x)))
  match(ans, CAUSAL_COUNTS)
}
display_method <- function(method, dimension) {
  method <- as.character(method)
  out <- method
  out[method == "Paintor"] <- "PAINTOR"
  out[method == "Paintor+annot"] <- if (dimension == "d = 10") "PAINTOR+" else "PAINTOR+H"
  out[method == "MAFPC"] <- if (dimension == "d = 10") "INFUSE" else "INFUSE-H"
  out
}

existing_sets <- list()
existing_curves <- list()
for (dimension in c("d = 10", "d = 187")) {
  e <- load_existing(dimension)
  s <- as.data.table(e$all_Set_data_dataframe)
  s[, h2_index := decode_h2(h2)]
  s[, causal_index := decode_causal(causal_num)]
  s[, `:=`(h2 = H2_VALUES[h2_index], causal_count = CAUSAL_COUNTS[causal_index])]
  s[, method := display_method(Method, dimension)]
  s[, locus := seq_len(.N), by = .(method, h2_index, causal_index)]
  if (any(s[, .N, by = .(method, h2_index, causal_index)]$N != N_LOCI)) {
    stop("Existing credible-set table does not contain 100 loci per condition", call. = FALSE)
  }
  existing_sets[[dimension]] <- s[, .(
    method, dimension, locus, causal_index, causal_count, h2_index, h2,
    set_power = as.numeric(Power), set_fdr = as.numeric(FDR), set_size = as.numeric(Size)
  )]

  add_existing_curve <- function(tab, endpoint) {
    x <- as.data.table(tab)
    ancestry <- rep(NA_character_, nrow(x))
    base_method <- as.character(x$Method)
    if (endpoint == "Ancestry-specific") {
      ancestry[grepl(" BB$", base_method)] <- "Black British"
      ancestry[grepl(" WB$", base_method)] <- "White British"
      base_method <- sub(" (BB|WB)$", "", base_method)
    }
    x[, h2_index := decode_h2(h2)]
    x[, causal_index := decode_causal(causal_num)]
    x[, `:=`(h2_value = H2_VALUES[h2_index], causal_count = CAUSAL_COUNTS[causal_index],
             method_display = display_method(base_method, dimension),
             endpoint_value = endpoint, ancestry_value = ancestry)]
    x[, .(method = method_display, dimension = dimension,
           causal_index, causal_count, h2_index, h2 = h2_value,
           endpoint = endpoint_value, ancestry = ancestry_value,
           Cutoff = as.numeric(as.character(Cutoff)),
           Power = as.numeric(Power), FDR = as.numeric(FDR))]
  }
  existing_curves[[paste0(dimension, "_either")]] <- add_existing_curve(e$either_all_ROC_data_dataframe, "Either")
  existing_curves[[paste0(dimension, "_shared")]] <- add_existing_curve(e$shared_all_ROC_data_dataframe, "Shared")
  existing_curves[[paste0(dimension, "_ancestry")]] <- add_existing_curve(e$ancestry_all_ROC_data_dataframe, "Ancestry-specific")
}

new_sets_both <- rbindlist(lapply(c("d = 10", "d = 187"), function(dimension) {
  x <- copy(new_cell_metrics)
  x[, dimension := dimension]
  x[, .(method, dimension, locus, causal_index, causal_count, h2_index, h2,
         set_power, set_fdr, set_size)]
}))
combined_sets <- rbindlist(c(existing_sets, list(new_sets_both)), use.names = TRUE)

new_curves_both <- rbindlist(lapply(c("d = 10", "d = 187"), function(dimension) {
  x <- copy(new_mean_curves)
  x[, dimension := dimension]
  x[, .(method, dimension, causal_index, causal_count, h2_index, h2,
         endpoint, ancestry, Cutoff, Power, FDR)]
}))
combined_curves <- rbindlist(c(existing_curves, list(new_curves_both)), use.names = TRUE)

target_fdr <- rbindlist(lapply(c(0.05, 0.10), function(target) {
  ans <- combined_curves[, {
    delta <- abs(FDR - target)
    candidate <- .SD[delta == min(delta)]
    candidate[which.max(Power)][1L]
  }, by = .(method, dimension, causal_index, causal_count, h2_index, h2,
            endpoint, ancestry)]
  ans[, target_FDR := target]
  ans
}))

set_summary <- combined_sets[, .(
  median_set_size = median(set_size), mean_set_size = mean(set_size),
  mean_set_power = mean(set_power), mean_set_fdr = mean(set_fdr), n_loci = .N
), by = .(method, dimension, causal_index, causal_count, h2_index, h2)]
overall_set_summary <- combined_sets[, .(
  median_set_size = median(set_size), mean_set_size = mean(set_size),
  mean_set_power = mean(set_power), mean_set_fdr = mean(set_fdr), n_evaluations = .N
), by = .(method, dimension)]
overall_target_fdr <- target_fdr[, .(
  mean_power = mean(Power), mean_realized_fdr = mean(FDR), n_conditions = .N
), by = .(method, dimension, target_FDR, endpoint, ancestry)]

fwrite(combined_sets, file.path(RESULT_DIR, "combined_credible_set_cell_metrics.tsv.gz"), sep = "\t")
fwrite(set_summary, file.path(RESULT_DIR, "combined_credible_set_summary_by_condition.tsv"), sep = "\t")
fwrite(overall_set_summary, file.path(RESULT_DIR, "combined_credible_set_summary_overall.tsv"), sep = "\t")
fwrite(combined_curves, file.path(RESULT_DIR, "combined_mean_threshold_curves.tsv.gz"), sep = "\t")
fwrite(target_fdr, file.path(RESULT_DIR, "combined_target_fdr_by_condition.tsv"), sep = "\t")
fwrite(target_fdr[target_FDR == 0.05],
       file.path(RESULT_DIR, "combined_target_fdr_0.05_by_condition.tsv"), sep = "\t")
fwrite(target_fdr[target_FDR == 0.10],
       file.path(RESULT_DIR, "combined_target_fdr_0.10_by_condition.tsv"), sep = "\t")
fwrite(overall_target_fdr[target_FDR == 0.05],
       file.path(RESULT_DIR, "combined_target_fdr_0.05_overall.tsv"), sep = "\t")
fwrite(overall_target_fdr[target_FDR == 0.10],
       file.path(RESULT_DIR, "combined_target_fdr_0.10_overall.tsv"), sep = "\t")
fwrite(new_cell_metrics[, .(
  median_minutes = median(runtime_minutes), mean_minutes = mean(runtime_minutes),
  max_minutes = max(runtime_minutes), n = .N
), by = method], file.path(RESULT_DIR, "new_method_runtime_summary.tsv"), sep = "\t")
saveRDS(list(combined_sets = combined_sets, set_summary = set_summary,
             combined_curves = combined_curves, target_fdr = target_fdr,
             new_cell_metrics = new_cell_metrics),
        file.path(RESULT_DIR, "baseline_susiex_multisusie_benchmark_summary.rds"),
        compress = "xz")
cat("Summary complete.\n")
print(overall_set_summary)
print(overall_target_fdr)
