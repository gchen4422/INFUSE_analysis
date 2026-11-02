suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

method_columns <- c(
  MESuSiE = "MESuSiE_PIP_Either",
  Classic_Annotations = "MESuSiE_PIP_Either_all_mlk"
)

method_colors <- c(
  MESuSiE = "#84c3b7",
  Classic_Annotations = "#71b8ed"
)

method_legend_labels <- c(
  MESuSiE = "MESuSiE",
  Classic_Annotations = "INFUSE"
)

load_raw_evo2_result <- function(path, trait) {
  if (!file.exists(path)) stop("Missing prepared result: ", path, call. = FALSE)
  env <- new.env(parent = emptyenv())
  loaded <- load(path, envir = env)
  if (!identical(loaded, "res_all_new")) {
    stop("Expected only res_all_new in ", path, call. = FALSE)
  }
  x <- as.data.table(env$res_all_new)
  needed <- c(
    "SNP", "Region", "CHR", "POS", "evo2_delta_score",
    unname(method_columns)
  )
  missing <- setdiff(needed, names(x))
  if (length(missing)) {
    stop("Prepared result is missing: ", paste(missing, collapse = ", "), call. = FALSE)
  }
  ans <- x[, c(
    list(
      Trait = trait,
      SNP = as.character(SNP),
      Region = as.character(Region),
      CHR = as.character(CHR),
      POS = as.integer(POS),
      evo2_delta_score = evo2_delta_score
    ),
    mget(unname(method_columns))
  )]
  rm(env, x)
  gc(verbose = FALSE)
  ans
}

# Assign one bootstrap cluster to every connected set of overlapping genomic
# regions. This is used after pooling related traits so that, for example, an
# SBP region and a DBP region covering the same locus are resampled together.
assign_overlapping_region_clusters <- function(data,
                                               source_trait_col = "Source_Trait") {
  if (!source_trait_col %in% names(data)) {
    stop("Missing source-trait column: ", source_trait_col, call. = FALSE)
  }
  needed <- c("Region", "CHR", "POS", source_trait_col)
  missing <- setdiff(needed, names(data))
  if (length(missing)) {
    stop("Cannot construct genomic clusters; missing: ",
         paste(missing, collapse = ", "), call. = FALSE)
  }

  d <- copy(as.data.table(data))
  d[, Source_Region := as.character(Region)]
  d[, Source_Trait := as.character(get(source_trait_col))]
  if (anyNA(d$CHR) || anyNA(d$POS)) {
    stop("CHR and POS must be observed for every variant.", call. = FALSE)
  }

  region_map <- d[, {
    chromosome <- unique(CHR)
    if (length(chromosome) != 1L) {
      stop("A source region spans multiple chromosomes.", call. = FALSE)
    }
    .(
      CHR = chromosome,
      region_start = min(POS),
      region_end = max(POS),
      n_variant_records = .N
    )
  }, by = .(Source_Trait, Source_Region)]

  setorder(region_map, CHR, region_start, region_end,
           Source_Trait, Source_Region)
  region_map[, cluster_index := {
    prior_running_end <- shift(cummax(as.numeric(region_end)), fill = -Inf)
    cumsum(region_start > prior_running_end)
  }, by = CHR]
  region_map[, `:=`(
    cluster_start = min(region_start),
    cluster_end = max(region_end)
  ), by = .(CHR, cluster_index)]
  region_map[, Bootstrap_Region := paste0(
    "chr", CHR, ":", cluster_start, "-", cluster_end
  )]

  d[
    region_map,
    on = .(Source_Trait, Source_Region),
    Bootstrap_Region := i.Bootstrap_Region
  ]
  if (anyNA(d$Bootstrap_Region)) {
    stop("Failed to assign a bootstrap region to every variant.", call. = FALSE)
  }
  d[, Region := Bootstrap_Region]
  d[, Bootstrap_Region := NULL]

  cluster_audit <- region_map[, .(
    cluster_start = unique(cluster_start),
    cluster_end = unique(cluster_end),
    n_source_regions = .N,
    source_traits = paste(sort(unique(Source_Trait)), collapse = ","),
    source_regions = paste(
      paste0(Source_Trait, "::", Source_Region),
      collapse = ","
    ),
    n_variant_records = sum(n_variant_records)
  ), by = .(Bootstrap_Region, CHR)]
  cluster_audit[, CHR_sort := suppressWarnings(as.integer(CHR))]
  setorder(cluster_audit, CHR_sort, CHR, cluster_start, cluster_end)
  cluster_audit[, CHR_sort := NULL]

  list(data = d, cluster_audit = cluster_audit)
}

cluster_bootstrap <- function(score, region, selected, n_boot = 1000L, seed = 1L) {
  z <- data.table(score = score, region = region, selected = selected)
  by_region <- z[, .(
    selected_sum = sum(score[selected]),
    selected_n = sum(selected),
    background_sum = sum(score[!selected]),
    background_n = sum(!selected)
  ), by = region]

  if (sum(by_region$selected_n) == 0L) {
    return(list(
      selected_mean = NA_real_, selected_ci_low = NA_real_, selected_ci_high = NA_real_,
      background_mean = NA_real_, difference = NA_real_,
      difference_ci_low = NA_real_, difference_ci_high = NA_real_
    ))
  }

  selected_mean <- sum(by_region$selected_sum) / sum(by_region$selected_n)
  background_mean <- if (sum(by_region$background_n) > 0L) {
    sum(by_region$background_sum) / sum(by_region$background_n)
  } else {
    NA_real_
  }

  set.seed(seed)
  n_region <- nrow(by_region)
  draw <- matrix(
    sample.int(n_region, n_region * n_boot, replace = TRUE),
    nrow = n_region,
    ncol = n_boot
  )
  boot_selected_n <- colSums(matrix(by_region$selected_n[draw], nrow = n_region))
  boot_selected <- colSums(matrix(by_region$selected_sum[draw], nrow = n_region)) /
    boot_selected_n

  boot_background_n <- colSums(matrix(by_region$background_n[draw], nrow = n_region))
  boot_background <- colSums(matrix(by_region$background_sum[draw], nrow = n_region)) /
    boot_background_n
  boot_difference <- boot_selected - boot_background

  list(
    selected_mean = selected_mean,
    selected_ci_low = unname(quantile(boot_selected, 0.025, na.rm = TRUE)),
    selected_ci_high = unname(quantile(boot_selected, 0.975, na.rm = TRUE)),
    background_mean = background_mean,
    difference = selected_mean - background_mean,
    difference_ci_low = unname(quantile(boot_difference, 0.025, na.rm = TRUE)),
    difference_ci_high = unname(quantile(boot_difference, 0.975, na.rm = TRUE))
  )
}

summarize_selection <- function(d, selected, trait, method, rule, cutoff,
                                n_boot = 1000L, seed = 1L) {
  selected <- selected & is.finite(d$evo2_delta_score)
  stat <- cluster_bootstrap(
    score = d$evo2_delta_score,
    region = d$Region,
    selected = selected,
    n_boot = n_boot,
    seed = seed
  )
  selected_score <- d$evo2_delta_score[selected]
  data.table(
    Trait = trait,
    method = method,
    selection_rule = rule,
    cutoff = cutoff,
    n_selected = sum(selected),
    n_regions_selected = uniqueN(d$Region[selected]),
    raw_evo2_mean = stat$selected_mean,
    raw_evo2_ci_low = stat$selected_ci_low,
    raw_evo2_ci_high = stat$selected_ci_high,
    raw_evo2_median = if (length(selected_score)) median(selected_score) else NA_real_,
    raw_evo2_q25 = if (length(selected_score)) unname(quantile(selected_score, 0.25)) else NA_real_,
    raw_evo2_q75 = if (length(selected_score)) unname(quantile(selected_score, 0.75)) else NA_real_,
    all_snp_mean = mean(d$evo2_delta_score),
    background_mean = stat$background_mean,
    mean_difference = stat$difference,
    difference_ci_low = stat$difference_ci_low,
    difference_ci_high = stat$difference_ci_high,
    deleterious_shift = -stat$difference,
    deleterious_shift_ci_low = -stat$difference_ci_high,
    deleterious_shift_ci_high = -stat$difference_ci_low
  )
}

summarize_thresholds <- function(data, thresholds = seq(0.50, 0.90, 0.05),
                                 n_boot = 1000L, seed = 20260828L) {
  rbindlist(lapply(unique(data$Trait), function(trait) {
    d <- data[Trait == trait]
    reference_pip <- d[[method_columns[["Classic_Annotations"]]]]
    rbindlist(lapply(seq_along(thresholds), function(i) {
      # Preserve the original script's comparison: the Classic-Annotations PIP
      # threshold determines N, and every method contributes its top N SNPs.
      # Thus the two methods are evaluated at identical list sizes.
      top_n <- sum(is.finite(reference_pip) & reference_pip > thresholds[i])
      if (top_n == 0L) return(NULL)
      rbindlist(lapply(seq_along(method_columns), function(m) {
        method <- names(method_columns)[m]
        pip <- d[[method_columns[[m]]]]
        ord <- order(-pip, d$SNP, na.last = NA)
        selected <- rep(FALSE, nrow(d))
        selected[ord[seq_len(min(top_n, length(ord)))]] <- TRUE
        summarize_selection(
          d = d,
          selected = selected,
          trait = trait,
          method = method,
          rule = "matched_top_N_from_Classic_PIP_threshold",
          cutoff = thresholds[i],
          n_boot = n_boot,
          seed = seed + 10000L * match(trait, unique(data$Trait)) + 100L * m + i
        )
      }))
    }))
  }))
}

summarize_top_n <- function(data, top_n = c(10L, 20L, 50L, 100L, 200L),
                            n_boot = 1000L, seed = 20260828L) {
  rbindlist(lapply(unique(data$Trait), function(trait) {
    d <- data[Trait == trait]
    rbindlist(lapply(seq_along(method_columns), function(m) {
      method <- names(method_columns)[m]
      pip <- d[[method_columns[[m]]]]
      ord <- order(-pip, d$SNP, na.last = NA)
      rbindlist(lapply(seq_along(top_n), function(i) {
        n_take <- min(top_n[i], length(ord))
        selected <- rep(FALSE, nrow(d))
        if (n_take > 0L) selected[ord[seq_len(n_take)]] <- TRUE
        summarize_selection(
          d = d,
          selected = selected,
          trait = trait,
          method = method,
          rule = "top_N",
          cutoff = top_n[i],
          n_boot = n_boot,
          seed = seed + 20000L * match(trait, unique(data$Trait)) + 100L * m + i
        )
      }))
    }))
  }))
}

publication_theme <- function() {
  theme_classic(base_size = 11) +
    theme(
      axis.title = element_text(face = "bold"),
      strip.text = element_text(face = "bold"),
      legend.position = "bottom",
      legend.title = element_blank(),
      plot.title = element_text(face = "bold", hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5, color = "grey30")
    )
}

make_mean_plot <- function(summary, x_var, x_label, title) {
  x_sym <- as.name(x_var)
  ggplot(summary, aes(x = !!x_sym, y = raw_evo2_mean,
                      color = method, group = method)) +
    geom_hline(
      data = unique(summary[, .(Trait, all_snp_mean)]),
      aes(yintercept = all_snp_mean),
      color = "grey45",
      linewidth = 0.45,
      linetype = "dashed"
    ) +
    geom_ribbon(
      aes(ymin = raw_evo2_ci_low, ymax = raw_evo2_ci_high, fill = method),
      color = NA,
      alpha = 0.14
    ) +
    geom_line(linewidth = 0.8) +
    geom_point(size = 1.8) +
    facet_wrap(~Trait) +
    scale_color_manual(
      values = method_colors,
      breaks = names(method_columns),
      labels = method_legend_labels,
      drop = FALSE
    ) +
    scale_fill_manual(values = method_colors, drop = FALSE, guide = "none") +
    scale_y_continuous(labels = scales::label_scientific(digits = 2)) +
    labs(
      title = title,
      subtitle = "Lower (more negative) Evo 2 delta scores indicate greater predicted disruption",
      x = x_label,
      y = "Mean raw Evo 2 delta score"
    ) +
    publication_theme()
}

make_difference_plot <- function(summary, x_var, x_label, title) {
  x_sym <- as.name(x_var)
  ggplot(summary, aes(x = !!x_sym, y = deleterious_shift,
                      color = method, group = method)) +
    geom_hline(yintercept = 0, color = "grey45", linewidth = 0.45,
               linetype = "dashed") +
    geom_ribbon(
      aes(ymin = deleterious_shift_ci_low, ymax = deleterious_shift_ci_high,
          fill = method),
      color = NA,
      alpha = 0.14
    ) +
    geom_line(linewidth = 0.8) +
    geom_point(size = 1.8) +
    facet_wrap(~Trait) +
    scale_color_manual(
      values = method_colors,
      breaks = names(method_columns),
      labels = method_legend_labels,
      drop = FALSE
    ) +
    scale_fill_manual(values = method_colors, drop = FALSE, guide = "none") +
    scale_y_continuous(labels = scales::label_scientific(digits = 2)) +
    labs(
      title = title,
      subtitle = "Positive values mean selected SNPs have more negative, more disruptive raw scores",
      x = x_label,
      y = "Deleterious shift\n(background mean - selected mean)"
    ) +
    publication_theme()
}

save_plot_pair <- function(plot, stem, figure_dir, width, height) {
  dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
  ggsave(file.path(figure_dir, paste0(stem, ".pdf")), plot,
         width = width, height = height, units = "in", device = cairo_pdf)
  ggsave(file.path(figure_dir, paste0(stem, ".png")), plot,
         width = width, height = height, units = "in", dpi = 500, bg = "white")
}

run_raw_evo2_validation <- function(data, analysis_label, module_root) {
  thresholds <- summarize_thresholds(data)
  top_n <- summarize_top_n(data)

  result_dir <- file.path(module_root, "results")
  figure_dir <- file.path(dirname(module_root), "figures", "evo2_raw_score_external_validation")
  dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)

  fwrite(
    thresholds,
    file.path(result_dir, paste0("raw_evo2_enrichment_by_pip_threshold_", analysis_label, ".tsv")),
    sep = "\t"
  )
  fwrite(
    top_n,
    file.path(result_dir, paste0("raw_evo2_enrichment_by_top_n_", analysis_label, ".tsv")),
    sep = "\t"
  )

  threshold_mean <- make_mean_plot(
    thresholds,
    "cutoff",
    "PIP cutoff",
    paste0("Raw Evo 2 score among prioritized SNPs (", toupper(analysis_label), ")")
  )
  threshold_difference <- make_difference_plot(
    thresholds,
    "cutoff",
    "PIP cutoff",
    paste0("Raw Evo 2 deleterious enrichment (", toupper(analysis_label), ")")
  )
  top_n_mean <- make_mean_plot(
    top_n,
    "cutoff",
    "Number of top-ranked SNPs",
    paste0("Raw Evo 2 score across matched top-N lists (", toupper(analysis_label), ")")
  )
  top_n_difference <- make_difference_plot(
    top_n,
    "cutoff",
    "Number of top-ranked SNPs",
    paste0("Raw Evo 2 top-N deleterious enrichment (", toupper(analysis_label), ")")
  )

  plot_height <- if (uniqueN(data$Trait) > 1L) 4.5 else 4.2
  save_plot_pair(threshold_mean, paste0("raw_evo2_mean_by_pip_threshold_", analysis_label),
                 figure_dir, 7.8, plot_height)
  save_plot_pair(threshold_difference, paste0("raw_evo2_difference_by_pip_threshold_", analysis_label),
                 figure_dir, 7.8, plot_height)
  save_plot_pair(top_n_mean, paste0("raw_evo2_mean_by_top_n_", analysis_label),
                 figure_dir, 7.8, plot_height)
  save_plot_pair(top_n_difference, paste0("raw_evo2_difference_by_top_n_", analysis_label),
                 figure_dir, 7.8, plot_height)

  invisible(list(thresholds = thresholds, top_n = top_n))
}
