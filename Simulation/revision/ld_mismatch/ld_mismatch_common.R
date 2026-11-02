.infuse_repo <- Sys.getenv("INFUSE_ANALYSIS_ROOT", "/path/to/INFUSE_analysis")
.infuse_project <- Sys.getenv("INFUSE_PROJECT", "/path/to/baseline_simulation")
.infuse_work <- Sys.getenv("INFUSE_REVISION", "/path/to/revision_output")
################################################################################
# Shared utilities for the INFUSE revision LD-mismatch sensitivity simulation.
#
# The summary statistics are NEVER changed in this module.  Each scenario changes
# only the LD matrix supplied to the existing single-ancestry and joint kriging
# diagnostics.
################################################################################

INF_DIR <- .infuse_project
BASE_DIR <- dirname(INF_DIR)
MODULE_DIR <- file.path(.infuse_work, "LD_mismatch_sensitivity")
RUN_ID <- "no_z_filter_ld_r0p1"
RESULT_DIR <- file.path(.infuse_work, "res", "LD_mismatch_sensitivity")
FIGURE_DIR <- file.path(.infuse_work, "figures", "LD_mismatch_sensitivity")
if (nzchar(RUN_ID)) {
  if (!grepl("^[A-Za-z0-9_.-]+$", RUN_ID)) {
    stop("LD_MISMATCH_RUN_ID contains unsupported characters.", call. = FALSE)
  }
  RESULT_DIR <- file.path(RESULT_DIR, RUN_ID)
  FIGURE_DIR <- file.path(FIGURE_DIR, RUN_ID)
}

Z_ROOT <- file.path(INF_DIR, "More_annotes")
LD_ROOT <- Sys.getenv("INFUSE_SHARED_LD", file.path(BASE_DIR, "shared_50"))
KRIGING_SOURCE <- file.path(Sys.getenv("INFUSE_LEGACY_MODEL_SOURCE", "/path/to/legacy_model_sources"), "modified_kriging_rss.R")
KRIGING_SEPARATE_SOURCE <- file.path(.infuse_repo, "Simulation/revision/ld_mismatch/modified_kriging_rss_separate_s.R")

scenario_grid <- function() {
  data.frame(
    scenario=c("matched","sparse_attenuation_025","sparse_attenuation_050",
               "sparse_negative_050","sparse_negative_100"),
    scenario_class=c("control",rep("sparse",4)),
    delta=c(0,.25,.50,NA,NA), lambda=c(0,rep(NA_real_,4)),
    anchor_scale=c(1,sqrt(.75),sqrt(.50),-.50,-1),
    mismatch_ancestry=c("none",rep("both",4)), stringsAsFactors=FALSE
  )
}

z_path <- function(causal_num, locus, h2_num) {
  file.path(
    Z_ROOT, paste0("causal_num_", causal_num), "summary_data", "w2_startover",
    sprintf("CAUSAL_%d_LOCI_%d_h2_%d", causal_num, locus, h2_num)
  )
}

ld_path <- function(ancestry_number, causal_num, locus, h2_num) {
  file.path(
    LD_ROOT, paste0("causal_num_", causal_num), "summary_data",
    sprintf("CAUSAL_%d_LOCI_%d_h2_%d.LD%d",
            causal_num, locus, h2_num, ancestry_number)
  )
}

load_zfile <- function(causal_num, locus, h2_num) {
  f <- z_path(causal_num, locus, h2_num)
  if (!file.exists(f)) stop("Missing summary-statistic input: ", f, call. = FALSE)
  z <- data.table::fread(f, data.table = FALSE)
  required <- c("RSID", "zscore_1", "zscore_2", "Signal", "N_1", "N_2")
  miss <- setdiff(required, names(z))
  if (length(miss)) stop("Missing z-file columns: ", paste(miss, collapse = ", "), call. = FALSE)
  z
}

load_true_ld <- function(causal_num, locus, h2_num, snp_ids) {
  read_one <- function(k) {
    f <- ld_path(k, causal_num, locus, h2_num)
    if (!file.exists(f)) stop("Missing LD input: ", f, call. = FALSE)
    R <- as.matrix(data.table::fread(f, data.table = FALSE))
    storage.mode(R) <- "double"
    if (nrow(R) != length(snp_ids) || ncol(R) != length(snp_ids)) {
      stop("LD dimension mismatch for ", f, call. = FALSE)
    }
    dimnames(R) <- list(snp_ids, snp_ids)
    R
  }
  list(EUR = read_one(1L), AFR = read_one(2L))
}

validate_correlation_matrix <- function(R, label, tol = 1e-7) {
  if (!is.matrix(R) || nrow(R) != ncol(R)) stop(label, " is not square.", call. = FALSE)
  if (any(!is.finite(R))) stop(label, " contains non-finite entries.", call. = FALSE)
  asym <- max(abs(R - t(R)))
  diag_err <- max(abs(diag(R) - 1))
  if (asym > tol) stop(label, " is not symmetric; max error = ", signif(asym, 4), call. = FALSE)
  if (diag_err > tol) stop(label, " does not have unit diagonal; max error = ",
                           signif(diag_err, 4), call. = FALSE)
  invisible(TRUE)
}

# Local PSD-preserving LD perturbation.
#
# For a selected anchor j, every correlation involving j is multiplied by
# sqrt(1-delta); correlations between two selected anchors are multiplied by
# 1-delta.  The diagonal is restored to one.  This is positive semidefinite:
# it is equivalent to shrinking each selected variable's original latent loading
# and adding delta independent variance unique to that variable.
scale_selected_ld <- function(R, selected, anchor_scale) {
  if (!is.finite(anchor_scale) || abs(anchor_scale) > 1) {
    stop("anchor_scale must be finite and lie in [-1,1].", call. = FALSE)
  }
  p <- nrow(R)
  selected <- sort(unique(as.integer(selected)))
  selected <- selected[selected >= 1L & selected <= p]
  scale_vec <- rep(1, p)
  scale_vec[selected] <- anchor_scale
  out <- R * tcrossprod(scale_vec)
  diag(out) <- 1
  dimnames(out) <- dimnames(R)
  out
}

attenuate_selected_ld <- function(R, selected, delta) {
  if (delta < 0 || delta >= 1) stop("delta must lie in [0,1).", call. = FALSE)
  scale_selected_ld(R, selected, sqrt(1 - delta))
}

# Select sparse perturbation anchors independently of z scores, causal status,
# and diagnostic output. An LD-connectivity requirement avoids perturbing a
# completely isolated variant for which changing the LD row/column would be
# essentially trivial. Sampling is random within this prespecified eligible set.
select_sparse_anchors <- function(R, fraction = 0.01, seed = 1L,
                                  min_abs_r = 0.10) {
  p <- nrow(R)
  k <- min(p, max(1L, ceiling(fraction * p)))
  abs_R <- abs(R)
  diag(abs_R) <- 0
  max_abs_r <- apply(abs_R, 1L, max)
  eligible <- which(is.finite(max_abs_r) & max_abs_r >= min_abs_r)

  set.seed(seed)
  if (length(eligible) >= k) return(sort(sample(eligible, k, replace = FALSE)))

  # A small locus may contain fewer than k variants above the connectivity
  # threshold. Preserve all eligible variants and fill deterministically by LD
  # connectivity alone. No z-score information enters this fallback.
  remaining <- setdiff(which(is.finite(max_abs_r)), eligible)
  fill <- remaining[order(max_abs_r[remaining], decreasing = TRUE)][
    seq_len(min(k - length(eligible), length(remaining)))
  ]
  ans <- sort(unique(c(eligible, fill)))
  sort(ans[seq_len(min(k, length(ans)))])
}

apply_ld_scenario <- function(scenario_row, R_true, zfile, seed_base,
                              sparse_fraction = 0.01) {
  cls <- scenario_row$scenario_class
  p <- nrow(R_true$EUR)
  empty_truth <- list(EUR = integer(0), AFR = integer(0))

  if (cls == "control") {
    return(list(R = R_true, truth = empty_truth, truth_scope = "none",
                target_rule = "no LD mismatch"))
  }

  if (cls == "sparse") {
    # Seeds do not depend on perturbation severity, so mild and strong scenarios
    # perturb the identical anchor variants within a simulation cell.
    target_eur <- select_sparse_anchors(R_true$EUR,
                                        fraction = sparse_fraction,
                                        seed = seed_base + 11L)
    target_afr <- select_sparse_anchors(R_true$AFR,
                                        fraction = sparse_fraction,
                                        seed = seed_base + 29L)
    anchor_scale <- scenario_row$anchor_scale
    if (!is.finite(anchor_scale)) {
      stop("Sparse scenario is missing anchor_scale.", call. = FALSE)
    }
    return(list(
      R = list(
        EUR = scale_selected_ld(R_true$EUR, target_eur, anchor_scale),
        AFR = scale_selected_ld(R_true$AFR, target_afr, anchor_scale)
      ),
      truth = list(EUR = target_eur, AFR = target_afr),
      truth_scope = "sparse_anchor",
      target_rule = sprintf(paste0("ceil(%.3f*p) anchors/ancestry sampled randomly ",
                                   "from max|r|>=0.10 variants, independent of Z ",
                                   "and causal status; anchor LD scale=%.3f"),
                            sparse_fraction, anchor_scale)
    ))
  }

  stop("Unknown scenario class: ", cls, call. = FALSE)
}

relative_frobenius <- function(A, B) sqrt(sum((A - B)^2) / sum(A^2))

ld_difference_summary <- function(R_true, R_used, ancestry) {
  d <- R_used - R_true
  off <- upper.tri(d)
  data.frame(
    ancestry = ancestry,
    relative_frobenius = relative_frobenius(R_true, R_used),
    mean_abs_offdiag_difference = mean(abs(d[off])),
    max_abs_offdiag_difference = max(abs(d[off])),
    stringsAsFactors = FALSE
  )
}

extract_plot_indices <- function(plot_obj, name) {
  env <- plot_obj$plot_env
  if (is.null(env) || !exists(name, envir = env, inherits = FALSE)) return(integer(0))
  as.integer(get(name, envir = env, inherits = FALSE))
}

run_kriging_diagnostics <- function(z1, z2, R1, R2, n1, n2) {
  # The published single-ancestry estimator provides s_EUR and s_AFR. The
  # project-specific joint extension estimates one shared s from the stacked
  # z-scores and block-diagonal LD. Saving all three estimates directly tests
  # whether the shared parameter compromises between a matched and mismatched
  # ancestry in the asymmetric reviewer scenarios.
  s_shared <- estimate_s_rss_joint(z1, z2, R1, R2, n1, n2)
  s_separate <- estimate_s_rss_joint_separate(z1, z2, R1, R2, n1, n2)

  eu <- susieR::kriging_rss(z1, R1, n = n1)
  afr <- susieR::kriging_rss(z2, R2, n = n2)
  # Pass the estimates explicitly. The original kriging_rss_joint() otherwise
  # shrinks z once and then calls estimate_s_rss_joint() with n again, applying
  # the finite-sample transformation a second time. Explicit parameters ensure
  # that both shared-s and separate-s models use one shrinkage step and that the
  # saved estimates are exactly those used in the covariance regularization.
  joint <- kriging_rss_joint(z1, z2, R1, R2, n1, n2, s = s_shared)
  joint_separate <- kriging_rss_joint_separate_s(
    z1, z2, R1, R2, n1, n2,
    s1 = unname(s_separate["s1"]),
    s2 = unname(s_separate["s2"])
  )
  s_eur <- unname(s_separate["s1"])
  s_afr <- unname(s_separate["s2"])

  p <- length(z1)
  eu_flags <- sort(unique(extract_plot_indices(eu$plot, "idx")))
  afr_flags <- sort(unique(extract_plot_indices(afr$plot, "idx")))
  jf1 <- extract_plot_indices(joint$plot, "idx_flag1")
  jf2 <- extract_plot_indices(joint$plot, "idx_flag2")
  jfj <- extract_plot_indices(joint$plot, "idx_flag_joint")
  jfj <- if (length(jfj)) ((jfj - 1L) %% p) + 1L else integer(0)
  joint_flags <- sort(unique(c(jf1, jf2, jfj)))
  joint_separate_flags <- joint_separate$flags$joint

  # The standardized-residual-only count is retained as an explicitly exploratory
  # broad inconsistency diagnostic.  The manuscript's allele-flip LR flag remains
  # the primary joint test.
  cd <- joint$conditional_dist
  cd_separate <- joint_separate$conditional_dist
  generic_eur <- which(abs(cd$z_std_diff[seq_len(p)]) > 4)
  generic_afr <- which(abs(cd$z_std_diff[p + seq_len(p)]) > 4)

  list(
    s_summary = data.frame(s_shared = s_shared, s_EUR = s_eur, s_AFR = s_afr),
    flags = list(EUR = eu_flags, AFR = afr_flags,
                 separate_union = sort(unique(c(eu_flags, afr_flags))),
                 joint = joint_flags,
                 joint_separate_s = joint_separate_flags,
                 residual_only = sort(unique(c(generic_eur, generic_afr)))),
    score_summary = data.frame(
      s_model = rep(c("shared_s", "separate_s"), each = 2),
      ancestry = rep(c("EUR", "AFR"), 2),
      mean_abs_zstd = c(mean(abs(cd$z_std_diff[seq_len(p)])),
                        mean(abs(cd$z_std_diff[p + seq_len(p)])),
                        mean(abs(cd_separate$z_std_diff[seq_len(p)])),
                        mean(abs(cd_separate$z_std_diff[p + seq_len(p)]))),
      q99_abs_zstd = c(stats::quantile(abs(cd$z_std_diff[seq_len(p)]), 0.99),
                       stats::quantile(abs(cd$z_std_diff[p + seq_len(p)]), 0.99),
                       stats::quantile(abs(cd_separate$z_std_diff[seq_len(p)]), 0.99),
                       stats::quantile(abs(cd_separate$z_std_diff[p + seq_len(p)]), 0.99)),
      max_abs_zstd = c(max(abs(cd$z_std_diff[seq_len(p)])),
                       max(abs(cd$z_std_diff[p + seq_len(p)])),
                       max(abs(cd_separate$z_std_diff[seq_len(p)])),
                       max(abs(cd_separate$z_std_diff[p + seq_len(p)]))),
      stringsAsFactors = FALSE
    )
  )
}

anchor_signal_bin <- function(abs_z) {
  cut(abs_z,
      breaks = c(-Inf, 1, 2, 4, Inf),
      labels = c("|Z|<1", "1<=|Z|<2", "2<=|Z|<4", "|Z|>=4"),
      right = FALSE)
}

# Per-anchor detection records for post-sampling signal-stratified evaluation.
# Signal strength is used only here, after anchor sampling, and never affects
# the truth set.
build_anchor_detection_rows <- function(diag_res, truth, z1, z2) {
  idx <- sort(unique(c(truth$EUR, truth$AFR)))
  if (!length(idx)) return(data.frame())

  in_eur <- idx %in% truth$EUR
  in_afr <- idx %in% truth$AFR
  abs_z_eur <- abs(z1[idx])
  abs_z_afr <- abs(z2[idx])
  abs_z_perturbed <- ifelse(
    in_eur & in_afr,
    pmax(abs_z_eur, abs_z_afr),
    ifelse(in_eur, abs_z_eur, abs_z_afr)
  )
  perturbation_ancestry <- ifelse(in_eur & in_afr, "both",
                                  ifelse(in_eur, "EUR", "AFR"))
  tests <- c("separate_union", "joint", "joint_separate_s")

  do.call(rbind, lapply(tests, function(test) {
    data.frame(
      test = test,
      SNP_index = idx,
      perturbation_ancestry = perturbation_ancestry,
      abs_z_EUR = abs_z_eur,
      abs_z_AFR = abs_z_afr,
      abs_z_perturbed = abs_z_perturbed,
      z_stratum = as.character(anchor_signal_bin(abs_z_perturbed)),
      detected = idx %in% diag_res$flags[[test]],
      stringsAsFactors = FALSE
    )
  }))
}

classification_metrics <- function(flags, truth, p, truth_scope) {
  flags <- sort(unique(flags[flags >= 1L & flags <= p]))
  truth <- sort(unique(truth[truth >= 1L & truth <= p]))
  n_flag <- length(flags)

  if (truth_scope != "sparse_anchor") {
    return(data.frame(
      n_flagged = n_flag, n_truth = if (truth_scope == "none") 0L else NA_integer_,
      TP = NA_integer_, FP = NA_integer_, FN = NA_integer_, TN = NA_integer_,
      sensitivity = NA_real_, FDR = NA_real_, FPR = NA_real_,
      flag_rate = n_flag / p
    ))
  }

  tp <- length(intersect(flags, truth))
  fp <- length(setdiff(flags, truth))
  fn <- length(setdiff(truth, flags))
  tn <- p - tp - fp - fn
  data.frame(
    n_flagged = n_flag, n_truth = length(truth), TP = tp, FP = fp, FN = fn, TN = tn,
    sensitivity = if (length(truth)) tp / length(truth) else NA_real_,
    FDR = if (n_flag) fp / n_flag else 0,
    FPR = if ((p - length(truth)) > 0) fp / (p - length(truth)) else NA_real_,
    flag_rate = n_flag / p
  )
}

build_metric_rows <- function(diag_res, truth, p, truth_scope) {
  tests <- c("EUR", "AFR", "separate_union", "joint", "joint_separate_s",
             "residual_only")
  truth_by_test <- list(
    EUR = truth$EUR,
    AFR = truth$AFR,
    separate_union = union(truth$EUR, truth$AFR),
    joint = union(truth$EUR, truth$AFR),
    joint_separate_s = union(truth$EUR, truth$AFR),
    residual_only = union(truth$EUR, truth$AFR)
  )
  do.call(rbind, lapply(tests, function(test) {
    data.frame(
      test = test,
      classification_metrics(diag_res$flags[[test]], truth_by_test[[test]],
                             p, truth_scope),
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  }))
}
