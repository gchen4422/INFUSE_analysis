.infuse_project <- Sys.getenv("INFUSE_PROJECT", "/path/to/baseline_simulation")
.infuse_work <- Sys.getenv("INFUSE_REVISION", "/path/to/revision_output")
# Shared simulation and inference helpers for variant-availability analyses.
# The native simulation workflow supplies its own genotype-based inputs.
# Masked-availability utilities in this file generate complete summary statistics
# before hiding selected ancestry entries; availability denotes input coverage.
# Their per-cell random seeds differ from the baseline generator. The causal
# sampling coefficient is 2, with five signal annotations and ten fitted inputs.
# Signal labels: 0 noncausal, 1 EUR-specific, 2 AFR-specific, 3 shared.
# Simulated configuration probabilities are (0.25, 0.25, 0.50); the fitted
# EUR-only, AFR-only and shared configuration priors are (3/7, 3/7, 1/7).

suppressMessages({
  library(data.table)
  library(dplyr)
  library(mvtnorm)
  library(snpStats)
})

INF_DIR   <- .infuse_project
BASE_DIR  <- Sys.getenv("INFUSE_SIMULATION", dirname(.infuse_project))
GENO_EUR  <- file.path(BASE_DIR, "risk_loci_ld_eur")
GENO_AFR  <- file.path(BASE_DIR, "risk_loci_ld_afr")
ANNOT_RAW <- file.path(BASE_DIR, "risk_loci_annot")

PKG_DIR   <- file.path(.infuse_work, "INFUSE_variant_availability")
PKG_LIB   <- file.path(PKG_DIR, ".lib")
V2_DIR    <- file.path(PKG_DIR, "simulation_v2")
V2_OUT    <- file.path(V2_DIR, "results")

# ancestry 1 = EUR (LD1 / zscore_1), ancestry 2 = AFR (LD2 / zscore_2).
# This ordering matches the existing pipeline, so INFUSE's column_config is
#   [[1]] = EUR, [[2]] = AFR, [[3]] = EUR+AFR (shared)
# and pip_config columns line up with Signal == 1, 2, 3 respectively.
ANCESTRY <- c("EUR", "AFR")

NUM_EU <- 300000
NUM_BB <- 300000

# Annotations.  The SIGNAL panel drives the simulation truth (causal variants are
# sampled with probability proportional to exp(A %*% enrichment_beta)); the FIT
# panel is what the annotation-informed methods actually receive, and adds five
# null distractors so the fitter has to identify the informative ones rather than
# being handed the exact truth basis.  The two lists match the low-dimensional
# (base_d = 10) panel of the ancestry-specific-prior module.
ANNOT_SIGNAL_COLS <- c("H3K27ac_Hnisz_common", "DHS_Trynka_common",
                       "H3K4me3_Trynka_common", "Conserved_LindbladToh_common",
                       "non_synonymous_common")
ANNOT_NULL_COLS <- c("Coding_UCSC_common", "Coding_UCSC.flanking.500_common",
                     "CTCF_Hoffman_common", "Enhancer_Andersson_common",
                     "DHS_peaks_Trynka_common")
ANNOT_FIT_COLS <- c(ANNOT_SIGNAL_COLS, ANNOT_NULL_COLS)
# Enrichment coefficient on each signal annotation used to build the truth.
ENRICHMENT_BETA <- 2


## ---------------------------------------------------------------- AS-SV grid ----

#' AS-SV scenario grid
#'
#' `frac` is the fraction of eligible variants made ancestry-specific-segregating.
#' `pool` says which variants are eligible:
#'   noncausal  -- only Signal == 0 variants
#'   ascv       -- only the AS-CV variants (Signal 1 or 2), availability ALIGNED with
#'                 the causal ancestry (see assign_assv())
#'   shared     -- only the shared causal variants (Signal == 3)
#'   any        -- any variant
assv_scenarios <- function() {
  rbind(
    data.frame(assv = "V0_none",         frac = 0.00, pool = "noncausal", stringsAsFactors = FALSE),
    data.frame(assv = "V1_noncausal_10", frac = 0.10, pool = "noncausal", stringsAsFactors = FALSE),
    data.frame(assv = "V2_ascv_aligned", frac = 1.00, pool = "ascv",      stringsAsFactors = FALSE),
    data.frame(assv = "V3_shared_causal",frac = 1.00, pool = "shared",    stringsAsFactors = FALSE),
    data.frame(assv = "V4_any_05",       frac = 0.05, pool = "any",       stringsAsFactors = FALSE),
    data.frame(assv = "V4_any_10",       frac = 0.10, pool = "any",       stringsAsFactors = FALSE),
    data.frame(assv = "V4_any_20",       frac = 0.20, pool = "any",       stringsAsFactors = FALSE)
  )
}


## Causal-variant sampling ----------------------------------------------------

process_geno <- function(geno_data) {
  apply(geno_data, 2, function(x) {
    x[is.na(x)] <- mean(x, na.rm = TRUE)
    x <- scale(x - mean(x))
    return(x)
  })
}


#' Load one locus exactly as 3_shared_50_simulation_mesusie_inf.R does
load_locus_v2 <- function(LD_BLOCK) {
  outfile_EU <- paste0(GENO_EUR, "/loci_", LD_BLOCK)
  EU_plink <- read.plink(paste0(outfile_EU, ".bed"))
  EU_plink_geno <- as(EU_plink$genotypes, "numeric")
  EU_plink_geno <- process_geno(EU_plink_geno)
  EU_cov <- as.matrix(fread(paste0(GENO_EUR, "/loci_", LD_BLOCK, ".ld")))

  outfile_BB <- paste0(GENO_AFR, "/loci_", LD_BLOCK)
  BB_plink <- read.plink(paste0(outfile_BB, ".bed"))
  BB_plink$map$CHR_POS <- paste0(BB_plink$map$chromosome, ":", BB_plink$map$position)
  BB_plink_geno <- as(BB_plink$genotypes, "numeric")
  colnames(BB_plink_geno) <- BB_plink$map$CHR_POS
  BB_plink_geno <- process_geno(BB_plink_geno)
  BB_cov <- as.matrix(fread(paste0(GENO_AFR, "/loci_", LD_BLOCK, ".ld")))
  colnames(BB_cov) <- BB_plink$map$CHR_POS

  snp <- as.character(EU_plink$map$snp.name)
  dimnames(EU_cov) <- list(snp, snp)
  dimnames(BB_cov) <- list(snp, snp)

  annot_file <- fread(file.path(ANNOT_RAW, paste0("loci_", LD_BLOCK, ".annot")),
                      data.table = FALSE)
  missing_cols <- setdiff(ANNOT_FIT_COLS, names(annot_file))
  if (length(missing_cols)) {
    stop("locus ", LD_BLOCK, ": annotation column(s) absent: ",
         paste(missing_cols, collapse = ", "), call. = FALSE)
  }
  # truth: the five signal annotations, each with coefficient ENRICHMENT_BETA
  annot_signal <- as.matrix(annot_file[, ANNOT_SIGNAL_COLS, drop = FALSE])
  rownames(annot_signal) <- snp
  enrichment_beta <- rep(ENRICHMENT_BETA, length(ANNOT_SIGNAL_COLS))
  lin <- annot_signal %*% enrichment_beta
  true_prior_prob <- exp(lin) / sum(exp(lin))
  # fitting: the ten-annotation panel (five signal + five null distractors)
  annot_fit <- as.matrix(annot_file[, ANNOT_FIT_COLS, drop = FALSE])
  rownames(annot_fit) <- snp

  list(locus = LD_BLOCK, snp = snp, p = length(snp),
       EU_geno = EU_plink_geno, BB_geno = BB_plink_geno,
       EU_cov = EU_cov, BB_cov = BB_cov,
       map = EU_plink$map, annot = annot_fit, annot_signal = annot_signal,
       true_prior_prob = as.numeric(true_prior_prob))
}


#' Draw the AS-CV causal configuration, reproducing the original seeding
#'
#' Returns `Signal` (0/1/2/3) plus the causal SNP list and their configs.
draw_ascv <- function(lo, num_causal) {
  set.seed(1128)
  seed_LOCIrator <- sample.int(1000, 100)[lo$locus]
  set.seed(seed_LOCIrator)

  num_causal_SNP <- c(2, 3, 5)[num_causal]
  causal_SNP_config <- sample(c(1, 2, 3), num_causal_SNP, replace = TRUE,
                              prob = c(0.25, 0.25, 0.5))
  causal_SNP_list <- lo$snp[sample(seq_len(lo$p), num_causal_SNP, replace = FALSE,
                                   lo$true_prior_prob)]

  list(num_causal_SNP = num_causal_SNP,
       causal_SNP_config = causal_SNP_config,
       causal_SNP_list = causal_SNP_list)
}


#' Generate z-scores for one locus, causal count and heritability setting
#'
#' Uses the baseline effect structure and marginal-beta z-score construction.
simulate_zfile_v2 <- function(lo, ascv, h2_num) {
  h2 <- c(1e-4, 2e-4)[h2_num]
  p <- lo$p
  h2_small_EU <- 0 / p
  h2_small_BB <- 0 / p

  signal <- beta_EUR_all <- beta_AFR_all <- rep(0, p)

  for (causal_index in seq_len(ascv$num_causal_SNP)) {
    cs <- ascv$causal_SNP_list[causal_index]
    j <- which(lo$snp %in% cs)

    if (ascv$causal_SNP_config[causal_index] == 1) {
      EU_causal <- matrix(lo$EU_geno[, colnames(lo$EU_geno) %in% cs], ncol = 1)
      b <- rnorm(1, 0, 1)
      beta_EUR_all[j] <- b * sqrt(h2 / drop(var(EU_causal * b)))
      signal[j] <- 1
    } else if (ascv$causal_SNP_config[causal_index] == 2) {
      BB_causal <- matrix(lo$BB_geno[, colnames(lo$BB_geno) %in% cs], ncol = 1)
      b <- rnorm(1, 0, 1)
      beta_AFR_all[j] <- b * sqrt(h2 / drop(var(BB_causal * b)))
      signal[j] <- 2
    } else {
      EU_causal <- matrix(lo$EU_geno[, colnames(lo$EU_geno) %in% cs], ncol = 1)
      BB_causal <- matrix(lo$BB_geno[, colnames(lo$BB_geno) %in% cs], ncol = 1)
      beta_shared <- rmvnorm(1, mean = rep(0, 2),
                             sigma = matrix(c(1, 0.8, 0.8, 1), ncol = 2, nrow = 2))
      beta_EUR_all[j] <- beta_shared[1] * sqrt(h2 / drop(var(EU_causal * beta_shared[1])))
      beta_AFR_all[j] <- beta_shared[2] * sqrt(h2 / drop(var(BB_causal * beta_shared[2])))
      signal[j] <- 3
    }
  }

  beta_EUR_small <- rnorm(p, 0, h2_small_EU)
  beta_EUR_marginal <- lo$EU_cov %*% (beta_EUR_all + beta_EUR_small)

  beta_AFR_small <- rnorm(p, 0, h2_small_BB)
  beta_AFR_marginal <- lo$BB_cov %*% (beta_AFR_all + beta_AFR_small)

  h2_causal_EUR <- drop(var(lo$EU_geno %*% beta_EUR_all))
  h2_causal_AFR <- drop(var(lo$BB_geno %*% beta_AFR_all))

  y_null_EU <- rnorm(nrow(lo$EU_geno), 0, sqrt(1 - h2_causal_EUR - h2_small_EU * p))
  y_null_EU <- y_null_EU - mean(y_null_EU)
  y_null_BB <- rnorm(nrow(lo$BB_geno), 0, sqrt(1 - h2_causal_AFR - h2_small_BB * p))
  y_null_BB <- y_null_BB - mean(y_null_BB)

  err_beta_EU <- t(lo$EU_geno) %*% y_null_EU / nrow(lo$EU_geno)
  err_beta_BB <- t(lo$BB_geno) %*% y_null_BB / nrow(lo$BB_geno)

  err_beta_EU_scale <- sqrt(nrow(lo$EU_geno) / NUM_EU) * err_beta_EU
  err_beta_BB_scale <- sqrt(nrow(lo$BB_geno) / NUM_BB) * err_beta_BB

  z_EU <- (beta_EUR_marginal + err_beta_EU_scale) * sqrt(NUM_EU)
  z_BB <- (beta_AFR_marginal + err_beta_BB_scale) * sqrt(NUM_BB)

  data.frame(CHR = as.numeric(lo$map$chromosome[1]), POS = lo$map$position,
             RSID = lo$snp, zscore_1 = as.numeric(z_EU), zscore_2 = as.numeric(z_BB),
             Signal = signal, N_1 = NUM_EU, N_2 = NUM_BB,
             stringsAsFactors = FALSE)
}


## ------------------------------------------------------------- AS-SV layer ------

#' Assign ancestry-specific SEGREGATING status
#'
#' @param zfile z-file with `Signal`.
#' @param scen One row of [assv_scenarios()].
#' @param seed RNG seed (kept separate from the AS-CV seed so the AS-CV truth is
#'   identical across AS-SV scenarios).
#' @return Logical `p x 2` availability matrix with columns EUR, AFR.
#'
#' Alignment rule for the `ascv` pool: an AS-CV variant whose effect lives in
#' ancestry k is made available in ancestry k ONLY.  That is the canonical
#' motivating case -- a variant that segregates only in one ancestry and is causal
#' only in that same ancestry.  The reverse combination (causal in EUR but
#' observable only in AFR) is deliberately never generated here: its effect is
#' unobservable in principle, so no method could recover it and including it would
#' only dilute the power estimates.  It CAN arise by chance in the `any` pool, and
#' is then reported as a separate `unobservable` stratum rather than mixed in.
assign_assv <- function(zfile, scen, seed = 1L) {
  set.seed(seed)
  p <- nrow(zfile)
  avail <- matrix(TRUE, nrow = p, ncol = 2,
                  dimnames = list(zfile$RSID, ANCESTRY))
  if (scen$frac <= 0) return(avail)

  elig <- switch(scen$pool,
    noncausal = which(zfile$Signal == 0),
    ascv      = which(zfile$Signal %in% c(1, 2)),
    shared    = which(zfile$Signal == 3),
    any       = seq_len(p),
    stop("unknown pool: ", scen$pool)
  )
  if (length(elig) == 0) return(avail)

  n_pick <- if (scen$frac >= 1) length(elig) else floor(scen$frac * length(elig))
  if (n_pick < 1) return(avail)
  pick <- if (n_pick == length(elig)) elig else sample(elig, n_pick)

  if (scen$pool == "ascv") {
    # ALIGNED: keep the variant only in the ancestry where its effect lives
    for (j in pick) {
      if (zfile$Signal[j] == 1) avail[j, "AFR"] <- FALSE   # EUR-causal -> EUR only
      if (zfile$Signal[j] == 2) avail[j, "EUR"] <- FALSE   # AFR-causal -> AFR only
    }
  } else {
    # split evenly: half lose AFR (EUR-only), half lose EUR (AFR-only)
    half <- floor(length(pick) / 2)
    if (half > 0)                   avail[pick[seq_len(half)], "AFR"] <- FALSE
    if (length(pick) > half)        avail[pick[setdiff(seq_along(pick), seq_len(half))], "EUR"] <- FALSE
  }

  # never let a variant vanish from both ancestries
  dead <- which(rowSums(avail) == 0)
  if (length(dead)) avail[dead, ] <- TRUE
  avail
}


#' Classify each variant's AS-CV x AS-SV cell, for stratified evaluation
#'
#' `unobservable` marks a causal variant whose effect lives ONLY in an ancestry
#' where the variant is not available -- power is 0 by construction for every
#' method, so it must be reported separately, never pooled.
effect_visible_mask <- function(signal, avail) {
  stopifnot(length(signal) == nrow(avail), all(c("EUR", "AFR") %in% colnames(avail)))
  effect_visible <- rep(TRUE, length(signal))
  effect_visible[signal == 1] <- avail[signal == 1, "EUR"]
  effect_visible[signal == 2] <- avail[signal == 2, "AFR"]
  # a shared effect stays visible as long as the variant survives anywhere
  effect_visible[signal == 3] <- rowSums(avail)[signal == 3] >= 1
  effect_visible
}

classify_variants <- function(zfile, avail) {
  sig <- zfile$Signal
  n_av <- rowSums(avail)
  assv <- n_av < 2
  effect_visible <- effect_visible_mask(sig, avail)

  data.frame(
    RSID = zfile$RSID,
    Signal = sig,
    signal_label = c("noncausal", "EUR_causal", "AFR_causal", "shared_causal")[sig + 1],
    avail_EUR = avail[, "EUR"],
    avail_AFR = avail[, "AFR"],
    is_assv = assv,
    effect_visible = effect_visible,
    stratum = ifelse(sig == 0, "noncausal",
              ifelse(!effect_visible, "unobservable",
              ifelse(assv, "causal_AS_SV", "causal_present_both"))),
    stringsAsFactors = FALSE
  )
}


## ------------------------------------------------------------------ methods -----

#' Restrict inputs to variants available in ALL ancestries (the status quo)
intersect_v2 <- function(R_list, ss_list, avail, annot = NULL) {
  keep <- which(rowSums(avail) == ncol(avail))
  list(idx = keep,
       R  = lapply(R_list, function(R) R[keep, keep, drop = FALSE]),
       ss = lapply(ss_list, function(s) { z <- s[keep, , drop = FALSE]; rownames(z) <- NULL; z }),
       annot = if (is.null(annot)) NULL else annot[keep, , drop = FALSE])
}

scatter <- function(v, idx, p, fill = 0) { out <- rep(fill, p); out[idx] <- v; out }

scatter_mat <- function(m, idx, p, ncol_out, fill = 0) {
  out <- matrix(fill, nrow = p, ncol = ncol_out)
  out[idx, ] <- m
  out
}


#' Build INFUSE-style summary statistics from a z-file
make_sumstats <- function(zfile) {
  ss <- list(
    EUR = data.frame(SNP = zfile$RSID, Beta = zfile$zscore_1 / sqrt(zfile$N_1),
                     Se = 1 / sqrt(zfile$N_1), Z = zfile$zscore_1, N = zfile$N_1,
                     stringsAsFactors = FALSE),
    AFR = data.frame(SNP = zfile$RSID, Beta = zfile$zscore_2 / sqrt(zfile$N_2),
                     Se = 1 / sqrt(zfile$N_2), Z = zfile$zscore_2, N = zfile$N_2,
                     stringsAsFactors = FALSE)
  )
  ss
}


#' Map a fitted object's pip_config onto (EUR, AFR, shared) by ancestry MEMBERSHIP
#'
#' Never by column position: the availability-aware build exposes `config_names`,
#' and for the upstream/MESuSiE objects the order follows column_config given the
#' EUR-then-AFR list order.
config_pips <- function(pip_config, config_names = NULL) {
  if (!is.null(config_names)) {
    key <- vapply(strsplit(config_names, "+", fixed = TRUE),
                  function(x) paste(sort(x), collapse = "+"), character(1))
    i_eur <- match("EUR", key); i_afr <- match("AFR", key)
    i_sh  <- match(paste(sort(c("EUR", "AFR")), collapse = "+"), key)
  } else {
    i_eur <- 1L; i_afr <- 2L; i_sh <- 3L
  }
  stopifnot(!is.na(i_eur), !is.na(i_afr), !is.na(i_sh))
  list(EUR = pip_config[, i_eur], AFR = pip_config[, i_afr], shared = pip_config[, i_sh])
}


#' Run every method on one (locus, num_causal, h2, AS-SV scenario)
run_methods_v2 <- function(zfile, avail, R_list, annot, L = 10,
                           annot_method = "mlk", verbose = FALSE,
                           map = NULL, locus = NULL, susiex_work = NULL,
                           run_susiex_flag = TRUE) {
  p <- nrow(zfile)
  ss <- make_sumstats(zfile)
  # Identity-pad unavailable variants before evaluating credible-set purity,
  # which reads the supplied LD matrices directly.
  R_list <- pad_ld_unavailable(R_list, avail)
  out <- list(); timing <- list()

  tick <- function(label, expr) {
    t0 <- Sys.time()
    val <- tryCatch(expr, error = function(e) {
      message("  [", label, "] FAILED: ", conditionMessage(e)); NULL })
    timing[[label]] <<- as.numeric(difftime(Sys.time(), t0, units = "mins"))
    val
  }

  ## 1. availability-aware INFUSE on the union -------------------------------
  f <- tick("infuse_union", {
    invisible(capture.output(
      x <- INFUSE::infuse_core_avail(R_list, ss, L = L, variant_available = avail,
                                          annot = annot, annot_method = annot_method,
                                          verbose = verbose)))
    x })
  if (!is.null(f)) {
    cp <- config_pips(f$pip_config, f$availability$config_names)
    out$infuse_union <- list(pip = as.numeric(f$pip), pip_EUR = cp$EUR,
                             pip_AFR = cp$AFR, pip_shared = cp$shared,
                             cs = f$cs$cs)
  }

  ## 1b. Annotation-free availability-aware fit for MESuSiE union ------------
  ##     annot = NULL omits annotation-weight updates.
  f <- tick("mesusie_union", {
    invisible(capture.output(
      x <- INFUSE::infuse_core_avail(R_list, ss, L = L, variant_available = avail,
                                          annot = NULL, annot_method = NULL,
                                          verbose = verbose)))
    x })
  if (!is.null(f)) {
    cp <- config_pips(f$pip_config, f$availability$config_names)
    out$mesusie_union <- list(pip = as.numeric(f$pip), pip_EUR = cp$EUR,
                                     pip_AFR = cp$AFR, pip_shared = cp$shared,
                                     cs = f$cs$cs)
  }

  ## 2. current INFUSE on the intersection ----------------------------------
  ii <- intersect_v2(R_list, ss, avail, annot)
  f <- tick("infuse_intersect", {
    invisible(capture.output(
      x <- INFUSE::infuse_core(ii$R, ii$ss, L = L, annot = ii$annot,
                               annot_method = annot_method)))
    x })
  if (!is.null(f)) {
    cp <- config_pips(f$pip_config)
    out$infuse_intersect <- list(
      pip = scatter(as.numeric(f$pip), ii$idx, p),
      pip_EUR = scatter(cp$EUR, ii$idx, p), pip_AFR = scatter(cp$AFR, ii$idx, p),
      pip_shared = scatter(cp$shared, ii$idx, p),
      cs = lapply(f$cs$cs, function(z) ii$idx[z]))
  }

  ## 3. MESuSiE on the intersection -----------------------------------------
  f <- tick("mesusie_intersect", {
    invisible(capture.output(x <- MESuSiE::meSuSie_core(ii$R, ii$ss, L = L)))
    x })
  if (!is.null(f)) {
    cp <- config_pips(f$pip_config)
    out$mesusie_intersect <- list(
      pip = scatter(as.numeric(f$pip), ii$idx, p),
      pip_EUR = scatter(cp$EUR, ii$idx, p), pip_AFR = scatter(cp$AFR, ii$idx, p),
      pip_shared = scatter(cp$shared, ii$idx, p),
      cs = lapply(f$cs$cs, function(z) ii$idx[z]))
  }

  ## 4. single-ancestry SuSiE ------------------------------------------------
  #    Run EXACTLY as 4_shared_50_mesusie_inf_susie_only.R does:
  #        susie_rss(z, R, check_prior = FALSE)
  #    i.e. no `n` (z-score-only likelihood) and susieR's default L.  Passing `n`
  #    would switch susie_rss to the full RSS likelihood and would no longer be the
  #    published baseline.
  #    Scored per 5_result_summarizing_pipcalibration.R:
  #      either  = pmax(pip_EUR, pip_AFR)
  #      shared  = pmin(pip_EUR, pip_AFR)
  #      EUR/AFR = own-ancestry pip, via compute_ROC_uni_susie
  sp <- list(); scs <- list()
  for (k in ANCESTRY) {
    idx_k <- which(avail[, k])
    zk <- if (k == "EUR") zfile$zscore_1[idx_k] else zfile$zscore_2[idx_k]
    fk <- tick(paste0("susie_", k), {
      suppressWarnings(susieR::susie_rss(zk, R_list[[k]][idx_k, idx_k, drop = FALSE],
                                         check_prior = FALSE)) })
    sp[[k]] <- if (is.null(fk)) rep(0, p) else scatter(as.numeric(fk$pip), idx_k, p)
    if (!is.null(fk) && !is.null(fk$sets$cs))
      scs <- c(scs, lapply(fk$sets$cs, function(z) idx_k[z]))
  }
  out$susie <- list(pip = pmax(sp$EUR, sp$AFR),
                    pip_EUR = sp$EUR, pip_AFR = sp$AFR,
                    pip_shared = pmin(sp$EUR, sp$AFR),
                    cs = scs)

  ## 5. SuSiEx, both input strategies ---------------------------------------
  if (run_susiex_flag && !is.null(map) && !is.null(locus) && !is.null(susiex_work)) {
    for (mode in c("intersect", "union")) {
      f <- tick(paste0("susiex_", mode), {
        run_susiex(zfile, avail, map, mode,
                   work = file.path(susiex_work, mode), locus = locus) })
      if (!is.null(f)) out[[paste0("susiex_", mode)]] <- f
    }
  }

  list(res = out, timing = unlist(timing))
}


## ------------------------------------------------------------------ SuSiEx -----
#
#  SuSiEx is the key external benchmark: it natively accepts DIFFERENT variant sets
#  per ancestry, so it already does for its own model what the availability-aware
#  build does for INFUSE.  Both of its input strategies are run, mirroring
#  <...>/susiex_variant_overlap_test/run_test.sh:
#
#    susiex_intersect     -- sst files first reduced to the ancestry intersection
#                            (what the current pipeline does)
#    susiex_union         -- each ancestry keeps its own variant set
#
#  Caveat on configuration-specific PIPs: SuSiEx reports per-SNP PIPs per credible
#  set, and ancestry specificity only as CS-level post-hoc probabilities
#  (POST-HOC_PROB_POP1/2).  It does NOT give per-SNP configuration posteriors, so
#  the shared / EUR / AFR quantities below are DERIVED approximations and are
#  labelled as such wherever they are reported.

SUSIEX_BIN <- Sys.getenv("INFUSE_SUSIEX_BIN", unname(Sys.which("SuSiEx")))
PLINK_BIN  <- Sys.getenv("INFUSE_PLINK_BIN", unname(Sys.which("plink")))

#' Write a SuSiEx sst file in the column order the reference test uses
#'
#' SNP(1) CHR(2) POS(3) PVAL(4) ALT(5) REF(6) N(7) Z(8) BETA(9) SE(10)
write_susiex_sst <- function(path, snp, chr, pos, z, n, alt, ref) {
  dt <- data.frame(SNP = snp, CHR = chr, POS = pos,
                   PVAL = 2 * pnorm(-abs(z)), ALT = alt, REF = ref, N = n,
                   Z = z, BETA = z / sqrt(n), SE = 1 / sqrt(n),
                   stringsAsFactors = FALSE)
  data.table::fwrite(dt, path, sep = "\t")
  invisible(dt)
}


#' Parse a SuSiEx .snp file into one overall PIP per SNP
#'
#' The file has one PIP column per credible set; the overall inclusion probability
#' is 1 - prod(1 - PIP_k) across sets.
parse_susiex_snp <- function(f, snp_universe) {
  if (!file.exists(f)) return(NULL)
  d <- tryCatch(data.table::fread(f, data.table = FALSE), error = function(e) NULL)
  if (is.null(d) || nrow(d) == 0) return(NULL)
  pip_cols <- grep("^PIP\\(CS", names(d), value = TRUE)
  if (length(pip_cols) == 0) return(NULL)
  P <- as.matrix(d[, pip_cols, drop = FALSE])
  P[!is.finite(P)] <- 0
  overall <- 1 - apply(1 - P, 1, prod)
  out <- rep(0, length(snp_universe))
  hit <- match(as.character(d$SNP), snp_universe)
  ok <- !is.na(hit)
  out[hit[ok]] <- overall[ok]
  list(pip = out, per_cs = P, snp = as.character(d$SNP))
}


#' Parse a SuSiEx .cs file into credible sets as indices into snp_universe
parse_susiex_cs <- function(f, snp_universe) {
  if (!file.exists(f)) return(list())
  d <- tryCatch(data.table::fread(f, data.table = FALSE), error = function(e) NULL)
  if (is.null(d) || nrow(d) == 0 || !"CS_ID" %in% names(d)) return(list())
  split(match(as.character(d$SNP), snp_universe), d$CS_ID) |>
    lapply(function(z) z[!is.na(z)])
}


#' Parse the CS-level post-hoc ancestry probabilities from .summary
parse_susiex_summary <- function(f) {
  if (!file.exists(f)) return(NULL)
  ln <- readLines(f, warn = FALSE)
  ln <- ln[!grepl("^#", ln)]
  if (length(ln) < 2) return(NULL)
  d <- tryCatch(data.table::fread(text = paste(ln, collapse = "\n"), data.table = FALSE),
                error = function(e) NULL)
  if (is.null(d) || nrow(d) == 0) return(NULL)
  p1 <- suppressWarnings(as.numeric(d[["POST-HOC_PROB_POP1"]]))
  p2 <- suppressWarnings(as.numeric(d[["POST-HOC_PROB_POP2"]]))
  data.frame(CS_ID = d$CS_ID, p_pop1 = p1, p_pop2 = p2)
}


#' Run SuSiEx for one locus in one input mode
#'
#' @param mode "intersect" (ancestry intersection) or "union" (per-ancestry sets).
#' @param work Scratch directory for sst / LD / output files.
#' @return list(pip, pip_EUR, pip_AFR, pip_shared, cs) or NULL on failure.
run_susiex <- function(zfile, avail, map, mode, work, locus, threads = 2) {
  if (!file.exists(SUSIEX_BIN) || !file.exists(PLINK_BIN)) {
    message("  [susiex_", mode, "] SKIPPED: binary not found")
    return(NULL)
  }
  dir.create(work, recursive = TRUE, showWarnings = FALSE)
  snp_universe <- zfile$RSID

  idx <- switch(mode,
    intersect = list(EUR = which(rowSums(avail) == 2), AFR = which(rowSums(avail) == 2)),
    union     = list(EUR = which(avail[, "EUR"]),      AFR = which(avail[, "AFR"])),
    stop("bad mode"))
  if (min(lengths(idx)) < 10) return(NULL)

  alt <- as.character(map$allele.1)
  ref <- as.character(map$allele.2)
  chr <- as.numeric(map$chromosome)
  pos <- as.numeric(map$position)

  f_eur <- file.path(work, "eur.tsv"); f_afr <- file.path(work, "afr.tsv")
  write_susiex_sst(f_eur, zfile$RSID[idx$EUR], chr[idx$EUR], pos[idx$EUR],
                   zfile$zscore_1[idx$EUR], NUM_EU, alt[idx$EUR], ref[idx$EUR])
  write_susiex_sst(f_afr, zfile$RSID[idx$AFR], chr[idx$AFR], pos[idx$AFR],
                   zfile$zscore_2[idx$AFR], NUM_BB, alt[idx$AFR], ref[idx$AFR])

  out_name <- paste0("susiex_", mode)
  args <- c(
    sprintf("--sst_file=%s,%s", f_eur, f_afr),
    sprintf("--ref_file=%s/loci_%d,%s/loci_%d", GENO_EUR, locus, GENO_AFR, locus),
    sprintf("--ld_file=%s/eur,%s/afr", work, work),
    sprintf("--n_gwas=%d,%d", NUM_EU, NUM_BB),
    sprintf("--out_dir=%s/", work),
    sprintf("--out_name=%s", out_name),
    sprintf("--chr=%d", chr[1]),
    sprintf("--bp=%d,%d", min(pos), max(pos)),
    "--chr_col=2,2", "--bp_col=3,3", "--snp_col=1,1",
    "--a1_col=5,5", "--a2_col=6,6", "--eff_col=9,9", "--se_col=10,10",
    "--pval_col=4,4",
    sprintf("--plink=%s", PLINK_BIN),
    "--keep-ambig=True", "--mult-step=True", "--maf=0.01", "--level=0.95",
    "--min_purity=0.5", "--pval_thresh=1e-05", "--tol=0.0001",
    "--n_sig=10", "--max_iter=200", sprintf("--threads=%d", threads))

  st <- suppressWarnings(system2(SUSIEX_BIN, args,
                                 stdout = file.path(work, "run.log"),
                                 stderr = file.path(work, "run.err")))

  snp_res <- parse_susiex_snp(file.path(work, paste0(out_name, ".snp")), snp_universe)
  if (is.null(snp_res)) return(NULL)
  cs <- parse_susiex_cs(file.path(work, paste0(out_name, ".cs")), snp_universe)
  smry <- parse_susiex_summary(file.path(work, paste0(out_name, ".summary")))

  ## derived (approximate) configuration PIPs -------------------------------
  # SuSiEx gives ancestry specificity only per credible set.  Distribute each
  # SNP's per-CS PIP over the three configurations using that CS's post-hoc
  # probabilities.  Clearly an approximation; reported as such.
  p <- length(snp_universe)
  pip_eur <- pip_afr <- pip_sh <- rep(0, p)
  if (!is.null(smry) && ncol(snp_res$per_cs) >= 1) {
    hit <- match(snp_res$snp, snp_universe)
    for (k in seq_len(ncol(snp_res$per_cs))) {
      cs_id <- as.integer(sub(".*CS([0-9]+).*", "\\1", colnames(snp_res$per_cs)[k]))
      row <- which(smry$CS_ID == cs_id)
      if (length(row) != 1) next
      q1 <- smry$p_pop1[row]; q2 <- smry$p_pop2[row]
      if (!is.finite(q1) || !is.finite(q2)) next
      v <- snp_res$per_cs[, k]
      ok <- !is.na(hit)
      pip_eur[hit[ok]] <- pmax(pip_eur[hit[ok]], v[ok] * q1 * (1 - q2))
      pip_afr[hit[ok]] <- pmax(pip_afr[hit[ok]], v[ok] * (1 - q1) * q2)
      pip_sh[hit[ok]]  <- pmax(pip_sh[hit[ok]],  v[ok] * q1 * q2)
    }
  }

  list(pip = snp_res$pip, pip_EUR = pip_eur, pip_AFR = pip_afr,
       pip_shared = pip_sh, cs = cs, exit_status = st)
}


## ------------------------------------------------------------------ metrics -----

#' Discovery COUNTS over PIP thresholds -- the utility.R cutoff grid
#'
#' Returns TP / FP / n_true rather than per-locus power and FDR rates, because
#' power at a target FDR must be computed from discoveries POOLED ACROSS LOCI.
#' Selecting the threshold within a single locus is not "power at 5% FDR": with
#' only a handful of causal variants per locus the realized FDR jumps between 0 and
#' 0.5, so "closest to 0.05" almost always lands on exactly 0 and the number
#' silently becomes power at a near-zero FDR.  The original pipeline pools first --
#' 5_result_summarizing_pipcalibration.R sums compute_ROC() across loci before
#' calling FDR_Power().
compute_ROC_v2 <- function(pip, true_causal) {
  n_true <- sum(true_causal != 0)
  do.call(rbind, lapply(seq(1, 0.01, -0.01), function(x) {
    called <- pip > x
    TP <- sum(called & true_causal != 0)
    FP <- sum(called) - TP
    data.frame(Cutoff = x, TP = TP, FP = FP, n_true = n_true, n_called = TP + FP)
  }))
}


#' Single-ancestry SuSiE scored the ancestry-specific way (utility.R:20)
#'
#' A variant counts as an ancestry-k call only when ancestry k's PIP clears the
#' threshold AND the other ancestry's does not.
compute_ROC_uni_susie_v2 <- function(pip_1, pip_2, true_causal) {
  n_true <- sum(true_causal != 0)
  do.call(rbind, lapply(seq(1, 0.01, -0.01), function(x) {
    called <- pip_1 > x & pip_2 < x
    TP <- sum(called & true_causal != 0)
    FP <- sum(called) - TP
    data.frame(Cutoff = x, TP = TP, FP = FP, n_true = n_true, n_called = TP + FP)
  }))
}


#' Exact PIP-calibration bin counts over ALL variants
#'
#' Computed in the worker over the full variant set, so the summariser only has to
#' add counts.  This removes the selection bias that arises from thinning the
#' per-variant output: retaining every causal and every high-PIP variant but only a
#' sample of the rest inflates the observed causal frequency in the low-PIP bins.
#'
#' @return data.frame with one row per decile bin: n, n_true, sum_pip.
calib_bins <- function(pip, is_true, breaks = seq(0, 1, by = 0.1)) {
  gr <- cut(pip, breaks = breaks, include.lowest = TRUE)
  data.frame(
    bin = levels(gr),
    n = as.integer(table(gr)),
    n_true = as.integer(tapply(as.integer(is_true), gr, sum)[levels(gr)]),
    sum_pip = as.numeric(tapply(pip, gr, sum)[levels(gr)]),
    stringsAsFactors = FALSE
  ) |> transform(n_true = ifelse(is.na(n_true), 0L, n_true),
                 sum_pip = ifelse(is.na(sum_pip), 0, sum_pip))
}


#' Wilson score interval for a binomial proportion
#'
#' The original pipeline uses a Wald interval (utility.R:434), which can fall below
#' 0 or exceed 1 -- and does, in the sparse high-PIP bins.  Wilson is reported here
#' and the Wald bounds are kept alongside for comparability with published tables.
wilson_ci <- function(k, n, conf = 0.95) {
  z <- stats::qnorm(1 - (1 - conf) / 2)
  phat <- ifelse(n > 0, k / n, NA_real_)
  denom <- 1 + z^2 / n
  centre <- (phat + z^2 / (2 * n)) / denom
  halfw <- z * sqrt(phat * (1 - phat) / n + z^2 / (4 * n^2)) / denom
  data.frame(obs_frq = phat,
             wilson_min = pmax(0, centre - halfw),
             wilson_max = pmin(1, centre + halfw),
             wald_min = phat - z * sqrt(phat * (1 - phat) / n),
             wald_max = phat + z * sqrt(phat * (1 - phat) / n))
}


#' Credible-set metrics, including the published combined-set definition
#'
#' The published pipeline (utility.R compute_set_susie / compute_set_pip) scores the
#' UNION of all reported credible sets at a locus:
#'
#'   Power = n_true_selected / n_causal
#'   FDR   = (n_selected - n_true_selected) / n_selected
#'   Size  = n_selected
#'
#' Those are the quantities behind the manuscript's set-size / set-power figure, so
#' they are computed here under exactly that definition.  The set-level statistics
#' kept alongside (`n_cs`, `coverage` = fraction of sets containing a causal
#' variant) are a different thing -- closer to set-level precision -- and must not
#' be reported as causal coverage.
cs_metrics_v2 <- function(cs, is_causal, eligible = rep(TRUE, length(is_causal))) {
  stopifnot(length(eligible) == length(is_causal), !anyNA(eligible))
  eligible <- as.logical(eligible)
  scoped_truth <- as.integer(is_causal == 1 & eligible)
  n_causal <- sum(scoped_truth)
  # Remove out-of-scope variants from both the selected set and the denominator.
  # Empty CSs are dropped; this makes the observable endpoint a true restricted
  # evaluation universe rather than counting structurally impossible positives.
  if (!is.null(cs) && length(cs)) {
    cs <- lapply(cs, function(z) {
      z <- as.integer(z)
      z[is.finite(z) & z >= 1L & z <= length(eligible) & eligible[z]]
    })
    cs <- cs[lengths(cs) > 0L]
  }
  if (is.null(cs) || length(cs) == 0)
    return(data.frame(n_cs = 0L, mean_size = NA_real_, median_size = NA_real_,
                      coverage = NA_real_, n_causal_captured = 0L,
                      # published combined-set definition
                      Power = 0, FDR = 0, Size = 0L, True_selected = 0L))
  hit <- vapply(cs, function(z) any(scoped_truth[z] == 1), logical(1))
  sizes <- vapply(cs, length, integer(1))
  selected <- unique(as.integer(unlist(cs, use.names = FALSE)))
  selected <- selected[is.finite(selected)]
  n_selected <- length(selected)
  n_true_selected <- if (n_selected) sum(scoped_truth[selected] == 1) else 0L
  data.frame(n_cs = length(cs), mean_size = mean(sizes),
             median_size = stats::median(sizes), coverage = mean(hit),
             n_causal_captured = n_true_selected,
             Power = if (n_causal) n_true_selected / n_causal else 0,
             FDR = if (n_selected) (n_selected - n_true_selected) / n_selected else 0,
             Size = n_selected, True_selected = n_true_selected)
}


#' Identity-pad the LD rows/columns of variants unavailable in that ancestry
#'
#' Unavailable ancestry entries have zero effect contributions. Identity padding
#' also excludes their LD correlations from cross-ancestry credible-set purity.
pad_ld_unavailable <- function(R_list, avail) {
  out <- R_list
  for (k in colnames(avail)) {
    hid <- which(!avail[, k])
    if (!length(hid)) next
    Rk <- out[[k]]
    Rk[hid, ] <- 0
    Rk[, hid] <- 0
    diag(Rk) <- 1
    out[[k]] <- Rk
  }
  out
}
