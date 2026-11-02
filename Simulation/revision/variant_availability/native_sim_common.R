.infuse_repo <- Sys.getenv("INFUSE_ANALYSIS_ROOT", "/path/to/INFUSE_analysis")
################################################################################
# Native AS-SV simulation: shared/AS-CV truth plus true ancestry-specific SNPs
################################################################################

this_native_dir <- file.path(.infuse_repo, "Simulation/revision/variant_availability")
source(file.path(.infuse_repo, "Simulation/revision/variant_availability/native_assv_config.R"))

# Shared inference wrappers and evaluation helpers.
source(file.path(.infuse_repo, "Simulation/revision/variant_availability/common/inference_helpers.R"))

# SuSiEx forwards this executable to its internal PLINK LD construction. Its own
# --threads flag is not forwarded to PLINK, so use a local wrapper to enforce the
# two CPUs requested from SLURM and avoid node-wide thread detection.
PLINK_BIN <- file.path(.infuse_repo, "Simulation/revision/variant_availability/plink_2threads.sh")

suppressMessages({
  library(data.table)
  library(mvtnorm)
  library(snpStats)
})

NATIVE_SIM_OUT <- file.path(NATIVE_DIR, "results_native_sim")
NATIVE_SIM_WORK <- file.path(NATIVE_DIR, ".work_native_sim")
NATIVE_MULTISUSIE_OUT <- file.path(NATIVE_SIM_OUT, "multisusie_intersect")
NATIVE_INFUSE_D187_OUT <- file.path(NATIVE_SIM_OUT, "infuse_union_d187")
NATIVE_MULTISUSIE_INPUT <- file.path(NATIVE_DIR, "multisusie_inputs")
NATIVE_ASSV_LEVELS <- 0:2
NATIVE_BASE_COUNTS <- c(2L, 3L, 5L)
NATIVE_H2 <- c(1e-4, 2e-4)
NATIVE_MAX_ASSV_EACH <- 2L

native_locus_dir <- function(locus) {
  file.path(PANEL_ROOT, sprintf("locus_%03d", as.integer(locus)))
}

native_annotation187_path <- function(locus) {
  file.path(native_locus_dir(locus), "annotation_187.rds")
}

standardize_native_genotype <- function(prefix, target_ids, flip_ids = character()) {
  pl <- read.plink(paste0(prefix, ".bed"))
  G <- as(pl$genotypes, "numeric")
  ids <- as.character(pl$map$snp.name)
  hit <- match(target_ids, ids)
  if (anyNA(hit)) stop("Genotype/BIM mismatch: ", prefix)
  G <- G[, hit, drop = FALSE]
  cm <- colMeans(G, na.rm = TRUE)
  miss <- which(is.na(G), arr.ind = TRUE)
  if (nrow(miss)) G[miss] <- cm[miss[, 2L]]
  G <- sweep(G, 2L, colMeans(G), "-")
  s <- sqrt(colSums(G^2) / (nrow(G) - 1))
  if (any(!is.finite(s) | s <= 0)) stop("Zero/non-finite genotype SD: ", prefix)
  G <- sweep(G, 2L, s, "/")
  colnames(G) <- target_ids
  j <- match(flip_ids, target_ids, nomatch = 0L)
  j <- j[j > 0L]
  if (length(j)) G[, j] <- -G[, j, drop = FALSE]
  G
}

load_native_locus <- function(locus, load_genotypes = TRUE) {
  d <- native_locus_dir(locus)
  required <- c(file.path(d, "COMPLETE"), file.path(d, "variant_map.tsv"),
                file.path(d, "annotation.rds"), file.path(d, "ld", "EUR.rds"),
                file.path(d, "ld", "AFR.rds"))
  missing <- required[!file.exists(required)]
  if (length(missing)) stop("Missing native-locus file(s):\n", paste(missing, collapse = "\n"))

  vm_all <- fread(file.path(d, "variant_map.tsv"))
  annot <- readRDS(file.path(d, "annotation.rds"))
  # annotation.rds contains only baselineUKB2.2-matched variants. Restrict every
  # method to this same analysis universe, matching the established simulation.
  keep <- match(rownames(annot), vm_all$SNP)
  if (anyNA(keep) || anyDuplicated(keep)) stop("Annotation/variant-map alignment failure")
  vm <- vm_all[keep]
  ids <- as.character(vm$SNP)
  avail <- as.matrix(vm[, .(available_EUR, available_AFR)])
  storage.mode(avail) <- "logical"
  colnames(avail) <- ANCESTRY
  rownames(avail) <- ids
  if (any(rowSums(avail) < 1L)) stop("Dead union variant in locus ", locus)

  annot <- annot[match(ids, rownames(annot)), , drop = FALSE]
  if (anyNA(annot) || !identical(rownames(annot), ids)) stop("Annotation alignment failure")

  idx <- list(EUR = which(avail[, "EUR"]), AFR = which(avail[, "AFR"]))
  R_obs <- list(EUR = readRDS(file.path(d, "ld", "EUR.rds")),
                AFR = readRDS(file.path(d, "ld", "AFR.rds")))
  for (anc in ANCESTRY) {
    expected <- ids[idx[[anc]]]
    hit <- match(expected, rownames(R_obs[[anc]]))
    if (anyNA(hit)) stop("LD does not contain annotation-complete SNPs: ", anc)
    R_obs[[anc]] <- R_obs[[anc]][hit, hit, drop = FALSE]
    if (!identical(rownames(R_obs[[anc]]), expected)) stop("LD order mismatch: ", anc)
  }

  # Full union matrices are neutral identity-padded. Inference also asserts the
  # availability mask, so padded cells never contribute evidence.
  p <- length(ids)
  R_union <- lapply(ANCESTRY, function(anc) {
    z <- diag(p)
    ii <- idx[[anc]]
    z[ii, ii] <- R_obs[[anc]]
    dimnames(z) <- list(ids, ids)
    z
  })
  names(R_union) <- ANCESTRY

  A1 <- ifelse(avail[, "EUR"], vm$A1_EUR, vm$A1_AFR)
  A2 <- ifelse(avail[, "EUR"], vm$A2_EUR, vm$A2_AFR)
  map <- data.frame(snp.name = ids, chromosome = vm$chromosome,
                    position = vm$position, allele.1 = A1, allele.2 = A2,
                    stringsAsFactors = FALSE)

  G <- NULL
  if (load_genotypes) {
    G <- list(
      EUR = standardize_native_genotype(file.path(d, "plink", "EUR"), ids[idx$EUR]),
      AFR = standardize_native_genotype(file.path(d, "plink", "AFR"), ids[idx$AFR],
                                        vm$SNP[vm$flip_AFR_to_EUR])
    )
  }

  # Truth enrichment, identical for all three causal candidate pools. Subtracting
  # the global maximum is numerical stabilization only and cancels on pool-wise
  # normalization: P(j | pool) is proportional to exp(A_j' * [2,2,2,2,2]).
  signal_A <- annot[, ANNOT_SIGNAL_COLS, drop = FALSE]
  lin <- drop(signal_A %*% rep(ENRICHMENT_BETA, ncol(signal_A)))
  prior_score <- exp(lin - max(lin))

  list(locus = as.integer(locus), dir = d, snp = ids, p = p, vm = vm,
       avail = avail, idx = idx, annot = annot, prior_score = prior_score,
       R_obs = R_obs, R_union = R_union, G = G, map = map)
}

weighted_sample_pool <- function(pool, n, weights) {
  if (length(pool) < n) stop("Causal pool has ", length(pool), " variants; needs ", n)
  sample(pool, n, replace = FALSE, prob = weights[pool])
}

# Draw the maximum nested truth once. The first n AS-SVs from each ancestry are
# used at level n, while the shared/AS-CV base truth never changes with AS-SV level.
draw_native_truth <- function(lo, base_count_index) {
  n_base <- NATIVE_BASE_COUNTS[base_count_index]
  if (is.na(n_base)) stop("base_count_index must be 1, 2, or 3")
  shared_pool <- which(rowSums(lo$avail) == 2L)
  eur_pool <- which(lo$avail[, "EUR"] & !lo$avail[, "AFR"])
  afr_pool <- which(lo$avail[, "AFR"] & !lo$avail[, "EUR"])

  set.seed(1128)
  locus_seed <- sample.int(1000, 100)[lo$locus]

  # Preserve the established base configuration probabilities and sampling order.
  set.seed(locus_seed)
  base_config <- sample(c(1L, 2L, 3L), n_base, replace = TRUE,
                        prob = c(0.25, 0.25, 0.5))
  base_idx <- weighted_sample_pool(shared_pool, n_base, lo$prior_score)

  # Independent deterministic streams prevent unequal pool sizes from coupling
  # EUR-only and AFR-only causal selection.
  set.seed(locus_seed + 30011L + 101L * base_count_index)
  eur_assv <- weighted_sample_pool(eur_pool, NATIVE_MAX_ASSV_EACH, lo$prior_score)
  set.seed(locus_seed + 60013L + 101L * base_count_index)
  afr_assv <- weighted_sample_pool(afr_pool, NATIVE_MAX_ASSV_EACH, lo$prior_score)

  list(locus_seed = locus_seed, n_base = n_base, base_idx = base_idx,
       base_config = base_config, eur_assv = eur_assv, afr_assv = afr_assv)
}

truth_at_assv_level <- function(lo, design, n_assv_each) {
  if (!n_assv_each %in% NATIVE_ASSV_LEVELS) stop("n_assv_each must be 0, 1, or 2")
  signal <- integer(lo$p)
  origin <- rep("noncausal", lo$p)
  signal[design$base_idx] <- design$base_config
  origin[design$base_idx] <- c("ASCV_EUR", "ASCV_AFR", "shared")[design$base_config]
  if (n_assv_each > 0L) {
    e <- design$eur_assv[seq_len(n_assv_each)]
    a <- design$afr_assv[seq_len(n_assv_each)]
    signal[e] <- 1L; origin[e] <- "native_ASSV_EUR"
    signal[a] <- 2L; origin[a] <- "native_ASSV_AFR"
  }
  data.frame(SNP = lo$snp, Signal = signal, causal_origin = origin,
             available_EUR = lo$avail[, "EUR"], available_AFR = lo$avail[, "AFR"],
             stringsAsFactors = FALSE)
}

scale_one_effect <- function(G, j, raw_b, h2) {
  v <- drop(var(G[, j] * raw_b))
  if (!is.finite(v) || v <= 0) stop("Non-finite causal genotype variance")
  raw_b * sqrt(h2 / v)
}

simulate_native_z <- function(lo, design, base_count_index, h2_index, n_assv_each) {
  h2 <- NATIVE_H2[h2_index]
  if (is.na(h2)) stop("h2_index must be 1 or 2")
  truth <- truth_at_assv_level(lo, design, n_assv_each)
  beta <- list(EUR = setNames(rep(0, length(lo$idx$EUR)), lo$snp[lo$idx$EUR]),
               AFR = setNames(rep(0, length(lo$idx$AFR)), lo$snp[lo$idx$AFR]))

  effect_seed <- design$locus_seed + 7919L * h2_index + 104729L * base_count_index
  set.seed(effect_seed)
  raw_base <- vector("list", design$n_base)
  for (i in seq_len(design$n_base)) {
    raw_base[[i]] <- if (design$base_config[i] == 3L)
      drop(rmvnorm(1L, c(0, 0), matrix(c(1, .8, .8, 1), 2L))) else rnorm(1L)
  }
  set.seed(effect_seed + 30011L); raw_eur_assv <- rnorm(NATIVE_MAX_ASSV_EACH)
  set.seed(effect_seed + 60013L); raw_afr_assv <- rnorm(NATIVE_MAX_ASSV_EACH)

  for (i in seq_len(design$n_base)) {
    snp <- lo$snp[design$base_idx[i]]
    cfg <- design$base_config[i]
    if (cfg %in% c(1L, 3L)) {
      j <- match(snp, colnames(lo$G$EUR))
      rb <- if (cfg == 3L) raw_base[[i]][1L] else raw_base[[i]]
      beta$EUR[snp] <- scale_one_effect(lo$G$EUR, j, rb, h2)
    }
    if (cfg %in% c(2L, 3L)) {
      j <- match(snp, colnames(lo$G$AFR))
      rb <- if (cfg == 3L) raw_base[[i]][2L] else raw_base[[i]]
      beta$AFR[snp] <- scale_one_effect(lo$G$AFR, j, rb, h2)
    }
  }
  if (n_assv_each > 0L) {
    for (i in seq_len(n_assv_each)) {
      se <- lo$snp[design$eur_assv[i]]
      sa <- lo$snp[design$afr_assv[i]]
      beta$EUR[se] <- scale_one_effect(lo$G$EUR, match(se, colnames(lo$G$EUR)),
                                      raw_eur_assv[i], h2)
      beta$AFR[sa] <- scale_one_effect(lo$G$AFR, match(sa, colnames(lo$G$AFR)),
                                      raw_afr_assv[i], h2)
    }
  }

  z_obs <- list()
  h2_realized <- numeric(2L); names(h2_realized) <- ANCESTRY
  for (kk in seq_along(ANCESTRY)) {
    anc <- ANCESTRY[kk]
    G <- lo$G[[anc]]; b <- beta[[anc]]
    marginal <- drop(lo$R_obs[[anc]] %*% b)
    g <- drop(G %*% b)
    h2_realized[anc] <- var(g)
    if (!is.finite(h2_realized[anc]) || h2_realized[anc] >= 1)
      stop("Invalid realized h2 in ", anc, ": ", h2_realized[anc])
    set.seed(effect_seed + 90001L + 1009L * kk)
    y0 <- rnorm(nrow(G), 0, sqrt(1 - h2_realized[anc]))
    y0 <- y0 - mean(y0)
    err <- drop(crossprod(G, y0)) / nrow(G)
    n_gwas <- if (anc == "EUR") NUM_EU else NUM_BB
    z_obs[[anc]] <- (marginal + sqrt(nrow(G) / n_gwas) * err) * sqrt(n_gwas)
  }

  z1 <- z2 <- rep(0, lo$p)
  z1[lo$idx$EUR] <- z_obs$EUR
  z2[lo$idx$AFR] <- z_obs$AFR
  zfile <- data.frame(CHR = lo$map$chromosome, POS = lo$map$position,
                      RSID = lo$snp, zscore_1 = z1, zscore_2 = z2,
                      Signal = truth$Signal, N_1 = NUM_EU, N_2 = NUM_BB,
                      stringsAsFactors = FALSE)
  list(zfile = zfile, truth = truth, beta = beta, h2_realized = h2_realized)
}

# SuSiEx writes a distinctive, complete output contract when all candidate CSs are
# filtered: exit 0, no PIP(CS*) columns in .snp, and literal NULL in both .cs and
# .summary. Its cal_pip() implementation zeros the native overall-PIP array and
# only recomputes it when ncs > 0, so this contract means overall PIP == 0 for all
# SNPs as well as zero retained credible sets.
is_susiex_explicit_no_cs <- function(exit_status, snp_file, cs_file, summary_file) {
  if (!identical(as.integer(exit_status), 0L) ||
      !all(file.exists(c(snp_file, cs_file, summary_file)))) return(FALSE)
  snp_lines <- trimws(readLines(snp_file, n = 2L, warn = FALSE))
  if (length(snp_lines) < 2L || !nzchar(snp_lines[1L]) ||
      grepl("PIP\\(CS", snp_lines[1L])) return(FALSE)
  cs_lines <- trimws(readLines(cs_file, warn = FALSE))
  summary_lines <- trimws(readLines(summary_file, warn = FALSE))
  any(cs_lines == "NULL") && any(summary_lines == "NULL")
}

# Native SuSiEx wrapper: actual ancestry-specific PLINK prefixes, and only native
# overall PIP/credible sets are returned. Its CS-level population marginals are not
# converted into invalid SNP-level configuration PIPs.
run_susiex <- function(zfile, avail, map, mode, work, locus, threads = 2) {
  if (!file.exists(SUSIEX_BIN) || !file.exists(PLINK_BIN)) return(NULL)
  dir.create(work, recursive = TRUE, showWarnings = FALSE)
  idx <- switch(mode,
    intersect = list(EUR = which(rowSums(avail) == 2L), AFR = which(rowSums(avail) == 2L)),
    union = list(EUR = which(avail[, "EUR"]), AFR = which(avail[, "AFR"])),
    stop("bad SuSiEx mode"))
  if (min(lengths(idx)) < 10L) return(NULL)

  f_eur <- file.path(work, "eur.tsv"); f_afr <- file.path(work, "afr.tsv")
  write_susiex_sst(f_eur, zfile$RSID[idx$EUR], map$chromosome[idx$EUR],
                   map$position[idx$EUR], zfile$zscore_1[idx$EUR], NUM_EU,
                   map$allele.1[idx$EUR], map$allele.2[idx$EUR])
  write_susiex_sst(f_afr, zfile$RSID[idx$AFR], map$chromosome[idx$AFR],
                   map$position[idx$AFR], zfile$zscore_2[idx$AFR], NUM_BB,
                   map$allele.1[idx$AFR], map$allele.2[idx$AFR])
  d <- native_locus_dir(locus)
  ref_e <- file.path(d, "plink", "EUR"); ref_a <- file.path(d, "plink", "AFR")
  out_name <- paste0("susiex_", mode)
  args <- c(sprintf("--sst_file=%s,%s", f_eur, f_afr),
            sprintf("--ref_file=%s,%s", ref_e, ref_a),
            sprintf("--ld_file=%s/eur,%s/afr", work, work),
            sprintf("--n_gwas=%d,%d", NUM_EU, NUM_BB),
            sprintf("--out_dir=%s/", work), sprintf("--out_name=%s", out_name),
            sprintf("--chr=%d", map$chromosome[1L]),
            sprintf("--bp=%d,%d", min(map$position), max(map$position)),
            "--chr_col=2,2", "--bp_col=3,3", "--snp_col=1,1",
            "--a1_col=5,5", "--a2_col=6,6", "--eff_col=9,9", "--se_col=10,10",
            "--pval_col=4,4", sprintf("--plink=%s", PLINK_BIN),
            "--keep-ambig=True", "--mult-step=True", "--maf=0.01", "--level=0.95",
            "--min_purity=0.5", "--pval_thresh=1e-05", "--tol=0.0001",
            "--n_sig=10", "--max_iter=200", sprintf("--threads=%d", threads))
  st <- suppressWarnings(system2(SUSIEX_BIN, args,
                                 stdout = file.path(work, "run.log"),
                                 stderr = file.path(work, "run.err")))
  snp_file <- file.path(work, paste0(out_name, ".snp"))
  cs_file <- file.path(work, paste0(out_name, ".cs"))
  summary_file <- file.path(work, paste0(out_name, ".summary"))
  parsed <- parse_susiex_snp(snp_file, zfile$RSID)
  if (is.null(parsed)) {
    # Accept zero PIPs only under the complete native no-CS output contract;
    # missing/malformed outputs and nonzero exits remain method failures.
    if (!is_susiex_explicit_no_cs(st, snp_file, cs_file, summary_file)) return(NULL)
    return(list(pip = rep(0, nrow(zfile)), pip_EUR = rep(NA_real_, nrow(zfile)),
                pip_AFR = rep(NA_real_, nrow(zfile)),
                pip_shared = rep(NA_real_, nrow(zfile)), cs = list(),
                exit_status = st, config_pip_available = FALSE,
                no_credible_sets = TRUE))
  }
  cs <- parse_susiex_cs(cs_file, zfile$RSID)
  list(pip = parsed$pip, pip_EUR = rep(NA_real_, nrow(zfile)),
       pip_AFR = rep(NA_real_, nrow(zfile)), pip_shared = rep(NA_real_, nrow(zfile)),
       cs = cs, exit_status = st, config_pip_available = FALSE)
}

native_cell_path <- function(locus, base_count_index, h2_index, n_assv_each) {
  file.path(NATIVE_SIM_OUT, "cells",
            sprintf("base_%d", NATIVE_BASE_COUNTS[base_count_index]),
            sprintf("h2_%d", h2_index), sprintf("assv_each_%d", n_assv_each),
            sprintf("locus_%03d.rds", locus))
}

native_multisusie_maf_path <- function(locus) {
  file.path(NATIVE_MULTISUSIE_INPUT,
            sprintf("locus_%03d_shared_maf.rds", as.integer(locus)))
}

native_multisusie_path <- function(locus, base_count_index, h2_index,
                                   n_assv_each) {
  file.path(NATIVE_MULTISUSIE_OUT,
            sprintf("base_%d", NATIVE_BASE_COUNTS[base_count_index]),
            sprintf("h2_%d", h2_index),
            sprintf("assv_each_%d", n_assv_each),
            sprintf("locus_%03d.rds", as.integer(locus)))
}

native_infuse_d187_path <- function(locus, base_count_index, h2_index,
                                    n_assv_each) {
  file.path(NATIVE_INFUSE_D187_OUT,
            sprintf("base_%d", NATIVE_BASE_COUNTS[base_count_index]),
            sprintf("h2_%d", h2_index),
            sprintf("assv_each_%d", n_assv_each),
            sprintf("locus_%03d.rds", as.integer(locus)))
}
