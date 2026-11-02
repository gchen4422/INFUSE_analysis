.infuse_project <- Sys.getenv("INFUSE_PROJECT", "/path/to/baseline_simulation")
.infuse_work <- Sys.getenv("INFUSE_REVISION", "/path/to/revision_output")
# Sensitivity simulation for ancestry-specific causal priors.
#
# This script intentionally retains the data-generating machinery in
# 3_shared_50_simulation_mesusie_inf.R. The changes are limited to:
#   1. appending EUR- and AFR-specific aligned effect-allele-frequency (EAF)
#      annotations;
#   2. constructing EUR, AFR, and shared causal priors;
#   3. sampling causal SNPs conditional on their ancestry configuration; and
#   4. saving the ancestry-specific truth for downstream evaluation.

suppressPackageStartupMessages({
  library(data.table)
  library(snpStats)
  library(mvtnorm)
  library(dplyr)
})

args <- commandArgs(trailingOnly = TRUE)

scenario <- if (length(args) >= 1L) args[[1L]] else "different_annotation"
base_d <- if (length(args) >= 2L) as.integer(args[[2L]]) else 10L
ld_block_arg <- if (length(args) >= 3L && nzchar(args[[3L]])) as.integer(args[[3L]]) else NA_integer_

valid_scenarios <- c(
  "baseline",
  "different_annotation",
  "different_weight",
  "combined"
)

if (!scenario %in% valid_scenarios) {
  stop(
    "scenario must be one of: ",
    paste(valid_scenarios, collapse = ", "),
    call. = FALSE
  )
}

if (!base_d %in% c(10L, 187L)) {
  stop("base_d must be 10 or 187.", call. = FALSE)
}

if (!is.na(ld_block_arg) && (ld_block_arg < 1L || ld_block_arg > 100L)) {
  stop("LD_BLOCK must be between 1 and 100.", call. = FALSE)
}

ld_blocks <- if (is.na(ld_block_arg)) seq_len(100L) else ld_block_arg

# Preserve the balanced baseline sample sizes.
num_EU <- 300000
num_BB <- 300000

project_dir <- .infuse_project
geno_dir_eur <- file.path(Sys.getenv("INFUSE_SIMULATION", dirname(.infuse_project)), "risk_loci_ld_eur")
LD_dir_eur <- geno_dir_eur
geno_dir_afr <- file.path(Sys.getenv("INFUSE_SIMULATION", dirname(.infuse_project)), "risk_loci_ld_afr")
LD_dir_afr <- geno_dir_afr
annotation_source_dir <- file.path(Sys.getenv("INFUSE_SIMULATION", dirname(.infuse_project)), "risk_loci_annot")

# EAF_EUR and EAF_AFR are appended to the existing panel.
input_d <- base_d + 2L
simulation_root <- file.path(.infuse_work, "ancestry_specific_eaf_prior_simulation", paste0("scenario_", scenario), paste0("base_d_", base_d, "_input_d_", input_d))
annot_dir <- file.path(simulation_root, "annotations")
truth_dir <- file.path(simulation_root, "truth")
dir.create(annot_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(truth_dir, recursive = TRUE, showWarnings = FALSE)

signal_annotation_names <- c(
  "H3K27ac_Hnisz_common",
  "DHS_Trynka_common",
  "H3K4me3_Trynka_common",
  "Conserved_LindbladToh_common",
  "non_synonymous_common"
)

low_dim_null_names <- c(
  "Coding_UCSC_common",
  "Coding_UCSC.flanking.500_common",
  "CTCF_Hoffman_common",
  "Enhancer_Andersson_common",
  "DHS_peaks_Trynka_common"
)

softmax_stable <- function(x) {
  x <- as.numeric(x)
  if (!length(x) || any(!is.finite(x))) {
    stop("Causal-prior score contains non-finite values.", call. = FALSE)
  }
  shifted <- x - max(x)
  ans <- exp(shifted)
  ans <- ans / sum(ans)
  if (any(!is.finite(ans)) || abs(sum(ans) - 1) > 1e-10) {
    stop("Failed to construct a normalized causal prior.", call. = FALSE)
  }
  ans
}

jensen_shannon_divergence <- function(p, q) {
  eps <- .Machine$double.eps
  p <- pmax(as.numeric(p), eps)
  q <- pmax(as.numeric(q), eps)
  p <- p / sum(p)
  q <- q / sum(q)
  m <- 0.5 * (p + q)
  0.5 * sum(p * log(p / m)) + 0.5 * sum(q * log(q / m))
}

calculate_eaf <- function(geno) {
  # For a SnpMatrix imported by read.plink(), numeric dosage is the count of
  # allele.2. Because both allele.1 and allele.2 are checked for exact EUR/AFR
  # agreement below, these frequencies refer to the same aligned allele.
  eaf <- colMeans(geno, na.rm = TRUE) / 2
  if (any(!is.finite(eaf)) || any(eaf < 0 | eaf > 1)) {
    stop("Invalid EAF values calculated from genotype dosages.", call. = FALSE)
  }
  eaf
}

process_geno <- function(geno_data) {
  processed <- apply(geno_data, 2, function(x) {
    x[is.na(x)] <- mean(x, na.rm = TRUE)
    x <- scale(x - mean(x))
    return(x)
  })
  processed <- as.matrix(processed)
  colnames(processed) <- colnames(geno_data)
  processed
}

sample_causal_by_configuration <- function(
    causal_config,
    snp_ids,
    prior_EUR,
    prior_AFR,
    prior_shared
) {
  selected_indices <- integer(0)

  for (causal_index in seq_along(causal_config)) {
    config <- causal_config[[causal_index]]
    sampling_prob <- switch(
      as.character(config),
      "1" = prior_EUR,
      "2" = prior_AFR,
      "3" = prior_shared,
      stop("Unknown causal configuration: ", config, call. = FALSE)
    )

    if (length(selected_indices)) {
      sampling_prob[selected_indices] <- 0
    }
    if (!is.finite(sum(sampling_prob)) || sum(sampling_prob) <= 0) {
      stop("No positive probability remains for causal-SNP sampling.", call. = FALSE)
    }
    sampling_prob <- sampling_prob / sum(sampling_prob)
    selected_indices <- c(
      selected_indices,
      sample.int(length(snp_ids), size = 1L, replace = FALSE, prob = sampling_prob)
    )
  }

  list(index = selected_indices, snp = snp_ids[selected_indices])
}

for (LD_BLOCK in ld_blocks) {
  # Preserve the original locus-specific seed construction so conditions can be
  # compared using matched random-number streams.
  set.seed(1128)
  seed_LOCIrator <- sample.int(1000, 100)[LD_BLOCK]
  set.seed(seed_LOCIrator)

  # Reading and processing EUR data.
  outfile_EU <- file.path(geno_dir_eur, paste0("loci_", LD_BLOCK))
  EU_plink <- read.plink(paste0(outfile_EU, ".bed"))
  EU_plink_geno_raw <- as(EU_plink$genotypes, "numeric")
  colnames(EU_plink_geno_raw) <- EU_plink$map$snp.name

  ld_EU_name <- file.path(LD_dir_eur, paste0("loci_", LD_BLOCK, ".ld"))
  EU_cov <- as.matrix(fread(ld_EU_name))

  # Reading and processing AFR data.
  outfile_BB <- file.path(geno_dir_afr, paste0("loci_", LD_BLOCK))
  BB_plink <- read.plink(paste0(outfile_BB, ".bed"))
  BB_plink$map$CHR_POS <- paste0(BB_plink$map$chromosome, ":", BB_plink$map$position)
  BB_plink_geno_raw <- as(BB_plink$genotypes, "numeric")
  colnames(BB_plink_geno_raw) <- BB_plink$map$snp.name

  ld_BB_name <- file.path(LD_dir_afr, paste0("loci_", LD_BLOCK, ".ld"))
  BB_cov <- as.matrix(fread(ld_BB_name))
  # Match the original simulation: the AFR LD file carries rsID headers, while
  # the genotype/annotation workflow uses chromosome:position identifiers.
  colnames(BB_cov) <- rownames(BB_cov) <- as.character(BB_plink$map$CHR_POS)

  snp_ids <- as.character(EU_plink$map$snp.name)
  afr_match <- match(snp_ids, as.character(BB_plink$map$snp.name))
  if (anyNA(afr_match)) {
    stop("EUR and AFR genotype panels do not contain the same SNP set at locus ", LD_BLOCK, call. = FALSE)
  }

  BB_plink_geno_raw <- BB_plink_geno_raw[, afr_match, drop = FALSE]
  BB_map_aligned <- BB_plink$map[afr_match, , drop = FALSE]
  colnames(BB_plink_geno_raw) <- snp_ids

  if (nrow(EU_cov) != length(snp_ids) || ncol(EU_cov) != length(snp_ids) ||
      nrow(BB_cov) != length(snp_ids) || ncol(BB_cov) != length(snp_ids)) {
    stop("LD and genotype dimensions differ at locus ", LD_BLOCK, call. = FALSE)
  }

  if (!identical(as.character(colnames(EU_cov)), snp_ids)) {
    stop("EUR LD SNP order differs from the EUR genotype order at locus ", LD_BLOCK, call. = FALSE)
  }

  bb_ld_match <- match(snp_ids, as.character(colnames(BB_cov)))
  if (anyNA(bb_ld_match)) {
    stop("AFR LD matrix does not contain the full SNP set at locus ", LD_BLOCK, call. = FALSE)
  }
  BB_cov <- BB_cov[bb_ld_match, bb_ld_match, drop = FALSE]
  colnames(BB_cov) <- rownames(BB_cov) <- snp_ids
  rownames(EU_cov) <- snp_ids

  allele_1_aligned <- as.character(EU_plink$map$allele.1) ==
    as.character(BB_map_aligned$allele.1)
  allele_2_aligned <- as.character(EU_plink$map$allele.2) ==
    as.character(BB_map_aligned$allele.2)
  if (any(!allele_1_aligned | !allele_2_aligned)) {
    stop("Allele mismatch between ancestries at locus ", LD_BLOCK, call. = FALSE)
  }

  eaf_EUR <- calculate_eaf(EU_plink_geno_raw)
  eaf_AFR <- calculate_eaf(BB_plink_geno_raw)

  EU_plink_geno <- process_geno(EU_plink_geno_raw)
  BB_plink_geno <- process_geno(BB_plink_geno_raw)

  # Functional annotations used in the original simulation.
  annotation_file_name <- file.path(annotation_source_dir, paste0("loci_", LD_BLOCK, ".annot"))
  annot_file <- fread(annotation_file_name)
  # The simulation inputs identify variants as CHR:BP, whereas the annotation
  # table's SNP column contains rsIDs. Match on genomic position so the added
  # EAF annotations and all existing annotations remain aligned to the LD and
  # genotype matrices.
  annot_file[, CHR_POS := paste0(CHR, ":", BP)]
  annot_match <- match(snp_ids, annot_file$CHR_POS)
  if (anyNA(annot_match)) {
    stop("Annotation table does not contain the full SNP set at locus ", LD_BLOCK, call. = FALSE)
  }
  annot_file <- annot_file[annot_match]

  all_annotation_names <- setdiff(
    colnames(annot_file),
    c("CHR", "BP", "SNP", "CM", "CHR_POS")
  )
  if (length(all_annotation_names) != 187L) {
    stop("Expected 187 baseline annotation columns but found ", length(all_annotation_names), call. = FALSE)
  }

  if (!all(signal_annotation_names %in% all_annotation_names)) {
    stop("One or more signal annotations are absent at locus ", LD_BLOCK, call. = FALSE)
  }

  base_annotation_names <- if (base_d == 10L) {
    c(signal_annotation_names, low_dim_null_names)
  } else {
    all_annotation_names
  }

  base_annotation_matrix <- as.matrix(annot_file[, ..base_annotation_names])
  storage.mode(base_annotation_matrix) <- "double"
  colnames(base_annotation_matrix) <- base_annotation_names

  candidate_annotation_matrix <- cbind(
    base_annotation_matrix,
    EAF_EUR = eaf_EUR,
    EAF_AFR = eaf_AFR
  )

  if (ncol(candidate_annotation_matrix) != input_d) {
    stop("Candidate annotation dimension does not equal input_d.", call. = FALSE)
  }

  annotation_out <- as.data.frame(candidate_annotation_matrix)
  write.table(
    annotation_out,
    file.path(annot_dir, paste0("Region_", LD_BLOCK, ".annot")),
    quote = FALSE,
    row.names = FALSE,
    col.names = TRUE,
    sep = " "
  )

  # Full coefficient vectors over the candidate annotation matrix. For the
  # different-annotation scenarios, the same conceptual EAF effect is applied
  # to the ancestry-matched EAF column.
  gamma_EUR <- setNames(rep(0, input_d), colnames(candidate_annotation_matrix))
  gamma_AFR <- setNames(rep(0, input_d), colnames(candidate_annotation_matrix))

  weight_EUR <- 2
  weight_AFR <- if (scenario %in% c("different_weight", "combined")) 1 else 2

  gamma_EUR[signal_annotation_names] <- weight_EUR
  gamma_AFR[signal_annotation_names] <- weight_AFR

  if (scenario %in% c("different_annotation", "combined")) {
    gamma_EUR["EAF_EUR"] <- weight_EUR
    gamma_AFR["EAF_AFR"] <- weight_AFR
  }

  eta_EUR <- drop(candidate_annotation_matrix %*% gamma_EUR)
  eta_AFR <- drop(candidate_annotation_matrix %*% gamma_AFR)

  true_prior_EUR <- softmax_stable(eta_EUR)
  true_prior_AFR <- softmax_stable(eta_AFR)
  true_prior_shared <- softmax_stable((eta_EUR + eta_AFR) / 2)

  prior_diagnostics <- data.frame(
    LD_BLOCK = LD_BLOCK,
    scenario = scenario,
    base_d = base_d,
    input_d = input_d,
    gamma_EUR = weight_EUR,
    gamma_AFR = weight_AFR,
    pearson_prior_cor = cor(true_prior_EUR, true_prior_AFR, method = "pearson"),
    spearman_prior_cor = cor(true_prior_EUR, true_prior_AFR, method = "spearman"),
    js_divergence = jensen_shannon_divergence(true_prior_EUR, true_prior_AFR),
    mean_abs_eaf_difference = mean(abs(eaf_EUR - eaf_AFR)),
    max_abs_eaf_difference = max(abs(eaf_EUR - eaf_AFR))
  )

  saveRDS(
    list(
      LD_BLOCK = LD_BLOCK,
      scenario = scenario,
      base_d = base_d,
      input_d = input_d,
      seed = seed_LOCIrator,
      snp_ids = snp_ids,
      allele_1 = as.character(EU_plink$map$allele.1),
      allele_2 = as.character(EU_plink$map$allele.2),
      eaf_counted_allele = "allele_2",
      candidate_annotation_names = colnames(candidate_annotation_matrix),
      signal_annotation_names = signal_annotation_names,
      gamma_EUR = gamma_EUR,
      gamma_AFR = gamma_AFR,
      eta_EUR = eta_EUR,
      eta_AFR = eta_AFR,
      true_prior_EUR = true_prior_EUR,
      true_prior_AFR = true_prior_AFR,
      true_prior_shared = true_prior_shared,
      eaf_EUR = eaf_EUR,
      eaf_AFR = eaf_AFR,
      prior_diagnostics = prior_diagnostics
    ),
    file.path(truth_dir, paste0("Region_", LD_BLOCK, "_prior_truth.rds"))
  )

  fwrite(
    prior_diagnostics,
    file.path(truth_dir, paste0("Region_", LD_BLOCK, "_prior_diagnostics.tsv")),
    sep = "\t"
  )

  # Non-zero infinitesimal effects remain disabled, as in the original code.
  h2_small_EU <- 0 / nrow(EU_cov)
  h2_small_BB <- 0 / nrow(BB_cov)

  for (num_causal in seq_len(3L)) {
    wrk_dir <- file.path(simulation_root, paste0("causal_num_", num_causal))
    data_dir <- file.path(wrk_dir, "summary_data")
    causal_truth_dir <- file.path(wrk_dir, "truth")
    result_dir <- file.path(wrk_dir, "result")
    out_dir <- file.path(wrk_dir, "out")
    dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)
    dir.create(causal_truth_dir, recursive = TRUE, showWarnings = FALSE)
    dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

    num_causal_SNP <- c(2L, 3L, 5L)[num_causal]
    causal_SNP_config <- sample(
      c(1L, 2L, 3L),
      num_causal_SNP,
      replace = TRUE,
      prob = c(0.25, 0.25, 0.5)
    )

    causal_draw <- sample_causal_by_configuration(
      causal_config = causal_SNP_config,
      snp_ids = snp_ids,
      prior_EUR = true_prior_EUR,
      prior_AFR = true_prior_AFR,
      prior_shared = true_prior_shared
    )
    causal_SNP_list <- causal_draw$snp

    message(
      "scenario=", scenario,
      " base_d=", base_d,
      " input_d=", input_d,
      " LD_BLOCK=", LD_BLOCK,
      " num_causal=", num_causal_SNP
    )

    for (h2_num in seq_len(2L)) {
      h2 <- c(1e-4, 2e-4)[h2_num]
      signal <- beta_EUR_all <- beta_AFR_all <- rep(0, nrow(EU_cov))

      for (causal_index in seq_len(num_causal_SNP)) {
        causal_snp_index <- causal_draw$index[causal_index]

        if (causal_SNP_config[causal_index] == 1L) {
          EU_plink_causal <- matrix(EU_plink_geno[, causal_snp_index], ncol = 1)
          beta_EUR_causal <- rnorm(1, 0, 1)
          beta_EUR_causal_scale <- beta_EUR_causal * sqrt(
            h2 / drop(var(EU_plink_causal * beta_EUR_causal))
          )
          beta_EUR_all[causal_snp_index] <- beta_EUR_causal_scale
          signal[causal_snp_index] <- 1

        } else if (causal_SNP_config[causal_index] == 2L) {
          BB_plink_causal <- matrix(BB_plink_geno[, causal_snp_index], ncol = 1)
          beta_AFR_causal <- rnorm(1, 0, 1)
          beta_AFR_causal_scale <- beta_AFR_causal * sqrt(
            h2 / drop(var(BB_plink_causal * beta_AFR_causal))
          )
          beta_AFR_all[causal_snp_index] <- beta_AFR_causal_scale
          signal[causal_snp_index] <- 2

        } else if (causal_SNP_config[causal_index] == 3L) {
          EU_plink_causal <- matrix(EU_plink_geno[, causal_snp_index], ncol = 1)
          BB_plink_causal <- matrix(BB_plink_geno[, causal_snp_index], ncol = 1)
          beta_shared <- rmvnorm(
            1,
            mean = rep(0, 2),
            sigma = matrix(c(1, 0.8, 0.8, 1), ncol = 2, nrow = 2)
          )
          beta_EUR_causal_scale <- beta_shared[1] * sqrt(
            h2 / drop(var(EU_plink_causal * beta_shared[1]))
          )
          beta_AFR_causal_scale <- beta_shared[2] * sqrt(
            h2 / drop(var(BB_plink_causal * beta_shared[2]))
          )
          beta_EUR_all[causal_snp_index] <- beta_EUR_causal_scale
          beta_AFR_all[causal_snp_index] <- beta_AFR_causal_scale
          signal[causal_snp_index] <- 3
        }
      }

      beta_EUR_small <- rnorm(nrow(EU_cov), 0, h2_small_EU)
      beta_EUR_marginal <- as.matrix(EU_cov) %*% (beta_EUR_all + beta_EUR_small)

      beta_AFR_small <- rnorm(nrow(BB_cov), 0, h2_small_BB)
      beta_AFR_marginal <- as.matrix(BB_cov) %*% (beta_AFR_all + beta_AFR_small)

      h2_causal_EUR <- drop(var(EU_plink_geno %*% beta_EUR_all))
      h2_causal_AFR <- drop(var(BB_plink_geno %*% beta_AFR_all))

      y_null_EU <- rnorm(
        nrow(EU_plink_geno),
        0,
        sqrt(1 - h2_causal_EUR - h2_small_EU * nrow(EU_cov))
      )
      y_null_EU <- y_null_EU - mean(y_null_EU)

      y_null_BB <- rnorm(
        nrow(BB_plink_geno),
        0,
        sqrt(1 - h2_causal_AFR - h2_small_BB * nrow(EU_cov))
      )
      y_null_BB <- y_null_BB - mean(y_null_BB)

      err_beta_EU <- t(EU_plink_geno) %*% y_null_EU / nrow(EU_plink_geno)
      err_beta_BB <- t(BB_plink_geno) %*% y_null_BB / nrow(BB_plink_geno)

      err_beta_EU_scale <- sqrt(nrow(EU_plink_geno) / num_EU) * err_beta_EU
      err_beta_BB_scale <- sqrt(nrow(BB_plink_geno) / num_BB) * err_beta_BB

      z_EU <- (beta_EUR_marginal + err_beta_EU_scale) * sqrt(num_EU)
      z_BB <- (beta_AFR_marginal + err_beta_BB_scale) * sqrt(num_BB)

      z_file <- data.frame(
        CHR = as.numeric(EU_plink$map$chr[1]),
        POS = EU_plink$map$position,
        RSID = snp_ids,
        zscore_1 = z_EU,
        zscore_2 = z_BB,
        Signal = signal,
        N_1 = num_EU,
        N_2 = num_BB
      )

      output_stem <- paste0(
        "CAUSAL_", num_causal,
        "_LOCI_", LD_BLOCK,
        "_h2_", h2_num
      )
      write.table(
        z_file,
        file.path(data_dir, output_stem),
        col.names = TRUE,
        row.names = FALSE,
        quote = FALSE,
        sep = " "
      )

      saveRDS(
        list(
          scenario = scenario,
          base_d = base_d,
          input_d = input_d,
          LD_BLOCK = LD_BLOCK,
          num_causal_index = num_causal,
          num_causal_SNP = num_causal_SNP,
          h2_index = h2_num,
          h2 = h2,
          causal_SNP_config = causal_SNP_config,
          causal_SNP_index = causal_draw$index,
          causal_SNP_list = causal_SNP_list,
          beta_EUR = beta_EUR_all[causal_draw$index],
          beta_AFR = beta_AFR_all[causal_draw$index]
        ),
        file.path(causal_truth_dir, paste0(output_stem, "_truth.rds"))
      )
    }
  }
}
