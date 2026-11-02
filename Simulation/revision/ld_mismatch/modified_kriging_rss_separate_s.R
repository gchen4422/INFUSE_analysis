################################################################################
# Minimal ancestry-specific-s extension of the original project code:
#   MESuSiE_original/modified_kriging_rss.R
#
# The calculation follows estimate_s_rss_joint() and kriging_rss_joint() as
# closely as possible. The substantive model change is only
#
#   shared:   blockdiag((1-s)  R1 + s I,  (1-s)  R2 + s I)
#   separate: blockdiag((1-s1) R1 + s1 I, (1-s2) R2 + s2 I).
#
# The finite-sample z shrinkage, conditional moments, mixture-of-normal residual
# fit, allele-flip likelihood ratio, joint LR sum and flag thresholds are retained.
################################################################################

estimate_s_rss_joint_separate <- function(z1, z2,
                                          R1, R2,
                                          n1 = NULL, n2 = NULL,
                                          r_tol = 1e-8) {
  #-- 1) finite-sample shrinkage on each z (same as original)
  shrink <- function(z, n) {
    if (!is.null(n)) {
      if (n <= 1) stop("n must be >1")
      sigma2 <- (n - 1) / (z^2 + n - 2)
      z <- sqrt(sigma2) * z
    }
    z[is.na(z)] <- 0
    z
  }
  z1 <- shrink(z1, n1)
  z2 <- shrink(z2, n2)

  #-- 2) fit the same null-MLE objective separately in each ancestry
  estimate_one <- function(z, R) {
    if (is.null(attr(R, "eigen")))
      attr(R, "eigen") <- eigen(R, symmetric = TRUE)
    eig <- attr(R, "eigen")
    if (any(eig$values < -r_tol))
      warning("R not PSD: flooring small negative eigenvalues")
    d <- eig$values
    d[d < r_tol] <- 0
    ztv <- crossprod(z, eig$vectors)

    neglog <- function(s) {
      tmp <- (1 - s) * d + s
      0.5 * sum(log(tmp)) + 0.5 * sum((ztv^2) / tmp)
    }
    stats::optim(0.5, fn = neglog, method = "Brent",
                 lower = 0, upper = 1)$par
  }

  c(s1 = estimate_one(z1, R1), s2 = estimate_one(z2, R2))
}


#' Joint kriging_rss for two ancestries with ancestry-specific s1 and s2
#'
#' Minimal modification of kriging_rss_joint(): the output fields and flagging
#' rules are retained and fitted s values plus explicit flag sets are returned.
kriging_rss_joint_separate_s <- function(z1, z2, R1, R2,
                                         n1 = NULL, n2 = NULL,
                                         r_tol = 1e-8,
                                         s1 = NULL, s2 = NULL) {
  if (length(z1) != length(z2))
    stop("z1 and z2 must contain the same ordered variants")
  if (nrow(R1) != length(z1) || ncol(R1) != length(z1) ||
      nrow(R2) != length(z2) || ncol(R2) != length(z2))
    stop("LD dimensions must match their ancestry-specific z-score vectors")

  #-- 1) finite-sample shrinkage & NA-impute (same as original)
  shrink <- function(z, n) {
    if (!is.null(n)) {
      if (n <= 1) stop("n must be > 1")
      sigma2 <- (n - 1) / (z^2 + n - 2)
      z <- sqrt(sigma2) * z
    }
    z[is.na(z)] <- 0
    z
  }
  z1 <- shrink(z1, n1)
  z2 <- shrink(z2, n2)

  #-- 2) stack z (same as original)
  z <- c(z1, z2)

  #-- 3) estimate s1 and s2 if missing. z has already received the same
  #      finite-sample shrinkage, so pass n=NULL to avoid shrinking twice.
  if (is.null(s1) || is.null(s2)) {
    s_hat <- estimate_s_rss_joint_separate(z1, z2, R1, R2,
                                            n1 = NULL, n2 = NULL,
                                            r_tol = r_tol)
    if (is.null(s1)) s1 <- unname(s_hat["s1"])
    if (is.null(s2)) s2 <- unname(s_hat["s2"])
  }
  if (s1 < 0 || s2 < 0) stop("s1 and s2 must be non-negative")
  if (s1 > 1 || s2 > 1) stop("s1 and s2 must not exceed 1")

  #-- 4/5) the only substantive model change: apply s1 and s2 to the
  #         eigenvalues of their corresponding unregularized LD blocks. This
  #         exactly preserves the original order of operations: floor the LD
  #         eigenvalues first, then add s. In particular, when 0 < s < r_tol,
  #         a null LD direction has precision 1/s; it must not be discarded by
  #         flooring the already-regularized covariance eigenvalue.
  precision_one <- function(R, s) {
    if (is.null(attr(R, "eigen")))
      attr(R, "eigen") <- eigen(R, symmetric = TRUE)
    eig <- attr(R, "eigen")
    d <- eig$values
    if (any(d < -r_tol))
      warning("R not PSD: flooring small negative eigenvalues")
    d[d < r_tol] <- 0
    dinv <- 1 / ((1 - s) * d + s)
    dinv[is.infinite(dinv)] <- 0
    eig$vectors %*% (t(eig$vectors) * dinv)
  }
  precision <- as.matrix(Matrix::bdiag(precision_one(R1, s1),
                                       precision_one(R2, s2)))

  #-- 6) conditional means and variances (same loop as original)
  p_total <- length(z)
  condmean <- numeric(p_total)
  condvar <- numeric(p_total)
  for (i in seq_len(p_total)) {
    condmean[i] <- -(1 / precision[i, i]) * precision[i, -i] %*% z[-i]
    condvar[i] <- 1 / precision[i, i]
  }
  z_std_diff <- (z - condmean) / sqrt(condvar)

  #-- 7) fit the same mixture of normal residual scales
  a_min <- 0.8
  a_max <- if (max(z_std_diff^2) < 1) 2 else 2 * sqrt(max(z_std_diff^2))
  npoint <- ceiling(log2(a_max / a_min) / log2(1.05))
  a_grid <- 1.05^(seq(-npoint, 0)) * a_max

  sd_mtx <- outer(sqrt(condvar), a_grid)
  matrix_llik <- stats::dnorm(z - condmean, sd = sd_mtx, log = TRUE)
  lfactors <- apply(matrix_llik, 1, max)
  matrix_llik <- matrix_llik - lfactors
  w <- mixsqp::mixsqp(matrix_llik, log = TRUE,
                      control = list(verbose = FALSE))$x

  #-- 8) compute the same allele-flip log likelihood ratios
  logl0mix <- drop(log(exp(matrix_llik) %*% (w + 1e-15))) + lfactors
  matrix_llik <- stats::dnorm(z + condmean, sd = sd_mtx, log = TRUE)
  lfactors <- apply(matrix_llik, 1, max)
  matrix_llik <- matrix_llik - lfactors
  logl1mix <- drop(log(exp(matrix_llik) %*% (w + 1e-15))) + lfactors
  logLRmix <- logl1mix - logl0mix

  #-- 9) split into LR1/LR2 and form the same summed joint LR
  p1 <- length(z1)
  logLR1 <- logLRmix[seq_len(p1)]
  logLR2 <- logLRmix[p1 + seq_len(p1)]
  logLR_joint <- logLR1 + logLR2

  ancestry <- rep(c("anc1", "anc2"), times = c(length(z1), length(z2)))
  res <- data.frame(
    ancestry = ancestry,
    z = z,
    condmean = condmean,
    condvar = condvar,
    z_std_diff = z_std_diff,
    logLR = logLRmix,
    logLR_joint = rep(logLR_joint, 2)
  )

  p <- ggplot2::ggplot(res, ggplot2::aes(x = condmean, y = z,
                                          color = ancestry)) +
    ggplot2::geom_point() +
    ggplot2::geom_abline(intercept = 0, slope = 1) +
    ggplot2::theme_bw() +
    ggplot2::labs(x = "Expected z", y = "Observed z")

  #-- 10) retain the original thresholds
  idx_flag1 <- which(logLR1 > 2 & abs(z_std_diff[seq_len(p1)]) > 4)
  idx_flag2 <- which(logLR2 > 2 &
                       abs(z_std_diff[p1 + seq_len(p1)]) > 4)
  # The original code implicitly recycles the length-p joint LR over the stacked
  # length-2p residual vector. Make the identical rule explicit.
  idx_flag_joint_stacked <- which(rep(logLR_joint, 2) > 2 &
                                    abs(z_std_diff) > 4)
  idx_flag_joint <- sort(unique(((idx_flag_joint_stacked - 1L) %% p1) + 1L))
  joint_flags <- sort(unique(c(idx_flag1, idx_flag2, idx_flag_joint)))

  if (length(idx_flag1) > 0) {
    p <- p + ggplot2::geom_point(data = res[idx_flag1, ],
                                 ggplot2::aes(x = condmean, y = z), color = "red")
  }
  if (length(idx_flag2) > 0) {
    p <- p + ggplot2::geom_point(data = res[p1 + idx_flag2, ],
                                 ggplot2::aes(x = condmean, y = z), color = "red")
  }
  if (length(idx_flag_joint_stacked) > 0) {
    p <- p + ggplot2::geom_point(data = res[idx_flag_joint_stacked, ],
                                 ggplot2::aes(x = condmean, y = z), color = "red")
  }

  list(
    plot = p,
    conditional_dist = res,
    s = c(s1 = s1, s2 = s2),
    flags = list(EUR = idx_flag1, AFR = idx_flag2,
                 joint_only = idx_flag_joint, joint = joint_flags)
  )
}
