#' RAFE-Stein-Corrected Mean
#'
#' Positive-part James--Stein shrinkage of a forecast mean toward a target, in
#' the \eqn{\hat\Sigma^{-1/2}}-whitened space.
#'
#' Writing \eqn{\hat m := \hat\Sigma^{-1/2}\hat\mu} and \eqn{m_0 :=
#' \hat\Sigma^{-1/2}\mu_0} for the whitened forecast and target, the James--Stein
#' map of Salcher, Stöckl & Hanke is
#' \deqn{JS(z) = m_0 + \left(1 - \frac{N-2}{T\,\|z - m_0\|_2^2}\right)_{+}(z - m_0),}
#' and the corrected mean is \eqn{\tilde\mu = \hat\Sigma^{1/2} JS(\hat m)}.
#' Because \eqn{\|\hat m - m_0\|_2^2 = (\hat\mu - \mu_0)^\top \hat\Sigma^{-1}
#' (\hat\mu - \mu_0)}, this is algebraically a linear shrinkage in the *raw*
#' space whose intensity is computed in the *risk* metric:
#' \eqn{\tilde\mu = (1 - \kappa)\hat\mu + \kappa\mu_0}.
#'
#' Supply `T_obs` to use the James--Stein intensity, or `kappa` to fix the
#' intensity directly (which is what [sep_tune()] and [joint_trafe_tune()]
#' search over). Exactly one of the two is required.
#'
#' @param mu_hat Forecast mean.
#' @param Sigma_hat Forecast covariance. Must be a covariance, never a
#'   precision matrix; the function inverts it internally.
#' @param kappa Optional shrinkage intensity in `[0, 1]`, where 0 leaves
#'   `mu_hat` untouched and 1 collapses it onto the target. If `NULL`, the
#'   James--Stein intensity is computed from `T_obs`.
#' @param target Shrinkage target. One of `"grandmean"` (default; shrinks
#'   toward the cross-sectional mean of `mu_hat`) or `"zero"`.
#' @param T_obs Number of observations the forecast was estimated from. Needed
#'   only when `kappa` is `NULL`.
#'
#' @return A numeric vector of corrected means (same length as `mu_hat`), with
#'   the applied intensity attached as attribute `kappa`.
#'
#' @examples
#' set.seed(1)
#' R <- matrix(rnorm(60 * 5, mean = 0.01, sd = 0.05), 60, 5)
#' mu_hat    <- colMeans(R)
#' Sigma_hat <- stats::cov(R)
#'
#' # James-Stein intensity, implied by the sample size:
#' mu_rafe_stein(mu_hat, Sigma_hat, T_obs = 60)
#'
#' # Fixed intensity:
#' mu_rafe_stein(mu_hat, Sigma_hat, kappa = 0.5)
#'
#' @references
#'   Stöckl, S., Salcher, L., & Hanke, M. Post-Optimal Moment Correction for
#'   Mean-Variance Portfolios. Working paper.
#' @export
mu_rafe_stein <- function(mu_hat, Sigma_hat, kappa = NULL,
                          target = c("grandmean", "zero"), T_obs = NULL) {
  target <- match.arg(target)
  mu_hat <- .check_vec(mu_hat, "mu_hat")
  n <- length(mu_hat)
  Sigma_hat <- .check_mat(Sigma_hat, n, "Sigma_hat")

  mu_0 <- switch(target,
                 grandmean = rep(mean(mu_hat), n),
                 zero      = rep(0, n))

  if (is.null(kappa)) {
    if (is.null(T_obs)) {
      stop("Supply either `kappa` or `T_obs`. To tune `kappa` on data, use ",
           "`sep_tune()` or `joint_trafe_tune()`.", call. = FALSE)
    }
    if (!is.numeric(T_obs) || length(T_obs) != 1L || !is.finite(T_obs) ||
        T_obs <= 0) {
      stop("`T_obs` must be a single positive number.", call. = FALSE)
    }
    if (n < 3L) {
      # The James-Stein dominance result needs N >= 3; below that, shrinking
      # has no theoretical justification, so leave the forecast alone.
      kappa <- 0
    } else {
      Sinv <- .pd_inverse(Sigma_hat)
      d <- mu_hat - mu_0
      dist_sq <- as.numeric(crossprod(d, Sinv %*% d))
      kappa <- if (dist_sq <= 0) 1 else (n - 2) / (T_obs * dist_sq)
      kappa <- max(0, min(1, kappa))
    }
  } else {
    if (!is.numeric(kappa) || length(kappa) != 1L || !is.finite(kappa)) {
      stop("`kappa` must be a single finite numeric, or NULL.", call. = FALSE)
    }
    kappa <- max(0, min(1, kappa))
  }

  out <- (1 - kappa) * mu_hat + kappa * mu_0
  attr(out, "kappa") <- kappa
  out
}

#' C-RAFE-Tuned Eigenvalue Floor for Covariance
#'
#' Cleans `Sigma_hat` by raising small eigenvalues to a floor set relative to
#' the mean eigenvalue: with \eqn{\hat\Sigma = Q\,\mathrm{diag}(\hat\lambda)\,Q^\top},
#' \deqn{\tilde\Sigma(\tau) = Q\,\mathrm{diag}\bigl(\max(\hat\lambda_i,\ \tau\,\bar\lambda)\bigr)\,Q^\top,}
#' where \eqn{\bar\lambda} is the mean of the sample eigenvalues. The relative
#' parameterisation makes `tau_rel` scale-free.
#'
#' `tau_rel` is tuned on inner-validation C-RAFE by [sep_tune()] and
#' [joint_trafe_tune()]; this function applies a given floor.
#'
#' @param Sigma_hat Forecast covariance.
#' @param tau_rel Relative floor, typically in `[0, 0.5]`. A value of 0 returns
#'   `Sigma_hat` symmetrised and unchanged.
#'
#' @return A cleaned covariance matrix (same dimension as `Sigma_hat`), with
#'   the applied absolute floor attached as attribute `tau_abs`.
#'
#' @examples
#' set.seed(1)
#' R <- matrix(rnorm(30 * 5), 30, 5)
#' Sigma_hat <- stats::cov(R)
#'
#' eigen(Sigma_hat, only.values = TRUE)$values
#' eigen(sigma_crafe_floor(Sigma_hat, 0.3), only.values = TRUE)$values
#'
#' @inherit mu_rafe_stein references
#' @export
sigma_crafe_floor <- function(Sigma_hat, tau_rel = NULL) {
  Sigma_hat <- .check_mat(Sigma_hat, NULL, "Sigma_hat")
  if (is.null(tau_rel)) {
    stop("Supply `tau_rel`. To tune it on data, use `sep_tune()` or ",
         "`joint_trafe_tune()`.", call. = FALSE)
  }
  if (!is.numeric(tau_rel) || length(tau_rel) != 1L || !is.finite(tau_rel) ||
      tau_rel < 0) {
    stop("`tau_rel` must be a single finite non-negative numeric.",
         call. = FALSE)
  }
  eg <- eigen(.as_sym(Sigma_hat), symmetric = TRUE)
  d <- eg$values
  tau_abs <- tau_rel * mean(d)
  d_floored <- pmax(d, tau_abs)
  out <- eg$vectors %*% (d_floored * t(eg$vectors))
  out <- .as_sym(out)
  dimnames(out) <- dimnames(Sigma_hat)
  attr(out, "tau_abs") <- tau_abs
  out
}

# Split a training window into inner-training and inner-validation blocks and
# return the sample moments of each. Shared by both tuners.
#' @noRd
.inner_moments <- function(R_train, inner_split) {
  R_train <- as.matrix(R_train)
  if (!is.numeric(R_train) || !all(is.finite(R_train))) {
    stop("`R_train` must be a numeric matrix of finite returns.", call. = FALSE)
  }
  if (length(inner_split) != 2L || any(inner_split < 2)) {
    stop("`inner_split` must be two integers, each at least 2.", call. = FALSE)
  }
  if (sum(inner_split) > nrow(R_train)) {
    stop("`inner_split` needs ", sum(inner_split), " rows but `R_train` has ",
         nrow(R_train), ".", call. = FALSE)
  }
  i_in  <- seq_len(inner_split[1])
  i_val <- inner_split[1] + seq_len(inner_split[2])
  list(
    mu_in     = colMeans(R_train[i_in, , drop = FALSE]),
    Sigma_in  = stats::cov(R_train[i_in, , drop = FALSE]),
    mu_val    = colMeans(R_train[i_val, , drop = FALSE]),
    Sigma_val = stats::cov(R_train[i_val, , drop = FALSE]),
    mu_full   = colMeans(R_train),
    Sigma_full = stats::cov(R_train)
  )
}

#' Sequentially Tuned Moment Correction
#'
#' Tunes \eqn{\kappa} on inner-validation RAFE and \eqn{\tau} on
#' inner-validation C-RAFE *independently* — the bound's two channels are
#' separable, so each parameter is chosen against its own channel.
#'
#' The training window is split into an inner-training block and an
#' inner-validation block (40 + 20 months by default). Moments are
#' estimated on the inner-training block, the realised moments of the
#' inner-validation block play the role of \eqn{(\mu, \Sigma)} in the metrics,
#' and the selected \eqn{(\kappa, \tau)} are then applied to the moments of the
#' *full* training window.
#'
#' @param R_train Numeric matrix of training returns, observations in rows and
#'   assets in columns.
#' @param kappa_grid Grid of shrinkage intensities to search. The default
#'   matches the working paper's grid, 101 points on `[0, 1]`.
#' @param tau_grid Grid of relative eigenvalue floors to search. The default
#'   matches the working paper's grid, 101 points on `[0, 0.5]`. Floors above
#'   0.5 flatten the spectrum so aggressively that the cleaned covariance
#'   carries little information about the original, and are not admitted.
#' @param inner_split Length-2 vector: rows used for inner training and for
#'   inner validation.
#' @param target Shrinkage target passed to [mu_rafe_stein()].
#'
#' @return A list with components `mu_tilde`, `Sigma_tilde`, `kappa`,
#'   `tau_rel`, and the inner-validation loss profiles `losses_kappa` and
#'   `losses_tau`.
#'
#' @examples
#' set.seed(1)
#' R <- matrix(rnorm(80 * 6, mean = 0.008, sd = 0.05), 80, 6)
#' fit <- sep_tune(R)
#' fit$kappa
#' fit$tau_rel
#'
#' @seealso [joint_trafe_tune()] for the joint search.
#' @inherit mu_rafe_stein references
#' @export
sep_tune <- function(R_train,
                     kappa_grid = seq(0, 1, length.out = 101),
                     tau_grid = seq(0, 0.5, length.out = 101),
                     inner_split = c(40, 20),
                     target = c("grandmean", "zero")) {
  target <- match.arg(target)
  m <- .inner_moments(R_train, inner_split)

  losses_kappa <- vapply(kappa_grid, function(k) {
    mu_t <- mu_rafe_stein(m$mu_in, m$Sigma_in, kappa = k, target = target)
    compute_rafe(mu_t, m$mu_val, Sigma = m$Sigma_val)
  }, numeric(1))
  kappa_star <- kappa_grid[which.min(losses_kappa)]

  losses_tau <- vapply(tau_grid, function(t) {
    compute_crafe(m$Sigma_val, sigma_crafe_floor(m$Sigma_in, t))
  }, numeric(1))
  tau_star <- tau_grid[which.min(losses_tau)]

  list(
    mu_tilde     = mu_rafe_stein(m$mu_full, m$Sigma_full, kappa = kappa_star,
                                 target = target),
    Sigma_tilde  = sigma_crafe_floor(m$Sigma_full, tau_star),
    kappa        = kappa_star,
    tau_rel      = tau_star,
    losses_kappa = stats::setNames(losses_kappa, kappa_grid),
    losses_tau   = stats::setNames(losses_tau, tau_grid)
  )
}

#' Jointly T-RAFE-Tuned Moment Correction
#'
#' Tunes \eqn{(\kappa, \tau)} jointly on the full grid to minimise
#' inner-validation T-RAFE. Empirically close to indistinguishable from the
#' sequential search of [sep_tune()], which is the practical evidence that the
#' bound's two channels are separable.
#'
#' @inheritParams sep_tune
#'
#' @return A list with components `mu_tilde`, `Sigma_tilde`, `kappa`,
#'   `tau_rel`, and the full inner-validation loss surface `loss_grid`
#'   (kappa in rows, tau in columns).
#'
#' @examples
#' set.seed(1)
#' R <- matrix(rnorm(80 * 6, mean = 0.008, sd = 0.05), 80, 6)
#' fit <- joint_trafe_tune(R, kappa_grid = seq(0, 1, by = 0.1),
#'                            tau_grid = seq(0, 1, by = 0.2))
#' fit$kappa
#' fit$tau_rel
#'
#' @seealso [sep_tune()] for the sequential search.
#' @inherit mu_rafe_stein references
#' @export
joint_trafe_tune <- function(R_train,
                             kappa_grid = seq(0, 1, length.out = 101),
                             tau_grid = seq(0, 0.5, length.out = 101),
                             inner_split = c(40, 20),
                             target = c("grandmean", "zero")) {
  target <- match.arg(target)
  m <- .inner_moments(R_train, inner_split)

  # Pre-compute the two correction families once; the grid search only
  # recombines them.
  mu_by_kappa <- lapply(kappa_grid, function(k) {
    mu_rafe_stein(m$mu_in, m$Sigma_in, kappa = k, target = target)
  })
  Sigma_by_tau <- lapply(tau_grid, function(t) sigma_crafe_floor(m$Sigma_in, t))

  loss_grid <- matrix(NA_real_, length(kappa_grid), length(tau_grid),
                      dimnames = list(kappa = kappa_grid, tau = tau_grid))
  for (i in seq_along(kappa_grid)) {
    for (j in seq_along(tau_grid)) {
      loss_grid[i, j] <- as.numeric(
        compute_trafe(mu_by_kappa[[i]], m$mu_val, m$Sigma_val, Sigma_by_tau[[j]])
      )
    }
  }

  # which.min() indexes column-major; recover the (kappa, tau) cell from it.
  # Avoid which(arr.ind = TRUE) here: the dimnames make its column labels
  # "kappa"/"tau" rather than "row"/"col".
  idx <- which.min(loss_grid)
  kappa_star <- kappa_grid[((idx - 1L) %% length(kappa_grid)) + 1L]
  tau_star   <- tau_grid[((idx - 1L) %/% length(kappa_grid)) + 1L]

  list(
    mu_tilde    = mu_rafe_stein(m$mu_full, m$Sigma_full, kappa = kappa_star,
                                target = target),
    Sigma_tilde = sigma_crafe_floor(m$Sigma_full, tau_star),
    kappa       = kappa_star,
    tau_rel     = tau_star,
    loss_grid   = loss_grid
  )
}
