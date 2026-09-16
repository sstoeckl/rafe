#' Risk-Adjusted Forecast Error (RAFE)
#'
#' Computes \eqn{\mathrm{RAFE} = \sqrt{(\hat\mu - \mu)^\top \Sigma^{-1} (\hat\mu - \mu)}},
#' the risk-adjusted distance between a forecast mean and the realised mean.
#'
#' The metric is a Mahalanobis distance in the population-risk metric, with no
#' \eqn{1/N} normalisation. Consequently, at \eqn{\Sigma = I} it reduces to the
#' Euclidean norm of the forecast error, which is \eqn{\sqrt{N}} times the
#' conventional RMSE — see the examples.
#'
#' @param mu_hat Numeric vector of forecast means.
#' @param mu Numeric vector of realised means (same length as `mu_hat`).
#' @param Sigma Realised covariance matrix. The argument order matches the
#'   paper's replication code, so `compute_rafe(mu_hat, mu, Sigma)` does the
#'   same thing here as it does there.
#' @param Sigma_inv Optional precision matrix, supplied instead of `Sigma` when
#'   it is already available (it avoids an inversion). Must be a *precision*,
#'   not a covariance.
#'
#' @return A length-1 numeric.
#'
#' @examples
#' set.seed(1)
#' mu     <- rnorm(5)
#' mu_hat <- mu + rnorm(5, sd = 0.1)
#'
#' compute_rafe(mu_hat, mu, Sigma = diag(5))
#'
#' # At Sigma = I the metric is sqrt(N) times RMSE:
#' sqrt(5) * sqrt(mean((mu_hat - mu)^2))
#'
#' @references
#'   Salcher, T., Stöckl, S., & Hanke, M. (2026). Lost in Translation?
#'   Risk-Adjusting RMSE for Economic Forecast Performance.
#'   *Journal of Forecasting*. doi:10.1002/for.70134
#' @export
compute_rafe <- function(mu_hat, mu, Sigma = NULL, Sigma_inv = NULL) {
  mu_hat <- .check_vec(mu_hat, "mu_hat")
  mu     <- .check_vec(mu, "mu")
  if (length(mu_hat) != length(mu)) {
    stop("`mu_hat` and `mu` must have the same length; got ",
         length(mu_hat), " and ", length(mu), ".", call. = FALSE)
  }
  Sinv <- .resolve_precision(Sigma_inv, Sigma, length(mu))
  d <- mu_hat - mu
  quad <- as.numeric(crossprod(d, Sinv %*% d))
  # A negative quadratic form can only arise from numerical noise at a
  # near-singular Sigma; clamp rather than return NaN.
  sqrt(max(quad, 0))
}

#' Covariance Forecast Error (C-RAFE)
#'
#' Computes the operator-norm precision distortion
#' \eqn{\mathrm{C\text{-}RAFE} = \|\Sigma^{1/2} \hat\Sigma^{-1} \Sigma^{1/2} - I\|_2}.
#'
#' The metric is zero if and only if \eqn{\hat\Sigma = \Sigma}, and is invariant
#' to the common scaling of both arguments.
#'
#' @param Sigma Realised covariance matrix.
#' @param Sigma_hat Forecast covariance matrix.
#'
#' @return A length-1 numeric.
#'
#' @section Argument order:
#' The realised covariance comes **first**, matching `compute_crafe()` in the
#' replication code of Salcher, Stöckl & Hanke (2026) so that paper code runs
#' unchanged against this package. The metric is not symmetric in its
#' arguments, so the order matters.
#'
#' @examples
#' Sigma <- diag(c(1, 2, 3))
#' compute_crafe(Sigma, Sigma)          # exactly 0
#' compute_crafe(Sigma, 2 * Sigma)      # forecast covariance twice too large
#'
#' # Under the paper's lambda-contamination scheme, Sigma_hat = (1-lambda)^2 Sigma
#' # and C-RAFE has a closed form:
#' lambda <- 0.05
#' c(computed  = compute_crafe(Sigma, (1 - lambda)^2 * Sigma),
#'   closedform = abs(1 / (1 - lambda)^2 - 1))
#'
#' @inherit compute_rafe references
#' @export
compute_crafe <- function(Sigma, Sigma_hat) {
  Sigma     <- .check_mat(Sigma, NULL, "Sigma")
  Sigma_hat <- .check_mat(Sigma_hat, nrow(Sigma), "Sigma_hat")
  Ssq  <- .sqrtm_sym(Sigma)
  Sinv <- .pd_inverse(Sigma_hat)
  M <- Ssq %*% Sinv %*% Ssq - diag(1, nrow(Sigma))
  .spec_norm(M)
}

#' Total RAFE (T-RAFE) — Sharpe-Gap Upper Bound
#'
#' Combines RAFE and C-RAFE into the upper bound on the Sharpe-ratio gap
#' \eqn{\Delta = SR^{*} - SR(\hat w)} of a plug-in mean-variance portfolio:
#' \eqn{\Delta \le c \cdot \mathrm{RAFE} + SR^{*} \cdot \mathrm{C\text{-}RAFE}}.
#'
#' Salcher, Stöckl & Hanke (2026) prove the bound with a *data-dependent*
#' constant: \eqn{c = 1} when \eqn{\mathrm{RAFE} \le SR^{*}}, and \eqn{c = 2}
#' otherwise. This is the default (`c = NULL`). Passing a fixed numeric `c`
#' overrides the rule and is intended for diagnostics only — in particular,
#' `c = 1` does **not** yield a valid upper bound in the regime
#' \eqn{\mathrm{RAFE} > SR^{*}}.
#'
#' @inheritParams compute_rafe
#' @inheritParams compute_crafe
#' @param c Scalar weighting for RAFE. If `NULL` (default), the theorem's rule
#'   is applied: 1 if `RAFE <= SR_star`, else 2. Note that the Part 1
#'   replication code computes `RAFE + |SR*| * C-RAFE`, i.e. a fixed `c = 1`;
#'   pass `c = 1` to reproduce those numbers exactly.
#' @param SR_star Optional oracle Sharpe ratio. If `NULL`, computed as
#'   \eqn{\sqrt{\mu^\top \Sigma^{-1} \mu}} from `mu` and `Sigma_inv` (or
#'   `Sigma`).
#'
#' @return A length-1 numeric, carrying the components as attributes `rafe`,
#'   `crafe`, `c` and `SR_star`.
#'
#' @examples
#' set.seed(1)
#' mu        <- rnorm(5) / 10
#' Sigma     <- diag(5)
#' mu_hat    <- mu + rnorm(5, sd = 0.02)
#' Sigma_hat <- Sigma + diag(5) * 0.05
#'
#' tr <- compute_trafe(mu_hat, mu, Sigma, Sigma_hat)
#' tr
#' attributes(tr)[c("rafe", "crafe", "c", "SR_star")]
#'
#' @inherit compute_rafe references
#' @export
compute_trafe <- function(mu_hat, mu, Sigma, Sigma_hat, c = NULL,
                          SR_star = NULL) {
  mu_hat <- .check_vec(mu_hat, "mu_hat")
  mu     <- .check_vec(mu, "mu")
  if (length(mu_hat) != length(mu)) {
    stop("`mu_hat` and `mu` must have the same length; got ",
         length(mu_hat), " and ", length(mu), ".", call. = FALSE)
  }
  n <- length(mu)
  Sigma     <- .check_mat(Sigma, n, "Sigma")
  Sigma_hat <- .check_mat(Sigma_hat, n, "Sigma_hat")

  Sinv <- .pd_inverse(Sigma)
  rafe  <- compute_rafe(mu_hat, mu, Sigma_inv = Sinv)
  crafe <- compute_crafe(Sigma, Sigma_hat)

  if (is.null(SR_star)) {
    SR_star <- sqrt(max(as.numeric(crossprod(mu, Sinv %*% mu)), 0))
  } else {
    if (!is.numeric(SR_star) || length(SR_star) != 1L || !is.finite(SR_star)) {
      stop("`SR_star` must be a single finite numeric, or NULL.", call. = FALSE)
    }
  }

  if (is.null(c)) {
    c <- if (rafe <= SR_star) 1 else 2
  } else if (!is.numeric(c) || length(c) != 1L || !is.finite(c)) {
    stop("`c` must be a single finite numeric, or NULL.", call. = FALSE)
  }

  out <- c * rafe + SR_star * crafe
  attr(out, "rafe")    <- rafe
  attr(out, "crafe")   <- crafe
  attr(out, "c")       <- c
  attr(out, "SR_star") <- SR_star
  out
}
