#' Risk-Adjusted Forecast Error (RAFE)
#'
#' Computes \eqn{\mathrm{RAFE} = \sqrt{(\hat\mu - \mu)^\top \Sigma^{-1} (\hat\mu - \mu)}},
#' the risk-adjusted distance between a forecast mean and the realised mean.
#'
#' @param mu_hat Numeric vector of forecast means.
#' @param mu Numeric vector of realised means (same length as `mu_hat`).
#' @param Sigma_inv Optional precision matrix (inverse covariance). If `NULL`,
#'   `Sigma` must be provided.
#' @param Sigma Optional covariance matrix (used if `Sigma_inv` is `NULL`).
#'
#' @return A length-1 numeric.
#'
#' @references
#'   Salcher, T., Stöckl, S., & Hanke, M. (2026). Lost in Translation?
#'   Risk-Adjusting RMSE for Economic Forecast Performance.
#'   *Journal of Forecasting*. doi:10.1002/for.70134
#' @export
compute_rafe <- function(mu_hat, mu, Sigma_inv = NULL, Sigma = NULL) {
  # TODO: implement.
  # Matches paper §3 definition; falls back to solve(Sigma) if Sigma_inv = NULL.
  stop("TODO: implement compute_rafe()")
}

#' Covariance Forecast Error (C-RAFE)
#'
#' Computes the operator-norm precision distortion
#' \eqn{\mathrm{C\text{-}RAFE} = \|\Sigma^{1/2} \hat\Sigma^{-1} \Sigma^{1/2} - I\|_2}.
#'
#' @param Sigma_hat Forecast covariance matrix.
#' @param Sigma Realised covariance matrix.
#'
#' @return A length-1 numeric.
#'
#' @inherit compute_rafe references
#' @export
compute_crafe <- function(Sigma_hat, Sigma) {
  # TODO: implement.
  stop("TODO: implement compute_crafe()")
}

#' Total RAFE (T-RAFE) — Sharpe-Gap Upper Bound
#'
#' Combines RAFE and C-RAFE into the upper bound on the Sharpe-ratio gap
#' of a plug-in mean-variance portfolio:
#' \eqn{\Delta \le c \cdot \mathrm{RAFE} + SR^{*} \cdot \mathrm{C\text{-}RAFE}}.
#'
#' @inheritParams compute_rafe
#' @inheritParams compute_crafe
#' @param c Scalar weighting for RAFE; defaults to 1. (The paper proves
#'   \eqn{c \in \{1, 2\}}.)
#' @param SR_star Optional oracle Sharpe ratio. If `NULL`, computed from
#'   `mu` and `Sigma_inv` (or `Sigma`).
#'
#' @return A length-1 numeric.
#'
#' @inherit compute_rafe references
#' @export
compute_trafe <- function(mu_hat, mu, Sigma_hat, Sigma, c = 1, SR_star = NULL) {
  # TODO: implement. Compose RAFE + SR_star * C-RAFE; compute SR_star if NULL.
  stop("TODO: implement compute_trafe()")
}
