#' RAFE-Stein-Corrected Mean
#'
#' Stein-style shrinkage of a forecast mean toward a target, in the
#' \eqn{\hat\Sigma^{-1/2}}-whitened space. The shrinkage intensity \eqn{\kappa}
#' can be passed explicitly or tuned on inner-validation RAFE.
#'
#' @param mu_hat Forecast mean.
#' @param Sigma_hat Forecast covariance.
#' @param kappa Optional shrinkage intensity in [0, 1]. If `NULL`, tuned on
#'   inner-validation RAFE.
#' @param target Shrinkage target. One of `"grandmean"` (default; shrinks
#'   toward the cross-sectional mean) or `"zero"`.
#'
#' @return A numeric vector of corrected means (same length as `mu_hat`).
#'
#' @references
#'   Stöckl, S., Salcher, T., & Hanke, M. (2026/27, working). Post-Optimal
#'   Moment Correction for Mean-Variance Portfolios.
#' @export
mu_rafe_stein <- function(mu_hat, Sigma_hat, kappa = NULL, target = c("grandmean", "zero")) {
  target <- match.arg(target)
  # TODO: implement. Matches paper §3 (jsr / mu_kappa).
  stop("TODO: implement mu_rafe_stein()")
}

#' C-RAFE-Tuned Eigenvalue Floor for Covariance
#'
#' Cleans `Sigma_hat` by raising small eigenvalues to a relative floor.
#' The floor \eqn{\tau} can be passed explicitly or tuned on inner-validation
#' C-RAFE.
#'
#' @param Sigma_hat Forecast covariance.
#' @param tau_rel Optional relative floor in [0, 0.5]. If `NULL`, tuned on
#'   inner-validation C-RAFE.
#'
#' @return A cleaned covariance matrix (same dimension as `Sigma_hat`).
#'
#' @inherit mu_rafe_stein references
#' @export
sigma_crafe_floor <- function(Sigma_hat, tau_rel = NULL) {
  # TODO: implement. Matches paper §3 (covtau / cc).
  stop("TODO: implement sigma_crafe_floor()")
}

#' Jointly T-RAFE-Tuned Moment Correction
#'
#' Tunes \eqn{(\kappa, \tau)} jointly to minimise inner-validation T-RAFE.
#' Empirically indistinguishable from sequential tuning (`sep_tune()`) — the
#' bound's two channels are separable.
#'
#' @param mu_hat Forecast mean.
#' @param Sigma_hat Forecast covariance.
#' @param ... Passed to [mu_rafe_stein()] and [sigma_crafe_floor()].
#'
#' @return A list with components `mu_tilde`, `Sigma_tilde`, `kappa`, `tau_rel`.
#'
#' @inherit mu_rafe_stein references
#' @export
joint_trafe_tune <- function(mu_hat, Sigma_hat, ...) {
  # TODO: implement. Matches paper §3 (joint).
  stop("TODO: implement joint_trafe_tune()")
}

#' Sequentially Tuned Moment Correction
#'
#' Tunes \eqn{\kappa} on inner-validation RAFE and \eqn{\tau} on inner-validation
#' C-RAFE independently. Empirically indistinguishable from joint tuning
#' (`joint_trafe_tune()`) at every cross-section size in the validation
#' empirics — the bound's two channels are separable.
#'
#' @param mu_hat Forecast mean.
#' @param Sigma_hat Forecast covariance.
#' @param ... Passed to [mu_rafe_stein()] and [sigma_crafe_floor()].
#'
#' @return A list with components `mu_tilde`, `Sigma_tilde`, `kappa`, `tau_rel`.
#'
#' @inherit mu_rafe_stein references
#' @export
sep_tune <- function(mu_hat, Sigma_hat, ...) {
  # TODO: implement. Matches paper §3 (sep).
  stop("TODO: implement sep_tune()")
}
