#' rafe: Decision-Aligned Forecast Evaluation and Moment Correction
#'
#' Evaluates return forecasts by the economic damage their errors cause
#' rather than by their size, and corrects the forecast moments that a
#' mean-variance optimiser is most sensitive to.
#'
#' Conventional accuracy measures weight every asset equally; a portfolio
#' optimiser does not. The risk-adjusted forecast error (RAFE) and the
#' operator-norm covariance forecast error (C-RAFE) reweight forecast errors
#' by the risk metric the decision actually uses, and combine into an upper
#' bound on the Sharpe-ratio gap of the plug-in portfolio (T-RAFE).
#'
#' Two groups of functions:
#'
#' - **Metrics**: [compute_rafe()], [compute_crafe()], [compute_trafe()].
#' - **Moment correction**: [mu_rafe_stein()] and [sigma_crafe_floor()] act on
#'   the two channels of the bound; [sep_tune()] and [joint_trafe_tune()]
#'   select their intensities on an inner-validation split.
#'
#' Two datasets, [ff12] and [ff49], make every example reproducible offline.
#'
#' @references
#'   Salcher, T., Stöckl, S., & Hanke, M. (2026). Lost in Translation?
#'   Risk-Adjusting RMSE for Economic Forecast Performance.
#'   *Journal of Forecasting*. \doi{10.1002/for.70134}
#'
#'   Stöckl, S., Salcher, T., & Hanke, M. Post-Optimal Moment Correction for
#'   Mean-Variance Portfolios. Working paper.
#'
#' @keywords internal
"_PACKAGE"
