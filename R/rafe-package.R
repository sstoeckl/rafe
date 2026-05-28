#' rafe: Decision-Aligned Forecast Evaluation and Moment Correction
#'
#' Computes risk-adjusted forecast errors (RAFE), covariance forecast
#' errors (C-RAFE), and the combined upper bound on the Sharpe-ratio gap
#' (T-RAFE) of plug-in mean-variance portfolios, and provides decision-
#' aligned moment-correction estimators and training-loss objectives.
#'
#' Companion implementation to Salcher, Stoeckl & Hanke (2026, *Journal
#' of Forecasting*) and the RAFE Program of papers.
#'
#' Three groups of functions:
#'
#' - Metrics: [compute_rafe()], [compute_crafe()], [compute_trafe()].
#' - Post-processing estimators (Paper A): [mu_rafe_stein()],
#'   [sigma_crafe_floor()], [joint_trafe_tune()], [sep_tune()].
#' - Training-loss objectives (Paper B): [xgb_trafe_objective()],
#'   [torch_trafe_loss()] (to be added).
#'
#' @keywords internal
"_PACKAGE"
