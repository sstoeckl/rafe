#' Fama-French 12 Industry Portfolios, Monthly Excess Returns
#'
#' Value-weighted monthly *excess* returns (industry return minus the
#' one-month risk-free rate) on Kenneth French's 12 industry portfolios,
#' January 1964 to December 2023. This is the universe used for the worked
#' examples in the package vignettes; it is included so that the vignettes
#' build reproducibly and without a network connection.
#'
#' Dates are stored at month *start* as a labelling convention; each row is the
#' return over that calendar month. Returns are simple (not log), expressed as
#' decimals rather than percentages, and are **excess of the one-month
#' risk-free rate** -- verified identical to `industry - RF` from the
#' Fama-French research factors over all 720 months. They can therefore be used
#' directly in Sharpe-ratio and tangency-portfolio calculations.
#'
#' @format A data frame with 720 rows and 13 columns:
#' \describe{
#'   \item{date}{Month, as a `Date` at the first of the month.}
#'   \item{NoDur}{Consumer non-durables, excess return.}
#'   \item{Durbl}{Consumer durables.}
#'   \item{Manuf}{Manufacturing.}
#'   \item{Enrgy}{Energy.}
#'   \item{Chems}{Chemicals.}
#'   \item{BusEq}{Business equipment.}
#'   \item{Telcm}{Telecoms.}
#'   \item{Utils}{Utilities.}
#'   \item{Shops}{Wholesale and retail.}
#'   \item{Hlth}{Healthcare.}
#'   \item{Money}{Finance.}
#'   \item{Other}{Other.}
#' }
#'
#' @seealso [ff49] for the 49-industry universe used in the evaluation vignette.
#'
#' @source Kenneth R. French's data library,
#'   \url{https://mba.tuck.dartmouth.edu/pages/faculty/ken.french/data_library.html}
#'
#' @examples
#' data(ff12)
#' R <- as.matrix(ff12[, -1])
#' dim(R)
#' round(colMeans(R) * 12, 3)   # annualised mean returns
"ff12"

#' Fama-French 49 Industry Portfolios, Monthly Excess Returns
#'
#' Value-weighted monthly excess returns (industry return minus the one-month
#' risk-free rate) on Kenneth French's 49 industry portfolios, July 1969 to
#' December 2023. Industries with any missing observation over the window were
#' dropped by the source preparation; none were, so all 49 survive.
#'
#' This is the universe of Salcher, Stöckl & Hanke (2026) and the sample start
#' matches that paper's replication configuration exactly. It is included so
#' `vignette("rafe-evaluation")` can reproduce the paper's experiment without a
#' network connection.
#'
#' Dates are stored at month *start* as a labelling convention. Returns are
#' simple (not log) and expressed as decimals rather than percentages.
#'
#' @format A data frame with 654 rows and 50 columns: `date`, plus one column
#'   of excess returns per industry (`Agric`, `Food`, `Soda`, ... `Other`).
#'
#'
#' @source Kenneth R. French's data library,
#'   \url{https://mba.tuck.dartmouth.edu/pages/faculty/ken.french/data_library.html}
#'
#' @seealso [ff12] for the 12-industry universe.
#'
#' @examples
#' data(ff49)
#' dim(ff49)
#' range(ff49$date)
"ff49"
