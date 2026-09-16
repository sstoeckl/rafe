#' Fama-French 12 Industry Portfolios, Monthly
#'
#' Value-weighted monthly returns on Kenneth French's 12 industry portfolios,
#' January 1964 to December 2023. This is the universe used for the worked
#' examples in the package vignettes; it is included so that the vignettes
#' build reproducibly and without a network connection.
#'
#' Dates are stored at month *start* as a labelling convention; each row is the
#' return over that calendar month. Returns are simple (not log) and are
#' expressed as decimals, not percentages.
#'
#' @format A data frame with 720 rows and 13 columns:
#' \describe{
#'   \item{date}{Month, as a `Date` at the first of the month.}
#'   \item{NoDur}{Consumer non-durables.}
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
#' @source Kenneth R. French's data library,
#'   \url{https://mba.tuck.dartmouth.edu/pages/faculty/ken.french/data_library.html}
#'
#' @examples
#' data(ff12)
#' R <- as.matrix(ff12[, -1])
#' dim(R)
#' round(colMeans(R) * 12, 3)   # annualised mean returns
"ff12"
