# Internal linear-algebra helpers. Not exported.
#
# Ported from 3212_recipe/R/00_helpers.R, with the
# MASS::ginv fallback replaced by an eigen-based Moore-Penrose inverse so the
# package keeps a base-only dependency footprint.

#' Symmetrise a matrix
#' @noRd
.as_sym <- function(S) 0.5 * (S + t(S))

#' Inverse of a symmetric positive-definite matrix
#'
#' Uses a Cholesky factorisation when `S` is numerically PD, and falls back to
#' an eigen-based Moore-Penrose pseudo-inverse otherwise.
#'
#' @noRd
.pd_inverse <- function(S, tol = 1e-12) {
  S <- .as_sym(S)
  ch <- tryCatch(chol(S), error = function(e) NULL)
  if (!is.null(ch)) return(chol2inv(ch))
  eg <- eigen(S, symmetric = TRUE)
  av <- abs(eg$values)
  keep <- av > tol * max(av)
  inv_vals <- numeric(length(eg$values))
  inv_vals[keep] <- 1 / eg$values[keep]
  eg$vectors %*% (inv_vals * t(eg$vectors))
}

#' Symmetric positive-semidefinite matrix square root
#' @noRd
.sqrtm_sym <- function(S, eps = 1e-12) {
  eg <- eigen(.as_sym(S), symmetric = TRUE)
  vals <- sqrt(pmax(eg$values, eps))
  eg$vectors %*% (vals * t(eg$vectors))
}

#' Spectral (operator-2) norm of a matrix
#' @noRd
.spec_norm <- function(M) base::norm(M, type = "2")

#' Validate a numeric vector argument
#' @noRd
.check_vec <- function(x, arg) {
  if (!is.numeric(x) || !is.null(dim(x)) && !(is.matrix(x) && min(dim(x)) == 1L)) {
    stop("`", arg, "` must be a numeric vector.", call. = FALSE)
  }
  x <- as.numeric(x)
  if (!all(is.finite(x))) {
    stop("`", arg, "` must contain only finite values.", call. = FALSE)
  }
  x
}

#' Validate a square matrix argument of a given dimension
#' @noRd
.check_mat <- function(S, n, arg) {
  if (!is.matrix(S) || !is.numeric(S)) {
    stop("`", arg, "` must be a numeric matrix.", call. = FALSE)
  }
  if (nrow(S) != ncol(S)) {
    stop("`", arg, "` must be square; got ", nrow(S), " x ", ncol(S), ".",
         call. = FALSE)
  }
  if (!is.null(n) && nrow(S) != n) {
    stop("`", arg, "` must be ", n, " x ", n, " to match the mean vector; got ",
         nrow(S), " x ", ncol(S), ".", call. = FALSE)
  }
  if (!all(is.finite(S))) {
    stop("`", arg, "` must contain only finite values.", call. = FALSE)
  }
  S
}

#' Resolve a precision matrix from either a precision or a covariance argument
#' @noRd
.resolve_precision <- function(Sigma_inv, Sigma, n,
                               inv_arg = "Sigma_inv", cov_arg = "Sigma") {
  if (!is.null(Sigma_inv)) {
    return(.as_sym(.check_mat(Sigma_inv, n, inv_arg)))
  }
  if (!is.null(Sigma)) {
    return(.pd_inverse(.check_mat(Sigma, n, cov_arg)))
  }
  stop("Supply either `", inv_arg, "` or `", cov_arg, "`.", call. = FALSE)
}
