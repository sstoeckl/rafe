# The two tests that check code against theory rather than against itself.
#
# 1. At Sigma = I the risk-adjusted metric must collapse onto the ordinary
#    Euclidean forecast error, i.e. sqrt(N) * RMSE. RAFE carries no 1/N
#    normalisation (Salcher, Stoeckl & Hanke 2026, eq. 1), so the sqrt(N)
#    factor is part of the claim, not slack in the test.
#
# 2. The Sharpe-gap bound Delta <= c * RAFE + SR* * C-RAFE must hold in every
#    Monte Carlo draw, with c = 1 if RAFE <= SR* and c = 2 otherwise.

rand_cov <- function(n) crossprod(matrix(rnorm(n * n), n)) / n + diag(0.1, n)

rmvn <- function(n_obs, mu, Sigma) {
  L <- chol(Sigma)
  matrix(rnorm(n_obs * length(mu)), n_obs) %*% L + rep(mu, each = n_obs)
}

# Realised Sharpe of the plug-in tangency portfolio built from (mu_hat, S_hat),
# evaluated under the true moments. Sharpe is scale-invariant in w, so the
# normalisation below does not affect the result.
plugin_sharpe <- function(mu_hat, S_hat, mu, Sigma) {
  iSh <- solve(S_hat)
  w   <- drop(iSh %*% mu_hat) / drop(t(mu_hat) %*% iSh %*% mu_hat)
  drop(t(mu) %*% w) / sqrt(drop(t(w) %*% Sigma %*% w))
}

test_that("at Sigma = I, RAFE collapses to sqrt(N) * RMSE", {
  set.seed(101)
  for (n in c(1L, 2L, 5L, 25L)) {
    mu     <- rnorm(n)
    mu_hat <- mu + rnorm(n, sd = 0.3)
    rmse   <- sqrt(mean((mu_hat - mu)^2))
    expect_equal(
      compute_rafe(mu_hat, mu, Sigma = diag(n)),
      sqrt(n) * rmse,
      tolerance = 1e-10,
      info = paste("N =", n)
    )
  }
})

test_that("at Sigma = sigma^2 I, RAFE is RMSE rescaled by sqrt(N)/sigma", {
  set.seed(102)
  n     <- 6
  sigma <- 0.4
  mu     <- rnorm(n)
  mu_hat <- mu + rnorm(n, sd = 0.2)
  expect_equal(
    compute_rafe(mu_hat, mu, Sigma = diag(sigma^2, n)),
    sqrt(n) * sqrt(mean((mu_hat - mu)^2)) / sigma,
    tolerance = 1e-10
  )
})

test_that("the Sharpe-gap bound holds in every Monte Carlo draw", {
  set.seed(2026)
  n_rep  <- 300L
  n      <- 5L
  n_obs  <- 60L

  violations <- 0L
  slack_min  <- Inf
  gap_max    <- -Inf

  for (i in seq_len(n_rep)) {
    Sigma <- rand_cov(n)
    mu    <- rnorm(n) * 0.05

    X       <- rmvn(n_obs, mu, Sigma)
    mu_hat  <- colMeans(X)
    S_hat   <- stats::cov(X)

    SR_star <- sqrt(drop(t(mu) %*% solve(Sigma) %*% mu))
    Delta   <- SR_star - plugin_sharpe(mu_hat, S_hat, mu, Sigma)

    bound <- compute_trafe(mu_hat, mu, Sigma, S_hat)

    # The gap is non-negative by construction: SR* maximises the Sharpe ratio.
    expect_gte(Delta, -1e-10)

    if (Delta > as.numeric(bound) + 1e-10) violations <- violations + 1L
    slack_min <- min(slack_min, as.numeric(bound) - Delta)
    gap_max   <- max(gap_max, Delta)
  }

  expect_identical(violations, 0L)
  # Sanity: the draws must actually exercise the bound, not sit at Delta ~ 0.
  expect_gt(gap_max, 0.01)
  expect_gte(slack_min, -1e-10)
})

test_that("both branches of the constant c are exercised by random draws", {
  set.seed(2027)
  cs <- replicate(200L, {
    n     <- 5L
    Sigma <- rand_cov(n)
    mu    <- rnorm(n) * 0.05
    X     <- rmvn(60L, mu, Sigma)
    attr(compute_trafe(colMeans(X), mu, Sigma, stats::cov(X)), "c")
  })
  expect_true(all(cs %in% c(1, 2)))
})

# ---------------------------------------------------------------------------
# Agreement with the published replication code
# (github.com/sstoeckl/Lost_in_Translation_Replication, code/functions/metrics.R)
#
# The helpers below are verbatim transcriptions of that file. If the package
# ever drifts from the published definitions, these fail.
# ---------------------------------------------------------------------------

ref_inv <- function(Sigma) chol2inv(chol(Sigma))

ref_const_corr_cov <- function(Sigma, rho) {
  sd_vec <- sqrt(diag(Sigma))
  R <- matrix(rho, ncol = length(sd_vec), nrow = length(sd_vec))
  diag(R) <- 1
  D <- diag(sd_vec)
  D %*% R %*% D
}

ref_rae <- function(e, W) sqrt(drop(t(e) %*% W %*% e))

ref_crafe <- function(Sigma_hat_inv, Sigma) {
  eig <- eigen(Sigma)
  Sigma_sqrt <- eig$vectors %*% diag(sqrt(eig$values)) %*% t(eig$vectors)
  A <- Sigma_sqrt %*% Sigma_hat_inv %*% Sigma_sqrt
  norm(A - diag(ncol(Sigma)), type = "2")
}

ref_forecast_errors <- function(mu_hat, Sigma_hat_inv, mu, Sigma) {
  e <- mu_hat - mu
  N <- length(e)
  Sigma_cc05_inv <- ref_inv(ref_const_corr_cov(Sigma, 0.5))
  Sigma_cc0_inv  <- ref_inv(ref_const_corr_cov(Sigma, 0))
  s2_avg <- mean(diag(Sigma))
  c(rafe       = ref_rae(e, solve(Sigma)),
    rafe_cc05  = ref_rae(e, Sigma_cc05_inv),
    rafe_cc0   = ref_rae(e, Sigma_cc0_inv),
    rafe_cv    = sqrt(sum(e^2) / s2_avg),
    rmse       = sqrt(sum(e^2)),
    crafe      = ref_crafe(Sigma_hat_inv, Sigma),
    crafe_cc05 = ref_crafe(Sigma_hat_inv, solve(Sigma_cc05_inv)),
    crafe_cc0  = ref_crafe(Sigma_hat_inv, solve(Sigma_cc0_inv)),
    crafe_cv   = ref_crafe(Sigma_hat_inv, diag(s2_avg, N)),
    crafe_i    = ref_crafe(Sigma_hat_inv, diag(1, N)))
}

test_that("every metric variant matches the published replication code", {
  set.seed(303)
  map_rafe  <- c(none = "rafe", cc05 = "rafe_cc05", cc0 = "rafe_cc0",
                 cv = "rafe_cv", i = "rmse")
  map_crafe <- c(none = "crafe", cc05 = "crafe_cc05", cc0 = "crafe_cc0",
                 cv = "crafe_cv", i = "crafe_i")

  for (i in 1:25) {
    n <- sample(3:10, 1)
    Sigma     <- crossprod(matrix(rnorm(n * n), n)) / n + diag(0.2, n)
    Sigma_hat <- crossprod(matrix(rnorm(n * n), n)) / n + diag(0.2, n)
    mu     <- rnorm(n) / 20
    mu_hat <- rnorm(n) / 20

    ref <- ref_forecast_errors(mu_hat, ref_inv(Sigma_hat), mu, Sigma)

    for (v in names(map_rafe)) {
      expect_equal(compute_rafe(mu_hat, mu, Sigma, variant = v),
                   ref[[map_rafe[[v]]]], tolerance = 1e-9,
                   info = paste("rafe variant", v))
      expect_equal(compute_crafe(Sigma, Sigma_hat, variant = v),
                   ref[[map_crafe[[v]]]], tolerance = 1e-9,
                   info = paste("crafe variant", v))
    }

    # Eq. (22): T-RAFE = RAFE + SR* * C-RAFE, with no constant on the mean
    # channel, and SR* taken from the unrestricted realised moments.
    SR_star <- sqrt(drop(t(mu) %*% solve(Sigma) %*% mu))
    for (v in names(map_rafe)) {
      expect_equal(
        as.numeric(compute_trafe(mu_hat, mu, Sigma, Sigma_hat, variant = v)),
        ref[[map_rafe[[v]]]] + SR_star * ref[[map_crafe[[v]]]],
        tolerance = 1e-9, info = paste("trafe variant", v)
      )
    }
  }
})

test_that("the -I restriction reproduces the paper's RMSE exactly", {
  # The paper's RMSE (Section 3.3, final step) is the Euclidean norm of the
  # forecast error, sqrt(sum(e^2)) -- it carries no 1/N.
  set.seed(304)
  for (n in c(2L, 5L, 12L)) {
    Sigma  <- crossprod(matrix(rnorm(n * n), n)) / n + diag(0.2, n)
    mu     <- rnorm(n)
    mu_hat <- mu + rnorm(n, sd = 0.3)
    expect_equal(compute_rafe(mu_hat, mu, Sigma, variant = "i"),
                 sqrt(sum((mu_hat - mu)^2)), tolerance = 1e-12)
    # and the unrestricted metric at Sigma = I coincides with it
    expect_equal(compute_rafe(mu_hat, mu, diag(n)),
                 compute_rafe(mu_hat, mu, Sigma, variant = "i"),
                 tolerance = 1e-12)
  }
})

test_that("restrict_cov imposes the documented structure", {
  set.seed(305)
  n <- 6L
  Sigma <- crossprod(matrix(rnorm(n * n), n)) / n + diag(0.3, n)

  expect_equal(restrict_cov(Sigma, "none"), Sigma, tolerance = 1e-12)
  expect_equal(restrict_cov(Sigma, "i"), diag(n))
  expect_equal(diag(restrict_cov(Sigma, "cv")), rep(mean(diag(Sigma)), n))

  for (v in c("cc05", "cc0")) {
    R <- stats::cov2cor(restrict_cov(Sigma, v))
    rho <- if (v == "cc05") 0.5 else 0
    expect_equal(R[upper.tri(R)], rep(rho, n * (n - 1) / 2), tolerance = 1e-10)
    # variances are preserved by the correlation-only restrictions
    expect_equal(diag(restrict_cov(Sigma, v)), diag(Sigma), tolerance = 1e-10)
  }
})

test_that("the package reproduces Table 4 Panel A of the published paper", {
  # Salcher, Stoeckl & Hanke (2026), Journal of Forecasting, Table 4 Panel A:
  # correlation between the Sharpe-ratio gap and each RAFE variant, on the
  # FF-12 industry portfolios with 60-month training and 36-month test windows.
  skip_on_cran()
  data(ff12, envir = environment())
  R <- as.matrix(ff12[, -1])
  train_len <- 60L; test_len <- 36L

  ridge_inverse <- function(Rw, alpha = 0.20) {
    S  <- stats::cov(Rw)
    Ss <- (1 - alpha) * S + alpha * diag(diag(S))
    ev <- eigen(Ss, symmetric = TRUE, only.values = TRUE)$values
    if (min(ev) < 1e-6) Ss <- Ss + diag(abs(min(ev)) + 1e-4, ncol(Ss))
    solve(Ss)
  }
  norm_w <- function(w) w / sum(w)
  rsharpe <- function(w, Rt) { r <- drop(Rt %*% w); mean(r) / stats::sd(r) }

  variants <- c("none", "cc05", "cc0", "cv", "i")
  starts <- seq_len(nrow(R) - train_len - test_len + 1L)

  panel <- do.call(rbind, lapply(starts, function(s) {
    tr <- R[s:(s + train_len - 1L), ]
    te <- R[(s + train_len):(s + train_len + test_len - 1L), ]
    mu_hat <- colMeans(tr)
    Sw <- ridge_inverse(tr)
    mu <- colMeans(te); Sigma <- stats::cov(te); Sinv <- solve(Sigma)
    gap_tan <- sqrt(drop(t(mu) %*% Sinv %*% mu)) -
      rsharpe(norm_w(drop(Sw %*% mu_hat)), te)
    gap_gmv <- rsharpe(norm_w(drop(Sinv %*% rep(1, ncol(R)))), te) -
      rsharpe(norm_w(drop(Sw %*% rep(1, ncol(R)))), te)
    c(vapply(variants, function(v) compute_rafe(mu_hat, mu, Sigma, variant = v),
             numeric(1)),
      gap_tan = gap_tan, gap_gmv = gap_gmv)
  }))
  panel <- as.data.frame(panel)
  expect_equal(nrow(panel), 625L)

  got_tan <- round(vapply(variants, function(v) cor(panel[[v]], panel$gap_tan),
                          numeric(1)), 3)
  got_gmv <- round(vapply(variants, function(v) cor(panel[[v]], panel$gap_gmv),
                          numeric(1)), 3)

  expect_equal(unname(got_tan), c(0.680, 0.440, 0.149, 0.190, 0.110))
  expect_equal(unname(got_gmv), c(0.065, 0.126, 0.049, 0.041, 0.020))
})
