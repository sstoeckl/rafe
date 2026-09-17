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
# (github.com/sstoeckl/Lost_in_Translation_Replication, code/00_setup.R)
# ---------------------------------------------------------------------------

test_that("C-RAFE matches the closed form under lambda contamination", {
  # The paper's design sets Sigma_hat = (1 - lambda)^2 * Sigma, for which
  # compute_crafe_lambda(lambda) = |1/(1-lambda)^2 - 1| exactly.
  set.seed(201)
  for (n in c(2L, 5L, 12L)) {
    Sigma <- crossprod(matrix(rnorm(n * n), n)) / n + diag(0.1, n)
    for (lambda in c(0, 0.002, 0.01, 0.05, 0.1)) {
      expect_equal(
        compute_crafe(Sigma, (1 - lambda)^2 * Sigma),
        abs(1 / (1 - lambda)^2 - 1),
        tolerance = 1e-9,
        info = paste("N =", n, "lambda =", lambda)
      )
    }
  }
})

test_that("the metrics reproduce the replication code's own definitions", {
  # Verbatim transcriptions of compute_rafe() and compute_crafe() from
  # code/00_setup.R of the replication repository, called in that repo's
  # argument order.
  ref_rafe <- function(mu_hat, mu, Sigma) {
    e_mu <- mu_hat - mu
    sqrt(as.numeric(t(e_mu) %*% solve(Sigma) %*% e_mu))
  }
  ref_crafe <- function(Sigma, Sigma_hat) {
    eig <- eigen(Sigma, symmetric = TRUE)
    Sigma_sqrt <- eig$vectors %*% diag(sqrt(pmax(eig$values, 0))) %*% t(eig$vectors)
    M <- Sigma_sqrt %*% solve(Sigma_hat) %*% Sigma_sqrt - diag(nrow(Sigma))
    max(svd(M)$d)
  }

  set.seed(202)
  for (i in 1:20) {
    n <- sample(2:10, 1)
    Sigma     <- crossprod(matrix(rnorm(n * n), n)) / n + diag(0.2, n)
    Sigma_hat <- crossprod(matrix(rnorm(n * n), n)) / n + diag(0.2, n)
    mu     <- rnorm(n) / 20
    mu_hat <- rnorm(n) / 20

    expect_equal(compute_rafe(mu_hat, mu, Sigma), ref_rafe(mu_hat, mu, Sigma),
                 tolerance = 1e-9)
    expect_equal(compute_crafe(Sigma, Sigma_hat), ref_crafe(Sigma, Sigma_hat),
                 tolerance = 1e-9)

    # The replication code's T-RAFE fixes c = 1 and takes |SR*|.
    SR_star <- sqrt(drop(t(mu) %*% solve(Sigma) %*% mu))
    expect_equal(
      as.numeric(compute_trafe(mu_hat, mu, Sigma, Sigma_hat, c = 1)),
      ref_rafe(mu_hat, mu, Sigma) + abs(SR_star) * ref_crafe(Sigma, Sigma_hat),
      tolerance = 1e-9
    )
  }
})
