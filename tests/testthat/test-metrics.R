# Unit tests for the three metrics. The two tests that check the metrics
# against the *theory* rather than against themselves live in test-bound.R.

test_that("compute_rafe is zero exactly at a perfect forecast", {
  mu <- c(0.01, -0.02, 0.03)
  expect_equal(compute_rafe(mu, mu, Sigma = diag(3)), 0)
})

test_that("compute_rafe agrees whether given Sigma or Sigma_inv", {
  set.seed(11)
  S <- crossprod(matrix(rnorm(25), 5)) / 5 + diag(0.1, 5)
  mu     <- rnorm(5)
  mu_hat <- mu + rnorm(5, sd = 0.1)
  expect_equal(
    compute_rafe(mu_hat, mu, Sigma = S),
    compute_rafe(mu_hat, mu, Sigma_inv = solve(S)),
    tolerance = 1e-10
  )
})

test_that("compute_rafe equals the Mahalanobis distance by direct computation", {
  set.seed(12)
  S <- crossprod(matrix(rnorm(16), 4)) / 4 + diag(0.1, 4)
  mu     <- rnorm(4)
  mu_hat <- rnorm(4)
  d <- mu_hat - mu
  expect_equal(
    compute_rafe(mu_hat, mu, Sigma = S),
    sqrt(drop(t(d) %*% solve(S) %*% d)),
    tolerance = 1e-10
  )
})

test_that("compute_rafe rejects malformed input", {
  expect_error(compute_rafe(rnorm(3), rnorm(4), Sigma = diag(3)),
               "same length")
  expect_error(compute_rafe(rnorm(3), rnorm(3)),
               "Supply either")
  expect_error(compute_rafe(rnorm(3), rnorm(3), Sigma = diag(4)),
               "to match the mean vector")
  expect_error(compute_rafe(c(1, NA, 3), rnorm(3), Sigma = diag(3)),
               "finite")
})

test_that("compute_crafe is zero exactly when the covariance forecast is right", {
  S <- diag(c(1, 2, 3))
  expect_equal(compute_crafe(S, S), 0, tolerance = 1e-12)
})

test_that("compute_crafe is invariant to common rescaling of both arguments", {
  set.seed(13)
  S     <- crossprod(matrix(rnorm(25), 5)) / 5 + diag(0.1, 5)
  S_hat <- crossprod(matrix(rnorm(25), 5)) / 5 + diag(0.1, 5)
  expect_equal(
    compute_crafe(S, S_hat),
    compute_crafe(7 * S, 7 * S_hat),
    tolerance = 1e-10
  )
})

test_that("compute_crafe matches the operator-norm definition directly", {
  set.seed(14)
  S     <- crossprod(matrix(rnorm(16), 4)) / 4 + diag(0.2, 4)
  S_hat <- crossprod(matrix(rnorm(16), 4)) / 4 + diag(0.2, 4)
  eg   <- eigen(S, symmetric = TRUE)
  Ssq  <- eg$vectors %*% diag(sqrt(eg$values)) %*% t(eg$vectors)
  M    <- Ssq %*% solve(S_hat) %*% Ssq - diag(4)
  expect_equal(compute_crafe(S, S_hat), max(svd(M)$d), tolerance = 1e-10)
})

test_that("compute_trafe applies the theorem's data-dependent c", {
  set.seed(15)
  n <- 4
  S <- diag(n)

  # Tiny mean error relative to SR*: the rule must select c = 1.
  mu_big  <- rep(1, n)
  mu_hat1 <- mu_big + rep(1e-4, n)
  t1 <- compute_trafe(mu_hat1, mu_big, S, S)
  expect_lte(attr(t1, "rafe"), attr(t1, "SR_star"))
  expect_equal(attr(t1, "c"), 1)

  # Large mean error relative to SR*: the rule must select c = 2.
  mu_small <- rep(1e-4, n)
  mu_hat2  <- mu_small + rep(1, n)
  t2 <- compute_trafe(mu_hat2, mu_small, S, S)
  expect_gt(attr(t2, "rafe"), attr(t2, "SR_star"))
  expect_equal(attr(t2, "c"), 2)
})

test_that("compute_trafe decomposes into its reported components", {
  set.seed(16)
  n <- 5
  S     <- crossprod(matrix(rnorm(n * n), n)) / n + diag(0.1, n)
  S_hat <- crossprod(matrix(rnorm(n * n), n)) / n + diag(0.1, n)
  mu     <- rnorm(n) / 20
  mu_hat <- mu + rnorm(n, sd = 0.01)

  tr <- compute_trafe(mu_hat, mu, S, S_hat)
  expect_equal(
    as.numeric(tr),
    attr(tr, "c") * attr(tr, "rafe") + attr(tr, "SR_star") * attr(tr, "crafe"),
    tolerance = 1e-12
  )
  expect_equal(attr(tr, "rafe"),  compute_rafe(mu_hat, mu, Sigma = S),
               tolerance = 1e-10)
  expect_equal(attr(tr, "crafe"), compute_crafe(S, S_hat), tolerance = 1e-10)
  expect_equal(attr(tr, "SR_star"),
               sqrt(drop(t(mu) %*% solve(S) %*% mu)), tolerance = 1e-10)
})

test_that("compute_trafe honours an explicit c and SR_star", {
  set.seed(17)
  n <- 3
  S <- diag(n)
  mu     <- rnorm(n)
  mu_hat <- mu + rnorm(n, sd = 0.1)
  tr <- compute_trafe(mu_hat, mu, S, S, c = 2, SR_star = 0.5)
  expect_equal(attr(tr, "c"), 2)
  expect_equal(attr(tr, "SR_star"), 0.5)
  expect_equal(as.numeric(tr), 2 * attr(tr, "rafe"), tolerance = 1e-12)

  expect_error(compute_trafe(mu_hat, mu, S, S, c = "a"), "single finite numeric")
  expect_error(compute_trafe(mu_hat, mu, S, S, SR_star = c(1, 2)),
               "single finite numeric")
})
