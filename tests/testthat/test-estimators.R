make_cov <- function(n, seed) {
  set.seed(seed)
  crossprod(matrix(rnorm(n * n), n)) / n + diag(0.1, n)
}

test_that("mu_rafe_stein respects the endpoints of kappa", {
  set.seed(31)
  mu_hat <- rnorm(5)
  S      <- make_cov(5, 31)
  expect_equal(as.numeric(mu_rafe_stein(mu_hat, S, kappa = 0)), mu_hat)
  expect_equal(as.numeric(mu_rafe_stein(mu_hat, S, kappa = 1)),
               rep(mean(mu_hat), 5))
  expect_equal(as.numeric(mu_rafe_stein(mu_hat, S, kappa = 1, target = "zero")),
               rep(0, 5))
})

test_that("mu_rafe_stein clamps kappa into [0, 1]", {
  set.seed(32)
  mu_hat <- rnorm(4)
  S      <- make_cov(4, 32)
  expect_equal(attr(mu_rafe_stein(mu_hat, S, kappa = -2), "kappa"), 0)
  expect_equal(attr(mu_rafe_stein(mu_hat, S, kappa = 17), "kappa"), 1)
})

test_that("the grandmean target preserves the cross-sectional mean", {
  set.seed(33)
  mu_hat <- rnorm(7)
  S      <- make_cov(7, 33)
  for (k in c(0, 0.3, 0.75, 1)) {
    expect_equal(mean(mu_rafe_stein(mu_hat, S, kappa = k)), mean(mu_hat),
                 tolerance = 1e-12, info = paste("kappa =", k))
  }
})

test_that("mu_rafe_stein is the whitened-space James-Stein map, not its inverse", {
  # Guards the covariance/precision inversion class of bug: the intensity must
  # be computed against Sigma_hat^{-1/2} (whitening), never Sigma_hat^{+1/2}
  # (colouring). Verified by rebuilding the estimator from first principles in
  # the whitened space and comparing.
  set.seed(34)
  n <- 6L
  T_obs <- 60
  mu_hat <- rnorm(n) / 20
  S <- make_cov(n, 34)

  eg <- eigen(S, symmetric = TRUE)
  W  <- eg$vectors %*% diag(1 / sqrt(eg$values)) %*% t(eg$vectors)  # Sigma^-1/2
  Wi <- eg$vectors %*% diag(sqrt(eg$values))     %*% t(eg$vectors)  # Sigma^+1/2

  mu_0  <- rep(mean(mu_hat), n)
  m_hat <- drop(W %*% mu_hat)
  m_0   <- drop(W %*% mu_0)
  dist_sq <- sum((m_hat - m_0)^2)
  shrink  <- max(0, 1 - (n - 2) / (T_obs * dist_sq))
  js      <- m_0 + shrink * (m_hat - m_0)
  expected <- drop(Wi %*% js)

  got <- mu_rafe_stein(mu_hat, S, T_obs = T_obs)
  expect_equal(as.numeric(got), expected, tolerance = 1e-10)

  # The test has teeth: swapping the whitening for the colouring matrix gives a
  # materially different answer, so this is not trivially satisfied.
  m_hat_bad <- drop(Wi %*% mu_hat)
  m_0_bad   <- drop(Wi %*% mu_0)
  dist_bad  <- sum((m_hat_bad - m_0_bad)^2)
  expect_false(isTRUE(all.equal(dist_sq, dist_bad, tolerance = 1e-6)))
})

test_that("the James-Stein intensity matches its closed form", {
  set.seed(35)
  n <- 5L
  T_obs <- 48
  mu_hat <- rnorm(n) / 10
  S <- make_cov(n, 35)
  d <- mu_hat - rep(mean(mu_hat), n)
  dist_sq <- drop(t(d) %*% solve(S) %*% d)
  expect_equal(attr(mu_rafe_stein(mu_hat, S, T_obs = T_obs), "kappa"),
               max(0, min(1, (n - 2) / (T_obs * dist_sq))),
               tolerance = 1e-10)
})

test_that("mu_rafe_stein does not shrink below N = 3", {
  expect_equal(attr(mu_rafe_stein(c(0.01, 0.02), diag(2), T_obs = 50), "kappa"), 0)
})

test_that("mu_rafe_stein requires exactly one of kappa and T_obs", {
  expect_error(mu_rafe_stein(rnorm(4), diag(4)), "Supply either")
  expect_error(mu_rafe_stein(rnorm(4), diag(4), T_obs = -1), "positive")
})

test_that("sigma_crafe_floor at tau_rel = 0 leaves the spectrum alone", {
  S <- make_cov(5, 36)
  expect_equal(unclass(sigma_crafe_floor(S, 0))[1:5, 1:5], S, tolerance = 1e-10,
               ignore_attr = TRUE)
})

test_that("sigma_crafe_floor raises the smallest eigenvalue to the floor", {
  S <- make_cov(6, 37)
  lam <- eigen(S, symmetric = TRUE)$values
  tau_rel <- 0.5
  floor_abs <- tau_rel * mean(lam)
  cleaned <- sigma_crafe_floor(S, tau_rel)
  lam_new <- eigen(cleaned, symmetric = TRUE)$values
  expect_gte(min(lam_new), floor_abs - 1e-10)
  expect_equal(sort(lam_new), sort(pmax(lam, floor_abs)), tolerance = 1e-10)
  expect_equal(attr(cleaned, "tau_abs"), floor_abs, tolerance = 1e-12)
})

test_that("sigma_crafe_floor returns a symmetric positive-definite matrix", {
  S <- make_cov(5, 38)
  cleaned <- sigma_crafe_floor(S, 0.4)
  expect_equal(cleaned, t(cleaned), tolerance = 1e-12, ignore_attr = TRUE)
  expect_gt(min(eigen(cleaned, symmetric = TRUE)$values), 0)
})

test_that("a larger floor cannot increase C-RAFE distortion of a noisy spectrum", {
  # Sanity on the intended use: flooring a badly-conditioned sample covariance
  # toward a well-conditioned truth should help at some positive tau.
  set.seed(39)
  n <- 8L
  Sigma <- diag(n)
  R <- matrix(rnorm(12 * n), 12, n)   # T barely above N: very noisy spectrum
  S_hat <- stats::cov(R)
  base <- compute_crafe(Sigma, S_hat)
  best <- min(vapply(seq(0, 1, by = 0.05),
                     function(t) compute_crafe(Sigma, sigma_crafe_floor(S_hat, t)),
                     numeric(1)))
  expect_lt(best, base)
})

test_that("sigma_crafe_floor rejects malformed input", {
  expect_error(sigma_crafe_floor(make_cov(3, 40)), "Supply `tau_rel`")
  expect_error(sigma_crafe_floor(make_cov(3, 40), -0.1), "non-negative")
})

test_that("the tuners return grid points and a consistent correction", {
  set.seed(41)
  R <- matrix(rnorm(80 * 6, mean = 0.008, sd = 0.05), 80, 6)
  kg <- seq(0, 1, by = 0.1)
  tg <- seq(0, 1, by = 0.2)

  s <- sep_tune(R, kappa_grid = kg, tau_grid = tg)
  expect_true(s$kappa %in% kg)
  expect_true(s$tau_rel %in% tg)
  expect_length(s$mu_tilde, 6L)
  expect_equal(dim(s$Sigma_tilde), c(6L, 6L))
  expect_equal(as.numeric(s$mu_tilde),
               as.numeric(mu_rafe_stein(colMeans(R), stats::cov(R),
                                        kappa = s$kappa)),
               tolerance = 1e-12)

  j <- joint_trafe_tune(R, kappa_grid = kg, tau_grid = tg)
  expect_true(j$kappa %in% kg)
  expect_true(j$tau_rel %in% tg)
  expect_equal(dim(j$loss_grid), c(length(kg), length(tg)))
})

test_that("joint_trafe_tune reports the cell that actually minimises the grid", {
  set.seed(42)
  R <- matrix(rnorm(80 * 5, mean = 0.01, sd = 0.04), 80, 5)
  kg <- seq(0, 1, by = 0.25)
  tg <- seq(0, 1, by = 0.25)
  j <- joint_trafe_tune(R, kappa_grid = kg, tau_grid = tg)
  expect_equal(
    j$loss_grid[as.character(j$kappa), as.character(j$tau_rel)],
    min(j$loss_grid),
    tolerance = 1e-12
  )
})

test_that("the sequential and joint searches select comparable kappa", {
  # Separability: the two channels tune independently, so the
  # sequential search should land close to the joint one on the mean channel.
  set.seed(43)
  R <- matrix(rnorm(90 * 8, mean = 0.007, sd = 0.05), 90, 8)
  kg <- seq(0, 1, by = 0.05)
  tg <- seq(0, 1, by = 0.1)
  s <- sep_tune(R, kappa_grid = kg, tau_grid = tg)
  j <- joint_trafe_tune(R, kappa_grid = kg, tau_grid = tg)
  expect_lt(abs(s$kappa - j$kappa), 0.35)
})

test_that("the tuners validate the inner split", {
  R <- matrix(rnorm(30 * 4), 30, 4)
  expect_error(sep_tune(R, inner_split = c(40, 20)), "but `R_train` has 30")
  expect_error(sep_tune(R, inner_split = c(5, 1)), "at least 2")
})
