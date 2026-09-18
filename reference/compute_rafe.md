# Risk-Adjusted Forecast Error (RAFE)

Computes \\\mathrm{RAFE} = \sqrt{(\hat\mu - \mu)^\top \Sigma^{-1}
(\hat\mu - \mu)}\\, the risk-adjusted distance between a forecast mean
and the realised mean.

## Usage

``` r
compute_rafe(
  mu_hat,
  mu,
  Sigma = NULL,
  Sigma_inv = NULL,
  variant = c("none", "cc05", "cc0", "cv", "i")
)
```

## Arguments

- mu_hat:

  Numeric vector of forecast means.

- mu:

  Numeric vector of realised means (same length as `mu_hat`).

- Sigma:

  Realised covariance matrix. The argument order matches the paper's
  replication code, so `compute_rafe(mu_hat, mu, Sigma)` does the same
  thing here as it does there.

- Sigma_inv:

  Optional precision matrix, supplied instead of `Sigma` when it is
  already available (it avoids an inversion). Must be a *precision*, not
  a covariance. Ignored unless `variant` is `"none"`.

- variant:

  Covariance restriction from the nested sequence of Section 3.3; see
  [`restrict_cov()`](https://sstoeckl.github.io/rafe/reference/restrict_cov.md).
  `"i"` returns the paper's RMSE.

## Value

A length-1 numeric.

## Details

The metric is a Mahalanobis distance in the population-risk metric, with
no \\1/N\\ normalisation. Consequently, at \\\Sigma = I\\ it reduces to
the Euclidean norm of the forecast error, which is \\\sqrt{N}\\ times
the conventional RMSE — see the examples.

## References

Salcher, L., Stöckl, S., & Hanke, M. (2026). Lost in Translation?
Risk-Adjusting RMSE for Economic Forecast Performance. *Journal of
Forecasting*. doi:10.1002/for.70134

## Examples

``` r
set.seed(1)
mu     <- rnorm(5)
mu_hat <- mu + rnorm(5, sd = 0.1)

compute_rafe(mu_hat, mu, Sigma = diag(5))
#> [1] 0.1371373

# At Sigma = I the metric is sqrt(N) times RMSE:
sqrt(5) * sqrt(mean((mu_hat - mu)^2))
#> [1] 0.1371373
```
