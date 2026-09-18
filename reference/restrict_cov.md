# Restricted Covariance Matrices for the Nested Metric Sequence

Applies the sequence of simplifying assumptions of Salcher, Stöckl &
Hanke (2026, Section 3.3) to a covariance matrix. Imposing them one
after another turns the risk-adjusted forecast error into the ordinary
root mean squared error, which is how the paper frames RMSE as a
severely restricted special case of RAFE.

## Usage

``` r
restrict_cov(Sigma, variant = c("none", "cc05", "cc0", "cv", "i"))
```

## Arguments

- Sigma:

  A covariance matrix.

- variant:

  One of `"none"`, `"cc05"`, `"cc0"`, `"cv"`, `"i"`.

## Value

A covariance matrix of the same dimension as `Sigma`.

## Details

- `"none"`:

  The covariance is left alone: the full RAFE.

- `"cc05"`:

  All pairwise correlations replaced by \\\rho = 0.5\\, variances kept —
  Equation (17).

- `"cc0"`:

  All correlations set to zero, variances kept; the covariance becomes
  diagonal.

- `"cv"`:

  Correlations zero and all variances set to the average variance of
  `Sigma`.

- `"i"`:

  Correlations zero and all variances set to one, giving the identity.
  With this restriction the mean error is the RMSE of the paper,
  \\\lVert\hat\mu - \mu\rVert_2\\.

## References

Salcher, L., Stöckl, S., & Hanke, M. (2026). Lost in Translation?
Risk-Adjusting RMSE for Economic Forecast Performance. *Journal of
Forecasting*. doi:10.1002/for.70134

## Examples

``` r
set.seed(1)
Sigma <- stats::cov(matrix(rnorm(200), 40, 5))
round(stats::cov2cor(restrict_cov(Sigma, "cc05")), 3)   # all rho = 0.5
#>      [,1] [,2] [,3] [,4] [,5]
#> [1,]  1.0  0.5  0.5  0.5  0.5
#> [2,]  0.5  1.0  0.5  0.5  0.5
#> [3,]  0.5  0.5  1.0  0.5  0.5
#> [4,]  0.5  0.5  0.5  1.0  0.5
#> [5,]  0.5  0.5  0.5  0.5  1.0
restrict_cov(Sigma, "i")                                # identity
#>      [,1] [,2] [,3] [,4] [,5]
#> [1,]    1    0    0    0    0
#> [2,]    0    1    0    0    0
#> [3,]    0    0    1    0    0
#> [4,]    0    0    0    1    0
#> [5,]    0    0    0    0    1
```
