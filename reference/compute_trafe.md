# Total RAFE (T-RAFE) — Sharpe-Gap Upper Bound

Combines RAFE and C-RAFE into the upper bound on the Sharpe-ratio gap
\\\Delta = SR^{\*} - SR(\hat w)\\ of a plug-in mean-variance portfolio:
\\\Delta \le c \cdot \mathrm{RAFE} + SR^{\*} \cdot
\mathrm{C\text{-}RAFE}\\.

## Usage

``` r
compute_trafe(
  mu_hat,
  mu,
  Sigma,
  Sigma_hat,
  c = 1,
  SR_star = NULL,
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

- Sigma_hat:

  Forecast covariance matrix.

- c:

  Scalar weighting for the mean channel. Defaults to 1, the published
  definition of T-RAFE in Equation (22). Set `c = NULL` to apply the
  data-dependent rule of Equation (4) (1 if `RAFE <= SR_star`, else 2),
  which yields a guaranteed upper bound on the Sharpe-ratio gap.

- SR_star:

  Optional oracle Sharpe ratio. If `NULL`, computed as \\\sqrt{\mu^\top
  \Sigma^{-1} \mu}\\ from `mu` and `Sigma_inv` (or `Sigma`).

- variant:

  Covariance restriction applied to both channels; see
  [`restrict_cov()`](https://sstoeckl.github.io/rafe/reference/restrict_cov.md).
  `SR_star` is always computed from the unrestricted `Sigma`, matching
  the published tables.

## Value

A length-1 numeric, carrying the components as attributes `rafe`,
`crafe`, `c` and `SR_star`.

## Details

Equation (22) of Salcher, Stöckl & Hanke (2026) defines T-RAFE with no
constant on the mean channel, i.e. \\c = 1\\, and that is the default
here so that `compute_trafe()` returns the published metric.

Equation (4) of the same paper gives the *bound*, which carries a
data-dependent constant: \\c = 1\\ when \\\mathrm{RAFE} \le SR^{\*}\\
and \\c = 2\\ otherwise. Pass `c = NULL` to apply that rule and obtain a
quantity guaranteed to dominate the Sharpe-ratio gap. The two coincide
whenever \\\mathrm{RAFE} \le SR^{\*}\\.

## References

Salcher, L., Stöckl, S., & Hanke, M. (2026). Lost in Translation?
Risk-Adjusting RMSE for Economic Forecast Performance. *Journal of
Forecasting*. doi:10.1002/for.70134

## Examples

``` r
set.seed(1)
mu        <- rnorm(5) / 10
Sigma     <- diag(5)
mu_hat    <- mu + rnorm(5, sd = 0.02)
Sigma_hat <- Sigma + diag(5) * 0.05

tr <- compute_trafe(mu_hat, mu, Sigma, Sigma_hat)
tr
#> [1] 0.03668315
#> attr(,"rafe")
#> [1] 0.02742747
#> attr(,"crafe")
#> [1] 0.04761905
#> attr(,"c")
#> [1] 1
#> attr(,"SR_star")
#> [1] 0.1943693
attributes(tr)[c("rafe", "crafe", "c", "SR_star")]
#> $rafe
#> [1] 0.02742747
#> 
#> $crafe
#> [1] 0.04761905
#> 
#> $c
#> [1] 1
#> 
#> $SR_star
#> [1] 0.1943693
#> 
```
