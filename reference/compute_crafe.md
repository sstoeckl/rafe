# Covariance Forecast Error (C-RAFE)

Computes the operator-norm precision distortion \\\mathrm{C\text{-}RAFE}
= \\\Sigma^{1/2} \hat\Sigma^{-1} \Sigma^{1/2} - I\\\_2\\.

## Usage

``` r
compute_crafe(Sigma, Sigma_hat, variant = c("none", "cc05", "cc0", "cv", "i"))
```

## Arguments

- Sigma:

  Realised covariance matrix.

- Sigma_hat:

  Forecast covariance matrix.

- variant:

  Covariance restriction applied to `Sigma`, from the nested sequence of
  Section 3.3; see
  [`restrict_cov()`](https://sstoeckl.github.io/rafe/reference/restrict_cov.md).

## Value

A length-1 numeric.

## Details

The metric is zero if and only if \\\hat\Sigma = \Sigma\\, and is
invariant to the common scaling of both arguments.

## Argument order

The realised covariance comes **first**, matching `compute_crafe()` in
the replication code of Salcher, Stöckl & Hanke (2026) so that paper
code runs unchanged against this package. The metric is not symmetric in
its arguments, so the order matters.

## References

Salcher, L., Stöckl, S., & Hanke, M. (2026). Lost in Translation?
Risk-Adjusting RMSE for Economic Forecast Performance. *Journal of
Forecasting*. doi:10.1002/for.70134

## Examples

``` r
Sigma <- diag(c(1, 2, 3))
compute_crafe(Sigma, Sigma)          # exactly 0
#> [1] 0
compute_crafe(Sigma, 2 * Sigma)      # forecast covariance twice too large
#> [1] 0.5

# Variants from the nested sequence of Section 3.3:
compute_crafe(Sigma, 2 * Sigma, variant = "cc0")
#> [1] 0.5
```
