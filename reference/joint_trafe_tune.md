# Jointly T-RAFE-Tuned Moment Correction

Tunes \\(\kappa, \tau)\\ jointly on the full grid to minimise
inner-validation T-RAFE. Empirically close to indistinguishable from the
sequential search of
[`sep_tune()`](https://sstoeckl.github.io/rafe/reference/sep_tune.md),
which is the practical evidence that the bound's two channels are
separable.

## Usage

``` r
joint_trafe_tune(
  R_train,
  kappa_grid = seq(0, 1, length.out = 101),
  tau_grid = seq(0, 0.5, length.out = 101),
  inner_split = c(40, 20),
  target = c("grandmean", "zero")
)
```

## Arguments

- R_train:

  Numeric matrix of training returns, observations in rows and assets in
  columns.

- kappa_grid:

  Grid of shrinkage intensities to search. The default matches the
  working paper's grid, 101 points on `[0, 1]`.

- tau_grid:

  Grid of relative eigenvalue floors to search. The default matches the
  working paper's grid, 101 points on `[0, 0.5]`. Floors above 0.5
  flatten the spectrum so aggressively that the cleaned covariance
  carries little information about the original, and are not admitted.

- inner_split:

  Length-2 vector: rows used for inner training and for inner
  validation.

- target:

  Shrinkage target passed to
  [`mu_rafe_stein()`](https://sstoeckl.github.io/rafe/reference/mu_rafe_stein.md).

## Value

A list with components `mu_tilde`, `Sigma_tilde`, `kappa`, `tau_rel`,
and the full inner-validation loss surface `loss_grid` (kappa in rows,
tau in columns).

## References

Stöckl, S., Salcher, L., & Hanke, M. Post-Optimal Moment Correction for
Mean-Variance Portfolios. Working paper.

## See also

[`sep_tune()`](https://sstoeckl.github.io/rafe/reference/sep_tune.md)
for the sequential search.

## Examples

``` r
set.seed(1)
R <- matrix(rnorm(80 * 6, mean = 0.008, sd = 0.05), 80, 6)
fit <- joint_trafe_tune(R, kappa_grid = seq(0, 1, by = 0.1),
                           tau_grid = seq(0, 1, by = 0.2))
fit$kappa
#> [1] 0
fit$tau_rel
#> [1] 1
```
