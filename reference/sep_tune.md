# Sequentially Tuned Moment Correction

Tunes \\\kappa\\ on inner-validation RAFE and \\\tau\\ on
inner-validation C-RAFE *independently* — the bound's two channels are
separable, so each parameter is chosen against its own channel.

## Usage

``` r
sep_tune(
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
  [`mu_rafe_stein()`](https://www.sebastianstoeckl.com/rafe/reference/mu_rafe_stein.md).

## Value

A list with components `mu_tilde`, `Sigma_tilde`, `kappa`, `tau_rel`,
and the inner-validation loss profiles `losses_kappa` and `losses_tau`.

## Details

The training window is split into an inner-training block and an
inner-validation block (40 + 20 months by default). Moments are
estimated on the inner-training block, the realised moments of the
inner-validation block play the role of \\(\mu, \Sigma)\\ in the
metrics, and the selected \\(\kappa, \tau)\\ are then applied to the
moments of the *full* training window.

## References

Stöckl, S., Salcher, L., & Hanke, M. Post-Optimal Moment Correction for
Mean-Variance Portfolios. Working paper.

## See also

[`joint_trafe_tune()`](https://www.sebastianstoeckl.com/rafe/reference/joint_trafe_tune.md)
for the joint search.

## Examples

``` r
set.seed(1)
R <- matrix(rnorm(80 * 6, mean = 0.008, sd = 0.05), 80, 6)
fit <- sep_tune(R)
fit$kappa
#> [1] 0
fit$tau_rel
#> [1] 0.5
```
