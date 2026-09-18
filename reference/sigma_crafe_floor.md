# C-RAFE-Tuned Eigenvalue Floor for Covariance

Cleans `Sigma_hat` by raising small eigenvalues to a floor set relative
to the mean eigenvalue: with \\\hat\Sigma =
Q\\\mathrm{diag}(\hat\lambda)\\Q^\top\\, \$\$\tilde\Sigma(\tau) =
Q\\\mathrm{diag}\bigl(\max(\hat\lambda_i,\\
\tau\\\bar\lambda)\bigr)\\Q^\top,\$\$ where \\\bar\lambda\\ is the mean
of the sample eigenvalues. The relative parameterisation makes `tau_rel`
scale-free.

## Usage

``` r
sigma_crafe_floor(Sigma_hat, tau_rel = NULL)
```

## Arguments

- Sigma_hat:

  Forecast covariance.

- tau_rel:

  Relative floor, typically in `[0, 0.5]`. A value of 0 returns
  `Sigma_hat` symmetrised and unchanged.

## Value

A cleaned covariance matrix (same dimension as `Sigma_hat`), with the
applied absolute floor attached as attribute `tau_abs`.

## Details

`tau_rel` is tuned on inner-validation C-RAFE by
[`sep_tune()`](https://www.sebastianstoeckl.com/rafe/reference/sep_tune.md)
and
[`joint_trafe_tune()`](https://www.sebastianstoeckl.com/rafe/reference/joint_trafe_tune.md);
this function applies a given floor.

## References

Stöckl, S., Salcher, L., & Hanke, M. Post-Optimal Moment Correction for
Mean-Variance Portfolios. Working paper.

## Examples

``` r
set.seed(1)
R <- matrix(rnorm(30 * 5), 30, 5)
Sigma_hat <- stats::cov(R)

eigen(Sigma_hat, only.values = TRUE)$values
#> [1] 1.1935890 0.9484871 0.8390561 0.7601640 0.2978117
eigen(sigma_crafe_floor(Sigma_hat, 0.3), only.values = TRUE)$values
#> [1] 1.1935890 0.9484871 0.8390561 0.7601640 0.2978117
```
