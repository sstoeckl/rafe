# RAFE-Stein-Corrected Mean

Positive-part James–Stein shrinkage of a forecast mean toward a target,
in the \\\hat\Sigma^{-1/2}\\-whitened space.

## Usage

``` r
mu_rafe_stein(
  mu_hat,
  Sigma_hat,
  kappa = NULL,
  target = c("grandmean", "zero"),
  T_obs = NULL
)
```

## Arguments

- mu_hat:

  Forecast mean.

- Sigma_hat:

  Forecast covariance. Must be a covariance, never a precision matrix;
  the function inverts it internally.

- kappa:

  Optional shrinkage intensity in `[0, 1]`, where 0 leaves `mu_hat`
  untouched and 1 collapses it onto the target. If `NULL`, the
  James–Stein intensity is computed from `T_obs`.

- target:

  Shrinkage target. One of `"grandmean"` (default; shrinks toward the
  cross-sectional mean of `mu_hat`) or `"zero"`.

- T_obs:

  Number of observations the forecast was estimated from. Needed only
  when `kappa` is `NULL`.

## Value

A numeric vector of corrected means (same length as `mu_hat`), with the
applied intensity attached as attribute `kappa`.

## Details

Writing \\\hat m := \hat\Sigma^{-1/2}\hat\mu\\ and \\m_0 :=
\hat\Sigma^{-1/2}\mu_0\\ for the whitened forecast and target, the
James–Stein map of Salcher, Stöckl & Hanke is \$\$JS(z) = m_0 +
\left(1 - \frac{N-2}{T\\\\z - m_0\\\_2^2}\right)\_{+}(z - m_0),\$\$ and
the corrected mean is \\\tilde\mu = \hat\Sigma^{1/2} JS(\hat m)\\.
Because \\\\\hat m - m_0\\\_2^2 = (\hat\mu - \mu_0)^\top \hat\Sigma^{-1}
(\hat\mu - \mu_0)\\, this is algebraically a linear shrinkage in the
*raw* space whose intensity is computed in the *risk* metric:
\\\tilde\mu = (1 - \kappa)\hat\mu + \kappa\mu_0\\.

Supply `T_obs` to use the James–Stein intensity, or `kappa` to fix the
intensity directly (which is what
[`sep_tune()`](https://www.sebastianstoeckl.com/rafe/reference/sep_tune.md)
and
[`joint_trafe_tune()`](https://www.sebastianstoeckl.com/rafe/reference/joint_trafe_tune.md)
search over). Exactly one of the two is required.

## References

Stöckl, S., Salcher, L., & Hanke, M. Post-Optimal Moment Correction for
Mean-Variance Portfolios. Working paper.

## Examples

``` r
set.seed(1)
R <- matrix(rnorm(60 * 5, mean = 0.01, sd = 0.05), 60, 5)
mu_hat    <- colMeans(R)
Sigma_hat <- stats::cov(R)

# James-Stein intensity, implied by the sample size:
mu_rafe_stein(mu_hat, Sigma_hat, T_obs = 60)
#> [1] 0.01167921 0.01167921 0.01167921 0.01167921 0.01167921
#> attr(,"kappa")
#> [1] 1

# Fixed intensity:
mu_rafe_stein(mu_hat, Sigma_hat, kappa = 0.5)
#> [1] 0.013530016 0.013634742 0.009676794 0.008699989 0.012854528
#> attr(,"kappa")
#> [1] 0.5
```
