# 0. rafe — getting started

``` r

library(rafe)
```

`rafe` implements decision-aligned forecast evaluation: metrics that
score a return forecast by the damage its errors do to a portfolio,
rather than by their size.

## The idea in one paragraph

Feed a forecast $`(\hat\mu, \hat\Sigma)`$ into a mean-variance optimiser
and you get a portfolio whose Sharpe ratio falls short of the oracle’s
by some gap $`\Delta`$. Salcher, Stöckl & Hanke (2026) show that gap is
bounded by two terms: a mean error measured in the risk metric, and a
covariance distortion measured in the operator norm.

``` math
\Delta \;\le\; c \cdot \underbrace{\lVert\hat\mu - \mu\rVert_{\Sigma^{-1}}}_{\mathrm{RAFE}}
 \;+\; SR^{*} \cdot \underbrace{\lVert\Sigma^{1/2}\hat\Sigma^{-1}\Sigma^{1/2} - I\rVert_2}_{\mathrm{C\text{-}RAFE}}
```

RMSE appears nowhere. That is the point: it weights every asset equally,
while the optimiser does not.

## The metrics

``` r

data(ff12)
R <- as.matrix(ff12[, -1])

train <- R[1:60, ]
eval  <- R[61:120, ]

mu_hat <- colMeans(train); Sigma_hat <- stats::cov(train)
mu     <- colMeans(eval);  Sigma     <- stats::cov(eval)

compute_rafe(mu_hat, mu, Sigma = Sigma)
#> [1] 0.9526
compute_crafe(Sigma, Sigma_hat)
#> [1] 4.557
```

[`compute_trafe()`](https://sstoeckl.github.io/rafe/reference/compute_trafe.md)
combines them into the paper’s total error and returns the decomposition
as attributes. Pass `c = NULL` for the guaranteed upper bound, which
carries a data-dependent constant ($`c = 1`$ if
$`\mathrm{RAFE} \le SR^{*}`$, else $`c = 2`$):

``` r

bound <- compute_trafe(mu_hat, mu, Sigma, Sigma_hat)
bound
#> [1] 3.594
#> attr(,"rafe")
#> [1] 0.9526
#> attr(,"crafe")
#> [1] 4.557
#> attr(,"c")
#> [1] 1
#> attr(,"SR_star")
#> [1] 0.5796
unlist(attributes(bound))
#>    rafe   crafe       c SR_star 
#>  0.9526  4.5570  1.0000  0.5796
```

## From RAFE down to RMSE

The paper’s framing is that RMSE is RAFE after the covariance has been
stripped of everything that makes it informative.
[`restrict_cov()`](https://sstoeckl.github.io/rafe/reference/restrict_cov.md)
performs that sequence and every metric takes a matching `variant`:

``` r

vapply(c("none", "cc05", "cc0", "cv", "i"),
       function(v) compute_rafe(mu_hat, mu, Sigma, variant = v), numeric(1))
#>    none    cc05     cc0      cv       i 
#> 0.95257 0.60464 0.77102 0.82509 0.04287
```

The last entry is the RMSE of the paper,
$`\lVert\hat\mu - \mu\rVert_2`$.

## The correctors

If you cannot retrain the forecaster, correct its output.
[`sep_tune()`](https://sstoeckl.github.io/rafe/reference/sep_tune.md)
picks a mean-shrinkage intensity $`\kappa`$ and a covariance eigenvalue
floor $`\tau`$ on an inner-validation split, then applies both:

``` r

fit <- sep_tune(train, kappa_grid = seq(0, 1, by = 0.1),
                       tau_grid   = seq(0, 0.5, by = 0.1))
c(kappa = fit$kappa, tau = fit$tau_rel)
#> kappa   tau 
#>   1.0   0.5

c(rafe_before  = compute_rafe(mu_hat, mu, Sigma = Sigma),
  rafe_after   = compute_rafe(fit$mu_tilde, mu, Sigma = Sigma),
  crafe_before = compute_crafe(Sigma, Sigma_hat),
  crafe_after  = compute_crafe(Sigma, fit$Sigma_tilde))
#>  rafe_before   rafe_after crafe_before  crafe_after 
#>       0.9526       0.5636       4.5570       2.2056
```

## Where to go next

- [`vignette("rafe-evaluation", package = "rafe")`](https://sstoeckl.github.io/rafe/articles/rafe-evaluation.md)
  — reproduces Table 4 of the published paper from the raw returns, and
  shows how much explanatory power each simplification of the covariance
  costs on the way down to RMSE.
- [`vignette("rafe-post-processing", package = "rafe")`](https://sstoeckl.github.io/rafe/articles/rafe-post-processing.md)
  — correcting moments you were handed: tuning $`(\kappa, \tau)`$, what
  the correction does to both channels of the bound, and what it does to
  portfolios.

## References

Salcher, L., Stöckl, S., & Hanke, M. (2026). Lost in Translation?
Risk-Adjusting RMSE for Economic Forecast Performance. *Journal of
Forecasting*.
