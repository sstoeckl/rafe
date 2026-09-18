# rafe

> Evaluate return forecasts by the economic damage their errors cause,
> not by their size.

A forecast of expected returns and their covariance is rarely an end in
itself. It is an input to a decision — most often a portfolio choice.
Yet forecasts are almost always scored with root mean squared error,
which weights every asset equally. A mean-variance optimiser does not:
it leans hardest on exactly the directions where the covariance estimate
is least reliable, so errors that RMSE treats as interchangeable can do
wildly different amounts of economic damage.

`rafe` implements metrics that close that gap, together with estimators
that act on them.

## The bound

Let $`\hat w`$ be the plug-in mean-variance portfolio built from a
forecast $`(\hat\mu, \hat\Sigma)`$, and let
$`\Delta = SR^{*} - SR(\hat w)`$ be the Sharpe ratio it gives up
relative to the oracle. Salcher, Stöckl and Hanke (2026) show that

``` math
\Delta \;\le\; c \cdot \underbrace{\lVert \hat\mu - \mu \rVert_{\Sigma^{-1}}}_{\textrm{RAFE}} \;+\; SR^{*} \cdot \underbrace{\lVert \Sigma^{1/2}\hat\Sigma^{-1}\Sigma^{1/2} - I \rVert_{2}}_{\textrm{C-RAFE}}
```

with $`c = 1`$ when $`\textrm{RAFE} \le SR^{*}`$ and $`c = 2`$
otherwise.
[`compute_trafe()`](https://www.sebastianstoeckl.com/rafe/reference/compute_trafe.md)
returns the published total error
$`\textrm{RAFE} + SR^{*}\cdot\textrm{C-RAFE}`$; pass `c = NULL` for the
guaranteed bound.

The right-hand side separates cleanly into a **mean channel** and a
**covariance channel**. Each is a quantity a practitioner can measure
and, as the second half of this package shows, act on.

## Installation

``` r

# install.packages("remotes")
remotes::install_github("sstoeckl/rafe")
```

## Usage

``` r

library(rafe)
data(ff12)

R     <- as.matrix(ff12[, -1])
train <- R[1:60, ]
eval  <- R[61:120, ]

mu_hat <- colMeans(train); Sigma_hat <- cov(train)
mu     <- colMeans(eval);  Sigma     <- cov(eval)

compute_rafe(mu_hat, mu, Sigma)     # mean channel
compute_crafe(Sigma, Sigma_hat)     # covariance channel

bound <- compute_trafe(mu_hat, mu, Sigma, Sigma_hat)
unlist(attributes(bound))           # bound plus its decomposition
```

If the forecaster cannot be retrained, correct its output before it
reaches the optimiser.
[`sep_tune()`](https://www.sebastianstoeckl.com/rafe/reference/sep_tune.md)
selects a mean-shrinkage intensity $`\kappa`$ and a covariance
eigenvalue floor $`\tau`$ on an inner-validation split:

``` r

fit <- sep_tune(train)
c(kappa = fit$kappa, tau = fit$tau_rel)

compute_crafe(Sigma, fit$Sigma_tilde)   # below the uncorrected distortion
```

The individual correctors are available directly as
[`mu_rafe_stein()`](https://www.sebastianstoeckl.com/rafe/reference/mu_rafe_stein.md)
and
[`sigma_crafe_floor()`](https://www.sebastianstoeckl.com/rafe/reference/sigma_crafe_floor.md).

Argument order follows the replication code of the published paper, so
scripts written against it run unchanged here. Note in particular that
[`compute_crafe()`](https://www.sebastianstoeckl.com/rafe/reference/compute_crafe.md)
takes the **realised** covariance first.

## Articles

- [Getting
  started](https://www.sebastianstoeckl.com/rafe/articles/rafe-getting-started.html)
  — a short tour of the metrics and correctors.
- [Evaluating
  forecasts](https://www.sebastianstoeckl.com/rafe/articles/rafe-evaluation.html)
  — reproduces Table 4 of the published paper, to all three reported
  decimals, using only package functions. As the covariance is
  progressively stripped out of the error measure, its correlation with
  realised economic loss falls from 0.68 to 0.11. The last step in that
  sequence is RMSE.
- [Post-processing
  forecasts](https://www.sebastianstoeckl.com/rafe/articles/rafe-post-processing.html)
  — tuning $`(\kappa, \tau)`$, what correction does to each channel of
  the bound, and what it does to realised portfolios.

## Data

`ff12` ships with the package: monthly excess returns on Kenneth
French’s 12 industry portfolios, January 1964 to December 2023 — the
sample of the published paper, so the vignettes reproduce its results
offline.

## Reproducibility

The package carries 443 unit tests. Beyond ordinary input and edge-case
coverage, three of them check the implementation against theory and
against published code rather than against itself:

- **Table 4 of the paper is reproduced to all three reported decimals**,
  from the raw returns, using only package functions.
- All ten metric variants agree to nine decimal places with verbatim
  transcriptions of the published replication code.
- The Sharpe-gap bound holds in every draw of a Monte Carlo study, and
  the fully restricted metric equals the paper’s RMSE exactly.

## Citation

For the evaluation framework:

> Salcher, L., Stöckl, S., & Hanke, M. (2026). Lost in Translation?
> Risk-Adjusting RMSE for Economic Forecast Performance. *Journal of
> Forecasting*.
> [doi:10.1002/for.70134](https://doi.org/10.1002/for.70134)

For the moment-correction estimators:

> Stöckl, S., Salcher, L., & Hanke, M. Post-Optimal Moment Correction
> for Mean-Variance Portfolios. Working paper.

## License

MIT © Sebastian Stöckl
