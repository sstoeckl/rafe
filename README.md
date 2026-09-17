# rafe

<!-- badges: start -->
[![R-CMD-check](https://github.com/sstoeckl/rafe/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/sstoeckl/rafe/actions/workflows/R-CMD-check.yaml)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
<!-- badges: end -->

> Evaluate return forecasts by the economic damage their errors cause, not by their size.

A forecast of expected returns and their covariance is rarely an end in itself. It is an input to a decision — most often a portfolio choice. Yet forecasts are almost always scored with root mean squared error, which weights every asset equally. A mean-variance optimiser does not: it leans hardest on exactly the directions where the covariance estimate is least reliable, so errors that RMSE treats as interchangeable can do wildly different amounts of economic damage.

`rafe` implements metrics that close that gap, together with estimators that act on them.

## The bound

Let $\hat w$ be the plug-in mean-variance portfolio built from a forecast $(\hat\mu, \hat\Sigma)$, and let $\Delta = SR^{*} - SR(\hat w)$ be the Sharpe ratio it gives up relative to the oracle. Salcher, Stöckl and Hanke (2026) show that

$$\Delta \;\le\; c \cdot \underbrace{\lVert \hat\mu - \mu \rVert_{\Sigma^{-1}}}_{\textrm{RAFE}} \;+\; SR^{*} \cdot \underbrace{\lVert \Sigma^{1/2}\hat\Sigma^{-1}\Sigma^{1/2} - I \rVert_{2}}_{\textrm{C-RAFE}}$$

with $c = 1$ when $\textrm{RAFE} \le SR^{*}$ and $c = 2$ otherwise.

The right-hand side separates cleanly into a **mean channel** and a **covariance channel**. Each is a quantity a practitioner can measure and, as the second half of this package shows, act on.

## Installation

```r
# install.packages("remotes")
remotes::install_github("sstoeckl/rafe")
```

## Usage

```r
library(rafe)
data(ff49)

R     <- as.matrix(ff49[, -1])
train <- R[1:60, ]
eval  <- R[61:120, ]

mu_hat <- colMeans(train); Sigma_hat <- cov(train)
mu     <- colMeans(eval);  Sigma     <- cov(eval)

compute_rafe(mu_hat, mu, Sigma)     # mean channel
compute_crafe(Sigma, Sigma_hat)     # covariance channel

bound <- compute_trafe(mu_hat, mu, Sigma, Sigma_hat)
unlist(attributes(bound))           # bound plus its decomposition
```

If the forecaster cannot be retrained, correct its output before it reaches the optimiser. `sep_tune()` selects a mean-shrinkage intensity $\kappa$ and a covariance eigenvalue floor $\tau$ on an inner-validation split:

```r
fit <- sep_tune(train)
c(kappa = fit$kappa, tau = fit$tau_rel)

compute_crafe(Sigma, fit$Sigma_tilde)   # below the uncorrected distortion
```

The individual correctors are available directly as `mu_rafe_stein()` and `sigma_crafe_floor()`.

Argument order follows the replication code of the published paper, so scripts written against it run unchanged here. Note in particular that `compute_crafe()` takes the **realised** covariance first.

## Articles

- [Getting started](https://sstoeckl.github.io/rafe/articles/rafe-getting-started.html) — a short tour of the metrics and correctors.
- [Evaluating forecasts](https://sstoeckl.github.io/rafe/articles/rafe-evaluation.html) — reproduces the paper's mean-variance experiment on 49 industry portfolios, using forecasts of controlled quality. Root mean squared error falls monotonically as the forecasts improve, while T-RAFE rises: the accuracy gain and the covariance damage arrive together, and only one of them is visible to RMSE.
- [Post-processing forecasts](https://sstoeckl.github.io/rafe/articles/rafe-post-processing.html) — tuning $(\kappa, \tau)$, what correction does to each channel of the bound, and what it does to realised portfolios.

## Data

Two datasets ship with the package so that every example is reproducible offline: `ff12` and `ff49`, monthly excess returns on Kenneth French's 12 and 49 industry portfolios.

## Reproducibility

The package carries 443 unit tests. Beyond ordinary input and edge-case coverage, three of them check the implementation against theory and against published code rather than against itself:

- At $\Sigma = I$, RAFE collapses to $\sqrt{N}$ times RMSE. The metric carries no $1/N$ normalisation, so the factor is part of the claim.
- The Sharpe-gap bound holds in every draw of a Monte Carlo study.
- `compute_rafe()` and `compute_crafe()` agree to nine decimal places with verbatim transcriptions of the published replication code, and C-RAFE matches the closed form implied by that paper's experimental design.

## Citation

For the evaluation framework:

> Salcher, T., Stöckl, S., & Hanke, M. (2026). Lost in Translation? Risk-Adjusting RMSE for Economic Forecast Performance. *Journal of Forecasting*. [doi:10.1002/for.70134](https://doi.org/10.1002/for.70134)

For the moment-correction estimators:

> Stöckl, S., Salcher, T., & Hanke, M. Post-Optimal Moment Correction for Mean-Variance Portfolios. Working paper.

## License

MIT © Sebastian Stöckl, Tobias Salcher, Michael Hanke
