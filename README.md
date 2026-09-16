# rafe

> Decision-aligned forecast evaluation and moment correction, in R.

**S1 of the RAFE Program.** Companion implementation of the methods described in:

- Salcher, Stöckl & Hanke (2026), *Journal of Forecasting* — the RAFE / C-RAFE / T-RAFE bound on the Sharpe-ratio gap of plug-in mean-variance portfolios.
- Stöckl, Salcher & Hanke (2026/27, Paper A) — post-optimal moment correction.
- Stöckl et al. (2027, Paper B) — decision-aware training of financial forecasters.

## Status

🚧 **Pre-alpha.** API is being scaffolded. Not yet on CRAN. GitHub preview release planned for 2026 Q3; CRAN submission targeted for 2027 Q2.

## Installation (when available)

```r
# CRAN (after 2027 Q2)
install.packages("rafe")

# GitHub preview (from 2026 Q3 onward)
remotes::install_github("sstoeckl/rafe")
```

## Quickstart

```r
library(rafe)
data(ff12)

R     <- as.matrix(ff12[, -1])
train <- R[1:60, ]
eval  <- R[61:120, ]

mu_hat <- colMeans(train); Sigma_hat <- cov(train)
mu     <- colMeans(eval);  Sigma     <- cov(eval)

# Metrics
compute_rafe(mu_hat, mu, Sigma = Sigma)
compute_crafe(Sigma_hat, Sigma)

# The bound, with its decomposition attached as attributes
bound <- compute_trafe(mu_hat, mu, Sigma_hat, Sigma)
unlist(attributes(bound))

# Post-processing: tune (kappa, tau) on an inner-validation split
fit <- sep_tune(train)
c(kappa = fit$kappa, tau = fit$tau_rel)

compute_crafe(fit$Sigma_tilde, Sigma)   # lower than the uncorrected value

# Or apply a correction directly
mu_rafe_stein(mu_hat, Sigma_hat, T_obs = 60)   # James-Stein intensity
sigma_crafe_floor(Sigma_hat, tau_rel = 0.3)
```

Three vignettes:

- `vignette("rafe-getting-started")` — five-minute tour.
- `vignette("rafe-evaluation")` — the Part 1 evaluation example on FF-12.
- `vignette("rafe-post-processing")` — the Paper A correction example on FF-12.

## Citation

If you use this package, please cite Salcher, Stöckl & Hanke (2026, *Journal of Forecasting*) for the framework, and the relevant cluster paper (A / B / C) for the specific method you use.

## License

MIT.
