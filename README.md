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

mu_hat    <- rnorm(10)
Sigma_hat <- diag(10)
mu        <- rnorm(10)
Sigma     <- diag(10)

# Metrics
compute_rafe(mu_hat, mu, Sigma = Sigma)
compute_crafe(Sigma_hat, Sigma)
compute_trafe(mu_hat, mu, Sigma_hat, Sigma)

# Post-processing
mu_tilde    <- mu_rafe_stein(mu_hat, Sigma_hat)
Sigma_tilde <- sigma_crafe_floor(Sigma_hat)

# Training objective (XGBoost) — Paper B
# obj <- xgb_trafe_objective(Sigma_hat)
# xgb.train(..., obj = obj, ...)
```

See `vignettes/rafe-getting-started.Rmd` (once populated) for a guided tour.

## Citation

If you use this package, please cite Salcher, Stöckl & Hanke (2026, *Journal of Forecasting*) for the framework, and the relevant cluster paper (A / B / C) for the specific method you use.

## License

MIT.
