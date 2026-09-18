# Package index

## Metrics

Score a forecast by the damage its errors do to a mean-variance
decision. Argument order matches the replication code of Salcher, Stöckl
& Hanke (2026).

- [`compute_rafe()`](https://sstoeckl.github.io/rafe/reference/compute_rafe.md)
  : Risk-Adjusted Forecast Error (RAFE)
- [`compute_crafe()`](https://sstoeckl.github.io/rafe/reference/compute_crafe.md)
  : Covariance Forecast Error (C-RAFE)
- [`compute_trafe()`](https://sstoeckl.github.io/rafe/reference/compute_trafe.md)
  : Total RAFE (T-RAFE) — Sharpe-Gap Upper Bound
- [`restrict_cov()`](https://sstoeckl.github.io/rafe/reference/restrict_cov.md)
  : Restricted Covariance Matrices for the Nested Metric Sequence

## Moment correction

Correct forecast moments you were handed, before they reach the
optimiser.

- [`mu_rafe_stein()`](https://sstoeckl.github.io/rafe/reference/mu_rafe_stein.md)
  : RAFE-Stein-Corrected Mean
- [`sigma_crafe_floor()`](https://sstoeckl.github.io/rafe/reference/sigma_crafe_floor.md)
  : C-RAFE-Tuned Eigenvalue Floor for Covariance

## Tuning

Choose the shrinkage intensity and eigenvalue floor on an
inner-validation split.

- [`sep_tune()`](https://sstoeckl.github.io/rafe/reference/sep_tune.md)
  : Sequentially Tuned Moment Correction
- [`joint_trafe_tune()`](https://sstoeckl.github.io/rafe/reference/joint_trafe_tune.md)
  : Jointly T-RAFE-Tuned Moment Correction

## Data

The Fama-French industry portfolios of the published paper.

- [`ff12`](https://sstoeckl.github.io/rafe/reference/ff12.md) :
  Fama-French 12 Industry Portfolios, Monthly Excess Returns

## Package

- [`rafe`](https://sstoeckl.github.io/rafe/reference/rafe-package.md)
  [`rafe-package`](https://sstoeckl.github.io/rafe/reference/rafe-package.md)
  : rafe: Decision-Aligned Forecast Evaluation and Moment Correction
