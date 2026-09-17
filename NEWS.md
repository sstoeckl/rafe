# rafe 0.1.0

First public release.

## Metrics

* `compute_rafe()`, `compute_crafe()` and `compute_trafe()` implement the
  risk-adjusted mean error, the operator-norm precision distortion, and their
  combination into the upper bound on the Sharpe-ratio gap.
* Argument order matches the replication code of Salcher, Stöckl & Hanke
  (2026), so scripts from that repository run unchanged: `compute_rafe()` takes
  the realised covariance as its third positional argument, and
  `compute_crafe()` takes the realised covariance **first**.
* `compute_trafe()` applies the theorem's data-dependent constant by default
  (`c = 1` when RAFE ≤ SR\*, else `c = 2`), and returns the decomposition as
  attributes `rafe`, `crafe`, `c` and `SR_star`. The published replication
  code fixes `c = 1`; pass `c = 1` to reproduce those numbers exactly.

## Moment correction

* `mu_rafe_stein()` applies positive-part James–Stein shrinkage in the
  risk-whitened space. Supply `T_obs` for the James–Stein intensity or `kappa`
  to fix it.
* `sigma_crafe_floor()` raises small eigenvalues to a floor relative to the
  mean eigenvalue.
* `sep_tune()` and `joint_trafe_tune()` select `(kappa, tau)` on an
  inner-validation split, sequentially or jointly. Default grids use 101
  points on `[0, 1]` for kappa and on `[0, 0.5]` for tau.

## Data

* `ff12`: monthly excess returns on the Fama-French 12 industry portfolios,
  January 1964 to December 2023 — the sample of the published paper.

## Vignettes

* `rafe-getting-started`, `rafe-evaluation` (reproduces the published
  mean-variance experiment on 49 industry portfolios) and
  `rafe-post-processing`.
