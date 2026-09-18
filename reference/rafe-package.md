# rafe: Decision-Aligned Forecast Evaluation and Moment Correction

Evaluates return forecasts by the economic damage their errors cause
rather than by their size, and corrects the forecast moments that a
mean-variance optimiser is most sensitive to.

## Details

Conventional accuracy measures weight every asset equally; a portfolio
optimiser does not. The risk-adjusted forecast error (RAFE) and the
operator-norm covariance forecast error (C-RAFE) reweight forecast
errors by the risk metric the decision actually uses, and combine into
an upper bound on the Sharpe-ratio gap of the plug-in portfolio
(T-RAFE).

Two groups of functions:

- **Metrics**:
  [`compute_rafe()`](https://www.sebastianstoeckl.com/rafe/reference/compute_rafe.md),
  [`compute_crafe()`](https://www.sebastianstoeckl.com/rafe/reference/compute_crafe.md),
  [`compute_trafe()`](https://www.sebastianstoeckl.com/rafe/reference/compute_trafe.md),
  and
  [`restrict_cov()`](https://www.sebastianstoeckl.com/rafe/reference/restrict_cov.md)
  for the nested sequence of covariance restrictions that reduces RAFE
  to RMSE.

- **Moment correction**:
  [`mu_rafe_stein()`](https://www.sebastianstoeckl.com/rafe/reference/mu_rafe_stein.md)
  and
  [`sigma_crafe_floor()`](https://www.sebastianstoeckl.com/rafe/reference/sigma_crafe_floor.md)
  act on the two channels of the bound;
  [`sep_tune()`](https://www.sebastianstoeckl.com/rafe/reference/sep_tune.md)
  and
  [`joint_trafe_tune()`](https://www.sebastianstoeckl.com/rafe/reference/joint_trafe_tune.md)
  select their intensities on an inner-validation split.

The [ff12](https://www.sebastianstoeckl.com/rafe/reference/ff12.md)
dataset, the universe of the published paper, makes every example
reproducible offline.

## References

Salcher, L., Stöckl, S., & Hanke, M. (2026). Lost in Translation?
Risk-Adjusting RMSE for Economic Forecast Performance. *Journal of
Forecasting*. [doi:10.1002/for.70134](https://doi.org/10.1002/for.70134)

Stöckl, S., Salcher, L., & Hanke, M. Post-Optimal Moment Correction for
Mean-Variance Portfolios. Working paper.

## See also

Useful links:

- <https://github.com/sstoeckl/rafe>

- <https://www.sebastianstoeckl.com/rafe/>

- Report bugs at <https://github.com/sstoeckl/rafe/issues>

## Author

**Maintainer**: Sebastian Stöckl <sebastian.stoeckl@uni.li>
([ORCID](https://orcid.org/0000-0002-4196-6093))
