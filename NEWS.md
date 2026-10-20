# imisr 0.1.0

First release: a rewrite of the archived 'IMIS' package (Raftery and Bao).

## Algorithm

* All computations are on the log scale: `imis()` takes `log_likelihood()`
  and `log_prior()`, and importance weights, mixture densities and the
  marginal likelihood are computed with log-sum-exp. Very small likelihoods
  (e.g. `exp(-5000)`) no longer underflow.
* `n_centers` adds several normal components per iteration, spread over
  separate high-weight regions.
* `sampling = "mixture"` draws each iteration from the whole current mixture
  (prior plus all normal components), as in Steele, Raftery and Emond (2006).
  The default, `"incremental"`, draws from the newest components, as in
  Raftery and Bao (2010). Draws are always weighted against the full mixture
  using the realized number of draws per component.
* The optimization stage (`n_optim`, "IMIS-opt") works on the log-posterior,
  copes with points where the posterior is zero, falls back gracefully when
  BFGS or the Hessian fails, and excludes draws near each mode using that
  mode's covariance, as described in Raftery and Bao (2010).
* The likelihood is not evaluated where the prior density is zero.
* `parallel = TRUE` evaluates the likelihood in parallel with 'future.apply',
  in blocks of draws; the backend is chosen with `future::plan()`.
* `catch_errors = TRUE` gives a log-likelihood of `-Inf` to parameter sets
  for which the model fails, instead of stopping.
* Default sample sizes follow Raftery and Bao (2010): `n_initial = 1000 * d`
  and `n_per_component = 100 * d`.
* The marginal likelihood comes with a standard error (Raftery and Bao, 2010,
  eq. 7).

## Interface

* `imis()` returns an `imis` object with `print()`, `summary()`
  (importance-weighted posterior summaries), `plot()` (diagnostics, prior vs
  posterior marginals, pairs), `as.matrix()` and `as.data.frame()` methods.
* `bayes_factor()` compares two models fitted with `imis()`.
* `imis_legacy()` is a drop-in replacement for `IMIS::IMIS()`.
* `vectorized` / `vectorized_prior` let the log-likelihood and log-prior be
  written either for a matrix of parameter sets or for one parameter vector.
* Example models `sick_sicker_markov()` and `crs_markov()`, with calibration
  targets `sick_sicker_targets` and `crs_targets`.
