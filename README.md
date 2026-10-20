# imisr <img src="man/figures/logo.png" align="right" height="139" alt="imisr hex logo" />

<!-- badges: start -->
[![R-CMD-check](https://github.com/FerAlLab-git/imisr/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/FerAlLab-git/imisr/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

**imisr** implements incremental mixture importance sampling (IMIS) for
Bayesian calibration of simulation models and for estimating marginal
likelihoods (model evidence). It is a rewrite of the archived
[IMIS](https://cran.r-project.org/package=IMIS) package (Raftery and Bao,
2010) with:

- **Log-scale computation throughout**, so likelihoods such as
  $e^{-2400}$ no longer underflow to zero.
- **Several components per iteration** (`n_centers`), and a choice between
  drawing from the newest components (`sampling = "incremental"`, Raftery and
  Bao, 2010) or from the **whole current mixture**
  (`sampling = "mixture"`, Steele, Raftery and Emond, 2006).
- **Parallel likelihood evaluation** through the
  [future](https://future.futureverse.org) framework: multicore,
  multisession, clusters and HPC schedulers.
- **Robust optimization stage** ("IMIS-opt", `n_optim`) that copes with
  parameter values where the posterior is zero, and with non-positive-definite
  Hessians.
- **No wasted model runs**: the likelihood is never evaluated where the prior
  is zero, and models that fail for some inputs can be handled
  (`catch_errors = TRUE`).
- **Marginal likelihood with a standard error**, and `bayes_factor()` to
  compare models.
- `print()`, `summary()`, `plot()` and `as.matrix()` methods, and a drop-in
  `imis_legacy()` for scripts written for `IMIS::IMIS()`.
- Example **state-transition models** with calibration targets
  (Sick-Sicker and cancer relative survival).

## Installation

```r
# install.packages("remotes")
remotes::install_github("FerAlLab-git/imisr")
```

## Example

Calibrate the three unknown inputs of the Sick-Sicker model to observed
survival, prevalence and the proportion of the sick in the "sick" state:

```r
library(imisr)

lst_targets <- sick_sicker_targets
v_target_names <- names(lst_targets)
v_param_names <- c("p_S1S2", "hr_S1", "hr_S2")
v_lb <- c(p_S1S2 = 0.01, hr_S1 = 1, hr_S2 = 1)  # lower bounds
v_ub <- c(p_S1S2 = 0.50, hr_S1 = 12, hr_S2 = 25) # upper bounds

# Prior: independent uniforms, one column per parameter, indexed by name
sample_prior <- function(n_samp) {
  m_param_samp <- matrix(NA_real_, nrow = n_samp,
                         ncol = length(v_param_names),
                         dimnames = list(NULL, v_param_names))
  for (s_par in v_param_names) {
    m_param_samp[, s_par] <- runif(n_samp, v_lb[[s_par]], v_ub[[s_par]])
  }
  m_param_samp
}
calc_log_prior <- function(v_params) {
  if (is.null(dim(v_params))) v_params <- t(v_params)
  v_lprior <- rep(0, nrow(v_params))
  for (s_par in v_param_names) {
    v_lprior <- v_lprior + dunif(v_params[, s_par], v_lb[[s_par]],
                                 v_ub[[s_par]], log = TRUE)
  }
  v_lprior
}

# Likelihood: targets normally distributed around the model outputs. The
# model runs one parameter set (row) at a time; parameter sets for which it
# fails get a log-likelihood of -Inf.
calc_log_lik <- function(v_params) {
  if (is.null(dim(v_params))) v_params <- t(v_params)
  m_llik <- matrix(-Inf, nrow = nrow(v_params), ncol = length(v_target_names),
                   dimnames = list(NULL, v_target_names))
  for (j in seq_len(nrow(v_params))) {
    l_out <- tryCatch(sick_sicker_markov(v_params[j, ]),
                      error = function(e) NULL)
    if (is.null(l_out)) next
    m_llik[j, "Surv"] <- sum(dnorm(lst_targets$Surv$value, l_out$Surv,
                                   lst_targets$Surv$se, log = TRUE))
    m_llik[j, "Prev"] <- sum(dnorm(lst_targets$Prev$value, l_out$Prev,
                                   lst_targets$Prev$se, log = TRUE))
    m_llik[j, "PropSick"] <- sum(dnorm(lst_targets$PropSick$value,
                                       l_out$PropSick,
                                       lst_targets$PropSick$se, log = TRUE))
  }
  rowSums(m_llik)
}

set.seed(2026)
fit <- imis(calc_log_lik, calc_log_prior, sample_prior)
fit
#> Incremental mixture importance sampling
#>   Converged after 16 iteration(s); 7800 draws, 16 normal component(s)
#>   7776 likelihood evaluations in 1.6 s
#>   ESS = 2961.7; expected unique points in resample = 1927 / 3000
#>   log marginal likelihood = 168.0653 (SE 0.0145)
#>
#> Posterior summary (importance-weighted):
#>  parameter     mean      sd    2.5%      50%    97.5%
#>     p_S1S2  0.07943 0.00483 0.07012  0.07937  0.08902
#>      hr_S1  3.92900 0.83590 2.26600  3.93900  5.51700
#>      hr_S2 11.22000 1.33700 8.73600 11.20000 13.88000

plot(fit)                      # convergence diagnostics
plot(fit, type = "pairs")      # posterior draws and mixture centers
m_post <- as.matrix(fit)       # equally weighted posterior draws for PSA
```

For a slow model, run the likelihood evaluations in parallel:

```r
future::plan(future::multisession, workers = 8)
fit <- imis(calc_log_lik, calc_log_prior, sample_prior, parallel = TRUE)
```

See `vignette("imisr")` to get started, and
`vignette("calibrating-simulation-models")` for the full calibration
workflow, including a non-identifiable model with a bimodal posterior.

## How IMIS works

1. Sample from the prior and weight each draw by its likelihood.
2. *(Optional)* Find local posterior modes by optimization and start the
   mixture with a normal component at each.
3. Repeatedly add a normal component centered at the highest-weight draw,
   with covariance estimated from its nearest neighbors; draw from it; and
   reweight **all** draws against the defensive mixture of the prior and every
   component.
4. Stop when the weights are close to uniform (the expected number of unique
   points in a resample reaches 63%, the value expected when the sampling
   function equals the posterior), then resample.

## References

- Raftery AE, Bao L (2010). Estimating and projecting trends in HIV/AIDS
  generalized epidemics using incremental mixture importance sampling.
  *Biometrics*, 66(4), 1162–1173.
  [doi:10.1111/j.1541-0420.2010.01399.x](https://doi.org/10.1111/j.1541-0420.2010.01399.x)
- Steele RJ, Raftery AE, Emond MJ (2006). Computing normalizing constants for
  finite mixture models via incremental mixture importance sampling (IMIS).
  *Journal of Computational and Graphical Statistics*, 15(3), 712–734.
  [doi:10.1198/106186006X132358](https://doi.org/10.1198/106186006X132358)
- Menzies NA, Soeteman DI, Pandya A, Kim JJ (2017). Bayesian methods for
  calibrating health policy models: a tutorial. *PharmacoEconomics*, 35(6),
  613–624.
  [doi:10.1007/s40273-017-0494-4](https://doi.org/10.1007/s40273-017-0494-4)
