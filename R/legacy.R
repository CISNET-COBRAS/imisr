#' Drop-in replacement for `IMIS::IMIS()`
#'
#' Runs [imis()] with the interface of the archived 'IMIS' package, so that
#' existing scripts keep working: `prior()` and `likelihood()` return
#' densities on the natural (not log) scale, and by default the three user
#' functions are looked up by name in the calling environment, as `IMIS()`
#' did.
#'
#' New code should call [imis()] with log-scale functions instead. Here the
#' likelihood is computed on the natural scale and then logged, so
#' likelihoods smaller than about `1e-308` underflow to zero, just as they did
#' in the original package.
#'
#' @param B Number of draws per new mixture component.
#' @param B.re Size of the posterior resample.
#' @param number_k Maximum number of stages (the initial stage counts as one).
#' @param D Number of optimizers (`0` skips the optimization stage).
#' @param sample.prior,prior,likelihood Functions that draw from the prior,
#'   and evaluate the prior density and the likelihood of a matrix of
#'   parameter sets.
#' @return A list with the components of the original `IMIS()` output,
#'   `stat` (matrix with one row per stage and columns `MargLike`,
#'   `UniquePoint`, `MaxWeight`, `ESS`, `Entropy` and `VarWeights`),
#'   `resample` and `center`, plus `fit`, the full [imis()] result.
#' @examples
#' sample.prior <- function(n) cbind(mu = rnorm(n))
#' prior <- function(x) dnorm(x[, 1])
#' likelihood <- function(x) dnorm(1, x[, 1], 0.5)
#' res <- imis_legacy(B = 100, B.re = 1000, number_k = 20)
#' colMeans(res$resample)
#' @export
imis_legacy <- function(B = 1000, B.re = 3000, number_k = 100, D = 0,
                        sample.prior = get("sample.prior",
                                           envir = parent.frame(),
                                           mode = "function"),
                        prior = get("prior", envir = parent.frame(),
                                    mode = "function"),
                        likelihood = get("likelihood", envir = parent.frame(),
                                         mode = "function")) {
  force(sample.prior)
  force(prior)
  force(likelihood)
  fit <- imis(
    log_likelihood = function(x) log(likelihood(x)),
    log_prior = function(x) log(prior(x)),
    sample_prior = sample.prior,
    n_initial = 10 * B, n_per_component = B, n_resample = B.re,
    max_iter = max(0, number_k - 1), n_optim = D
  )
  df <- fit$diagnostics
  m_stat <- cbind(MargLike = df$log_ml, UniquePoint = df$expected_unique,
                  MaxWeight = df$max_weight, ESS = df$ess,
                  Entropy = df$entropy, VarWeights = df$var_weights)
  list(stat = m_stat, resample = fit$resample, center = fit$centers,
       fit = fit)
}
