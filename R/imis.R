#' Incremental mixture importance sampling
#'
#' Samples from the posterior distribution \eqn{p(\theta \mid y) \propto
#' p(\theta) L(\theta)} and estimates the marginal likelihood \eqn{p(y) = \int
#' p(\theta) L(\theta) d\theta} with incremental mixture importance sampling
#' (IMIS; Steele, Raftery and Emond, 2006; Raftery and Bao, 2010).
#'
#' @section Algorithm:
#' 1. **Initial stage.** Draw `n_initial` points from the prior and compute
#'    their importance weights, which are proportional to the likelihood.
#' 2. **Optimization stage** (only when `n_optim > 0`, "IMIS-opt"). Run
#'    `n_optim` local optimizations of the log-posterior, each started from the
#'    highest-weight draw not already near a mode found earlier, and add a
#'    normal component centered at each mode with covariance equal to the
#'    inverse Hessian. This lets the sampler reach modes the prior sample
#'    missed.
#' 3. **Importance sampling stage.** Repeat until the stopping rule holds or
#'    `max_iter` is reached: add `n_centers` normal components centered at the
#'    highest-weight draws, each with the weighted covariance of the
#'    `n_per_component` draws nearest to it; draw new points; and recompute
#'    every importance weight against the *whole* mixture,
#'    \deqn{w_i \propto \frac{p(\theta_i) L(\theta_i)}{q(\theta_i)}, \qquad
#'    q(\theta) = \frac{N_0}{N} p(\theta) + \sum_c \frac{N_c}{N}
#'    \phi(\theta; \mu_c, \Sigma_c),}
#'    where \eqn{N_c} is the number of draws taken so far from component
#'    \eqn{c}.
#' 4. **Resample.** Draw `n_resample` points with replacement, with
#'    probabilities equal to the final importance weights.
#'
#' Because the prior is always part of the mixture (a "defensive" mixture),
#' the importance weights are bounded and the marginal-likelihood estimator is
#' unbiased, consistent and asymptotically normal (Raftery and Bao, 2010,
#' Theorem 1).
#'
#' @section Stopping rule:
#' IMIS stops when the expected number of unique points in a resample of size
#' `n_resample` reaches `target_unique * n_resample`. The default target,
#' \eqn{1 - e^{-1} \approx 0.632}, is the value expected if the importance
#' sampling function were exactly the posterior.
#'
#' @section Choosing how new draws are generated (`sampling`):
#' * `"incremental"` (default; Raftery and Bao, 2010): each iteration takes
#'   `n_per_component` draws from each *newly added* component only.
#' * `"mixture"` (Steele, Raftery and Emond, 2006): each iteration takes
#'   `n_centers * n_per_component` draws from the *whole current mixture*:
#'   a fraction `prior_weight` from the prior and the rest split equally
#'   across all normal components added so far. Earlier components keep
#'   receiving draws, which helps when the posterior has several separated
#'   regions, at the cost of re-sampling regions that are already well
#'   covered.
#'
#' In both cases every draw is weighted against the full mixture using the
#' realized number of draws per component, so the estimator stays unbiased.
#'
#' @section Parallel computing:
#' The likelihood is evaluated in parallel when `parallel = TRUE`, using
#' \pkg{future.apply}. Choose the backend before calling `imis()`, e.g.
#' `future::plan(future::multisession, workers = 4)`. The draws to evaluate
#' are split into `n_chunks` blocks (one per worker by default), so overhead
#' is paid once per block, not once per draw. Parallelism helps when a single
#' model run is slow (milliseconds or more); for fast models, the overhead of
#' sending work to other processes makes it slower than running serially.
#' Random numbers on the workers are reproducible (`future.seed = TRUE`).
#'
#' @param log_likelihood Function returning the log-likelihood. With
#'   `vectorized = TRUE` (default) it receives a matrix with one parameter set
#'   per row (and column names) and must return a numeric vector with one
#'   value per row. With `vectorized = FALSE` it receives one named numeric
#'   vector and returns a single number. Return `-Inf` for parameter sets that
#'   are impossible.
#' @param log_prior Function returning the log prior density, with the same
#'   calling convention as `log_likelihood`. It must be a *normalized* density
#'   (e.g. `dunif(..., log = TRUE)`, not just a value proportional to it) for
#'   the marginal likelihood to be correct, and it must be consistent with
#'   `sample_prior`. The likelihood is not evaluated where the prior is zero.
#' @param sample_prior Function of one argument `n` returning an `n`-row
#'   matrix of draws from the prior (a vector is accepted when there is a
#'   single parameter). Column names are used as parameter names.
#' @param n_initial Number of prior draws at the initial stage (\eqn{N_0}).
#'   Defaults to `1000 * d`, where `d` is the number of parameters.
#' @param n_per_component Number of draws per new mixture component (\eqn{B}).
#'   It is also the number of nearest neighbors used to estimate each new
#'   component's covariance. Defaults to `100 * d`.
#' @param n_resample Size of the posterior resample (\eqn{J}).
#' @param max_iter Maximum number of importance sampling iterations.
#' @param n_optim Number of optimizers in the optimization stage
#'   (\eqn{D}). `0` (default) skips the stage; Raftery and Bao (2010)
#'   recommend 10 when the posterior may be multimodal.
#' @param n_centers Number of normal components added per iteration.
#' @param sampling How new draws are generated, `"incremental"` or
#'   `"mixture"`. See *Choosing how new draws are generated*.
#' @param prior_weight Share of each iteration's draws that come from the
#'   prior when `sampling = "mixture"`.
#' @param target_unique Stopping threshold, as a fraction of `n_resample`.
#'   See *Stopping rule*.
#' @param vectorized Whether `log_likelihood` takes a matrix of parameter
#'   sets (`TRUE`) or a single named parameter vector (`FALSE`). Simulation
#'   models usually run one parameter set at a time: either loop over the
#'   rows inside `log_likelihood` (see the vignettes), or write it for one
#'   parameter vector and set `vectorized = FALSE`.
#' @param vectorized_prior Same as `vectorized`, for `log_prior`. Defaults to
#'   the value of `vectorized`.
#' @param parallel Whether to evaluate the likelihood in parallel with
#'   \pkg{future.apply}. See *Parallel computing*.
#' @param n_chunks Number of blocks the draws are split into when
#'   `parallel = TRUE`. Defaults to [future::nbrOfWorkers()].
#' @param catch_errors If `TRUE`, parameter sets for which `log_likelihood`
#'   throws an error get a log-likelihood of `-Inf` (with a warning reporting
#'   how many) instead of stopping IMIS. Useful when some parameter
#'   combinations make a model invalid.
#' @param optim_control List with the maximum number of iterations for the
#'   optimization stage: `maxit_nm` (Nelder-Mead, default 1000) and
#'   `maxit_bfgs` (BFGS, default 1000).
#' @param verbose Whether to print progress.
#'
#' @return An object of class `imis`, a list with:
#'   \item{resample}{`n_resample` x `d` matrix of approximate posterior draws.}
#'   \item{draws}{All draws, one row per draw.}
#'   \item{weights}{Normalized importance weights of `draws`. Weighted
#'     summaries over `draws` are more precise than summaries of `resample`.}
#'   \item{log_weights, log_prior, log_lik}{Unnormalized log importance
#'     weights, log prior and log-likelihood of each draw.}
#'   \item{log_ml, log_ml_se}{Estimated log marginal likelihood and its
#'     approximate standard error.}
#'   \item{diagnostics}{Data frame with one row per stage: number of draws
#'     and components, log marginal likelihood and standard error, effective
#'     sample size (ESS), maximum weight, expected number of unique points
#'     in the resample, relative entropy and variance of the weights,
#'     cumulative likelihood evaluations and elapsed seconds.}
#'   \item{components}{List of the normal components (mean, covariance,
#'     number of draws, origin `"optim"` or `"weights"`, iteration).}
#'   \item{centers}{Matrix of component means.}
#'   \item{n_prior_draws}{Number of draws taken from the prior.}
#'   \item{n_initial}{Number of draws in the initial stage.}
#'   \item{converged}{Whether the stopping rule was met.}
#'   \item{n_iter}{Number of importance sampling iterations run.}
#'   \item{n_lik_evals}{Total number of likelihood evaluations, including
#'     those made by the optimizers.}
#'   \item{settings, call, elapsed}{Settings used, the call and run time
#'     in seconds.}
#'
#' @references
#' Raftery AE, Bao L (2010). Estimating and projecting trends in HIV/AIDS
#' generalized epidemics using incremental mixture importance sampling.
#' *Biometrics*, 66(4), 1162-1173. \doi{10.1111/j.1541-0420.2010.01399.x}
#'
#' Steele RJ, Raftery AE, Emond MJ (2006). Computing normalizing constants
#' for finite mixture models via incremental mixture importance sampling
#' (IMIS). *Journal of Computational and Graphical Statistics*, 15(3),
#' 712-734. \doi{10.1198/106186006X132358}
#'
#' Menzies NA, Soeteman DI, Pandya A, Kim JJ (2017). Bayesian methods for
#' calibrating health policy models: a tutorial. *PharmacoEconomics*, 35(6),
#' 613-624. \doi{10.1007/s40273-017-0494-4}
#'
#' @seealso [summary.imis()], [plot.imis()], [bayes_factor()],
#'   [imis_legacy()] for the interface of the original 'IMIS' package.
#'
#' @examples
#' # Normal likelihood with a normal prior: the posterior and the marginal
#' # likelihood are known exactly.
#' y <- c(1.2, 0.8, 1.9, 1.4)
#' fit <- imis(
#'   # outer() gives y_i - mu with one column per parameter set, so the
#'   # log-likelihood of every parameter set is computed at once
#'   log_likelihood = function(theta) {
#'     colSums(dnorm(outer(y, theta[, "mu"], "-"), log = TRUE))
#'   },
#'   log_prior = function(theta) dnorm(theta[, "mu"], 0, 2, log = TRUE),
#'   sample_prior = function(n) cbind(mu = rnorm(n, 0, 2)),
#'   verbose = FALSE
#' )
#' fit
#' summary(fit)
#'
#' # Exact posterior mean for comparison
#' sum(y) / (length(y) + 1 / 4)
#' @export
imis <- function(log_likelihood, log_prior, sample_prior,
                 n_initial = NULL, n_per_component = NULL,
                 n_resample = 3000L, max_iter = 100L, n_optim = 0L,
                 n_centers = 1L, sampling = c("incremental", "mixture"),
                 prior_weight = 0.1, target_unique = 1 - exp(-1),
                 vectorized = TRUE, vectorized_prior = vectorized,
                 parallel = FALSE, n_chunks = NULL,
                 catch_errors = FALSE, optim_control = list(),
                 verbose = TRUE) {
  cl <- match.call()
  t_start <- proc.time()[[3]]
  sampling <- match.arg(sampling)
  check_fun(log_likelihood, "log_likelihood")
  check_fun(log_prior, "log_prior")
  check_fun(sample_prior, "sample_prior")
  check_flag(vectorized, "vectorized")
  check_flag(vectorized_prior, "vectorized_prior")
  check_flag(parallel, "parallel")
  check_flag(catch_errors, "catch_errors")
  check_count(n_resample, "n_resample")
  check_count(max_iter, "max_iter", min = 0)
  check_count(n_optim, "n_optim", min = 0)
  check_count(n_centers, "n_centers")
  check_prob(prior_weight, "prior_weight", open_right = TRUE)
  check_prob(target_unique, "target_unique", open_left = TRUE)
  if (!is.null(n_chunks)) check_count(n_chunks, "n_chunks")
  optim_control <- utils::modifyList(list(maxit_nm = 1000L, maxit_bfgs = 1000L),
                                     optim_control)

  # Probe the prior sampler to learn the dimension and parameter names.
  m_probe <- as_draws(sample_prior(2L), 2L)
  n_d <- ncol(m_probe)
  v_names <- colnames(m_probe)
  n_initial <- n_initial %||% (1000L * n_d)
  n_per_component <- n_per_component %||% (100L * n_d)
  check_count(n_initial, "n_initial", min = 2)
  check_count(n_per_component, "n_per_component", min = 2)

  eval_prior <- function(m_x) {
    eval_log_fun(m_x, log_prior, vectorized_prior, what = "log_prior")
  }
  n_lik_evals <- 0L
  eval_lik <- function(m_x, v_lp) {
    v_ll <- rep(-Inf, nrow(m_x))
    v_ok <- is.finite(v_lp)
    if (any(v_ok)) {
      v_ll[v_ok] <- eval_log_fun(m_x[v_ok, , drop = FALSE], log_likelihood,
                                 vectorized, parallel, n_chunks, catch_errors,
                                 what = "log_likelihood")
      n_lik_evals <<- n_lik_evals + sum(v_ok)
    }
    v_ll
  }

  # ---- Initial stage -------------------------------------------------------
  m_x <- as_draws(sample_prior(n_initial), n_initial, v_names)
  v_lp <- eval_prior(m_x)
  if (any(!is.finite(v_lp))) {
    warning(sum(!is.finite(v_lp)), " draws from `sample_prior()` have zero ",
            "density under `log_prior()`; check that the two functions ",
            "describe the same prior.", call. = FALSE)
  }
  v_prior_var <- apply(m_x, 2L, stats::var)
  if (any(!is.finite(v_prior_var) | v_prior_var <= 0)) {
    stop("Every parameter must vary in the prior sample.", call. = FALSE)
  }
  if (verbose) {
    message(sprintf("IMIS: %d parameter(s), %d initial draws from the prior",
                    n_d, n_initial))
  }
  t_lik <- proc.time()[[3]]
  v_ll <- eval_lik(m_x, v_lp)
  if (verbose) {
    message(sprintf("  %d likelihoods evaluated in %.2f s", n_lik_evals,
                    proc.time()[[3]] - t_lik))
  }
  # Log of sum_c N_c * density_c(x): the mixture density times N. The prior is
  # the first component and the normals are added as IMIS proceeds.
  v_log_env <- log(n_initial) + v_lp
  n_prior <- n_initial
  l_comp <- list()
  l_diag <- list()
  converged <- FALSE
  l_stats <- NULL
  log_post1 <- function(theta) {
    m_theta <- matrix(theta, nrow = 1L, dimnames = list(NULL, v_names))
    n_lp <- eval_prior(m_theta)
    if (!is.finite(n_lp)) return(-Inf)
    n_lp + eval_lik(m_theta, n_lp)
  }

  for (k in 0:max_iter) {
    n <- nrow(m_x)
    v_lw <- v_lp + v_ll - (v_log_env - log(n))
    l_stats <- weight_stats(v_lw, n_resample)
    if (is.null(l_stats)) {
      stop("All importance weights are zero: the likelihood is -Inf (or the ",
           "prior density is 0) at every draw. Check `log_likelihood()`, or ",
           "widen the prior / increase `n_initial`.", call. = FALSE)
    }
    l_diag[[k + 1L]] <- data.frame(
      iter = k, n_draws = n, n_components = length(l_comp),
      log_ml = l_stats$log_ml, log_ml_se = l_stats$log_ml_se,
      ess = l_stats$ess, max_weight = l_stats$max_weight,
      expected_unique = l_stats$expected_unique,
      entropy = l_stats$entropy, var_weights = l_stats$var_weights,
      n_lik_evals = n_lik_evals, elapsed = proc.time()[[3]] - t_start
    )
    if (verbose) print_progress(l_diag[[k + 1L]], k == 0L)
    if (l_stats$expected_unique >= target_unique * n_resample) {
      converged <- TRUE
      break
    }
    if (k == max_iter) {
      break
    }

    # ---- Add components ----------------------------------------------------
    if (k == 0L && n_optim > 0L) {
      if (verbose) message("Optimization stage")
      l_new <- optim_components(m_x, l_stats$w, n_optim, log_post1,
                                v_prior_var, optim_control, verbose)
    } else {
      l_new <- weight_components(m_x, l_stats$w, n_centers, n_per_component,
                                 v_prior_var, iter = k + 1L)
    }
    l_comp <- c(l_comp, l_new)
    n_comp <- length(l_comp)

    # ---- Decide how many draws each component gets ---------------------------
    if (sampling == "incremental") {
      v_alloc <- c(0L, rep(0L, n_comp - length(l_new)),
                   rep(as.integer(n_per_component), length(l_new)))
    } else {
      v_alloc <- allocate_draws(
        n_per_component * length(l_new),
        c(prior_weight, rep((1 - prior_weight) / n_comp, n_comp))
      )
    }

    # ---- Draw ----------------------------------------------------------------
    l_new_x <- list()
    if (v_alloc[1L] > 0L) {
      l_new_x[[1L]] <- as_draws(sample_prior(v_alloc[1L]), v_alloc[1L],
                                v_names)
    }
    for (j in which(v_alloc[-1L] > 0L)) {
      l_new_x[[length(l_new_x) + 1L]] <- rmvnorm_chol(
        v_alloc[j + 1L], l_comp[[j]]$mean, l_comp[[j]]$chol)
    }
    m_new <- do.call(rbind, l_new_x)
    colnames(m_new) <- v_names

    # ---- Update the mixture density ----------------------------------------
    # Old draws: add the terms of the components that just received draws.
    m_add <- NULL
    if (v_alloc[1L] > 0L) {
      m_add <- cbind(m_add, log(v_alloc[1L]) + v_lp)
    }
    for (j in which(v_alloc[-1L] > 0L)) {
      m_add <- cbind(m_add, log(v_alloc[j + 1L]) +
                       dmvnorm_log(m_x, l_comp[[j]]$mean, l_comp[[j]]$chol))
    }
    v_log_env <- log_add(v_log_env, row_logsumexp(m_add))

    n_prior <- n_prior + v_alloc[1L]
    for (j in seq_len(n_comp)) {
      l_comp[[j]]$n <- l_comp[[j]]$n + v_alloc[j + 1L]
    }

    # New draws: every component with at least one draw so far.
    v_lp_new <- eval_prior(m_new)
    m_env_new <- matrix(log(n_prior) + v_lp_new, ncol = 1L)
    for (j in seq_len(n_comp)) {
      if (l_comp[[j]]$n > 0L) {
        m_env_new <- cbind(m_env_new, log(l_comp[[j]]$n) +
                             dmvnorm_log(m_new, l_comp[[j]]$mean,
                                         l_comp[[j]]$chol))
      }
    }
    v_ll_new <- eval_lik(m_new, v_lp_new)

    m_x <- rbind(m_x, m_new)
    v_lp <- c(v_lp, v_lp_new)
    v_ll <- c(v_ll, v_ll_new)
    v_log_env <- c(v_log_env, row_logsumexp(m_env_new))
  }

  if (!converged) {
    warning(sprintf(paste0(
      "IMIS did not meet the stopping rule within max_iter = %d iterations ",
      "(expected unique points: %.0f of the %.0f targeted). The posterior ",
      "sample may be poor; consider increasing `max_iter` or ",
      "`n_per_component`, or setting `n_optim`."),
      max_iter, l_stats$expected_unique, target_unique * n_resample),
      call. = FALSE)
  }

  # ---- Resample ------------------------------------------------------------
  v_w <- l_stats$w
  v_idx <- sample.int(length(v_w), n_resample, replace = TRUE, prob = v_w)
  m_resample <- m_x[v_idx, , drop = FALSE]
  rownames(m_resample) <- NULL

  m_centers <- if (length(l_comp) > 0L) {
    do.call(rbind, lapply(l_comp, `[[`, "mean"))
  } else {
    matrix(numeric(0), 0L, n_d)
  }
  colnames(m_centers) <- v_names
  df_diag <- do.call(rbind, l_diag)
  rownames(df_diag) <- NULL

  structure(list(
    resample = m_resample,
    draws = m_x,
    weights = v_w,
    log_weights = v_lp + v_ll - (v_log_env - log(nrow(m_x))),
    log_prior = v_lp,
    log_lik = v_ll,
    log_ml = l_stats$log_ml,
    log_ml_se = l_stats$log_ml_se,
    diagnostics = df_diag,
    components = lapply(l_comp, function(x) x[c("mean", "cov", "n", "origin",
                                                 "iter")]),
    centers = m_centers,
    n_prior_draws = n_prior,
    n_initial = as.integer(n_initial),
    converged = converged,
    n_iter = max(df_diag$iter),
    n_lik_evals = n_lik_evals,
    settings = list(n_initial = n_initial, n_per_component = n_per_component,
                    n_resample = n_resample, max_iter = max_iter,
                    n_optim = n_optim, n_centers = n_centers,
                    sampling = sampling, prior_weight = prior_weight,
                    target_unique = target_unique, parallel = parallel),
    call = cl,
    elapsed = proc.time()[[3]] - t_start
  ), class = "imis")
}

print_progress <- function(df_row, header) {
  if (header) {
    message(sprintf("%6s %9s %5s %12s %9s %10s %10s %9s",
                    "stage", "draws", "comps", "log_ml", "se", "ESS",
                    "unique", "max_w"))
  }
  message(sprintf("%6d %9d %5d %12.3f %9.4f %10.1f %10.1f %9.4f",
                  df_row$iter, df_row$n_draws, df_row$n_components,
                  df_row$log_ml, df_row$log_ml_se, df_row$ess,
                  df_row$expected_unique, df_row$max_weight))
}

check_fun <- function(x, name) {
  if (!is.function(x)) {
    stop("`", name, "` must be a function.", call. = FALSE)
  }
}

check_flag <- function(x, name) {
  if (!is.logical(x) || length(x) != 1L || is.na(x)) {
    stop("`", name, "` must be TRUE or FALSE.", call. = FALSE)
  }
}

check_count <- function(x, name, min = 1) {
  if (!is.numeric(x) || length(x) != 1L || is.na(x) || x < min ||
      x != round(x)) {
    stop("`", name, "` must be a single whole number >= ", min, ".",
         call. = FALSE)
  }
}

check_prob <- function(x, name, open_left = FALSE, open_right = FALSE) {
  ok <- is.numeric(x) && length(x) == 1L && !is.na(x) &&
    (if (open_left) x > 0 else x >= 0) && (if (open_right) x < 1 else x <= 1)
  if (!ok) {
    stop("`", name, "` must be a single number in ",
         if (open_left) "(" else "[", "0, 1", if (open_right) ")" else "]",
         ".", call. = FALSE)
  }
}
