#' Summarize an IMIS fit
#'
#' Posterior summaries computed from *all* IMIS draws and their importance
#' weights (more precise than summarizing the resample), together with
#' sampling diagnostics.
#'
#' @param object,x An object of class `imis` (or `summary.imis`).
#' @param probs Quantiles to report.
#' @param digits Number of significant digits to print.
#' @param ... Unused.
#' @return `summary()` returns an object of class `summary.imis`: a list with
#'   `estimates` (data frame of posterior mean, sd and quantiles per
#'   parameter), `log_ml`, `log_ml_se`, `ess`, `expected_unique`,
#'   `converged`, `n_iter`, `n_draws` and `n_lik_evals`.
#' @examples
#' fit <- imis(
#'   log_likelihood = function(theta) dnorm(1, theta[, 1], 0.5, log = TRUE),
#'   log_prior = function(theta) dnorm(theta[, 1], 0, 1, log = TRUE),
#'   sample_prior = function(n) cbind(mu = rnorm(n)),
#'   verbose = FALSE
#' )
#' summary(fit)
#' @export
summary.imis <- function(object, probs = c(0.025, 0.5, 0.975), ...) {
  v_w <- object$weights
  m_x <- object$draws
  df_est <- data.frame(
    parameter = colnames(m_x),
    mean = colSums(m_x * v_w),
    sd = apply(m_x, 2L, function(v) {
      n_mu <- sum(v * v_w)
      sqrt(sum(v_w * (v - n_mu)^2))
    }),
    row.names = NULL
  )
  m_q <- t(apply(m_x, 2L, weighted_quantile, w = v_w, probs = probs))
  if (length(probs) == 1L) m_q <- t(m_q)
  colnames(m_q) <- paste0(100 * probs, "%")
  df_est <- cbind(df_est, m_q)
  df_last <- object$diagnostics[nrow(object$diagnostics), ]
  structure(list(
    estimates = df_est,
    log_ml = object$log_ml,
    log_ml_se = object$log_ml_se,
    ess = df_last$ess,
    expected_unique = df_last$expected_unique,
    n_resample = object$settings$n_resample,
    converged = object$converged,
    n_iter = object$n_iter,
    n_draws = nrow(m_x),
    n_components = length(object$components),
    n_lik_evals = object$n_lik_evals,
    elapsed = object$elapsed
  ), class = "summary.imis")
}

#' @rdname summary.imis
#' @export
print.summary.imis <- function(x, digits = 4, ...) {
  cat("Incremental mixture importance sampling\n")
  cat(sprintf("  %s after %d iteration(s); %d draws, %d normal component(s)\n",
              if (x$converged) "Converged" else "NOT converged", x$n_iter,
              x$n_draws, x$n_components))
  cat(sprintf("  %d likelihood evaluations in %.1f s\n", x$n_lik_evals,
              x$elapsed))
  cat(sprintf("  ESS = %.1f; expected unique points in resample = %.0f / %d\n",
              x$ess, x$expected_unique, x$n_resample))
  cat(sprintf("  log marginal likelihood = %.4f (SE %.4f)\n\n", x$log_ml,
              x$log_ml_se))
  cat("Posterior summary (importance-weighted):\n")
  df <- x$estimates
  v_num <- vapply(df, is.numeric, logical(1))
  df[v_num] <- lapply(df[v_num], signif, digits = digits)
  print(df, row.names = FALSE)
  invisible(x)
}

#' @rdname summary.imis
#' @export
print.imis <- function(x, digits = 4, ...) {
  print(summary(x), digits = digits)
  invisible(x)
}

#' Plot an IMIS fit
#'
#' @param x An object of class `imis`.
#' @param type `"diagnostics"` plots the log marginal likelihood (with
#'   \eqn{\pm 2} standard errors), the effective sample size, the expected
#'   number of unique resampled points (with the stopping threshold) and the
#'   maximum weight against the IMIS stage. `"marginals"` compares the
#'   marginal prior (from the initial prior draws) and posterior densities.
#'   `"pairs"` shows a pairs plot of the posterior resample with the mixture
#'   component centers in red.
#' @param ... Further arguments passed to the underlying plotting function.
#' @return `x`, invisibly.
#' @examples
#' fit <- imis(
#'   log_likelihood = function(theta) dnorm(1, theta[, 1], 0.5, log = TRUE),
#'   log_prior = function(theta) dnorm(theta[, 1], 0, 1, log = TRUE),
#'   sample_prior = function(n) cbind(mu = rnorm(n)),
#'   verbose = FALSE
#' )
#' plot(fit)
#' plot(fit, type = "marginals")
#' @export
plot.imis <- function(x, type = c("diagnostics", "marginals", "pairs"), ...) {
  type <- match.arg(type)
  old_par <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old_par))
  if (type == "diagnostics") {
    df <- x$diagnostics
    graphics::par(mfrow = c(2, 2), mar = c(4, 4, 2, 1))
    v_lo <- df$log_ml - 2 * df$log_ml_se
    v_hi <- df$log_ml + 2 * df$log_ml_se
    graphics::plot(df$iter, df$log_ml, type = "b", pch = 20,
                   ylim = range(c(v_lo, v_hi), finite = TRUE),
                   xlab = "Stage", ylab = "log marginal likelihood",
                   main = "Marginal likelihood", ...)
    graphics::lines(df$iter, v_lo, lty = 2)
    graphics::lines(df$iter, v_hi, lty = 2)
    graphics::plot(df$iter, df$ess, type = "b", pch = 20, xlab = "Stage",
                   ylab = "ESS", main = "Effective sample size", ...)
    n_target <- x$settings$target_unique * x$settings$n_resample
    graphics::plot(df$iter, df$expected_unique, type = "b", pch = 20,
                   ylim = c(0, max(c(df$expected_unique, n_target))),
                   xlab = "Stage", ylab = "Expected unique points",
                   main = "Stopping rule", ...)
    graphics::abline(h = n_target, col = "red", lty = 2)
    graphics::plot(df$iter, df$max_weight, type = "b", pch = 20, log = "y",
                   xlab = "Stage", ylab = "Maximum weight",
                   main = "Largest importance weight", ...)
  } else if (type == "marginals") {
    n_d <- ncol(x$draws)
    n_col <- ceiling(sqrt(n_d))
    graphics::par(mfrow = c(ceiling(n_d / n_col), n_col),
                  mar = c(4, 4, 2, 1))
    v_prior_idx <- seq_len(x$n_initial)
    for (j in seq_len(n_d)) {
      v_post <- stats::density(x$draws[, j], weights = x$weights,
                               warnWbw = FALSE)
      v_prior <- stats::density(x$draws[v_prior_idx, j])
      graphics::plot(v_post, main = colnames(x$draws)[j], xlab = "",
                     ylim = c(0, max(v_post$y, v_prior$y)),
                     col = "black", lwd = 2, ...)
      graphics::lines(v_prior, col = "grey60", lty = 2, lwd = 2)
      if (j == 1L) {
        graphics::legend("topright", c("posterior", "prior"),
                         col = c("black", "grey60"), lty = c(1, 2), lwd = 2,
                         bty = "n")
      }
    }
  } else {
    if (ncol(x$resample) < 2L) {
      stop("A pairs plot needs at least two parameters.", call. = FALSE)
    }
    m_all <- rbind(x$resample, x$centers)
    v_col <- c(rep(grDevices::adjustcolor("black", 0.3), nrow(x$resample)),
               rep("red", nrow(x$centers)))
    v_pch <- c(rep(20, nrow(x$resample)), rep(8, nrow(x$centers)))
    graphics::pairs(m_all, col = v_col, pch = v_pch, ...)
  }
  invisible(x)
}

#' Extract posterior draws from an IMIS fit
#'
#' @param x An object of class `imis`.
#' @param which `"resample"` (default) returns the equally weighted posterior
#'   resample. `"weighted"` returns all IMIS draws with a `weight` column
#'   holding their normalized importance weights.
#' @param ... Unused.
#' @return A matrix (`as.matrix()`) or data frame (`as.data.frame()`).
#' @examples
#' fit <- imis(
#'   log_likelihood = function(theta) dnorm(1, theta[, 1], 0.5, log = TRUE),
#'   log_prior = function(theta) dnorm(theta[, 1], 0, 1, log = TRUE),
#'   sample_prior = function(n) cbind(mu = rnorm(n)),
#'   verbose = FALSE
#' )
#' head(as.matrix(fit))
#' head(as.data.frame(fit, which = "weighted"))
#' @export
as.matrix.imis <- function(x, which = c("resample", "weighted"), ...) {
  which <- match.arg(which)
  if (which == "resample") {
    return(x$resample)
  }
  cbind(x$draws, weight = x$weights)
}

#' @rdname as.matrix.imis
#' @param row.names,optional Unused; for compatibility with the generic.
#' @export
as.data.frame.imis <- function(x, row.names = NULL, optional = FALSE,
                               which = c("resample", "weighted"), ...) {
  as.data.frame(as.matrix(x, which = which))
}

#' Bayes factor between two models fitted with IMIS
#'
#' Computes the log Bayes factor of model 1 against model 2 from their IMIS
#' marginal likelihood estimates, with a standard error obtained by adding
#' the two (independent) Monte Carlo variances, as in Raftery and Bao (2010,
#' Section 5.2).
#'
#' @param fit1,fit2 Objects of class `imis` for the two models, fitted to the
#'   same data.
#' @return A named numeric vector: `log_bf`, `log_bf_se` and `bf`.
#' @examples
#' y <- c(0.4, 1.1, 0.9)
#' m_free <- imis(
#'   function(t) colSums(dnorm(outer(y, t[, "mu"], "-"), log = TRUE)),
#'   function(t) dnorm(t[, "mu"], 0, 1, log = TRUE),
#'   function(n) cbind(mu = rnorm(n)),
#'   verbose = FALSE
#' )
#' m_shift <- imis(
#'   function(t) colSums(dnorm(outer(y, t[, "mu"], "-"), log = TRUE)),
#'   function(t) dnorm(t[, "mu"], 3, 1, log = TRUE),
#'   function(n) cbind(mu = rnorm(n, 3)),
#'   verbose = FALSE
#' )
#' bayes_factor(m_free, m_shift)
#' @export
bayes_factor <- function(fit1, fit2) {
  if (!inherits(fit1, "imis") || !inherits(fit2, "imis")) {
    stop("`fit1` and `fit2` must both be `imis` objects.", call. = FALSE)
  }
  n_log_bf <- fit1$log_ml - fit2$log_ml
  c(log_bf = n_log_bf,
    log_bf_se = sqrt(fit1$log_ml_se^2 + fit2$log_ml_se^2),
    bf = exp(n_log_bf))
}

weighted_quantile <- function(x, w, probs) {
  v_ord <- order(x)
  x <- x[v_ord]
  v_cw <- cumsum(w[v_ord])
  v_cw <- v_cw / v_cw[length(v_cw)]
  vapply(probs, function(p) x[which(v_cw >= p - 1e-12)[1L]], numeric(1))
}
