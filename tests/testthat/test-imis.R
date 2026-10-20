test_that("IMIS recovers a conjugate normal posterior and marginal likelihood", {
  set.seed(101)
  fit <- imis(conj_ll, conj_lp, conj_sp, verbose = FALSE)
  expect_s3_class(fit, "imis")
  expect_true(fit$converged)
  s <- summary(fit)$estimates
  expect_equal(s$mean, conj_post_mean, tolerance = 0.03)
  expect_equal(s$sd, conj_post_sd, tolerance = 0.05)
  expect_lt(abs(fit$log_ml - conj_log_ml), 4 * fit$log_ml_se + 0.01)
  expect_equal(mean(fit$resample[, "mu"]), conj_post_mean, tolerance = 0.05)
  expect_equal(colnames(fit$resample), "mu")
  expect_equal(nrow(fit$resample), 3000)
  expect_equal(sum(fit$weights), 1)
})

test_that("both ways of drawing new points are unbiased", {
  for (s in c("incremental", "mixture")) {
    set.seed(102)
    fit <- imis(conj_ll, conj_lp, conj_sp, sampling = s, verbose = FALSE)
    expect_lt(abs(fit$log_ml - conj_log_ml), 4 * fit$log_ml_se + 0.01)
    expect_equal(summary(fit)$estimates$mean, conj_post_mean,
                 tolerance = 0.03)
  }
})

test_that("mixture sampling keeps drawing from the prior", {
  set.seed(103)
  fit <- no_conv_warning(imis(conj_ll, conj_lp, conj_sp, sampling = "mixture",
                              prior_weight = 0.2, max_iter = 3,
                              target_unique = 1, n_per_component = 100,
                              verbose = FALSE))
  expect_equal(fit$n_prior_draws, 1000 + 3 * 20)
  expect_equal(sum(vapply(fit$components, `[[`, 0, "n")) +
                 fit$n_prior_draws, nrow(fit$draws))
})

test_that("a correlated 2-d posterior is recovered", {
  # y ~ N(theta, I) observed once per coordinate, theta ~ N(0, S): the
  # posterior is N(S (S + I)^-1 y, S - S (S + I)^-1 S).
  s0 <- matrix(c(1, 0.8, 0.8, 1), 2)
  y <- c(1, 2)
  ch <- chol(s0)
  fit <- withr_seed(104, imis(
    log_likelihood = function(th) {
      dnorm(y[1], th[, 1], 1, log = TRUE) + dnorm(y[2], th[, 2], 1, log = TRUE)
    },
    log_prior = function(th) imisr:::dmvnorm_log(th, c(0, 0), ch),
    sample_prior = function(n) {
      x <- imisr:::rmvnorm_chol(n, c(0, 0), ch)
      colnames(x) <- c("a", "b")
      x
    },
    verbose = FALSE
  ))
  k <- s0 %*% solve(s0 + diag(2))
  expect_equal(summary(fit)$estimates$mean, drop(k %*% y), tolerance = 0.03)
  expect_equal(cov.wt(fit$draws, fit$weights)$cov, s0 - k %*% s0,
               tolerance = 0.05, ignore_attr = TRUE)
  log_ml <- imisr:::dmvnorm_log(t(y), c(0, 0), chol(s0 + diag(2)))
  expect_lt(abs(fit$log_ml - log_ml), 4 * fit$log_ml_se + 0.01)
})

test_that("IMIS-opt and multiple centers find both modes of a bimodal target", {
  configs <- list(list(n_optim = 5), list(n_centers = 3),
                  list(n_optim = 5, sampling = "mixture"))
  for (cfg in configs) {
    set.seed(105)
    fit <- do.call(imis, c(list(bim_ll, bim_lp, bim_sp, verbose = FALSE),
                           cfg))
    expect_true(fit$converged)
    expect_equal(mean(fit$resample[, "a"] > 4.5), 0.5, tolerance = 0.1)
    expect_lt(abs(fit$log_ml - bim_log_ml), 4 * fit$log_ml_se + 0.02)
  }
  expect_true(any(vapply(fit$components, `[[`, "", "origin") == "optim"))
})

test_that("tiny likelihoods do not underflow", {
  # Shifting the log-likelihood by a constant must shift log_ml by exactly
  # that constant and leave the posterior unchanged. exp(-5000) is 0 in
  # double precision, which made the original implementation fail.
  set.seed(106)
  a <- imis(conj_ll, conj_lp, conj_sp, verbose = FALSE)
  set.seed(106)
  b <- imis(function(th) conj_ll(th) - 5000, conj_lp, conj_sp,
            verbose = FALSE)
  expect_equal(b$log_ml, a$log_ml - 5000)
  expect_equal(b$weights, a$weights)
  expect_equal(summary(b)$estimates, summary(a)$estimates)
})

test_that("non-vectorized functions give the same answer", {
  set.seed(107)
  a <- imis(conj_ll, conj_lp, conj_sp, verbose = FALSE)
  set.seed(107)
  b <- imis(function(th) sum(dnorm(conj_y, th[["mu"]], 1, log = TRUE)),
            function(th) dnorm(th[["mu"]], 0, conj_s0, log = TRUE),
            conj_sp, vectorized = FALSE, verbose = FALSE)
  expect_equal(b$log_ml, a$log_ml)
  expect_equal(b$resample, a$resample)
})

test_that("the likelihood is never evaluated outside the prior support", {
  seen <- numeric(0)
  ll <- function(th) {
    seen <<- c(seen, th[, 1])
    dnorm(0.95, th[, 1], 0.05, log = TRUE)
  }
  set.seed(108)
  fit <- imis(ll, function(th) dunif(th[, 1], 0, 1, log = TRUE),
              function(n) cbind(p = runif(n)), verbose = FALSE)
  expect_true(all(seen >= 0 & seen <= 1))
  expect_equal(length(seen), fit$n_lik_evals)
  expect_lt(fit$n_lik_evals, nrow(fit$draws)) # some draws fell outside
  expect_true(all(fit$resample >= 0 & fit$resample <= 1))
})

test_that("NA log-likelihoods and errors are handled", {
  ll_na <- function(th) ifelse(th[, 1] > 3, NA, conj_ll(th))
  set.seed(109)
  expect_warning(no_conv_warning(imis(ll_na, conj_lp, conj_sp, max_iter = 2,
                                      target_unique = 1, verbose = FALSE)),
                 "NA/NaN")

  ll_err <- function(th) {
    if (th[["mu"]] > 3) stop("invalid model")
    sum(dnorm(conj_y, th[["mu"]], 1, log = TRUE))
  }
  set.seed(110)
  lp_scalar <- function(th) dnorm(th[["mu"]], 0, conj_s0, log = TRUE)
  expect_error(imis(ll_err, lp_scalar, conj_sp, vectorized = FALSE,
                    verbose = FALSE), "invalid model")
  set.seed(110)
  expect_warning(
    fit <- no_conv_warning(imis(ll_err, lp_scalar, conj_sp,
                                vectorized = FALSE, catch_errors = TRUE,
                                max_iter = 2, target_unique = 1,
                                verbose = FALSE)),
    "failed for .* invalid model")
  expect_true(all(fit$log_lik[fit$draws[, 1] > 3] == -Inf))

  # A vectorized function that errors on a batch is retried row by row.
  ll_vec_err <- function(th) {
    if (any(th[, 1] > 3)) stop("invalid model")
    conj_ll(th)
  }
  set.seed(111)
  expect_warning(
    fit2 <- no_conv_warning(imis(ll_vec_err, conj_lp, conj_sp,
                                 catch_errors = TRUE, max_iter = 2,
                                 target_unique = 1, verbose = FALSE)),
    "failed")
  expect_true(any(is.finite(fit2$log_lik[fit2$draws[, 1] <= 3])))
})

test_that("results are reproducible with set.seed()", {
  set.seed(112)
  a <- imis(conj_ll, conj_lp, conj_sp, verbose = FALSE)
  set.seed(112)
  b <- imis(conj_ll, conj_lp, conj_sp, verbose = FALSE)
  expect_identical(a$resample, b$resample)
  expect_identical(a$log_ml, b$log_ml)
})

test_that("non-convergence and impossible likelihoods are reported", {
  set.seed(113)
  expect_warning(fit <- imis(conj_ll, conj_lp, conj_sp, max_iter = 1,
                             verbose = FALSE),
                 "did not meet the stopping rule")
  expect_false(fit$converged)
  expect_equal(fit$n_iter, 1)
  expect_error(imis(function(th) rep(-Inf, nrow(th)), conj_lp, conj_sp,
                    verbose = FALSE),
               "All importance weights are zero")
  expect_error(imis(function(th) rep(Inf, nrow(th)), conj_lp, conj_sp,
                    verbose = FALSE),
               "\\+Inf")
})

test_that("invalid arguments are rejected", {
  expect_error(imis(1, conj_lp, conj_sp), "must be a function")
  expect_error(imis(conj_ll, conj_lp, conj_sp, n_resample = 0), "n_resample")
  expect_error(imis(conj_ll, conj_lp, conj_sp, n_centers = 1.5), "n_centers")
  expect_error(imis(conj_ll, conj_lp, conj_sp, prior_weight = 1),
               "prior_weight")
  expect_error(imis(conj_ll, conj_lp, conj_sp, target_unique = 0),
               "target_unique")
  expect_error(imis(conj_ll, conj_lp, conj_sp, sampling = "foo"))
  expect_error(imis(function(th) 0, conj_lp, conj_sp, verbose = FALSE),
               "one value per row")
  expect_error(imis(conj_ll, conj_lp, function(n) cbind(mu = rep(1, n)),
                    verbose = FALSE),
               "must vary")
})

test_that("progress is reported when verbose", {
  set.seed(114)
  msgs <- capture_messages(imis(conj_ll, conj_lp, conj_sp, max_iter = 1,
                                target_unique = 0.01))
  expect_true(any(grepl("log_ml", msgs)))
  expect_true(any(grepl("initial draws", msgs)))
})
