test_that("supplied initial draws and log-likelihoods are reused, not recomputed", {
  m_init <- withr_seed(401, conj_sp(4000))
  v_llik <- conj_ll(m_init)
  n_calls <- 0L
  ll <- function(th) {
    n_calls <<- n_calls + nrow(th)
    conj_ll(th)
  }
  set.seed(402)
  fit <- imis(ll, conj_lp, conj_sp, initial_draws = m_init,
              initial_log_lik = v_llik, verbose = FALSE)
  expect_true(fit$converged)
  expect_equal(fit$n_initial, 4000L)
  expect_equal(fit$draws[1:4000, , drop = FALSE], m_init)
  expect_equal(fit$log_lik[1:4000], v_llik)
  # The likelihood was only evaluated on the draws added after the initial
  # stage
  expect_equal(n_calls, fit$n_lik_evals)
  expect_equal(fit$n_lik_evals, nrow(fit$draws) - 4000L)
  s <- summary(fit)$estimates
  expect_equal(s$mean, conj_post_mean, tolerance = 0.03)
  expect_lt(abs(fit$log_ml - conj_log_ml), 4 * fit$log_ml_se + 0.01)
})

test_that("initial draws without log-likelihoods are evaluated as usual", {
  m_init <- withr_seed(403, conj_sp(3000))
  fit_eval <- withr_seed(404, imis(conj_ll, conj_lp, conj_sp,
                                   initial_draws = m_init, verbose = FALSE))
  fit_given <- withr_seed(404, imis(conj_ll, conj_lp, conj_sp,
                                    initial_draws = m_init,
                                    initial_log_lik = conj_ll(m_init),
                                    verbose = FALSE))
  expect_equal(fit_eval$resample, fit_given$resample)
  expect_equal(fit_eval$n_lik_evals - fit_given$n_lik_evals, 3000L)
})

test_that("initial draws are matched to the parameters by name", {
  m_init <- withr_seed(405, bim_sp(4000))
  v_llik <- bim_ll(m_init)
  fit_with <- function(m) {
    withr_seed(406, no_conv_warning(
      imis(bim_ll, bim_lp, bim_sp, initial_draws = m, initial_log_lik = v_llik,
           max_iter = 2, verbose = FALSE)
    ))
  }
  fit_ref <- fit_with(m_init)
  expect_equal(fit_with(m_init[, c("b", "a")])$resample, fit_ref$resample)
  expect_equal(fit_with(as.data.frame(m_init))$resample, fit_ref$resample)
  expect_equal(fit_with(unname(m_init))$resample, fit_ref$resample)
})

test_that("initial draws outside the prior support are flagged and skipped", {
  ll <- function(th) dnorm(0.95, th[, 1], 0.05, log = TRUE)
  lp <- function(th) dunif(th[, 1], 0, 1, log = TRUE)
  sp <- function(n) cbind(p = runif(n))
  m_init <- cbind(p = c(-0.5, withr_seed(407, runif(999))))
  set.seed(408)
  expect_warning(
    fit <- imis(ll, lp, sp, initial_draws = m_init,
                initial_log_lik = ll(m_init), verbose = FALSE),
    "1 rows of `initial_draws` have zero density"
  )
  expect_equal(fit$log_lik[1], -Inf)
  expect_true(all(fit$resample >= 0 & fit$resample <= 1))
})

test_that("invalid initial draws and log-likelihoods are rejected", {
  m_init <- withr_seed(409, conj_sp(500))
  run <- function(...) {
    imis(conj_ll, conj_lp, conj_sp, verbose = FALSE, ...)
  }
  expect_error(run(initial_log_lik = conj_ll(m_init)), "needs the draws")
  expect_error(run(initial_draws = m_init, initial_log_lik = 1:3),
               "one value per row")
  expect_error(run(initial_draws = m_init, n_initial = 100),
               "differs from the number of rows")
  expect_error(run(initial_draws = cbind(sigma = m_init[, 1])),
               "do not match the parameters")
  expect_error(run(initial_draws = cbind(m_init, m_init)),
               "one column per parameter")
  expect_error(run(initial_draws = rbind(m_init, NA)), "missing values")
  expect_error(run(initial_draws = m_init,
                   initial_log_lik = c(Inf, conj_ll(m_init)[-1])), "\\+Inf")
  set.seed(410)
  expect_warning(run(initial_draws = m_init,
                     initial_log_lik = c(NA, conj_ll(m_init)[-1])),
                 "NA/NaN for 1 of 500")
})

test_that("progress reports where the initial draws came from", {
  m_init <- withr_seed(411, conj_sp(1000))
  set.seed(412)
  expect_message(
    no_conv_warning(imis(conj_ll, conj_lp, conj_sp, initial_draws = m_init,
                         initial_log_lik = conj_ll(m_init), max_iter = 0)),
    "supplied in `initial_draws`"
  )
})
