test_that("parallel evaluation with future gives the serial answer", {
  skip_if_not_installed("future")
  skip_if_not_installed("future.apply")
  old_plan <- future::plan(future::sequential)
  on.exit(future::plan(old_plan), add = TRUE)

  m <- cbind(mu = seq(-2, 2, length.out = 17))
  expect_equal(
    imisr:::eval_log_fun(m, conj_ll, TRUE, parallel = TRUE, n_chunks = 4),
    conj_ll(m)
  )
  set.seed(301)
  fit <- imis(conj_ll, conj_lp, conj_sp, parallel = TRUE, n_chunks = 3,
              verbose = FALSE)
  expect_true(fit$converged)
  expect_equal(summary(fit)$estimates$mean, conj_post_mean, tolerance = 0.03)
})

test_that("parallel evaluation works on separate R sessions", {
  skip_on_cran()
  skip_if_not_installed("future")
  skip_if_not_installed("future.apply")
  old_plan <- future::plan(future::multisession, workers = 2)
  on.exit(future::plan(old_plan), add = TRUE)

  y <- conj_y
  ll <- function(th) {
    colSums(dnorm(outer(y, th[, 1], "-"), log = TRUE))
  }
  set.seed(302)
  fit <- imis(ll, conj_lp, conj_sp, parallel = TRUE, verbose = FALSE)
  expect_equal(summary(fit)$estimates$mean, conj_post_mean, tolerance = 0.03)

  # Errors raised on the workers are caught row by row, too.
  set.seed(303)
  expect_warning(no_conv_warning(
    imis(function(th) sum(dnorm(y, th[["mu"]], 1, log = TRUE)) +
           if (th[["mu"]] > 3) stop("boom") else 0,
         function(th) dnorm(th[["mu"]], 0, 2, log = TRUE), conj_sp,
         vectorized = FALSE, parallel = TRUE, catch_errors = TRUE,
         max_iter = 1, target_unique = 1, verbose = FALSE)),
    "boom")
})
