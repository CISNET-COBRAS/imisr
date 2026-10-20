test_that("log-scale helpers are exact and handle -Inf", {
  x <- c(-1000, -1001, -1002)
  expect_equal(imisr:::logsumexp(x), -1000 + log(1 + exp(-1) + exp(-2)))
  expect_equal(imisr:::logsumexp(c(-Inf, -Inf)), -Inf)
  m <- rbind(c(0, 0), c(-Inf, -Inf), c(-800, -Inf))
  expect_equal(imisr:::row_logsumexp(m), c(log(2), -Inf, -800))
  expect_equal(imisr:::log_add(c(-Inf, 0, -Inf), c(-Inf, 0, 3)),
               c(-Inf, log(2), 3))
})

test_that("multivariate normal helpers match mvtnorm", {
  skip_if_not_installed("mvtnorm")
  set.seed(1)
  s <- matrix(c(2, 0.6, 0.3, 0.6, 1, -0.2, 0.3, -0.2, 0.5), 3)
  mu <- c(1, -1, 0.5)
  x <- matrix(rnorm(30), 10, 3)
  expect_equal(imisr:::dmvnorm_log(x, mu, chol(s)),
               mvtnorm::dmvnorm(x, mu, s, log = TRUE))
  z <- imisr:::rmvnorm_chol(2e5, mu, chol(s))
  expect_equal(colMeans(z), mu, tolerance = 0.02)
  expect_equal(cov(z), s, tolerance = 0.02)
})

test_that("make_pd repairs non-positive-definite matrices", {
  s <- matrix(c(1, 1, 1, 1), 2) # singular
  out <- imisr:::make_pd(s)
  expect_true(all(eigen(out$cov)$values > 0))
  expect_equal(crossprod(out$chol), out$cov)
  fb <- diag(2) * 3
  expect_equal(imisr:::make_pd(matrix(0, 2, 2), fallback = fb)$cov, fb)
})

test_that("weighted covariance matches stats::cov.wt", {
  set.seed(2)
  x <- matrix(rnorm(40), 20, 2)
  w <- runif(20)
  ctr <- c(0.1, -0.2)
  expect_equal(imisr:::weighted_cov(x, w, ctr),
               cov.wt(x, w, center = ctr, method = "unbiased")$cov,
               ignore_attr = TRUE)
})

test_that("allocate_draws hands out exactly n draws", {
  expect_equal(imisr:::allocate_draws(10, c(1, 1, 1)), c(4L, 3L, 3L))
  a <- imisr:::allocate_draws(1001, c(0.1, 0.45, 0.45))
  expect_equal(sum(a), 1001)
  expect_equal(imisr:::allocate_draws(5, c(0, 1)), c(0L, 5L))
})

test_that("as_draws coerces and validates prior samples", {
  x <- imisr:::as_draws(1:3, 3)
  expect_equal(dim(x), c(3, 1))
  expect_equal(colnames(x), "theta1")
  expect_equal(colnames(imisr:::as_draws(data.frame(a = 1:2, b = 3:4), 2)),
               c("a", "b"))
  one <- imisr:::as_draws(c(a = 1, b = 2), 1)
  expect_equal(dim(one), c(1, 2))
  expect_equal(colnames(one), c("a", "b"))
  expect_error(imisr:::as_draws(matrix(1:4, 2), 3), "rows")
  expect_error(imisr:::as_draws(matrix(1:4, 2), 2, c("a", "b", "c")),
               "columns")
})
