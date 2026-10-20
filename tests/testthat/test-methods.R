fit_conj <- local({
  set.seed(201)
  imis(conj_ll, conj_lp, conj_sp, verbose = FALSE)
})

test_that("print and summary report the fit", {
  s <- summary(fit_conj)
  expect_s3_class(s, "summary.imis")
  expect_named(s$estimates, c("parameter", "mean", "sd", "2.5%", "50%",
                              "97.5%"))
  expect_equal(s$estimates$`50%`, conj_post_mean, tolerance = 0.05)
  expect_output(print(fit_conj), "log marginal likelihood")
  expect_output(print(fit_conj), "Converged")
  s1 <- summary(fit_conj, probs = 0.5)
  expect_named(s1$estimates, c("parameter", "mean", "sd", "50%"))
})

test_that("plot methods run", {
  pdf(NULL)
  on.exit(dev.off())
  expect_invisible(plot(fit_conj))
  expect_invisible(plot(fit_conj, type = "marginals"))
  expect_error(plot(fit_conj, type = "pairs"), "two parameters")
  set.seed(202)
  fit2 <- imis(bim_ll, bim_lp, bim_sp, n_optim = 3, verbose = FALSE)
  expect_invisible(plot(fit2, type = "pairs"))
})

test_that("draws can be extracted", {
  expect_identical(as.matrix(fit_conj), fit_conj$resample)
  w <- as.matrix(fit_conj, which = "weighted")
  expect_equal(colnames(w), c("mu", "weight"))
  expect_equal(sum(w[, "weight"]), 1)
  expect_s3_class(as.data.frame(fit_conj), "data.frame")
})

test_that("bayes_factor combines two fits", {
  set.seed(203)
  other <- imis(conj_ll, function(th) dnorm(th[, 1], 3, 1, log = TRUE),
                function(n) cbind(mu = rnorm(n, 3)), verbose = FALSE)
  bf <- bayes_factor(fit_conj, other)
  expect_named(bf, c("log_bf", "log_bf_se", "bf"))
  expect_equal(bf[["log_bf"]], fit_conj$log_ml - other$log_ml)
  expect_equal(bf[["bf"]], exp(bf[["log_bf"]]))
  expect_error(bayes_factor(fit_conj, list()), "imis")
})

test_that("imis_legacy() reproduces the IMIS::IMIS() interface", {
  sample.prior <- function(n) cbind(mu = rnorm(n, 0, conj_s0))
  prior <- function(x) dnorm(x[, 1], 0, conj_s0)
  likelihood <- function(x) exp(conj_ll(x))
  set.seed(204)
  res <- suppressMessages(imis_legacy(B = 100, B.re = 1000, number_k = 50))
  expect_named(res, c("stat", "resample", "center", "fit"))
  expect_equal(colnames(res$stat), c("MargLike", "UniquePoint", "MaxWeight",
                                     "ESS", "Entropy", "VarWeights"))
  expect_equal(mean(res$resample), conj_post_mean, tolerance = 0.05)
  expect_equal(res$fit$n_initial, 1000)
})
