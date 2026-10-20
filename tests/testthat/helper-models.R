# Conjugate normal model with known posterior and marginal likelihood:
# y_i ~ N(mu, 1), mu ~ N(0, s0^2).
conj_y <- c(1.2, 0.8, 1.9, 1.4)
conj_s0 <- 2
conj_ll <- function(theta) {
  colSums(dnorm(outer(conj_y, theta[, 1], "-"), log = TRUE))
}
conj_lp <- function(theta) dnorm(theta[, 1], 0, conj_s0, log = TRUE)
conj_sp <- function(n) cbind(mu = rnorm(n, 0, conj_s0))
conj_post_prec <- length(conj_y) + 1 / conj_s0^2
conj_post_mean <- sum(conj_y) / conj_post_prec
conj_post_sd <- sqrt(1 / conj_post_prec)
conj_log_ml <- {
  n <- length(conj_y)
  m_s <- diag(n) + conj_s0^2
  -0.5 * (n * log(2 * pi) + as.numeric(determinant(m_s)$modulus) +
            drop(t(conj_y) %*% solve(m_s, conj_y)))
}

# Bimodal target of Raftery and Bao (2010, Section 4.3) in 2 dimensions:
# equal mixture of N(0, AR(-0.95)) and N(9, AR(0.95)) on a U(-3, 12)^2 prior.
bim_d <- 2
bim_chol <- function(r) chol(matrix(c(1, r, r, 1), 2))
bim_ll <- function(theta) {
  l1 <- imisr:::dmvnorm_log(theta, c(0, 0), bim_chol(-0.95))
  l2 <- imisr:::dmvnorm_log(theta, c(9, 9), bim_chol(0.95))
  log(0.5) + imisr:::log_add(l1, l2)
}
bim_lp <- function(theta) rowSums(dunif(theta, -3, 12, log = TRUE))
bim_sp <- function(n) {
  matrix(runif(2 * n, -3, 12), n, 2, dimnames = list(NULL, c("a", "b")))
}
bim_log_ml <- -2 * log(15) # likelihood mass inside the prior box is ~1

# Evaluate `expr` with a fixed seed, restoring the RNG state afterwards.
withr_seed <- function(seed, expr) {
  old <- if (exists(".Random.seed", globalenv())) get(".Random.seed", globalenv())
  on.exit(if (!is.null(old)) assign(".Random.seed", old, globalenv()))
  set.seed(seed)
  expr
}

# Run IMIS for a fixed, small number of iterations without the (expected)
# non-convergence warning.
no_conv_warning <- function(expr) {
  withCallingHandlers(expr, warning = function(w) {
    if (grepl("stopping rule", conditionMessage(w))) {
      invokeRestart("muffleWarning")
    }
  })
}
