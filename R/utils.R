# Internal numerical helpers. Everything here works on the log scale so that
# likelihoods such as exp(-2400) (common when calibrating to many targets) do
# not underflow to zero, which is what broke the original implementation.

# log(sum(exp(x))) without overflow/underflow.
logsumexp <- function(x) {
  m <- max(x)
  if (!is.finite(m)) {
    return(m)
  }
  m + log(sum(exp(x - m)))
}

# Row-wise log-sum-exp of a matrix.
row_logsumexp <- function(m_log) {
  if (ncol(m_log) == 1L) {
    return(m_log[, 1L])
  }
  v_max <- m_log[, 1L]
  for (j in 2:ncol(m_log)) {
    v_max <- pmax(v_max, m_log[, j])
  }
  v_out <- v_max + log(rowSums(exp(m_log - v_max)))
  v_out[is.infinite(v_max) & v_max < 0] <- -Inf
  v_out
}

# Element-wise log(exp(a) + exp(b)).
log_add <- function(a, b) {
  v_max <- pmax(a, b)
  v_out <- v_max + log1p(exp(-abs(a - b)))
  v_out[is.infinite(v_max) & v_max < 0] <- -Inf
  v_out
}

# Symmetrize a covariance matrix and make it positive definite by flooring its
# eigenvalues. Returns the covariance and its upper Cholesky factor. When the
# matrix is degenerate (all eigenvalues ~0) `fallback` is used instead.
make_pd <- function(m_cov, fallback = NULL, tol = 1e-10) {
  m_cov <- (m_cov + t(m_cov)) / 2
  m_chol <- tryCatch(chol(m_cov), error = function(e) NULL)
  if (!is.null(m_chol) && all(is.finite(m_chol))) {
    return(list(cov = m_cov, chol = m_chol))
  }
  l_eig <- eigen(m_cov, symmetric = TRUE)
  n_max <- max(l_eig$values)
  if ((!is.finite(n_max) || n_max <= 0) && !is.null(fallback)) {
    return(make_pd(fallback))
  }
  v_val <- pmax(l_eig$values, max(n_max * tol, .Machine$double.eps))
  m_cov <- l_eig$vectors %*% (v_val * t(l_eig$vectors))
  m_cov <- (m_cov + t(m_cov)) / 2
  list(cov = m_cov, chol = chol(m_cov))
}

# Log-density of a multivariate normal at the rows of `m_x`, given the upper
# Cholesky factor of the covariance.
dmvnorm_log <- function(m_x, v_mean, m_chol) {
  n_d <- length(v_mean)
  m_z <- backsolve(m_chol, t(m_x) - v_mean, transpose = TRUE)
  -0.5 * colSums(m_z^2) - sum(log(diag(m_chol))) - 0.5 * n_d * log(2 * pi)
}

# Draws from a multivariate normal given the upper Cholesky factor.
rmvnorm_chol <- function(n, v_mean, m_chol) {
  n_d <- length(v_mean)
  m_z <- matrix(stats::rnorm(n * n_d), nrow = n, ncol = n_d)
  m_z %*% m_chol + rep(v_mean, each = n)
}

# Squared Mahalanobis distance with a diagonal scale (vector of variances).
sq_dist_diag <- function(m_x, v_center, v_var) {
  colSums((t(m_x) - v_center)^2 / v_var)
}

# Weighted covariance around a fixed center, with the unbiased correction used
# by stats::cov.wt(method = "unbiased").
weighted_cov <- function(m_x, v_w, v_center) {
  v_w <- v_w / sum(v_w)
  m_xc <- sweep(m_x, 2L, v_center)
  m_cov <- crossprod(m_xc * sqrt(v_w))
  n_denom <- 1 - sum(v_w^2)
  if (n_denom > 0) {
    m_cov <- m_cov / n_denom
  }
  m_cov
}

# Split `n` draws across components with probabilities `v_prob`, deterministically
# (largest-remainder rounding). Using the realized counts in the mixture
# denominator keeps the estimator unbiased whatever the allocation.
allocate_draws <- function(n, v_prob) {
  v_prob <- v_prob / sum(v_prob)
  v_exact <- n * v_prob
  v_n <- floor(v_exact)
  n_left <- n - sum(v_n)
  if (n_left > 0) {
    v_idx <- order(v_exact - v_n, decreasing = TRUE)[seq_len(n_left)]
    v_n[v_idx] <- v_n[v_idx] + 1
  }
  as.integer(v_n)
}

# Coerce whatever `sample_prior()` returned into a numeric matrix with names.
as_draws <- function(x, n, v_names = NULL) {
  if (is.data.frame(x)) {
    x <- as.matrix(x)
  }
  if (is.null(dim(x))) {
    # A single draw of several parameters often comes back as a vector
    # (e.g. from sapply()); otherwise a vector is one parameter.
    x <- if (n == 1L && length(x) > 1L) {
      matrix(x, nrow = 1L, dimnames = list(NULL, names(x)))
    } else {
      matrix(x, ncol = 1L)
    }
  }
  if (!is.numeric(x) || nrow(x) != n) {
    stop("`sample_prior(n)` must return a numeric matrix (or vector, for a ",
         "single parameter) with `n` rows; got ", nrow(x), " rows for n = ", n,
         ".", call. = FALSE)
  }
  storage.mode(x) <- "double"
  if (!is.null(v_names)) {
    if (ncol(x) != length(v_names)) {
      stop("`sample_prior()` returned a different number of columns than ",
           "before.", call. = FALSE)
    }
    colnames(x) <- v_names
  } else if (is.null(colnames(x))) {
    colnames(x) <- paste0("theta", seq_len(ncol(x)))
  }
  x
}

# Importance-weight diagnostics (Raftery and Bao, 2010, Section 3.1) plus the
# marginal likelihood estimate and its delta-method standard error on the log
# scale (their equations 4-7).
weight_stats <- function(v_lw, n_resample) {
  n <- length(v_lw)
  n_lse <- logsumexp(v_lw)
  if (!is.finite(n_lse)) {
    return(NULL)
  }
  v_w <- exp(v_lw - n_lse)
  v_u <- exp(v_lw - max(v_lw))
  v_pos <- v_w > 0
  list(
    w = v_w,
    log_ml = n_lse - log(n),
    log_ml_se = sqrt(stats::var(v_u) / n) / mean(v_u),
    ess = 1 / sum(v_w^2),
    max_weight = max(v_w),
    expected_unique = sum(-expm1(n_resample * log1p(-v_w))),
    entropy = -sum(v_w[v_pos] * log(v_w[v_pos])) / log(n),
    var_weights = mean((n * v_w - 1)^2)
  )
}

`%||%` <- function(x, y) if (is.null(x)) y else x
