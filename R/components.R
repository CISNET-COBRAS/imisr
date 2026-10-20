# Building the multivariate normal components of the IMIS mixture.

new_component <- function(v_mean, l_pd, origin, iter) {
  list(mean = v_mean, cov = l_pd$cov, chol = l_pd$chol, n = 0L,
       origin = origin, iter = iter)
}

# Standard IMIS step (Raftery and Bao, 2010, step 2a), generalized to add
# `n_centers` components at once. Each center is the highest-weight draw not
# already used; its covariance is the weighted covariance of the `n_nearest`
# draws closest to it (Mahalanobis distance scaled by the prior variances),
# with weights proportional to the importance weights plus 1/N. Draws used for
# one component are excluded as centers for the next, which spreads several
# components across different regions of high weight.
weight_components <- function(m_x, v_w, n_centers, n_nearest, v_prior_var,
                              iter) {
  n <- nrow(m_x)
  v_w_cov <- v_w + 1 / n
  v_avail <- v_w > 0
  m_fallback <- diag(v_prior_var / 100, nrow = length(v_prior_var))
  l_out <- list()
  for (j in seq_len(n_centers)) {
    v_cand <- which(v_avail)
    if (length(v_cand) == 0L) {
      break
    }
    i_star <- v_cand[which.max(v_w[v_cand])]
    v_center <- m_x[i_star, ]
    v_near <- order(sq_dist_diag(m_x, v_center, v_prior_var))[
      seq_len(min(n_nearest, n))]
    m_cov <- weighted_cov(m_x[v_near, , drop = FALSE], v_w_cov[v_near],
                          v_center)
    l_out[[j]] <- new_component(v_center, make_pd(m_cov, m_fallback),
                                origin = "weights", iter = iter)
    v_avail[c(i_star, v_near)] <- FALSE
  }
  l_out
}

# Optimization stage of IMIS-opt (Raftery and Bao, 2010, Section 3.2). Finds up
# to `n_optim` local posterior modes, starting each search from the
# highest-weight draw not yet "explained" by a previous mode, and returns a
# normal component per mode with covariance equal to the inverse Hessian of the
# negative log-posterior.
optim_components <- function(m_x, v_w, n_optim, log_post1, v_prior_var,
                             control, verbose) {
  n <- nrow(m_x)
  n_d <- ncol(m_x)
  v_names <- colnames(m_x)
  v_scale <- sqrt(v_prior_var)
  m_prior_prec <- diag(1 / v_prior_var, nrow = n_d)
  m_fallback <- diag(v_prior_var / 100, nrow = n_d)
  n_exclude <- floor(n / n_optim)

  neg_log_post <- function(theta) {
    names(theta) <- v_names
    v <- -log_post1(theta)
    if (is.na(v)) Inf else v
  }

  v_remain <- which(v_w > 0)
  l_out <- list()
  for (i in seq_len(n_optim)) {
    if (length(v_remain) == 0L) {
      break
    }
    i_start <- v_remain[which.max(v_w[v_remain])]
    v_start <- m_x[i_start, ]
    t_opt <- proc.time()[[3]]

    # Nelder-Mead copes with -Inf log-posteriors (e.g. outside the prior's
    # support); BFGS then polishes the solution. Nelder-Mead is unreliable in
    # one dimension, so there we go straight to BFGS.
    l_opt <- list(par = v_start, value = neg_log_post(v_start), convergence = NA)
    if (n_d > 1L) {
      l_opt <- tryCatch(
        stats::optim(v_start, neg_log_post, method = "Nelder-Mead",
                     control = list(maxit = control$maxit_nm,
                                    parscale = v_scale)),
        error = function(e) l_opt)
    }
    l_bfgs <- tryCatch(
      stats::optim(l_opt$par, neg_log_post, method = "BFGS", hessian = TRUE,
                   control = list(maxit = control$maxit_bfgs,
                                  parscale = v_scale)),
      error = function(e) NULL)
    if (!is.null(l_bfgs) && is.finite(l_bfgs$value) &&
        l_bfgs$value <= l_opt$value) {
      l_opt <- l_bfgs
      m_hess <- l_bfgs$hessian
    } else {
      m_hess <- tryCatch(stats::optimHess(l_opt$par, neg_log_post,
                                          control = list(parscale = v_scale)),
                         error = function(e) NULL)
    }
    v_mode <- l_opt$par
    names(v_mode) <- v_names
    l_pd <- hessian_to_cov(m_hess, m_prior_prec, m_fallback)

    if (verbose) {
      message(sprintf(
        "  optimizer %d: log-posterior = %.3f, convergence = %s, %.2f s",
        i, -l_opt$value, format(l_opt$convergence),
        proc.time()[[3]] - t_opt))
    }
    l_out[[i]] <- new_component(v_mode, l_pd, origin = "optim", iter = 0L)

    # Exclude the starting point and the n/D draws closest to the new mode,
    # measured in the metric of that mode's covariance.
    v_dist <- dmvnorm_log(m_x[v_remain, , drop = FALSE], v_mode, l_pd$chol)
    v_near <- v_remain[order(v_dist, decreasing = TRUE)[
      seq_len(min(n_exclude, length(v_remain)))]]
    v_remain <- setdiff(v_remain, c(i_start, v_near))
  }
  l_out
}

# Covariance from the Hessian of the negative log-posterior at a mode. When the
# Hessian is not positive definite (a ridge, or a boundary optimum) its
# negative eigenvalues are set to zero and the prior precision is added, as in
# the reference implementation.
hessian_to_cov <- function(m_hess, m_prior_prec, m_fallback) {
  if (is.null(m_hess) || any(!is.finite(m_hess))) {
    return(make_pd(m_fallback))
  }
  m_hess <- (m_hess + t(m_hess)) / 2
  m_chol <- tryCatch(chol(m_hess), error = function(e) NULL)
  if (is.null(m_chol)) {
    l_eig <- eigen(m_hess, symmetric = TRUE)
    m_hess <- l_eig$vectors %*% (pmax(l_eig$values, 0) * t(l_eig$vectors)) +
      m_prior_prec
    m_chol <- chol(m_hess)
  }
  make_pd(chol2inv(m_chol), m_fallback)
}
