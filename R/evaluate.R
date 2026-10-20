# Evaluating the user's log-prior and log-likelihood functions.
#
# The likelihood is almost always the expensive part of a calibration (it runs
# the simulation model), so this is where parallelism happens. Parallel
# evaluation goes through future.apply, which means the user picks the backend
# with future::plan() (multisession, multicore, a cluster, HPC schedulers via
# future.batchtools, ...) and imisr never creates or tears down workers itself.

# Evaluate `fn` on the rows of `m_x` in the current R process.
eval_chunk <- function(m_x, fn, vectorized, catch_errors) {
  n <- nrow(m_x)
  v_names <- colnames(m_x)
  # A vectorized function still gets a (one-row) matrix when rows are
  # evaluated one at a time after a failed batch.
  eval_row <- function(i) {
    if (vectorized) {
      x <- m_x[i, , drop = FALSE]
    } else {
      x <- m_x[i, ]
      names(x) <- v_names
    }
    out <- fn(x)
    if (length(out) != 1L) {
      stop("The function must return a single number for a single ",
           "parameter set; it returned ", length(out), " values.",
           call. = FALSE)
    }
    as.numeric(out)
  }
  eval_rows <- function() {
    if (!catch_errors) {
      return(list(value = vapply(seq_len(n), eval_row, numeric(1)),
                  n_error = 0L, message = NULL))
    }
    v_out <- numeric(n)
    n_error <- 0L
    s_msg <- NULL
    for (i in seq_len(n)) {
      v_out[i] <- tryCatch(eval_row(i), error = function(e) {
        n_error <<- n_error + 1L
        if (is.null(s_msg)) s_msg <<- conditionMessage(e)
        -Inf
      })
    }
    list(value = v_out, n_error = n_error, message = s_msg)
  }

  if (!vectorized) {
    return(eval_rows())
  }
  out <- if (catch_errors) {
    tryCatch(fn(m_x), error = function(e) NULL)
  } else {
    fn(m_x)
  }
  if (is.null(out)) {
    # The batch failed; find the offending rows one at a time.
    return(eval_rows())
  }
  if (length(out) != n) {
    stop("With `vectorized = TRUE` the function must return one value per ",
         "row of its input: got ", length(out), " values for ", n, " rows.",
         call. = FALSE)
  }
  list(value = as.numeric(out), n_error = 0L, message = NULL)
}

# Evaluate `fn` on all rows of `m_x`, optionally splitting the rows into chunks
# that are evaluated in parallel with future.apply.
eval_log_fun <- function(m_x, fn, vectorized, parallel = FALSE,
                         n_chunks = NULL, catch_errors = FALSE,
                         what = "log_likelihood") {
  n <- nrow(m_x)
  if (n == 0L) {
    return(numeric(0))
  }
  if (isTRUE(parallel) && n > 1L) {
    check_future()
    n_chunks <- min(n, n_chunks %||% max(1L, future::nbrOfWorkers()))
    l_idx <- split(seq_len(n), cut(seq_len(n), n_chunks, labels = FALSE))
    l_res <- future.apply::future_lapply(
      l_idx,
      function(v_idx) {
        eval_chunk(m_x[v_idx, , drop = FALSE], fn, vectorized, catch_errors)
      },
      future.seed = TRUE
    )
    l_res <- list(
      value = unlist(lapply(l_res, `[[`, "value"), use.names = FALSE),
      n_error = sum(vapply(l_res, `[[`, integer(1), "n_error")),
      message = unlist(lapply(l_res, `[[`, "message"))[1]
    )
  } else {
    l_res <- eval_chunk(m_x, fn, vectorized, catch_errors)
  }
  v_out <- l_res$value
  if (l_res$n_error > 0L) {
    warning(sprintf(
      "`%s` failed for %d of %d parameter sets; they were given a value of -Inf. First error: %s",
      what, l_res$n_error, n, l_res$message
    ), call. = FALSE)
  }
  v_bad <- is.na(v_out)
  if (any(v_bad)) {
    warning(sprintf(
      "`%s` returned NA/NaN for %d of %d parameter sets; they were given a value of -Inf.",
      what, sum(v_bad), n
    ), call. = FALSE)
    v_out[v_bad] <- -Inf
  }
  if (any(v_out == Inf)) {
    stop("`", what, "` returned +Inf. It must return log-densities, which ",
         "are finite or -Inf.", call. = FALSE)
  }
  v_out
}

check_future <- function() {
  if (!requireNamespace("future.apply", quietly = TRUE) ||
      !requireNamespace("future", quietly = TRUE)) {
    stop("`parallel = TRUE` requires the 'future' and 'future.apply' ",
         "packages. Install them with ",
         "install.packages(c(\"future\", \"future.apply\")).", call. = FALSE)
  }
  invisible(TRUE)
}
