# Hex sticker for imisr, "eye of Sauron" version.
#
# The image is the imisr posterior of the two parameters of a Gompertz
# survival model (rate and shape) calibrated to a simulated cohort. The two
# parameters are strongly negatively correlated, so the posterior is a long,
# thin cloud. Each draw is colored by the quantile rank of its posterior
# density (magma palette), so the bright core is the highest-density region
# and the bands are highest-density credible regions. The cloud is rotated so
# its main axis is horizontal.
#
# Run from the package root: source("data-raw/hex-sticker-eye.R")

pkgload::load_all(quiet = TRUE)

# ---- Data: Gompertz event times, censored at the end of follow-up ------------
set.seed(263)
n_obs <- 1000
n_cens <- 80
v_t_event <- log(1 - 0.0135 / 0.0045 * log(runif(n_obs))) / 0.0135
v_time <- pmin(v_t_event, n_cens)
v_event <- as.numeric(v_t_event <= n_cens)

# ---- Fit ---------------------------------------------------------------------
# Hazard h(t) = rate * exp(shape * t). The log-likelihood is the sum over
# people of event * log h(t) + log S(t), computed for all parameter sets at
# once with one column per parameter set.
calc_log_lik <- function(m_params) {
  v_rate <- exp(m_params[, "log_rate"])
  v_shape <- m_params[, "shape"]
  m_shape_t <- outer(v_time, v_shape)
  colSums(v_event * (rep(log(v_rate), each = n_obs) + m_shape_t)) -
    colSums(sweep(expm1(m_shape_t), 2, v_rate / v_shape, "*"))
}
calc_log_prior <- function(m_params) {
  dnorm(m_params[, "log_rate"], log(0.005), 1, log = TRUE) +
    dunif(m_params[, "shape"], 0.0001, 0.05, log = TRUE)
}
sample_prior <- function(n_samp) {
  cbind(log_rate = rnorm(n_samp, log(0.005), 1),
        shape = runif(n_samp, 0.0001, 0.05))
}
fit <- imis(calc_log_lik, calc_log_prior, sample_prior, n_resample = 10000,
            verbose = FALSE)

# Quantile rank of each draw's posterior density
m_post <- fit$resample
v_lpost <- calc_log_lik(m_post) + calc_log_prior(m_post)
v_qrank <- rank(v_lpost) / length(v_lpost)

# ---- Geometry ----------------------------------------------------------------
# Pointy-top hexagon with circumradius 1, centered at the origin.
v_ang <- pi / 2 + (0:5) * pi / 3
m_hex <- cbind(cos(v_ang), sin(v_ang))

# Plot rate (not log-rate) against shape, as in the original figure.
# Standardize both, rotate so the main axis is horizontal, scale and shift.
m_xy <- scale(cbind(exp(m_post[, "log_rate"]), m_post[, "shape"]))
v_axis <- eigen(cov(m_xy))$vectors[, 1]
n_theta <- -atan2(v_axis[2], v_axis[1])
m_rot <- matrix(c(cos(n_theta), sin(n_theta), -sin(n_theta), cos(n_theta)), 2)
m_xy <- m_xy %*% t(m_rot)
# Scale so 99.5% of draws fit across the sticker, open the eye vertically
# (the axes are standardized, so their relative scale is arbitrary), and drop
# the few draws that would touch the border.
n_scale <- 0.7 / quantile(abs(m_xy[, 1]), 0.995)
m_xy <- sweep(m_xy * n_scale, 2, c(1, 1.6), "*")
v_keep <- abs(m_xy[, 1]) < 0.8
m_xy <- sweep(m_xy, 2, c(0, 0.32), "+")

# ---- Palette -----------------------------------------------------------------
col_bg <- "#0E0B16"
col_border <- "#FCA311"
col_text <- "#FFFFFF"
color_by_rank <- function(v_q, deciles) {
  if (deciles) {
    v_col <- viridisLite::magma(10)[cut(v_q, seq(0, 1, by = 0.1),
                                        include.lowest = TRUE)]
  } else {
    v_col <- viridisLite::magma(256)[pmax(1, ceiling(v_q * 256))]
  }
  # Transparency grows with rank, from 0.1 to 1 (ggplot2's default range)
  v_alpha <- 0.1 + 0.9 * (v_q - min(v_q)) / diff(range(v_q))
  mapply(adjustcolor, v_col, alpha.f = v_alpha, USE.NAMES = FALSE)
}

# ---- Draw --------------------------------------------------------------------
draw_sticker <- function(deciles = FALSE) {
  par(mar = c(0, 0, 0, 0))
  plot.new()
  plot.window(xlim = c(-0.866, 0.866), ylim = c(-1, 1), asp = 1, xaxs = "i",
              yaxs = "i")
  polygon(m_hex * 0.985, col = col_bg, border = NA)

  # Highest-density draws last, so the bright core sits on top
  v_ord <- order(v_qrank)
  v_ord <- v_ord[v_keep[v_ord]]
  points(m_xy[v_ord, ], pch = 16, cex = 0.75,
         col = color_by_rank(v_qrank[v_ord], deciles))

  text(0, -0.06, "calibrating models that might not identify", cex = 1.45,
       col = adjustcolor(col_text, 0.75), font = 3, family = "Avenir Next")
  text(0, -0.42, "imisr", col = col_text, cex = 7.5, font = 2,
       family = "Avenir Next")
  text(0.40, -0.70, "FerAlLab-git/imisr", srt = 30, cex = 1.15,
       col = adjustcolor(col_text, 0.7), family = "Avenir Next")

  polygon(m_hex * 0.985, col = NA, border = col_border, lwd = 14)
}

for (deciles in c(FALSE, TRUE)) {
  s_file <- if (deciles) "data-raw/logo-eye-deciles.png" else
    "data-raw/logo-eye.png"
  png(s_file, width = 1036, height = 1200, res = 144, type = "quartz",
      bg = "transparent")
  draw_sticker(deciles)
  dev.off()
}
