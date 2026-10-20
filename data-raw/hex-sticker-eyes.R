# Hex sticker for imisr, "two eyes" versions: bimodal imisr posteriors with
# each draw colored by the quantile rank of its posterior density (magma
# palette, deciles), so each mode glows like an eye.
#
# 1. data-raw/logo-eyes-crs.png: the cancer relative survival (CRS) model
#    calibrated to survival alone, as in the "calibrating-simulation-models"
#    vignette. It is non-identifiable, so the posterior is a ridge with two
#    mirror-image modes. Rotated 45 degrees and flipped so the ridge arches
#    over the name.
# 2. data-raw/logo-eyes-mixture.png: the bimodal mixture of two correlated
#    normals in the "imisr" vignette (Raftery and Bao, 2010, Section 4.3).
#
# Run from the package root: source("data-raw/hex-sticker-eyes.R")

pkgload::load_all(quiet = TRUE)

# ---- CRS posterior -----------------------------------------------------------
# Only every other survival time point is used, so the targets carry less
# information and the two modes are wider.
v_lb <- c(p_Mets = 0.04, p_DieMets = 0.04)
v_ub <- c(p_Mets = 0.16, p_DieMets = 0.16)
v_idx <- seq(1, nrow(crs_targets$Surv), by = 2)
calc_log_lik_crs <- function(v_params) {
  sum(dnorm(crs_targets$Surv$value[v_idx], crs_markov(v_params)$Surv[v_idx],
            crs_targets$Surv$se[v_idx], log = TRUE))
}
calc_log_prior_crs <- function(v_params) {
  sum(dunif(v_params, v_lb, v_ub, log = TRUE))
}
sample_prior_crs <- function(n_samp) {
  cbind(p_Mets = runif(n_samp, v_lb[["p_Mets"]], v_ub[["p_Mets"]]),
        p_DieMets = runif(n_samp, v_lb[["p_DieMets"]], v_ub[["p_DieMets"]]))
}
set.seed(2018)
fit_crs <- imis(calc_log_lik_crs, calc_log_prior_crs, sample_prior_crs,
                vectorized = FALSE, n_optim = 5, n_resample = 10000,
                verbose = FALSE)
m_crs <- fit_crs$resample
v_lpost_crs <- apply(m_crs, 1, calc_log_lik_crs)

# ---- Mixture posterior -------------------------------------------------------
# Equal mixture of N(0, AR(-0.95)) and N(9, AR(0.95)) on a U(-3, 12)^2 prior
calc_log_lik_mix <- function(m_params) {
  log_dens <- function(v_mean, n_rho) {
    m_chol <- chol(matrix(c(1, n_rho, n_rho, 1), 2))
    dmvnorm_log(m_params, v_mean, m_chol)
  }
  log(0.5) + log_add(log_dens(c(0, 0), -0.95), log_dens(c(9, 9), 0.95))
}
calc_log_prior_mix <- function(m_params) {
  rowSums(dunif(m_params, -3, 12, log = TRUE))
}
sample_prior_mix <- function(n_samp) {
  matrix(runif(2 * n_samp, -3, 12), n_samp, 2,
         dimnames = list(NULL, c("a", "b")))
}
set.seed(2)
fit_mix <- imis(calc_log_lik_mix, calc_log_prior_mix, sample_prior_mix,
                n_optim = 5, n_resample = 10000, verbose = FALSE)
m_mix <- fit_mix$resample
v_lpost_mix <- calc_log_lik_mix(m_mix)

# ---- Geometry ----------------------------------------------------------------
# Pointy-top hexagon with circumradius 1, centered at the origin.
v_ang <- pi / 2 + (0:5) * pi / 3
m_hex <- cbind(cos(v_ang), sin(v_ang))

rotate <- function(m, n_theta) {
  m %*% t(matrix(c(cos(n_theta), sin(n_theta), -sin(n_theta), cos(n_theta)),
                 2))
}
# CRS: center on the prior box, rotate so the diagonal is vertical, flip so
# the ridge arches, scale and shift up (as in data-raw/hex-sticker.R).
m_xy_crs <- rotate(sweep(m_crs, 2, c(0.1, 0.1)), pi / 4) * 7.3
m_xy_crs[, 2] <- -m_xy_crs[, 2] + 0.22
# The ridge is very thin, so widen it for display: fit its center line as a
# polynomial in x and stretch each draw's distance from that line by
# `n_widen`. Positions along the ridge, and the density ranks, are unchanged.
fit_ridge <- lm(y ~ poly(x, 4), data.frame(x = m_xy_crs[, 1], y = m_xy_crs[, 2]))
v_center <- fitted(fit_ridge)
widen_crs <- function(n_widen) {
  cbind(m_xy_crs[, 1], v_center + n_widen * (m_xy_crs[, 2] - v_center))
}
# Upside-down version: rotate the image 180 degrees in place, so the ridge
# dips into a "U" in the same area of the sticker.
flip_crs <- function(m_xy) {
  v_mid <- colMeans(apply(m_xy, 2, quantile, probs = c(0.02, 0.98)))
  sweep(-sweep(m_xy, 2, v_mid), 2, v_mid, "+")
}
# Mixture: center between the modes, rotate so they sit side by side, scale
# and shift up.
m_xy_mix <- rotate(sweep(m_mix, 2, c(4.5, 4.5)), -pi / 4) * 0.072
m_xy_mix[, 2] <- m_xy_mix[, 2] + 0.36

# ---- Palette -----------------------------------------------------------------
col_bg <- "#0E0B16"
col_border <- "#FCA311"
col_text <- "#FFFFFF"
color_by_rank <- function(v_q) {
  v_col <- viridisLite::magma(10)[cut(v_q, seq(0, 1, by = 0.1),
                                      include.lowest = TRUE)]
  # Transparency grows with rank, from 0.1 to 1 (ggplot2's default range)
  v_alpha <- 0.1 + 0.9 * (v_q - min(v_q)) / diff(range(v_q))
  mapply(adjustcolor, v_col, alpha.f = v_alpha, USE.NAMES = FALSE)
}

# ---- Draw --------------------------------------------------------------------
draw_sticker <- function(m_xy, v_lpost, n_cex) {
  par(mar = c(0, 0, 0, 0))
  plot.new()
  plot.window(xlim = c(-0.866, 0.866), ylim = c(-1, 1), asp = 1, xaxs = "i",
              yaxs = "i")
  polygon(m_hex * 0.985, col = col_bg, border = NA)

  # Highest-density draws last, so the bright cores sit on top. Leave out the
  # few draws that would touch the border.
  v_qrank <- rank(v_lpost) / length(v_lpost)
  v_inside <- abs(m_xy[, 1]) < 0.8 & abs(m_xy[, 2]) < 0.92 - abs(m_xy[, 1]) / sqrt(3)
  v_ord <- order(v_qrank)
  v_ord <- v_ord[v_inside[v_ord]]
  points(m_xy[v_ord, ], pch = 16, cex = n_cex,
         col = color_by_rank(v_qrank[v_ord]))

  text(0, -0.06, "calibrating models that might not identify", cex = 1.45,
       col = adjustcolor(col_text, 0.75), font = 3, family = "Avenir Next")
  text(0, -0.42, "imisr", col = col_text, cex = 7.5, font = 2,
       family = "Avenir Next")
  text(0.40, -0.70, "FerAlLab-git/imisr", srt = 30, cex = 1.15,
       col = adjustcolor(col_text, 0.7), family = "Avenir Next")

  polygon(m_hex * 0.985, col = NA, border = col_border, lwd = 14)
}

save_sticker <- function(s_file, ...) {
  png(s_file, width = 1036, height = 1200, res = 144, type = "quartz",
      bg = "transparent")
  draw_sticker(...)
  dev.off()
}
save_sticker("data-raw/logo-eyes-crs.png", widen_crs(8), v_lpost_crs, 1.1)
save_sticker("data-raw/logo-eyes-crs-flipped.png", flip_crs(widen_crs(8)),
             v_lpost_crs, 1.1)
# Narrower, angrier eyes
for (n_widen in c(4, 3)) {
  save_sticker(sprintf("data-raw/logo-eyes-crs-flipped-%d.png", n_widen),
               flip_crs(widen_crs(n_widen)), v_lpost_crs, 1.1)
}
save_sticker("data-raw/logo-eyes-mixture.png", m_xy_mix, v_lpost_mix, 0.8)
