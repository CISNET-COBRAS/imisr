# Hex sticker for imisr: man/figures/logo.png
#
# The image is the imisr posterior of the cancer relative survival (CRS) model
# calibrated to survival alone, as in the "calibrating-simulation-models"
# vignette. The model is non-identifiable: survival is symmetric in p_Mets and
# p_DieMets, so the posterior is a ridge with two mirror-image modes. Each draw
# is colored by the decile of its posterior density rank (magma palette), so
# each mode glows like an eye.
#
# Run from the package root: source("data-raw/hex-sticker.R")

pkgload::load_all(quiet = TRUE)

# ---- Posterior ---------------------------------------------------------------
# Only every other survival time point is used, so the targets carry less
# information and the two modes are wider.
v_lb <- c(p_Mets = 0.04, p_DieMets = 0.04)
v_ub <- c(p_Mets = 0.16, p_DieMets = 0.16)
v_idx <- seq(1, nrow(crs_targets$Surv), by = 2)
calc_log_lik <- function(v_params) {
  sum(dnorm(crs_targets$Surv$value[v_idx], crs_markov(v_params)$Surv[v_idx],
            crs_targets$Surv$se[v_idx], log = TRUE))
}
calc_log_prior <- function(v_params) {
  sum(dunif(v_params, v_lb, v_ub, log = TRUE))
}
sample_prior <- function(n_samp) {
  cbind(p_Mets = runif(n_samp, v_lb[["p_Mets"]], v_ub[["p_Mets"]]),
        p_DieMets = runif(n_samp, v_lb[["p_DieMets"]], v_ub[["p_DieMets"]]))
}
set.seed(2018)
fit <- imis(calc_log_lik, calc_log_prior, sample_prior, vectorized = FALSE,
            n_optim = 5, n_resample = 10000, verbose = FALSE)
m_post <- fit$resample
# The prior is uniform, so the log-likelihood ranks the posterior density
v_lpost <- apply(m_post, 1, calc_log_lik)

# ---- Geometry ----------------------------------------------------------------
# Pointy-top hexagon with circumradius 1, centered at the origin.
v_ang <- pi / 2 + (0:5) * pi / 3
m_hex <- cbind(cos(v_ang), sin(v_ang))

# Center on the prior box, rotate 45 degrees so the line p_Mets = p_DieMets is
# vertical, flip so the ridge arches, then scale and shift up.
n_theta <- pi / 4
m_rot <- matrix(c(cos(n_theta), sin(n_theta), -sin(n_theta), cos(n_theta)), 2)
m_xy <- sweep(m_post, 2, c(0.1, 0.1)) %*% t(m_rot) * 7.3
m_xy[, 2] <- -m_xy[, 2] + 0.22
# The ridge is very thin, so widen it for display: fit its center line as a
# polynomial in x and stretch each draw's distance from that line by
# `n_widen`. Positions along the ridge, and the density ranks, are unchanged.
n_widen <- 4
fit_ridge <- lm(y ~ poly(x, 4), data.frame(x = m_xy[, 1], y = m_xy[, 2]))
v_center <- fitted(fit_ridge)
m_xy[, 2] <- v_center + n_widen * (m_xy[, 2] - v_center)
# Rotate 180 degrees in place, so the ridge dips into a "U" and the two modes
# slant down toward the center like angry eyes.
v_mid <- colMeans(apply(m_xy, 2, quantile, probs = c(0.02, 0.98)))
m_xy <- sweep(-sweep(m_xy, 2, v_mid), 2, v_mid, "+")

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
draw_sticker <- function() {
  par(mar = c(0, 0, 0, 0))
  plot.new()
  plot.window(xlim = c(-0.866, 0.866), ylim = c(-1, 1), asp = 1, xaxs = "i",
              yaxs = "i")
  polygon(m_hex * 0.985, col = col_bg, border = NA)

  # Highest-density draws last, so the bright cores sit on top. Leave out the
  # few draws that would touch the border.
  v_qrank <- rank(v_lpost) / length(v_lpost)
  v_inside <- abs(m_xy[, 1]) < 0.8 &
    abs(m_xy[, 2]) < 0.92 - abs(m_xy[, 1]) / sqrt(3)
  v_ord <- order(v_qrank)
  v_ord <- v_ord[v_inside[v_ord]]
  points(m_xy[v_ord, ], pch = 16, cex = 1.1,
         col = color_by_rank(v_qrank[v_ord]))

  text(0, -0.06, "calibrating models that might not identify", cex = 1.45,
       col = adjustcolor(col_text, 0.75), font = 3, family = "Avenir Next")
  text(0, -0.42, "imisr", col = col_text, cex = 7.5, font = 2,
       family = "Avenir Next")
  text(0.40, -0.70, "FerAlLab-git/imisr", srt = 30, cex = 1.15,
       col = adjustcolor(col_text, 0.7), family = "Avenir Next")

  polygon(m_hex * 0.985, col = NA, border = col_border, lwd = 14)
}

# 518 x 600 px is the size pkgdown and usethis::use_logo() expect; draw at
# twice that for high-density screens.
n_width <- 1036 / 144
n_height <- 1200 / 144
dir.create("man/figures", showWarnings = FALSE, recursive = TRUE)
png("man/figures/logo.png", width = 1036, height = 1200, res = 144,
    type = "quartz", bg = "transparent")
draw_sticker()
dev.off()

# Vector files for printing, kept out of the package (data-raw/ is in
# .Rbuildignore). The Quartz PDF embeds the Avenir Next font, so it displays
# the same everywhere; the SVG names the font, so it needs Avenir Next
# installed to display as designed.
quartz(type = "pdf", file = "data-raw/logo.pdf", width = n_width,
       height = n_height, bg = "transparent")
draw_sticker()
dev.off()
svglite::svglite("data-raw/logo.svg", width = n_width, height = n_height,
                 bg = "transparent")
draw_sticker()
dev.off()
