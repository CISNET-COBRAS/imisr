# Hex sticker for imisr.
#
# The image is the posterior of the cancer relative survival (CRS) model
# calibrated to survival alone (see the "calibrating-simulation-models"
# vignette). The model is non-identifiable: survival is symmetric in p_Mets
# and p_DieMets, so the posterior has two mirror-image modes. The plot is
# rotated 45 degrees so that the line of symmetry (p_Mets = p_DieMets) is
# vertical. Dots are posterior draws and the ellipses are the normal mixture
# components IMIS added, colored from the first to the last.
#
# Run from the package root: source("data-raw/hex-sticker.R")

pkgload::load_all(quiet = TRUE)

# ---- Fit ---------------------------------------------------------------------
v_lb <- c(p_Mets = 0.04, p_DieMets = 0.04)
v_ub <- c(p_Mets = 0.16, p_DieMets = 0.16)
calc_log_lik <- function(v_params) {
  sum(dnorm(crs_targets$Surv$value, crs_markov(v_params)$Surv,
            crs_targets$Surv$se, log = TRUE))
}
calc_log_prior <- function(v_params) sum(dunif(v_params, v_lb, v_ub, log = TRUE))
sample_prior <- function(n_samp) {
  cbind(p_Mets = runif(n_samp, v_lb[["p_Mets"]], v_ub[["p_Mets"]]),
        p_DieMets = runif(n_samp, v_lb[["p_DieMets"]], v_ub[["p_DieMets"]]))
}
set.seed(2018)
fit <- imis(calc_log_lik, calc_log_prior, sample_prior, vectorized = FALSE,
            n_optim = 5, verbose = FALSE)

# Posterior density on a grid, for the contour lines
v_grid <- seq(0.04, 0.16, length.out = 160)
m_grid <- as.matrix(expand.grid(p_Mets = v_grid, p_DieMets = v_grid))
v_llik <- apply(m_grid, 1, calc_log_lik)
m_dens <- matrix(exp(v_llik - max(v_llik)), length(v_grid))

# ---- Geometry ----------------------------------------------------------------
# Pointy-top hexagon with circumradius 1, centered at the origin.
v_ang <- pi / 2 + (0:5) * pi / 3
m_hex <- cbind(cos(v_ang), sin(v_ang))

# Map parameter space to sticker space: center on the prior box, rotate 45
# degrees so the diagonal is vertical, flip so the ridge arches over the name,
# scale, and shift up.
v_mid <- c(0.1, 0.1)
n_scale <- 7.3
v_shift <- c(0, 0.22)
m_rot <- matrix(c(cos(pi / 4), sin(pi / 4), -sin(pi / 4), cos(pi / 4)), 2)
to_hex <- function(m) {
  m <- sweep(as.matrix(m), 2, v_mid) %*% t(m_rot) * n_scale
  m[, 2] <- -m[, 2]
  sweep(m, 2, v_shift, "+")
}
ellipse <- function(v_mean, m_cov, n_sd = 2, n_pts = 100) {
  v_t <- seq(0, 2 * pi, length.out = n_pts)
  m_circ <- cbind(cos(v_t), sin(v_t)) * n_sd
  sweep(m_circ %*% chol(m_cov), 2, v_mean, "+")
}

# ---- Palette -----------------------------------------------------------------
col_bg <- "#14213D"
col_border <- "#FCA311"
col_draws <- adjustcolor("#E5E5E5", 0.3)
col_text <- "#FFFFFF"
v_col_comp <- hcl.colors(length(fit$components), "Plasma")

# ---- Draw --------------------------------------------------------------------
draw_sticker <- function() {
  par(mar = c(0, 0, 0, 0))
  plot.new()
  plot.window(xlim = c(-0.866, 0.866), ylim = c(-1, 1), asp = 1, xaxs = "i",
              yaxs = "i")
  polygon(m_hex * 0.985, col = col_bg, border = NA)

  # Line of symmetry: p_Mets = p_DieMets
  segments(0, 0.1, 0, 0.82, col = adjustcolor(col_text, 0.3), lty = 2,
           lwd = 1.5)
  # Posterior density contours
  for (l_cl in contourLines(v_grid, v_grid, m_dens,
                            levels = exp(-c(0.5, 2, 4.5, 8, 12.5)))) {
    lines(to_hex(cbind(l_cl$x, l_cl$y)), col = adjustcolor("#8ECAE6", 0.45),
          lwd = 1.3)
  }

  points(to_hex(fit$resample), pch = 16, cex = 0.45, col = col_draws)
  for (i in seq_along(fit$components)) {
    l_comp <- fit$components[[i]]
    lines(to_hex(ellipse(l_comp$mean, l_comp$cov)), col = v_col_comp[i],
          lwd = 2.2)
  }

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
dir.create("man/figures", showWarnings = FALSE, recursive = TRUE)
png("man/figures/logo.png", width = 1036, height = 1200, res = 144,
    type = "quartz", bg = "transparent")
draw_sticker()
dev.off()
svglite::svglite("data-raw/logo.svg", width = 1036 / 144, height = 1200 / 144,
                 bg = "transparent")
draw_sticker()
dev.off()
