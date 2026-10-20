#' Example simulation models for calibration
#'
#' Two cohort state-transition (Markov) models used in the package vignettes
#' and in health decision science courses to teach model calibration. Both
#' take a named vector of the parameters to be calibrated and return the
#' model outputs that correspond to the calibration targets
#' [crs_targets] and [sick_sicker_targets].
#'
#' **Cancer relative survival (CRS) model** (`crs_markov()`), first described
#' by Alarid-Escudero et al. (2018) to illustrate non-identifiability in
#' calibration. Three states: no evidence of disease (NED), metastatic
#' disease (Mets) and death. Each annual cycle, patients in NED develop
#' metastases with probability `p_Mets`, and patients with metastases die
#' with probability `p_DieMets`.
#' The model runs for 60 cycles and returns overall survival at cycles 2 to
#' 60. Survival depends on the two parameters symmetrically, so survival data
#' alone cannot tell them apart: the model is *non-identifiable* and its
#' posterior has two modes.
#'
#' **Sick-Sicker model** (`sick_sicker_markov()`), used throughout the DARTH
#' tutorials on decision models in R (Krijkamp et al., 2018; Alarid-Escudero
#' et al., 2023a, 2023b). Four states: healthy (H), sick (S1), sicker (S2)
#' and dead (D), for a cohort followed from age 25 to 55 in annual cycles.
#' Calibrated parameters are the probability of progressing from S1 to S2
#' (`p_S1S2`), and the hazard ratios of death in S1 (`hr_S1`) and S2
#' (`hr_S2`) relative to healthy. Other parameters are fixed:
#' probability of death when healthy 0.005, of becoming sick 0.15, and of
#' recovering from S1 0.5. Returns survival and disease prevalence for cycles
#' 1 to 30, and the proportion of the sick who are in S1 at cycles 10, 20 and
#' 30.
#'
#' Both functions stop with an error when a parameter is missing or the
#' implied transition matrix is invalid (e.g. `p_S1S2` large enough to make
#' the probability of staying in S1 negative).
#'
#' @param params Named numeric vector of model parameters.
#' @return A named list of numeric vectors: `Surv` for `crs_markov()`;
#'   `Surv`, `Prev` and `PropSick` for `sick_sicker_markov()`.
#' @references
#' Alarid-Escudero F, MacLehose RF, Peralta Y, Kuntz KM, Enns EA (2018).
#' Nonidentifiability in model calibration and implications for medical
#' decision making. *Medical Decision Making*, 38(7), 810-821.
#' \doi{10.1177/0272989X18792283}
#'
#' Alarid-Escudero F, Krijkamp E, Enns EA, Yang A, Hunink MGM, Pechlivanoglou
#' P, Jalal H (2023a). An introductory tutorial on cohort state-transition
#' models in R using a cost-effectiveness analysis example. *Medical Decision
#' Making*, 43(1), 3-20. \doi{10.1177/0272989X221103163}
#'
#' Alarid-Escudero F, Krijkamp E, Enns EA, Yang A, Hunink MGM, Pechlivanoglou
#' P, Jalal H (2023b). A tutorial on time-dependent cohort state-transition
#' models in R using a cost-effectiveness analysis example. *Medical Decision
#' Making*, 43(1), 21-41. \doi{10.1177/0272989X221121747}
#'
#' Krijkamp EM, Alarid-Escudero F, Enns EA, Jalal HJ, Hunink MGM,
#' Pechlivanoglou P (2018). Microsimulation modeling for health decision
#' sciences using R: a tutorial. *Medical Decision Making*, 38(3), 400-422.
#' \doi{10.1177/0272989X18754513}
#' @examples
#' crs_markov(c(p_Mets = 0.10, p_DieMets = 0.05))$Surv[1:5]
#' str(sick_sicker_markov(c(p_S1S2 = 0.105, hr_S1 = 3, hr_S2 = 10)))
#' @name example_models
NULL

#' @rdname example_models
#' @export
crs_markov <- function(params) {
  check_params(params, c("p_Mets", "p_DieMets"))
  n_t <- 60L
  v_n <- c("NED", "Mets", "Death")
  m_p <- matrix(0, 3, 3, dimnames = list(v_n, v_n))
  m_p["NED", "NED"] <- 1 - params[["p_Mets"]]
  m_p["NED", "Mets"] <- params[["p_Mets"]]
  m_p["Mets", "Mets"] <- 1 - params[["p_DieMets"]]
  m_p["Mets", "Death"] <- params[["p_DieMets"]]
  m_p["Death", "Death"] <- 1
  check_transition_matrix(m_p)

  m_trace <- run_cohort(m_p, c(1, 0, 0), n_t)
  v_os <- 1 - m_trace[, "Death"]
  list(Surv = unname(v_os[-(1:2)]))
}

#' @rdname example_models
#' @export
sick_sicker_markov <- function(params) {
  check_params(params, c("p_S1S2", "hr_S1", "hr_S2"))
  n_t <- 55L - 25L
  v_n <- c("H", "S1", "S2", "D")
  p_HD <- 0.005
  p_HS1 <- 0.15
  p_S1H <- 0.5
  r_HD <- -log(1 - p_HD)
  p_S1D <- 1 - exp(-params[["hr_S1"]] * r_HD)
  p_S2D <- 1 - exp(-params[["hr_S2"]] * r_HD)
  p_S1S2 <- params[["p_S1S2"]]

  m_p <- matrix(0, 4, 4, dimnames = list(v_n, v_n))
  m_p["H", ] <- c(1 - (p_HS1 + p_HD), p_HS1, 0, p_HD)
  m_p["S1", ] <- c(p_S1H, 1 - (p_S1H + p_S1S2 + p_S1D), p_S1S2, p_S1D)
  m_p["S2", ] <- c(0, 0, 1 - p_S2D, p_S2D)
  m_p["D", "D"] <- 1
  check_transition_matrix(m_p)

  m_trace <- run_cohort(m_p, c(1, 0, 0, 0), n_t)
  v_os <- 1 - m_trace[, "D"]
  v_sick <- m_trace[, "S1"] + m_trace[, "S2"]
  list(Surv = unname(v_os[-1]),
       Prev = unname((v_sick / v_os)[-1]),
       PropSick = unname((m_trace[, "S1"] / v_sick)[c(11, 21, 31)]))
}

run_cohort <- function(m_p, v_init, n_t) {
  m_trace <- matrix(NA_real_, n_t + 1L, length(v_init),
                    dimnames = list(0:n_t, colnames(m_p)))
  m_trace[1, ] <- v_init
  for (t in seq_len(n_t)) {
    m_trace[t + 1L, ] <- m_trace[t, ] %*% m_p
  }
  m_trace
}

check_params <- function(params, v_required) {
  v_missing <- setdiff(v_required, names(params))
  if (length(v_missing) > 0L) {
    stop("Missing model parameter(s): ", paste(v_missing, collapse = ", "),
         call. = FALSE)
  }
}

# A valid transition matrix has every entry in [0, 1] and rows summing to 1.
# Checking the row sums alone is not enough: c(1.2, -0.2, 0) sums to 1.
check_transition_matrix <- function(m_p) {
  if (any(!is.finite(m_p)) || any(m_p < 0 | m_p > 1)) {
    stop("Invalid transition matrix: entries outside [0, 1].", call. = FALSE)
  }
  if (!isTRUE(all.equal(unname(rowSums(m_p)), rep(1, nrow(m_p))))) {
    stop("Invalid transition matrix: rows do not sum to 1.", call. = FALSE)
  }
}
