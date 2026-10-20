#' Calibration targets for the cancer relative survival model
#'
#' Simulated observed survival of a cohort of 200 patients, used to calibrate
#' [crs_markov()], the cancer relative survival model of Alarid-Escudero et
#' al. (2018).
#'
#' @format A list with one element, `Surv`, a data frame with 59 rows (cycles
#'   2 to 60) and columns:
#' \describe{
#'   \item{time}{Cycle (years).}
#'   \item{n_total}{Number of patients at the start of follow-up.}
#'   \item{n_alive}{Number of patients alive.}
#'   \item{value}{Observed proportion alive.}
#'   \item{se}{Standard error of `value`.}
#'   \item{lb, ub}{Lower and upper bounds of the 95% confidence interval.}
#' }
#' @source Simulated for teaching calibration of decision models (DARTH
#'   workgroup; Stanford HRP 263).
#' @references
#' Alarid-Escudero F, MacLehose RF, Peralta Y, Kuntz KM, Enns EA (2018).
#' Nonidentifiability in model calibration and implications for medical
#' decision making. *Medical Decision Making*, 38(7), 810-821.
#' \doi{10.1177/0272989X18792283}
"crs_targets"

#' Calibration targets for the Sick-Sicker model
#'
#' Simulated observed survival, disease prevalence and proportion of the sick
#' in the "sick" (S1) state, for a cohort of 500 people, used to calibrate
#' [sick_sicker_markov()], the Sick-Sicker model of the DARTH tutorials
#' (Krijkamp et al., 2018; Alarid-Escudero et al., 2023a, 2023b).
#'
#' @format A list of three data frames, each with columns `time`, `value`
#'   (observed proportion), `se` (its standard error) and `lb`, `ub` (95%
#'   confidence bounds), plus the counts the proportion was computed from:
#' \describe{
#'   \item{Surv}{Survival at cycles 1 to 30 (`n_total`, `n_alive`).}
#'   \item{Prev}{Prevalence of disease (S1 or S2) among the living at cycles
#'     1 to 30 (`n_alive`, `n_sick`).}
#'   \item{PropSick}{Proportion of the sick who are in S1 at cycles 10, 20 and
#'     30 (`n_total_sick`, `n_sick`).}
#' }
#' @source Simulated for teaching calibration of decision models (DARTH
#'   workgroup; Stanford HRP 263).
#' @references
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
"sick_sicker_targets"
