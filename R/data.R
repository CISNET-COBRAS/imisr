#' Calibration targets for the cancer relative survival model
#'
#' Simulated observed survival of a cohort of 200 patients, used to calibrate
#' [crs_markov()].
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
"crs_targets"

#' Calibration targets for the Sick-Sicker model
#'
#' Simulated observed survival, disease prevalence and proportion of the sick
#' in the "sick" (S1) state, for a cohort of 500 people, used to calibrate
#' [sick_sicker_markov()].
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
"sick_sicker_targets"
