test_that("the CRS model returns survival and is symmetric in its parameters", {
  out <- crs_markov(c(p_Mets = 0.1, p_DieMets = 0.05))
  expect_named(out, "Surv")
  expect_length(out$Surv, nrow(crs_targets$Surv))
  expect_true(all(diff(out$Surv) <= 0))
  # Non-identifiability: swapping the parameters gives the same survival.
  expect_equal(out$Surv, crs_markov(c(p_Mets = 0.05, p_DieMets = 0.1))$Surv)
})

test_that("the Sick-Sicker model returns outputs matching its targets", {
  out <- sick_sicker_markov(c(p_S1S2 = 0.105, hr_S1 = 3, hr_S2 = 10))
  expect_named(out, c("Surv", "Prev", "PropSick"))
  for (nm in names(out)) {
    expect_length(out[[nm]], nrow(sick_sicker_targets[[nm]]))
    expect_true(all(out[[nm]] >= 0 & out[[nm]] <= 1))
  }
})

test_that("models validate their inputs", {
  expect_error(crs_markov(c(p_mets = 0.1, p_DieMets = 0.05)), "p_Mets")
  expect_error(crs_markov(c(p_Mets = 1.2, p_DieMets = 0.05)), "Invalid")
  expect_error(sick_sicker_markov(c(p_S1S2 = 0.49, hr_S1 = 12, hr_S2 = 5)),
               "Invalid")
})

test_that("target data sets are well formed", {
  expect_named(sick_sicker_targets, c("Surv", "Prev", "PropSick"))
  for (df in c(crs_targets, sick_sicker_targets)) {
    expect_true(all(c("time", "value", "se", "lb", "ub") %in% names(df)))
    expect_true(all(df$se > 0))
  }
})
