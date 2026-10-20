# Calibration targets for the example models, taken from the HRP 263
# (Advanced Decision Science) calibration lab, R-labs/04-calibration.
src <- "~/GitHub/hrp-263-adv-dec-sci/R-labs/04-calibration"

e <- new.env()
load(file.path(src, "CRS_CalibTargets.RData"), envir = e)
load(file.path(src, "SickSicker_CalibTargets.RData"), envir = e)

crs_targets <- e$CRS_targets
sick_sicker_targets <- e$SickSicker_targets

save(crs_targets, file = "data/crs_targets.rda", compress = "xz")
save(sick_sicker_targets, file = "data/sick_sicker_targets.rda",
     compress = "xz")
