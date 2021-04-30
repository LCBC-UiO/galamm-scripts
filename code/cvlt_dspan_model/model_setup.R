dat <- readRDS("data/cvlt_dspan_model/simulated_data.rds")
formula <- cbind(successes, failures) ~ 0 + a1 + a2 + a3 + a4 + a5 + 
  bwd + fwd + retestcvlt + retestdigitspan +
  t2(age_z, itemgroup, by = weight, bs = c("cr", "re"), k = 10, full = TRUE)

random <- ~ (0 + weight:itemgroup | id / timepoint) 
