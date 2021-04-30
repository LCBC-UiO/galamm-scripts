dat <- readRDS("data/cvlt_model/simulated_data.rds")

lambda_init <- c(1, runif(4, min = 1, max = 2))
names(lambda_init) <- paste0("a", 1:5)

formula <- as.formula(
  cbind(successes, failures) ~ 0 + a2 + a3 + a4 + a5 + retest +
    s(age_z, by = weight, k = 15, bs = "cr")
)
random <- ~ (0 + weight | id / timepoint)
