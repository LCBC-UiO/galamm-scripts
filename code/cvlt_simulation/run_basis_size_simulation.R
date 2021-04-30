# Additional simulations to check the effect of basis dimension on coverage
library(tidyverse)
library(numDeriv)
library(furrr)
library(gamm4)

source("code/cvlt_simulation/sample_data_fun.R")
source("code/cvlt_simulation/model_setup.R")
sim_params <- readRDS("data/cvlt_simulation/sim_params.rds")

# Overrun the formula from the model setup
formula <- as.formula(
  cbind(successes, failures) ~ 0 + a2 + a3 + a4 + a5 + retest +
    s(age_z, by = weight, k = 30, bs = "cr")
)

dimensions <- crossing(
  iter = seq_len(500),
  n = 1000
)

plan(multisession)
sims <- future_pwalk(dimensions, function(iter, n){
  outfile <- file.path("results/cvlt_simulation/model_fits", 
                       paste0("mod_n", n, "_iter", iter, "_dimension30.rds"))
  if(file.exists(outfile)) return(0)
  
  res <- run_simulation(n, iter, sim_params)
  
  saveRDS(res, outfile)
  }, .options = furrr_options(seed = 4433L))

plan("default")
