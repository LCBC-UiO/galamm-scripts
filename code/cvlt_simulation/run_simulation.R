library(tidyverse)
library(numDeriv)
library(furrr)
library(gamm4)

source("code/cvlt_simulation/sample_data_fun.R")
source("code/cvlt_simulation/model_setup.R")
sim_params <- readRDS("data/cvlt_simulation/sim_params.rds")

dimensions <- crossing(
  iter = seq_len(500),
  n = c(200, 500, 1000)
)

plan(multisession)
sims <- future_pwalk(dimensions, function(iter, n){
  outfile <- file.path("results/cvlt_simulation/model_fits", 
                       paste0("mod_n", n, "_iter", iter, ".rds"))
  if(file.exists(outfile)) return(0)
  
  res <- run_simulation(n, iter, sim_params)
  
  saveRDS(res, outfile)
  }, .options = furrr_options(seed = 4433L))

plan("default")
