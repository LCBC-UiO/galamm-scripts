library(tidyr)
library(furrr)
plan(multisession)

params <- crossing(
  loadings = c("parents_equal", "itemgroups_equal", "free_loadings"),
  interactions = c("yes", "no", "nono")
)

future_pwalk(params, function(loadings, interactions){
    source("code/ses_model/model_setup.R")
    dat <- readRDS("data/ses_model/simulated_data.rds")
    outfile <- if(interactions == "yes") {
      paste0("results/ses_model/model_fits/opt_", loadings, ".rds")
    } else if(interactions == "no") {
      paste0("results/ses_model/model_fits/opt_", loadings, "_no_intrx.rds")
    } else if(interactions == "nono") {
      paste0("results/ses_model/model_fits/opt_", loadings, "_nono_intrx.rds")
    }    
    if(file.exists(outfile)) return(0)
    
    lambda_init <- c(runif(6, min = .9, max = 1.5), .1, .1)
    names(lambda_init) <- lambda_names
    
    inds_list <- indsfun(loadings)
    fixed_inds <- inds_list$fixed_inds
    equal_inds <- inds_list$equal_inds
    if(interactions == "no") {
      zero_inds <- 8L
    } else if(interactions == "nono"){
      zero_inds <- 7L:8L
    } else {
      zero_inds <- integer()
    }
    
    lambda_init[fixed_inds] <- 1
    lambda_init <- constrain_loadings(lambda_init, equal_inds, zero_inds)
    
    x0 <- lambda_init[-c(fixed_inds, equal_inds[-1], zero_inds)]
    
    lambda_init[setdiff(seq_along(lambda_init), fixed_inds)] <- NA_real_
    opt <- optim(x0, function(x){
      ll <- lambda_init
      ll[names(x)] <- x
      ll <- constrain_loadings(ll, equal_inds, zero_inds)

      mod <- modcall(formula, random, dat, weights, ll)
      -as.numeric(logLik(mod$lme))
    }, method = "L-BFGS-B", control = list(trace = 4), hessian = TRUE)
    
    saveRDS(opt, outfile)
    
  }, .options = furrr_options(packages = "mgcv", seed = 1234L))

plan("default")
