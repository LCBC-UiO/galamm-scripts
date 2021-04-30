# Usage: Rscript code/ses_simulation/run_simulation.R "iter"
iter <- as.integer(commandArgs(trailingOnly = TRUE)[[1]])
library(tidyverse)
library(numDeriv)
library(furrr)
library(gamm4)
options(bitmapType = "cairo")

source("code/ses_simulation/model_setup.R")
sim_params <- readRDS("data/ses_simulation/sim_params.rds")

dimensions <- crossing(
  lambda_interaction = sim_params$lambda_interaction,
  n = sim_params$n
) %>% 
  mutate(
    lambda_order = ntile(lambda_interaction, length(unique(lambda_interaction)))
    )

plan(multisession)
sims <- future_pwalk(
  dimensions, function(lambda_interaction, n, lambda_order){

    outfile <- file.path(
      "results/ses_simulation/model_fits", 
      paste0("mod_n", n, "_iter", iter, "_lambda", lambda_order, ".rds"))
  
  if(file.exists(outfile)) return(0)
  
  dat <- tibble(
    id = seq_len(n),
    age_at_baseline = sim_params$age_at_baseline_dist(runif(n)),
    item = map(seq_len(n), function(i){
      v1 <- c("edu_father", "edu_mother", "edu_self", 
              "income_father", "income_mother", "income_self")
      v2 <- rep("hippocampus", sample(seq_along(sim_params$timepoints_dist), 1, 
                                      prob = sim_params$timepoints_dist))
        c(v1, v2)}),
      eta1 = rnorm(n, sd = sim_params$psi1),
      eta2 = rnorm(n, sd = sim_params$psi2)
    ) %>% 
    unnest(cols = item) %>% 
    filter(
      (age_at_baseline < 20 & !str_detect(item, "self")) |
        (age_at_baseline >= 20 & str_detect(item, "self")) |
        item == "hippocampus"
    ) %>% 
    mutate(time_interval = sim_params$interval_dist(runif(nrow(.)))) %>% 
    group_by(id, item) %>% 
    mutate(
      age = age_at_baseline + cumsum(c(0, time_interval[-1]))
    ) %>% 
    ungroup() %>% 
    filter(between(age, sim_params$age_min, sim_params$age_max)) %>% 
    select(-age_at_baseline, -time_interval) %>% 
    mutate(
      residual = case_when(
        str_detect(item, "edu") ~ 
          rnorm(nrow(.), sd = sim_params$sigma_edu),
        str_detect(item, "income") ~ 
          rnorm(nrow(.), sd = sim_params$sigma_income),
        item == "hippocampus" ~
          rnorm(nrow(.), sd = sim_params$sigma_hippocampus)
      ),
      age_z = (age - sim_params$age_mean) / sim_params$age_sd,
      smooth_age = (item == "hippocampus") * sim_params$smooth_age(age_z),
      intercept = if_else(item == "hippocampus", 0, sim_params$beta[paste0("item", item)]),
      loading = sim_params$lambda[item] * eta1 + 
        (item == "hippocampus") * (!!lambda_interaction * age_z * eta1 + eta2),
      value_z = smooth_age + intercept + loading + residual,
      itemgroup = str_extract(item, "edu|income|hippocampus")
    ) %>% 
      bind_cols(
        as_tibble(model.matrix(~ 0 + item + itemgroup, data = .))
      )
  
  # Fit model with and without interaction
  
  models <- map(c("yes", "no"), function(interactions){
    
    lambda_init <- c(runif(3, min = .9, max = 1.1), runif(3, min = .4, max = .5),
                     .1, .05)
    names(lambda_init) <- lambda_names
    
    inds_list <- indsfun("itemgroups_equal")
    fixed_inds <- inds_list$fixed_inds
    equal_inds <- inds_list$equal_inds
    if(interactions == "no") {
      zero_inds <- 8L
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
      
    }, method = "L-BFGS-B", control = list(trace = 4, maxit = 5), hessian = TRUE)
    
    
    
    
    lambda_est <- numeric(length(lambda_names))
    names(lambda_est) <- lambda_names
    lambda_est[fixed_inds] <- 1
    lambda_est[names(opt$par)] <- opt$par
    lambda_est <- constrain_loadings(lambda_est, equal_inds, zero_inds)
    
    final_gamm <- modcall(formula, random, dat, weights, lambda_est)
    
    beta_hat <- coef(final_gamm$gam)
    cov_lambda <- solve(opt$hessian)
    cov_beta_naive <- vcov(final_gamm$gam)
    
    beta_deriv <- jacobian(function(x){
      ll <- lambda_est
      ll[names(x)] <- x
      ll <- constrain_loadings(ll, equal_inds, zero_inds)
      
      mod <- modcall(formula, random, dat, weights, ll)
      coef(mod$gam)
    }, x = opt$par, method = "simple")
    
    cov_beta <- cov_beta_naive + beta_deriv %*% cov_lambda %*% t(beta_deriv)
    
    res <- list(
      final_gamm = final_gamm,
      beta_deriv = beta_deriv,
      cov_beta_naive = cov_beta_naive,
      cov_beta = cov_beta,
      cov_lambda = cov_lambda,
      lambda_est = lambda_est
    )
  })
  
  
  saveRDS(models, outfile)
  }, .options = furrr_options(seed = iter))

plan("default")
