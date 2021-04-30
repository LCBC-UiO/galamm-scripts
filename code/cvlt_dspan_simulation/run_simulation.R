# Simulations are run with factor loadings fixed to their true level,
# as the goal is to study the correlation
library(tidyverse)
options(bitmapType = "cairo")
library(mvtnorm)
library(numDeriv)
library(furrr)
library(gamm4)

sim_params <- readRDS("data/cvlt_dspan_simulation/sim_params.rds")

# 
dimensions <- crossing(
  iter = seq_len(100),
  rho = .5,
  timepoints = c("original", "twoormore", "threeormore"),
  items = c("partial", "full")
)

formula <- cbind(successes, failures) ~ 0 + a1 + a2 + a3 + a4 + a5 + 
  bwd + fwd + retestcvlt + retestdigitspan +
  t2(age_z, itemgroup, by = weight, bs = c("cr", "re"), k = 20, full = TRUE)

random <- ~ (0 + weight:itemgroup | id / timepoint) 

plan(multisession)
sims <- future_pwalk(dimensions, function(iter, rho, timepoints, items){
  
  outfile <- file.path(
    "results/cvlt_dspan_simulation/model_fits", 
    paste0("mod_items_", items, "_timepoints", timepoints, "_iter", 
           iter, "_rho", rho, ".rds"))
  
  if(file.exists(outfile)) return(0)
  
  if(items == "partial"){
    d0 <- sim_params$data_skeleton
  } else if(items == "full"){
    d0 <- sim_params$data_skeleton %>%
      distinct(id, timepoint) %>%
      mutate(item = list(c(paste0("a", 1:5), c("bwd", "fwd")))) %>%
      unnest(cols = item) %>%
      mutate(
        itemgroup = factor(if_else(item %in% paste0("a", 1:5), "cvlt", "digitspan"))
        ) %>% 
      select(id, item, itemgroup, timepoint)
  }
  
  if(timepoints == "original"){
    d0 <- d0
  } else if(timepoints == "twoormore"){
    d0 <- d0 %>% 
      nest_by(id, .keep = TRUE) %>% 
      pmap_dfr(function(id, data){
        if(max(data$timepoint) == 1){
          data %>% 
            bind_rows(
              data %>% 
                mutate(timepoint = 2)
            )
        } else {
          data
        }
      })
  } else if(timepoints == "threeormore"){
  d0 <- d0 %>% 
    nest_by(id, .keep = TRUE) %>% 
    pmap_dfr(function(id, data){
      if(max(data$timepoint) <= 2){
        data %>% 
          bind_rows(
            data %>% 
              mutate(timepoint = 2),
            data %>% 
              mutate(timepoint = 3)
          )
      } else {
        data
      }
    })
}  
  
  zeta3 <- rmvnorm(length(unique(d0$id)), sigma = sim_params$Psi3)
  
  zeta2 <- rmvnorm(length(unique(paste0(d0$id, d0$timepoint))),
                   sigma = matrix(c(
                     sim_params$psi12, 
                     rep(rho * sqrt(sim_params$psi12 * sim_params$psi22), 2),
                     sim_params$psi22
                   ), nrow = 2))
  
  bl_age <- sim_params$age_at_baseline_dist(runif(length(unique(d0$id))))
  dt <- sim_params$interval_dist(runif(length(unique(paste0(d0$id, d0$timepoint)))))
  
  dat <- d0 %>% 
    mutate(
      timepoint_number = as.integer(factor(
        paste0(id, timepoint), levels = unique(paste0(id, timepoint)))),
      age_at_baseline = bl_age[id],
      zeta3 = if_else(itemgroup == "cvlt", zeta3[id, 1], zeta3[id, 2]),
      zeta2 = if_else(itemgroup == "cvlt", zeta2[timepoint_number, 1], 
                      zeta2[timepoint_number, 2]),
      time_interval = dt[timepoint_number]
    ) %>% 
    group_by(id) %>% 
    mutate(
      age = age_at_baseline + cumsum(c(0, unique(time_interval)))[timepoint]
    ) %>% 
    ungroup() %>% 
    filter(between(age, sim_params$age_min, sim_params$age_max)) %>% 
    select(-age_at_baseline, -time_interval, -timepoint_number) %>% 
    group_by(id, itemgroup) %>% 
    mutate(retest_group = as.integer(dense_rank(age) > 1)) %>% 
    ungroup() %>% 
    mutate(
      age_z = (age - sim_params$age_mean) / sim_params$age_sd,
      item_effect = sim_params$beta[item],
      retest_effect = sim_params$beta[paste0("retest", itemgroup)] * 
        retest_group,
      weight = sim_params$lambda[item],
      smooth_age = if_else(itemgroup == "cvlt",
                           sim_params$smooth_age_z_cvlt(age_z),
                           sim_params$smooth_age_z_digitspan(age_z)),
      linear_predictor = item_effect + retest_effect + 
        weight * smooth_age + weight * zeta2 + weight * zeta3,
      successes = rbinom(nrow(.), size = sim_params$trials, 
                         prob = plogis(linear_predictor)), 
      failures = sim_params$trials - successes
    ) %>% 
    bind_cols(as.data.frame(model.matrix(
      ~ 0 + item + retest_group:itemgroup, data = .))) %>% 
    rename_with(~ str_remove(., "^item"), matches("itema[1-5]")) %>% 
    rename_with(~ str_remove(., "^item"), matches("item.wd")) %>% 
    rename_with(~ str_remove(., "\\_group\\:itemgroup"), 
                matches("retest\\_group\\:itemgroup")) %>% 
    select(id, item, itemgroup, timepoint, age_z, successes, failures, 
           a1, a2, a3, a4, a5, bwd, fwd, retestcvlt, retestdigitspan)
  
  # Set lambda to true value
  lambda_est <- sim_params$lambda
  dat$weight <- lambda_est[dat$item]
  
  mod <- gamm4(
    formula = formula, random = random,
    family = binomial(), data = dat, verbose = 2L)
  
  saveRDS(mod, outfile)
  }, .options = furrr_options(seed = 4433L))

plan("default")
