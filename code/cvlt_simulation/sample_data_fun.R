
sample_data <- function(n, sim_params){
  with(
    sim_params,
    tibble(
      id = seq_len(n),
      age_at_baseline = age_at_baseline_dist((runif(n))),
      timepoints = sample(seq_along(timepoints_dist), size = n, 
                          replace = TRUE, prob = timepoints_dist),
      item = list(paste0("a", 1:5)),
      zeta3 = rnorm(n, sd = sqrt(vc["psi3"]))
    ) %>% 
      uncount(timepoints, .id = "timepoint") %>% 
      mutate(
        time_interval = interval_dist(runif(nrow(.))),
        zeta2 = rnorm(nrow(.), sd = sqrt(vc["psi2"]))
      ) %>% 
      group_by(id) %>% 
      mutate(
        age = age_at_baseline + cumsum(c(0, time_interval[-1]))
      ) %>% 
      ungroup() %>% 
      filter(between(age, age_min, age_max)) %>% 
      unnest(cols = item) %>% 
      mutate(retest = as.numeric(timepoint > 1)) %>% 
      select(-age_at_baseline, -time_interval) %>% 
      mutate(
        age_z = (age - age_mean) / age_sd,
        
        item_effect = beta[item],
        retest_effect = retest * beta["retest"],
        weight = lambda[item],
        h = smooth_age_z(age_z),
        linear_predictor = item_effect + retest_effect + 
          weight * (h + zeta2 + zeta3),
        successes = rbinom(nrow(.), size = trials, 
                           prob = plogis(linear_predictor)), 
        failures = trials - successes
      ) %>% 
      select(id, age, age_z, timepoint, retest, item, successes, failures) %>% 
      bind_cols(as.data.frame(model.matrix(~ 0 + item, data = .))) %>% 
      rename_with(~ str_remove(., "item"), matches("itema[1-5]"))
  )
}
