library(tidyverse)
theme_set(theme_bw())
theme_update(
  panel.grid = element_blank(),
  strip.background = element_blank()
)
library(patchwork)
library(ggrepel)

dat <- readRDS("data/cognition_model/dat.rds")
files <- list.files("results/cognition_model/simulation/summaries/", 
                    pattern = "theta0.rds", full.names = TRUE)

simres <- map(files, readRDS)
params <- readRDS("results/cognition_model/params.rds")

beta_df <- map_dfr(simres, function(x){
  as_tibble(x$beta_df) %>% 
    mutate(
      parameter = rownames(x$beta_df),
      true_value = params$beta_true,
      asym_se = params$beta_se
      )
}) %>% 
  mutate(
    bias = (beta_est - true_value) / true_value
  ) %>% 
  na.omit()

beta_agg <- beta_df  %>% 
  group_by(parameter, true_value, asym_se) %>% 
  summarise(
    mean_est = mean(beta_est),
    mean_bias = mean(bias),
    bootstrap_se = sd(beta_est),
    mean_se = mean(beta_se),
    rmse = sqrt(sum((true_value - beta_est)^2)),
    bootstrap_ci_lower = quantile(beta_est, probs = .025, names = FALSE),
    bootstrap_ci_upper = quantile(beta_est, probs = .975, names = FALSE),
    .groups = "drop"
  ) %>% 
  mutate(
    asymptotic_ci_lower = true_value + qnorm(.025) * asym_se,
    asymptotic_ci_upper = true_value + qnorm(.975) * asym_se
  )

lambda_df <- map_dfr(simres, function(x){
  as_tibble(x$lambda_df) %>% 
    mutate(
      parameter = rownames(x$lambda_df),
      true_value = params$lambda_true,
      asym_se = params$lambda_se
    )
}) %>% 
  mutate(
    bias = (lambda_est - true_value) / true_value
  )

lambda_agg <- lambda_df %>% 
  group_by(parameter, true_value, asym_se) %>% 
  summarise(
    mean_est = mean(lambda_est),
    mean_bias = mean(bias),
    bootstrap_se = sd(lambda_est),
    mean_se = mean(lambda_se),
    rmse = sqrt(sum((true_value - lambda_est)^2)),
    bootstrap_ci_lower = quantile(lambda_est, probs = .025, names = FALSE),
    bootstrap_ci_upper = quantile(lambda_est, probs = .975, names = FALSE),
    .groups = "drop"
  ) %>% 
  mutate(
    asymptotic_ci_lower = true_value + qnorm(.025) * asym_se,
    asymptotic_ci_upper = true_value + qnorm(.975) * asym_se
  )

agg_df <- bind_rows(
  "Regression coefficients" = beta_agg,
  "Factor loadings" = lambda_agg,
  .id = "type"
) %>% 
  mutate(
    parameter = recode(
      parameter,
      "stroop_time_1" = "Stroop 1",
      "stroop_time_2" = "Stroop 2",
      "stroop_time_3" = "Stroop 3",
      "stroop_time_4" = "Stroop 4",
      "teststroop_time_1" = "Stroop 1",
      "teststroop_time_2" = "Stroop 2",
      "teststroop_time_3" = "Stroop 3",
      "teststroop_time_4" = "Stroop 4"
      )
  )

p1 <- ggplot(agg_df, aes(x = mean_est, y = true_value)) + 
  geom_abline(slope = 1, intercept = 0, color = "gray") + 
  geom_point() + 
  xlab("Average estimate") + 
  ylab("True value") + 
  facet_wrap(vars(type), scales = "free")

p2 <- ggplot(agg_df, aes(x = mean_se, y = bootstrap_se)) +
  geom_abline(slope = 1, intercept = 0, color = "gray") + 
  geom_point() + 
  geom_text_repel(data = filter(agg_df, abs(mean_se - bootstrap_se) > .05), 
                  size = 2.5,
                  aes(label = parameter)) +
  xlab("Average SE") + 
  ylab("Bootstrap SE") + 
  facet_wrap(vars(type), scales = "free")

p1 / p2
ggsave("results/cognition_model/simulation/figures/cog_sim_parametric.pdf",
       plot = p1 / p2, width = 12, height = 12, units = "cm")


transf_df <- readRDS("results/cognition_model/transf_df.rds")

fff <- function(x){
  x %>% 
    transmute(
      parameter = parameter,
      estimate = round(true_value, 2),
      asymptotic_ci = paste0("(", round(asymptotic_ci_lower, 2), ", ", round(asymptotic_ci_upper, 2), ")"),
      bootstrap_estimate = round(mean_est, 2),
      bootstrap_ci = paste0("(", round(bootstrap_ci_lower, 2), ", ", round(bootstrap_ci_upper, 2), ")")
    )
}

beta_agg %>% 
  select(parameter, true_value, mean_est, contains("ci")) %>% 
  fff()

# Deal with Stroop conditions
beta_agg %>% 
  filter(parameter == "retest:itemgroup_reteststroop_12") %>% 
  select(parameter, true_value, mean_est, contains("ci")) %>% 
  mutate(
    across(c(mean_est, true_value, contains("ci")), ~ -.x * transf_df$scale[[1]])
  ) %>% 
  fff()

beta_agg %>% 
  filter(parameter == "retest:itemgroup_reteststroop_34") %>% 
  select(parameter, true_value, mean_est, contains("ci")) %>% 
  mutate(
    across(c(mean_est, true_value, contains("ci")), ~ -.x * transf_df$scale[[2]])
  ) %>% 
  fff()

lambda_agg %>% 
  select(parameter, true_value, mean_est, contains("ci")) %>% 
  transmute(
    parameter = parameter,
    estimate = round(true_value, 2),
    asymptotic_ci = paste0("(", round(asymptotic_ci_lower, 2), ", ", round(asymptotic_ci_upper, 2), ")"),
    bootstrap_estimate = round(mean_est, 2),
    bootstrap_ci = paste0("(", round(bootstrap_ci_lower, 2), ", ", round(bootstrap_ci_upper, 2), ")")
  )

true_plot_dat <- readRDS("results/cognition_model/plot_dat.rds")

smooth_df <- imap_dfr(simres, function(x, i){
  x$plot_dat %>% 
    mutate(
      iteration = i,
      true_fit = true_plot_dat$fit
    )
}) %>% 
  na.omit() %>% 
  mutate(
    age = age_z * sd(dat$age) + mean(dat$age)
  )

mean_fit_df <- smooth_df %>% 
  group_by(itemgroup, age, true_fit) %>% 
  summarise(mean_fit = mean(fit), .groups = "drop") %>% 
  pivot_longer(cols = c(mean_fit, true_fit)) %>% 
  mutate(
    name = recode(name, "true_fit" = "True", "mean_fit" = "Average")
  )

# This plot goes in supplement
p <-  
  ggplot(smooth_df) + 
  geom_line(aes(x = age, y = fit, group = iteration), color = "gray", size = .1) + 
  geom_line(data = mean_fit_df, aes(x = age, y = value, group = name, color = name)) +
  facet_wrap(vars(itemgroup), scales = "free_y", nrow = 1) +
  xlab("Age") +
  ggthemes::scale_color_colorblind() +
  theme(
    axis.title.y = element_blank(),
    legend.position = c(.15, .3)
    ) + 
  labs(color = NULL)

ggsave("results/cognition_model/simulation/figures/smooth_curves_repeated_sampling.pdf",
       plot = p, width = 16, height = 6, units = "cm")

covr_df <- smooth_df %>% 
  mutate(
    lower_ptw = fit + qnorm(.025) * se,
    upper_ptw = fit + qnorm(.975) * se,
    lower_sim = fit - crit * se,
    upper_sim = fit + crit * se
  )


p1 <- covr_df %>% 
  filter(itemgroup == "Executive function", age < 15) %>% 
  pivot_longer(cols = c(lower_sim, upper_sim)) %>% 
  mutate(
    name = recode(name, "lower_sim" = "lower", "upper_sim" = "upper")
  ) %>% 
  ggplot(aes(x = age, y = value, group = interaction(iteration, name), color = name)) + 
  geom_line() + 
  geom_line(aes(y = true_fit), color = "red", size = 1) +
  ggthemes::scale_color_colorblind() +
  theme(
    legend.position = c(.7, .3),
    axis.title.y = element_blank()
  ) + 
  labs(color = NULL) + 
  xlab("Age") 


ggsave("results/cognition_model/simulation/figures/smooth_coverage_executive_function.pdf",
       plot = p1, width = 8, height = 6, units = "cm")

smooth_df %>% 
  group_by(itemgroup) %>% 
  summarise(
    rmse = sqrt(mean((fit - true_fit)^2)),
    .groups = "drop"
  )

smooth_df %>% 
  group_by(itemgroup) %>% 
  summarise(
    round(diff(range(fit)), 2),
    .groups = "drop"
  )

covr_df %>% 
  group_by(itemgroup, iteration) %>% 
  summarise(
    simultaneous_coverage = all(true_fit > lower_sim & true_fit < upper_sim),
    pointwise_coverage = mean(true_fit > lower_ptw & true_fit < upper_ptw),
    .groups = "drop_last"
  ) %>% 
  summarise(
    simultaneous_est = mean(simultaneous_coverage),
    pointwise_est = mean(pointwise_coverage),
    .groups = "drop"
  )

edf <- map_dfr(simres, ~ tibble(edf = .x$edf))
actual_edf <- readRDS("results/cognition_model/edf.rds")


ggplot(edf, aes(x = edf)) + 
  geom_histogram(binwidth = .1) + 
  geom_vline(xintercept = actual_edf, linetype = "dashed") +
  xlab("Effective degrees of freedom")

ggsave("results/cognition_model/simulation/figures/edf.pdf",
       width = 8, height = 6, units = "cm")

edf %>% 
  summarise(mean(edf))


