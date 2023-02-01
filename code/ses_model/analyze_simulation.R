library(tidyverse)
library(furrr)
library(patchwork)
library(latex2exp)
library(gghalves)
options(bitmapType = "cairo")
theme_set(theme_bw())
library(mgcv)

sim_params <- readRDS("data/ses_model/sim_params.rds")
files <- list.files("results/ses_model/simulation/summaries/",
                    pattern = ".rds", full.names = TRUE)

simres <- map(files, function(f){
    res <- readRDS(f)
    tryCatch({
      list(
        comp_df = map_df(res, ~ .x$comp_df),
        lambda_df = map_df(res, ~ .x$lambda_df),
        smooth_df = map2_df(res, sim_params$lambda_interaction, function(r, l){
          r$smooth_df %>% 
            mutate(true_interaction = l)
        }),
        vc_df = map_dfr(res, ~ .x$vc_df)
      )
    }, 
    error = function(e) list())
  })


comp_df <- map_dfr(simres, ~ .x$comp_df) %>% 
  nest_by(true_interaction) %>% 
  pmap_dfr(function(true_interaction, data){
    list(
      power = prop.test(sum(data$sig), length(data$sig)),
      aic = prop.test(sum(data$chosen), length(data$chosen))
    ) %>% 
      imap_dfr(~ tibble(name = .y,
                        estimate = .x$estimate,
                        lower = .x$conf.int[[1]],
                        upper = .x$conf.int[[2]])) %>% 
      mutate(true_interaction = !!true_interaction)
  })

pd <- position_dodge(.1)
p <- ggplot(comp_df, aes(x = factor(true_interaction), y = estimate, 
                         ymin = lower, ymax = upper, 
                         group = name, color = name)) + 
  geom_point(position = pd) + 
  geom_line(position = pd) + 
  geom_errorbar(position = pd, width = .3) + 
  ggthemes::scale_color_colorblind(labels = c("AIC", "p<.05")) +
  theme(
    panel.grid = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(),
    legend.position = c(.7, .3)
  ) +
  ylab("Probability") + 
  xlab(TeX("True interaction, $\\lambda_{8}$")) +
  labs(color = NULL) + 
  geom_hline(yintercept = .05, color = "gray")

ggsave(filename = "results/ses_model/simulation/figures/power_curve.pdf",
       plot = p, width = 8, height = 6, units = "cm")

lambda_df <- map_dfr(simres, ~ .x$lambda_df) %>% 
  filter(parameter != "interaction") %>% 
  mutate(
    covered = lower < true_lambda & upper > true_lambda
  ) %>% 
  nest_by(parameter, true_interaction, has_interaction, true_lambda) %>% 
  pmap_dfr(function(parameter, true_interaction, 
                    has_interaction, true_lambda, data){
    cvr <- prop.test(sum(data$covered), length(data$covered))
    tibble(
      parameter = parameter,
      true_interaction = true_interaction,
      has_interaction = has_interaction,
      true_lambda = true_lambda,
      coverage_est = cvr$estimate,
      coverage_lower = cvr$conf.int[[1]],
      coverage_upper = cvr$conf.int[[2]],
      est = mean(data$estimate),
      est_lower = est + qnorm(.025) * sd(data$estimate) / sqrt(nrow(.)),
      est_upper = est + qnorm(.975) * sd(data$estimate) / sqrt(nrow(.))
    )
  }) %>%
  mutate(parameter = fct_rev(factor(parameter)), 
         has_interaction = factor(has_interaction))

levels(lambda_df$parameter) <- c(
  income_father = c(income_father = TeX("Income, $\\lambda_{4}= \\lambda_{5}= \\lambda_{6}$"),
                    hippocampus = TeX("Hippocampus, $\\lambda_{7}$"))
)
levels(lambda_df$has_interaction) <- c(no = "Model without interaction",
                                       yes = "Model with interaction")

pd <- position_dodge(.2)

p2 <- ggplot(lambda_df, aes(x = factor(true_interaction), group = has_interaction,
                            y = est, ymin = est_lower,
                            ymax = est_upper, color = has_interaction)) + 
  geom_point(position = pd) + 
  geom_line(position = pd) + 
  geom_errorbar(width = .3, position = pd) +
  geom_hline(aes(yintercept = true_lambda), color = "gray") +
  facet_wrap(vars(parameter), scales = "free_y", labeller = label_parsed) +
  ggthemes::scale_color_colorblind() +
  theme(
    panel.grid = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(),
    strip.background = element_blank(),
    legend.position = c(.3, .2)
  ) +
  xlab(TeX("Interaction, $\\lambda_{8}$")) +
  ylab("Estimate") +
  labs(color = NULL)

ggsave(
  "results/ses_model/simulation/figures/lambda_estimates.pdf",
  plot = p2, width = 16, height = 8, units = "cm"
)

p3 <- map_dfr(simres, ~ .x$lambda_df) %>% 
  filter(parameter == "interaction") %>% 
  ggplot(aes(x = true_interaction, y = estimate)) + 
  geom_half_point(aes(color = factor(true_interaction)), width = .01, size = .3) +
  geom_half_violin(aes(color = factor(true_interaction))) +
  geom_point(aes(y = true_interaction), alpha = .3, size = .4) +
  geom_abline(slope = 1, intercept = 0, alpha = .3) +
  theme(
    panel.grid = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(),
    strip.background = element_blank(),
    legend.position = "none"
  ) +
  scale_x_continuous(breaks = c(0, 0.024, 0.048, 0.072, 0.096, 0.12)) +
  xlab(TeX("True interaction, $\\lambda_{8}$")) +
  ylab("Estimate") +
  ggthemes::scale_color_colorblind()

ggsave(
  "results/ses_model/simulation/figures/lambda_interaction.pdf",
  plot = p3, width = 8, height = 6, units = "cm"
)

# Find across-the-function coverage of smooth term
covr_df <- map_dfr(simres, ~ .x$smooth_df, .id = "iteration") %>% 
  group_by(true_interaction, ses_sd, iteration) %>% 
  summarise(
    coverage_cross = mean(true_fit > fit_lower_pt & true_fit < fit_upper_pt), 
    covered_whole = all(true_fit > fit_lower_sim & true_fit < fit_upper_sim),
    .groups = "drop_last") %>% 
  summarise(
    coverage_cross_est = mean(coverage_cross),
    coverage_cross_se = sd(coverage_cross) / sqrt(n()),
    coverage_cross_lower = coverage_cross_est + qnorm(.025) * coverage_cross_se,
    coverage_cross_upper = coverage_cross_est + qnorm(.975) * coverage_cross_se,
    coverage_sim = map2(sum(covered_whole), n(), ~ prop.test(.x, n = .y)),
    .groups = "drop"
  ) %>% 
  mutate(
    coverage_sim_est = map_dbl(coverage_sim, ~ .x$estimate),
    coverage_sim_lower = map_dbl(coverage_sim, ~ .x$conf.int[[1]]),
    coverage_sim_upper = map_dbl(coverage_sim, ~ .x$conf.int[[2]]),
    ses_sd = factor(if_else(ses_sd == 0, "Mean", sprintf("Mean %+d SD", ses_sd)))
  )

levels(covr_df$ses_sd) <- levels(covr_df$ses_sd)[c(3, 2, 1, 4, 5)]

pd <- position_dodge(.4)
p4 <- ggplot(covr_df, aes(x = factor(true_interaction), y = coverage_cross_est, 
                          ymin = coverage_cross_lower, ymax = coverage_cross_upper,
                          group = ses_sd, color = ses_sd)) + 
  geom_hline(yintercept = .95, color = "gray") +
  geom_point(position = pd) +
  geom_errorbar(width = .5, position = pd) +
  ggthemes::scale_color_colorblind() +
  theme(
    panel.grid = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(),
    strip.background = element_blank(),
    legend.position = c(.3, .9),
    legend.text = element_text(size = 6),
    legend.title = element_text(size = 6),
    legend.key.height = unit(.2, 'cm'),
    axis.title = element_text(size = 8)
  ) +
  xlab(TeX("True interaction, $\\lambda_{8}$")) +
  ylab("Across-the-function coverage") +
  labs(color = NULL) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1))

p5 <- ggplot(covr_df, aes(x = factor(true_interaction), y = coverage_sim_est, 
                          ymin = coverage_sim_lower, ymax = coverage_sim_upper,
                          group = ses_sd, color = ses_sd)) + 
  geom_hline(yintercept = .95, color = "gray") +
  geom_point(position = pd) +
  geom_errorbar(width = .5, position = pd) +
  ggthemes::scale_color_colorblind() +
  theme(
    panel.grid = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(),
    axis.title = element_text(size = 8),
    strip.background = element_blank(),
    legend.position = "none"
  ) +
  xlab(TeX("True interaction, $\\lambda_{8}$")) +
  ylab("Simultaneous coverage") +
  labs(color = NULL) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1))

ggsave(
  "results/ses_model/simulation/figures/ses_sim_smooth_coverage.pdf",
  plot = p4 + p5, width = 16, height = 6, units = "cm")

