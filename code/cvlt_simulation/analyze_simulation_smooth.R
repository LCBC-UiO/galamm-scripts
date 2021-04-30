library(tidyverse)
library(furrr)
library(latex2exp)
library(patchwork)
options(bitmapType = "cairo")
theme_set(theme_bw())
library(gamm4)

sim_params <- readRDS("data/cvlt_simulation/sim_params.rds")
res_dir <- "results/cvlt_simulation/model_fits/"

grid <- tibble(
  a2 = 0, a3 = 0, a4 = 0, a5 = 0, retest = 0, weight = 1, 
  age = seq(from = sim_params$age_min, to = sim_params$age_max, 
            length.out = 100) 
  ) %>% 
  mutate(
    age_z = (age - sim_params$age_mean) / sim_params$age_sd,
    f_true = sim_params$smooth_age_z(age_z)
    )

nmc <- 10000

resfile <- "results/cvlt_simulation/smooth_simres.rds"
if(file.exists(resfile)){
  simres <- readRDS(resfile)
} else {
  plan(multisession)
  simres <- future_map_dfr(c(200, 500, 1000), function(n) {
    list.files(res_dir, pattern = paste0("n", n)) %>% 
      map_dfr(function(f) {
        res <- readRDS(file.path(res_dir, f))
        if(str_detect(f, "dimension")) {
          n <- paste0(n, "_30")
        } else if(str_detect(f, "\\_a")){
          n <- paste0(n, "_", str_extract(f, "a[:digit:]+\\.[:digit:]?"))
        } else {
          n <- as.character(n)
        }
        
        Xp <- predict(res$final_gamm$gam, newdata = grid, type = "lpmatrix")
        bhat <- coef(res$final_gamm$gam)
        S <- res$cov_beta
        
        f_est <- as.numeric(Xp %*% bhat)
        f_se <- as.numeric(sqrt(rowSums((Xp %*% S) * Xp)))
        
        betas <- rmvn(nmc, bhat, S)
        fits <- Xp %*% t(betas)
        simDev <- fits - f_est
        absDev <- apply(simDev, 2L, function(x) abs(x / f_se))
        masd <- apply(absDev, 2L, max)
        crit <- quantile(masd, prob = .95, type = 8)
        
        grid %>% 
          mutate(
            f_est = f_est,
            f_lower_pt = f_est + qnorm(.025) * f_se,
            f_upper_pt = f_est + qnorm(.975) * f_se,
            f_lower_sim = f_est - crit * f_se,
            f_upper_sim = f_est + crit * f_se,
            n = n
          )
      }, .id = "iteration")
  }, .options = furrr_options(seed = 2233L))
  plan("default")
  saveRDS(simres, file = resfile)
}

alpha_level <- .5
pd <- position_dodge(.2)
plot_df1 <- simres %>% 
  group_by(n, iteration) %>% 
  summarise(
    coverage_cross = mean(f_true >= f_lower_pt & f_true <= f_upper_pt),
    covered_whole = all(f_true >= f_lower_sim & f_true <= f_upper_sim),
    .groups = "drop_last"
  ) %>% 
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
    coverage_sim_upper = map_dbl(coverage_sim, ~ .x$conf.int[[2]])
  ) %>% 
  select(-coverage_sim, -coverage_cross_se) %>% 
  pivot_longer(cols = starts_with("coverage"), 
               names_to = c("type", ".value"),
               names_pattern = c("coverage_(.*)_(.*)")) %>% 
  mutate(
    type = if_else(n == "1000_30", "sim_30", 
                   if_else(str_detect(n, "a"), "a", type)),
    n = as.integer(str_extract(n, "^[:digit:]+"))
  )

p1 <- plot_df1 %>% 
  filter(type != "sim_30", type != "a") %>% 
  ggplot(aes(x = factor(n), y = est, group = type,
             color = type, ymin = lower, ymax = upper)) +
  geom_line(position = pd, alpha = alpha_level) + 
  geom_point(position = pd) +
  geom_errorbar(width = .1, position = pd, alpha = alpha_level) +
  geom_hline(aes(yintercept = .95), color = "gray") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  theme(
    panel.grid = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(),
    legend.position = c(.4, .15),
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 8),
    legend.key.height = unit(.2, 'cm')
  ) +
  ggthemes::scale_color_colorblind(
    labels = c("Pointwise", "Simultaneous")
  ) +
  labs(color = NULL) +
  xlab("Subjects") +
  ylab("Coverage")


p2 <- simres %>% 
  filter(n != "1000_30", !str_detect(n, "a")) %>% 
  mutate(n = as.integer(n)) %>% 
  mutate(covered = f_true >= f_lower_pt & f_true <= f_upper_pt) %>% 
  group_by(age, n) %>% 
  summarise(
    coverage_cross = map2(sum(covered), n(), ~ prop.test(.x, n = .y)),
    .groups = "drop"
  ) %>% 
  mutate(
    coverage_cross_est = map_dbl(coverage_cross, ~ .x$estimate),
    coverage_cross_lower = map_dbl(coverage_cross, ~ .x$conf.int[[1]]),
    coverage_cross_upper = map_dbl(coverage_cross, ~ .x$conf.int[[2]])
  ) %>% 
  ggplot(aes(x = age, y = coverage_cross_est, group = factor(n),
             ymin = coverage_cross_lower,
             ymax = coverage_cross_upper)) +
  geom_hline(yintercept = .95, alpha = alpha_level) +
  geom_line(aes(color = factor(n))) +
  geom_ribbon(aes(fill = factor(n)), alpha = .2) +
  ggthemes::scale_color_colorblind(
    labels = paste("n =", c(200, 500, 1000))
  ) +
  ggthemes::scale_fill_colorblind(
    labels = paste("n =", c(200, 500, 1000))
  ) +
  theme(
    panel.grid = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(),
    legend.position = c(.5, .2),
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 8),
    legend.key.height = unit(.2, 'cm')
  ) +
  labs(color = NULL, fill = NULL) +
  ylab("Coverage") +
  xlab("Age") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1))

mse_df <- simres %>% 
  filter(n != "1000_30", !str_detect(n, "a")) %>% 
  mutate(n = as.integer(n)) %>% 
  group_by(n, age) %>% 
  summarise(
    squared_bias = (mean(f_est - f_true))^2,
    variance = var(f_est - f_true) * (n() - 1) / n(),
    mse = squared_bias + variance,
    .groups = "drop_last"
  ) %>% 
  summarise(
    across(c(squared_bias, variance, mse), mean), .groups = "drop"
  ) %>% 
  pivot_longer(cols = c("squared_bias", "variance", "mse")) %>% 
  mutate(name = factor(name, levels = c("mse", "squared_bias", "variance")))

p3 <- ggplot(mse_df, aes(x = factor(n), y = value, group = name, 
             color = name)) +
  geom_line(alpha = .5) +
  geom_point() +
  theme(
    panel.grid = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(),
    legend.position = c(.75, .85),
    legend.text.align = 0,
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 8),
    legend.key.height = unit(.2, 'cm')
  ) +
  labs(color = NULL) +
  ggthemes::scale_color_colorblind(
    labels = unname(TeX(c("MSE", "$Bias^{2}$", "$Var.$")))
  ) +
  xlab("Subjects") +
  ylab("Value")

plt <- p1 + p2 + p3

ggsave("results/cvlt_simulation/figures/cvlt_sim_smooth.pdf",
       plot = plt, width = 18, height = 6, units = "cm")



# Find percentage bias
mse_df %>% 
  pivot_wider(names_from = name, values_from = value) %>% 
  mutate(
    pct_bias = squared_bias / mse
  )

mse_df %>% 
  filter(name == "mse") %>% 
  arrange(n) %>% 
  mutate(factor1 = last(value) / first(value))

# Read off simultaneous coverage with 30 basis functions
plot_df1 %>% 
  filter(type == "sim_30")
