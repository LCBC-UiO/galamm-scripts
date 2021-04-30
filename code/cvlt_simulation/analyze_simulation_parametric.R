library(tidyverse)
library(furrr)
library(latex2exp)
library(patchwork)
options(bitmapType = "cairo")
theme_set(theme_bw())
library(gamm4)

sim_params <- readRDS("data/cvlt_simulation/sim_params.rds")
res_dir <- "results/cvlt_simulation/model_fits/"

resfile <- "results/cvlt_simulation/param_simres.rds"

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
                
        betas_est <- coef(res$final_gamm$gam)[c(paste0("a", 2:5), "retest")]
        betas_se <- sqrt(diag(res$cov_beta))[names(betas_est)]
        names(betas_est) <- names(betas_se) <- 
          str_replace(names(betas_est), "a", "beta")
        
        lambdas_est <- res$lambda_est[-1]
        lambdas_se <- sqrt(diag(res$cov_lambda))
        names(lambdas_est) <- names(lambdas_se) <- 
          str_replace(names(lambdas_est), "a", "lambda")
        
        varcomps <- getME(res$final_gamm$mer, "theta")[c("timepoint:id.weight",
                                                         "id.weight")]^2
        
        estimates <- c(betas_est, lambdas_est, varcomps)
        std_errors <- c(betas_se, lambdas_se, varcomps * NA)
        
        tibble(
          param = names(estimates),
          estimate = estimates,
          std_error = std_errors,
          n = n
        )
      })
  })
  plan("default")
  saveRDS(simres, resfile)
}


# Parametric effects
param_df <- simres %>% 
  filter(n != "1000_30", !str_detect(n, "a")) %>% 
  filter(!param %in% c("timepoint:id.weight", "id.weight")) %>% 
  mutate(
    n = as.integer(n),
    true_val = case_when(
      str_detect(param, "beta[:digit:]") ~ 
        sim_params$beta[str_replace(param, "beta", "a")],
      param == "retest" ~ sim_params$beta[["retest"]],
      str_detect(param, "lambda[:digit:]") ~
        sim_params$lambda[str_replace(param, "lambda", "a")],
      TRUE ~ NA_real_
    ),
    pct_bias = (estimate - true_val) / true_val,
    covered = true_val > (estimate + qnorm(.025) * std_error) &
      true_val < (estimate + qnorm(.975) * std_error)
  ) %>% 
  group_by(param, n) %>% 
  summarise(
    pct_bias_est = mean(pct_bias),
    pct_bias_se = sd(pct_bias) / sqrt(n()),
    coverage = map2(sum(covered), n(),  ~ prop.test(.x, n = .y)),
    rmse = sqrt(mean((estimate - true_val)^2)),
    .groups = "drop"
  ) %>% 
  mutate(
    pct_bias_lower = pct_bias_est + qnorm(.025) * pct_bias_se,
    pct_bias_upper = pct_bias_est + qnorm(.975) * pct_bias_se,
    coverage_est = map_dbl(coverage, ~ .x$estimate),
    coverage_lower = map_dbl(coverage, ~ .x$conf.int[[1]]),
    coverage_upper = map_dbl(coverage, ~ .x$conf.int[[2]])
  ) %>% 
  select(-coverage) %>% 
  mutate(
    param_group = case_when(
      str_detect(param, "beta") ~ "Item effects",
      param == "retest" ~ "Retest effect",
      str_detect(param, "lambda") ~ "Factor loadings",
      TRUE ~ NA_character_
    ),
    param_group = factor(
      param_group, 
      levels = c("Item effects", "Retest effect", "Factor loadings"))
  )

alpha_level <- .5
plots1 <- param_df %>% 
  nest_by(param_group) %>% 
  pmap(function(param_group, data){
    pd <- position_dodge(.2)
    
    p <- data %>% 
      ggplot(aes(x = factor(n), y = pct_bias_est, group = param, color = param,
               ymin = pct_bias_lower, ymax = pct_bias_upper)) +
      geom_hline(yintercept = 0, color = "gray") +
      geom_line(position = pd, alpha = alpha_level) +
      geom_point(position = pd) +
      ylab("Bias") +
      xlab("Subjects")
    
    if(param_group == "Retest effect") {
      p <- p + geom_errorbar(width = .08, alpha = alpha_level)
    } else {
      p <- p + geom_errorbar(width = .3, position = pd, alpha = alpha_level)
    }
      
    p <- p + 
      scale_y_continuous(labels = scales::percent_format(accuracy = .1)) +
      labs(color = NULL) +
      ggthemes::scale_color_colorblind(
        labels = function(lb) {
          case_when(
            str_detect(lb, "beta[:digit:]") ~
              unname(TeX(paste0("$\\beta_{", str_extract(lb, "[:digit:]+"), "}"))),
            lb == "retest" ~ unname(TeX("$\\beta_{r}$")),
            str_detect(lb, "lambda[:digit:]") ~ 
              unname(TeX(paste0("$\\lambda_{", str_extract(lb, "[:digit:]+"), "}"))),
            TRUE ~ unname(TeX("NA_character_"))
          )}) +
      theme(
        panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8),
        legend.key.height = unit(.2, 'cm'),
        legend.position = c(.7, .85)
      )
  })
  

plots2 <- param_df %>% 
  nest_by(param_group) %>% 
  pmap(function(param_group, data){
    pd <- position_dodge(.3)
    
    p <- data %>% 
      ggplot(aes(x = factor(n), y = coverage_est, group = param, color = param,
                 ymin = coverage_lower, ymax = coverage_upper)) +
      geom_hline(yintercept = 0.95, color = "gray") +
      geom_line(position = pd, alpha = alpha_level) +
      geom_point(position = pd) +
      ylab("Coverage")
    
    if(param_group == "Retest effect") {
      p <- p + geom_errorbar(width = .08, alpha = alpha_level)
    } else {
      p <- p + geom_errorbar(width = .3, position = pd, alpha = alpha_level) 
    }
    
    p <- p + 
      scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
      theme(legend.position = "none") +
      ggthemes::scale_color_colorblind() +
      theme(
        panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line()
      ) + 
      xlab("Subjects")
    
  })

plots3 <- param_df %>% 
  nest_by(param_group) %>% 
  pmap(function(param_group, data){
    p <- data %>% 
      ggplot(aes(x = factor(n), y = rmse^2, group = param, color = param)) +
      geom_line(alpha = alpha_level) +
      geom_point() + 
      ylab("MSE") +
      xlab("Subject")
    
    p <- p + 
      theme(legend.position = "none") +
      ggthemes::scale_color_colorblind() +
      theme(
        panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line()
      )
    
  })


plt1 <- wrap_plots(plots1[[1]], plots2[[1]], plots3[[1]], ncol = 3) +
  plot_annotation(title = "Item effects")
ggsave("results/cvlt_simulation/figures/cvlt_sim_item_effects.pdf",
       plot = plt1, width = 16, height = 6, units = "cm")

plt2 <- wrap_plots(plots1[[2]], plots2[[2]], plots3[[2]], ncol = 3) +
  plot_annotation(title = "Retest effect")
ggsave("results/cvlt_simulation/figures/cvlt_sim_retest_effect.pdf",
       plot = plt2, width = 16, height = 6, units = "cm")

plt3 <- wrap_plots(plots1[[3]], plots2[[3]], plots3[[3]], ncol = 3) +
  plot_annotation(title = "Factor loadings")
ggsave("results/cvlt_simulation/figures/cvlt_sim_factor_loadings.pdf",
       plot = plt3, width = 16, height = 6, units = "cm")


for(f in list.files("results/cvlt_simulation/figures/")){
    to = file.path("/tsd/p23/data/durable/file-export/", f),
    overwrite = TRUE
  )
}

# root-n consistency
param_df %>% 
  select(param, n, rmse) %>% 
  mutate(mse = rmse^2) %>% 
  arrange(param, n) %>%
  group_by(param) %>% 
  summarise(
    factor1 = last(mse) / first(mse),
    factor2 = last(mse) / nth(mse, 2),
    .groups = "drop"
  ) %>% 
  pull(factor1) %>% 
  range()

# variance components
bfun <- function(a) (sum(sim_params$vc) - a * sim_params$vc[["psi2"]]) / sim_params$vc[["psi3"]]

# Table of variance components
vc_df <- simres %>% 
  filter(param %in% c("timepoint:id.weight", "id.weight"),
         !str_detect(n, "\\_30"), !str_detect(n, "a")) %>% 
  mutate(
    param = recode(param, `timepoint:id.weight` = "psi2",
                   `id.weight` = "psi3"),
    true_value = case_when(
      n %in% c("200", "500", "1000") ~ sim_params$vc[param],
      str_detect(n, "a0\\.5") & param == "psi2" ~ 0.5 * sim_params$vc[param],
      str_detect(n, "a0\\.5") & param == "psi3" ~ bfun(.5) * sim_params$vc[param],
      str_detect(n, "a2") & param == "psi2" ~ 2 * sim_params$vc[param],
      str_detect(n, "a2") & param == "psi3" ~ bfun(2) * sim_params$vc[param],
      TRUE ~ NA_real_
    ),
    param = factor(param)
  ) %>% 
  group_by(n, param, true_value) %>% 
  summarise(
    mean_estimate = mean(estimate),
    mean_bias = mean((estimate - true_value) / true_value),
    se_bias = sd((estimate - true_value) / true_value) / sqrt(n()),
    rmse = sqrt(mean((estimate - true_value)^2)),
    .groups = "drop"
  ) %>% 
  mutate(
    type = factor(case_when(
      !str_detect(n, "[:alpha:]") ~ "a",
      str_detect(n, "a0\\.5") ~ "b",
      str_detect(n, "a2\\.") ~ "c",
      TRUE ~ NA_character_
    )),
    n = factor(as.integer(str_extract(n, "[:digit:]+")))
  )

levels(vc_df$param) <- c(psi2 = TeX("$\\psi^{(2)}$"), psi3 = TeX("$\\psi^{(3)}$"))

pd <- position_dodge(.2)

ggplot(vc_df, aes(x = n, y = mean_bias, group = interaction(param, type),
             ymin = mean_bias + qnorm(.025) * se_bias, 
             ymax = mean_bias + qnorm(.975) * se_bias)) + 
  geom_line(position = pd) + 
  geom_point(position = pd) +
  geom_errorbar(position = pd, width = .3) +
  facet_wrap(vars(param), labeller = label_parsed, scales = "free_y") +
  geom_hline(yintercept = 0, color = "gray") +
  ylab("Bias") +
  scale_y_continuous(labels = scales::percent_format(accuracy = .1)) +
  theme(
    panel.grid = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(),
    strip.background = element_blank()
  )

ggsave("results/cvlt_simulation/figures/cvlt_sim_variance_components.pdf",
       width = 20, height = 8, units = "cm")


