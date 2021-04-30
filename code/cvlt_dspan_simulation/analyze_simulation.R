library(tidyverse)
library(gghalves)
library(furrr)
options(bitmapType = "cairo")
theme_set(theme_bw())
library(mgcv)

sim_params <- readRDS("data/cvlt_dspan_simulation/sim_params.rds")
res_dir <- "results/cvlt_dspan_simulation/model_fits/"
save_path <- "results/cvlt_dspan_simulation/simres.rds"

if(file.exists(save_path)){
  simres <- readRDS(save_path)
} else {
  plan(multisession)
  simres <- future_map(list.files(res_dir), function(f){
    mod <- readRDS(file.path(res_dir, f))
    iteration <- str_extract(f, "[:digit:]+")
    true_cor <- as.numeric(str_extract_all(f, "[:digit:]+")[[1]][[3]]) / 10
    
    setting <- if_else(str_detect(f, "full"), "full", "partial")
    timepoints <- if_else(str_detect(f, "original"), "original", 
                          paste0(str_extract(f, "two|three"), "ormore"))
    vc <- as_tibble(VarCorr(mod$mer)) %>% 
      filter(!is.na(var2)) %>% 
      pull(sdcor) %>% 
      tibble(
        parameter = c("corr2", "corr3"),
        estimate = .,
        true_value = c(true_cor, cov2cor(sim_params$Psi3)[[1,2]])
      )
    
    bhat <- tibble(
      parameter = names(sim_params$beta),
      estimate = coef(mod$gam)[names(sim_params$beta)],
      true_value = sim_params$beta,
      ci_lower = as.numeric(confint(mod$gam)[names(sim_params$beta), 1]),
      ci_upper = as.numeric(confint(mod$gam)[names(sim_params$beta), 2])
    )
    
    grid <- crossing(
      age = seq(from = sim_params$age_min,
                to = sim_params$age_max, length.out = 100),
      a2 = 0, a3 = 0, a4 = 0, a5 = 0, fwd = 0,
      retestcvlt = 0, retestdigitspan = 0, weight = 1,
      itemgroup = c("cvlt", "digitspan")
    ) %>% 
      mutate(
        a1 = as.numeric(itemgroup == "cvlt"),
        bwd = as.numeric(itemgroup == "digitspan"),
        age_z = (age - sim_params$age_mean) / sim_params$age_sd,
        true_value = case_when(
          itemgroup == "cvlt" ~ sim_params$smooth_age_z_cvlt(age_z) +
            sim_params$beta[["a1"]],
          itemgroup == "digitspan" ~ sim_params$smooth_age_z_digitspan(age_z) +
            sim_params$beta[["bwd"]]
        )
      )
    
    fit <- predict(mod$gam, newdata = grid, se.fit = TRUE)
    smooth <- grid %>% 
      mutate(
        setting = setting,
        timepoints = timepoints,
        iteration = iteration,
        pred = fit$fit,
        se = fit$se.fit,
        covered = (true_value > (pred + qnorm(.025) * se)) &
          (true_value < (pred + qnorm(.975) * se))
      )
    
    list(
      parametric = vc %>% 
        bind_rows(bhat) %>% 
        mutate(
          setting = setting,
          timepoints = timepoints,
          iteration = iteration,
        ),
      smooth = smooth)
  })
  saveRDS(simres, file = save_path)
  plan("default")
}

setting_factor <- function(s){
  factor(recode(
    s, 
    "full" = "Both tests at each timepoint",
    "partial" = "Original missingness"),
    levels = c("Original missingness", "Both tests at each timepoint"))
}

tp_factor <- function(tp){
  factor(recode(
    tp,
    "original" = "Original timepoints",
    "twoormore" = "At least two timepoints",
    "threeormore" = "At least three timepoints",
  ), levels = c("Original timepoints",
                "At least two timepoints",
                "At least three timepoints"))
}

item_factor <- function(ig){
  factor(
    recode(ig, "cvlt" = "Episodic memory", 
           "digitspan" = "Working memory"),
    levels = c("Episodic memory", "Working memory"))
}

pd <- position_dodge(.2)
map_dfr(simres, ~ .x$smooth) %>% 
  select(setting, timepoints, iteration, itemgroup, 
         pred, se, covered, true_value, age) %>% 
  group_by(setting, timepoints, itemgroup, iteration) %>% 
  summarise(
    coverage = mean(covered),
    .groups = "drop_last"
  ) %>% 
  summarise(
    coverage_mean = mean(coverage),
    coverage_se = sd(coverage) / sqrt(n()),
    .groups = "drop"
  ) %>% 
  mutate(
    coverage_low = coverage_mean + qnorm(.025) * coverage_se,
    coverage_high = coverage_mean + qnorm(.975) * coverage_se,
    setting = setting_factor(setting),
    timepoints = tp_factor(timepoints),
    itemgroup = item_factor(itemgroup)
  ) %>% 
  ggplot(aes(x = timepoints, y = coverage_mean,
             ymin = coverage_low, ymax = coverage_high,
             color = itemgroup)) + 
  geom_point(position = pd) + 
  geom_errorbar(width = .2, position = pd) + 
  facet_grid(rows = vars(setting)) +
  geom_hline(yintercept = .95, color = "gray") + 
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) + 
  theme(
    panel.grid = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(),
    strip.background = element_blank(),
    legend.position = c(.15, .6)
  ) +
  ggthemes::scale_color_colorblind() + 
  labs(color = NULL) +
  xlab("") + 
  ylab("Coverage")

ggsave("results/cvlt_dspan_simulation/figures/cvlt_dspan_sim_covr.pdf",
       width = 16, height = 12, units = "cm")


corr_simres <- map_dfr(simres, ~ .x$parametric) %>% 
  filter(parameter %in% c("corr2", "corr3")) %>% 
  mutate(
    setting = setting_factor(setting),
    timepoints = tp_factor(timepoints)
  )


create_corrplot <- function(dat, true_corr){
  ggplot(dat, aes(x = timepoints, y = estimate,
             group = interaction(setting, timepoints), color = setting)) + 
    geom_half_violin() + 
    geom_half_point(size = .2) +
    geom_hline(yintercept = true_corr, color = "gray") +
    xlab("") +
    ylab("Estimated correlation") +
    theme(
      panel.grid = element_blank(),
      panel.border = element_blank(),
      axis.line = element_line(),
      strip.background = element_blank(),
      legend.position = c(.75, .25)
    ) + 
    ggthemes::scale_color_colorblind() +
    labs(color = NULL)
}
  
corr_simres %>% 
  filter(parameter == "corr2") %>% 
  create_corrplot(true_corr = .5)

ggsave("results/cvlt_dspan_simulation/figures/cvlt_dspan_level2_corr.pdf",
       width = 12, height = 6, units = "cm")


corr_simres %>% 
  filter(parameter == "corr3") %>% 
  create_corrplot(true_corr = cov2cor(sim_params$Psi3)[[2,1]]) +
  theme(legend.position = "none")

ggsave("results/cvlt_dspan_simulation/figures/cvlt_dspan_level3_corr.pdf",
       width = 12, height = 6, units = "cm")



map_dfr(simres, ~ .x$smooth) %>% 
  group_by(setting, timepoints, iteration, itemgroup) %>% 
  summarise(rmse = sqrt(mean((pred - true_value)^2)), .groups = "drop") %>% 
  inner_join(
    map_dfr(simres, ~ .x$parametric) %>% 
      filter(parameter == "corr2") %>% 
      select(estimate, setting, timepoints, iteration),
    by = c("setting", "timepoints", "iteration")
  ) %>% 
  mutate(
    setting = setting_factor(setting),
    timepoints = tp_factor(timepoints),
    itemgroup = item_factor(itemgroup)
  ) %>% 
  ggplot(aes(x = estimate, y = rmse, color = itemgroup)) + 
  geom_point(size = .8) + 
  facet_grid(
    rows = vars(setting),
    cols = vars(timepoints)
  ) + 
  ggthemes::scale_color_colorblind() +
  xlab("Estimated level-2 correlation") + 
  ylab("Root-mean-square error of smooth function") +
  theme(
    panel.grid = element_blank(),
    #panel.border = element_blank(),
    axis.line = element_line(),
    strip.background = element_blank(),
    legend.position = c(.85, .9)
  ) +
  labs(color = NULL)

ggsave(
  "results/cvlt_dspan_simulation/figures/cvlt_dspan_sim_rmse.pdf",
  width = 16, height = 12, units = "cm"
)

