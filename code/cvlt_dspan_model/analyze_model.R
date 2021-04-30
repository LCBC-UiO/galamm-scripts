library(tidyverse)
library(patchwork)
library(latex2exp)
options(bitmapType = "cairo")
theme_set(theme_bw())
library(gamm4)
library(HDInterval)

mod <- readRDS("results/cvlt_dspan_model/res.rds")
dat <- readRDS("data/cvlt_dspan_model/simulated_data.rds")

dat %>% 
    group_by(id) %>% 
    filter(!any(itemgroup == "cvlt")) %>% 
    ungroup() %>% 
    summarise(n_distinct(id), n())

dat %>% 
    group_by(id) %>% 
    filter(!any(itemgroup == "digitspan")) %>% 
    ungroup() %>% 
    summarise(n_distinct(id), n())

dat %>% 
    group_by(id) %>% 
    filter(any(itemgroup == "digitspan"),
           any(itemgroup == "cvlt")) %>% 
    ungroup() %>% 
    summarise(n_distinct(id), n())


dat %>% summarise(n_distinct(id), n(), min(age), max(age))
dat %>% distinct(id, timepoint) %>% summarise(n_distinct(id), n())

dat %>% 
    filter(itemgroup == "digitspan") %>% 
    mutate(item = factor(item, levels = c("fwd", "bwd"))) %>% 
    ggplot(aes(x = age, y = successes, group = id)) +
    geom_point(size = .3) +
    geom_line(alpha = .2) +
    facet_grid(
        cols = vars(item), labeller = 
            as_labeller(function(x) if_else(x == "bwd", "Backward", "Forward"))) +
    xlab("Age") + ylab("Length of list") +
    theme(
        panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(),
        strip.background = element_blank()
    )

ggsave("results/cvlt_dspan_model/figures/dspan_sample_plot.pdf", 
       width = 12, height = 6, units = "cm")



grid <- crossing(
    a1 = 0, a2 = 0, a3 = 0, a4 = 0, a5 = 0,
    bwd = 0, fwd = 0, retestcvlt = 0, retestdigitspan = 0,
    weight = 1, 
    age = seq(from = min(dat$age), to = max(dat$age), length.out = 300),
    itemgroup = levels(dat$itemgroup)
    ) %>% 
    mutate(age_z = (age - mean(dat$age)) / sd(dat$age)) %>% 
    mutate(fit = as.numeric(predict(mod$final_gamm$gam, newdata = .)))

# Use posterior sampling to investigate the difference in peak
# age of episodic and working memory
nmc <- 10000
grid_epmem <- filter(grid, itemgroup == "cvlt", age < 60, age > 10)
Xp_epmem <- predict(
    mod$final_gamm$gam, newdata = grid_epmem, 
    type = "lpmatrix")


grid_wmem <- filter(grid, itemgroup == "digitspan", age < 60, age > 10)
Xp_wmem <- predict(
    mod$final_gamm$gam, newdata = grid_wmem,
    type = "lpmatrix"
)

bhat <- coef(mod$final_gamm$gam)
S <- mod$cov_beta
betas <- rmvn(nmc, bhat, S)

posterior <- tibble(
    epmem_age = grid_epmem$age[apply(betas %*% t(Xp_epmem), 1, which.max)],
    wmem_age = grid_wmem$age[apply(betas %*% t(Xp_wmem), 1, which.max)]
    ) %>% 
    mutate(diff = wmem_age - epmem_age)

hdi_df <- tibble(
    l = hdi(posterior$diff)[["lower"]],
    u = hdi(posterior$diff)[["upper"]])


p1 <- ggplot(grid, aes(x = age, y = fit, group = itemgroup,
                 color = itemgroup)) +
    geom_line() +
    ggthemes::scale_color_colorblind(
        labels = c("Episodic memory", "Working memory")
    ) +
    labs(color = NULL) +
    theme(
        legend.position = c(.5, .3),
        panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line()        
    ) +
    ylab("Level") +
    xlab("Age")

p2 <- ggplot(posterior, aes(x = diff)) + 
    geom_density() +
    geom_segment(data = hdi_df, aes(x = l, xend = u, y = 0, yend = 0),
                 size = 2) +
    theme(
        panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line()             
    ) +
    xlab("Difference in age at maximum") +
    ylab("Posterior density")

p1 + p2

ggsave("results/cvlt_dspan_model/figures/cvlt_dspan_figure.pdf",
       width = 16, height = 6, units = "cm")

