library(tidyverse)
theme_set(theme_bw())
theme_update(
  panel.grid = element_blank(),
  strip.background = element_blank()
)
library(latex2exp)

files <- list.files("results/cognition_model/simulation/model_fits/",
                    pattern = "opt", full.names = TRUE)

mod <- readRDS("results/cognition_model/cognition_model.rds")
Lambdat <- mod$lmod$reTrms$Lambdat
Lambdat@x <- mod$opt$par[mod$theta_inds][mod$lmod$reTrms$Lind]
psi_3 <- diag((t(Lambdat[11116:11118, 11116:11118])) %*% Lambdat[11116:11118, 11116:11118] * mod$final_mod$phi[[1]])[2:3]
rm(mod)

simres <- map_dfr(files, function(f){
  iter <- as.integer(str_extract(f, "(?<=_)[:digit:]+(?=_)"))
  theta <- as.numeric(str_extract(f, "(?<=_)[0-9.]+(?=\\.)"))
  
  mod <- readRDS(f)
  tibble(
    iter = iter,
    true_value = theta,
    est_zero = mod$par[1:2] == 0,
    parameter = c("Executive function", "Working memory")
  )
})

names(psi_3) <- c("Working memory", "Executive function")
plot_dat <- simres %>% 
  group_by(
    parameter, true_value
  ) %>% 
  summarise(
    prop_nonzero = mean(!est_zero),
    .groups = "drop"
  ) %>% 
  mutate(
    psi_3 = psi_3[parameter],
    relative = (true_value^2) / (true_value^2 + psi_3)
  )

ggplot(plot_dat, aes(x = relative, y = prop_nonzero, group = parameter, color = parameter)) +
  geom_point() + 
  geom_line() + 
  ggthemes::scale_color_colorblind() + 
  xlab(TeX("$\\psi_{m}^{(2)} / (\\psi_{m}^{(2)} + \\psi_{m}^{(3)})$ ")) +
  ylab("Nonzero estimates") +
  theme(legend.position = c(.6, .4)) + 
  labs(color = NULL)

ggsave("results/cognition_model/figures/theta_simulation.pdf",
       width = 8, height = 6, units = "cm")

