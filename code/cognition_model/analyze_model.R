library(tidyverse)
library(patchwork)
library(galamm)
library(mgcv)
options(bitmapType = "cairo")
theme_set(theme_bw())
theme_update(panel.grid = element_blank(),
             strip.background = element_blank())
library(Matrix)

mod <- readRDS("results/cognition_model/cognition_model.rds")
dat <- readRDS("data/cognition_model/dat.rds")
S <- solve(-mod$final_mod$hessian)

# Update matrices with final parameters
Lambdat <- mod$lmod$reTrms$Lambdat
Lambdat@x <- mod$opt$par[mod$theta_inds][mod$lmod$reTrms$Lind]
Zt <- mod$lmod$reTrms$Zt
Zt@x <- Zt@x * c(1, mod$opt$par[mod$lambda_inds])[mod$lambda_mapping_Zt + 2L]

# Plot the splines
mod$lmod$reTrms$cnms
mod$lmod$reTrms$Gp

span_stroop <- seq(from = mod$lmod$reTrms$Gp[[5]] + 1L, to = mod$lmod$reTrms$Gp[[6]])
span_digitspan <- seq(from = mod$lmod$reTrms$Gp[[6]] + 1L, to = mod$lmod$reTrms$Gp[[7]])
span_cvlt <- seq(from = mod$lmod$reTrms$Gp[[7]] + 1L, to = mod$lmod$reTrms$Gp[[8]])

grid <- crossing(
  age_z = seq(from = min(mod$mdat$age_z), to = max(mod$mdat$age_z), length.out = 1000),
  itemgroup = factor(levels(mod$mdat$itemgroup))
) %>% 
  mutate(age = age_z * sd(dat$age) + mean(dat$age))

Xp <- do.call(cbind, map(rev(mod$sm), ~ PredictMat(.x, data = grid)))
transD <- c(mod$re$stroop$trans.D, mod$re$digitspan$trans.D, 
            mod$re$cvlt$trans.D)
transU <- matrix(0, nrow = length(transD), ncol = length(transD))
transU[1:9, 1:9] <- mod$re$stroop$trans.U
transU[10:18, 10:18] <- mod$re$digitspan$trans.U
transU[19:27, 19:27] <- mod$re$cvlt$trans.U

b <- c(as.numeric(Lambdat[span_stroop, span_stroop] %*% mod$final_mod$u[span_stroop]), 
       mod$opt$par[mod$beta_inds[20]],
       as.numeric(Lambdat[span_digitspan, span_digitspan] %*% mod$final_mod$u[span_digitspan]), 
       mod$opt$par[mod$beta_inds[19]],
       as.numeric(Lambdat[span_cvlt, span_cvlt] %*% mod$final_mod$u[span_cvlt]), 
       mod$opt$par[mod$beta_inds[18]])

beta_spline <- transU %*% (transD * b)

Xfp <- cbind(as(mod$lmod$X[, str_detect(colnames(mod$lmod$X), "^Xf", negate = TRUE)], "dgCMatrix"), 
             mod$re$stroop$rand$Xr, mod$re$stroop$Xf,
             mod$re$digitspan$rand$Xr, mod$re$digitspan$Xf,
             mod$re$cvlt$rand$Xr, mod$re$cvlt$Xf)

GXf <- do.call(cbind, map(rev(mod$sm), ~ PredictMat(.x, data = dat)))

B <- Matrix(0, ncol(Xfp), ncol(Xfp))
diag(B) <- 1
tmp_ind <- seq(from = ncol(Xfp) - length(beta_spline) + 1L, to = ncol(Xfp))
B[tmp_ind, tmp_ind] <- t(transD * t(transU))

span <- c(span_stroop, span_digitspan, span_cvlt)
V <- Matrix::Diagonal(length(mod$final_mod$V), mod$final_mod$phi[[1]] / mod$final_mod$V) +
  crossprod(Lambdat[-span, -span] %*% Zt[-span, ]) * mod$final_mod$phi[[1]]

R <- Matrix::chol(V, pivot = TRUE)
piv <- attr(R, "pivot")

WX <- as(solve(t(R), Xfp[piv, ]), "matrix")
XVX <- as(solve(t(R), GXf[piv, ]), "matrix")
qrz <- qr(XVX, LAPACK = TRUE)
R <- qr.R(qrz)
R[, qrz$pivot] <- R
XVX <- crossprod(R)

Sp <- matrix(0, ncol(Xfp), ncol(Xfp))
diag(Sp)[tmp_ind] <- c(1 / rep(mod$opt$par[mod$theta_inds][[10]], 8), 0,
                       1 / rep(mod$opt$par[mod$theta_inds][[11]], 8), 0,
                       1 / rep(mod$opt$par[mod$theta_inds][[12]], 8), 0)

qrx <- qr(rbind(WX, Sp), LAPACK = TRUE)
Ri <- backsolve(qr.R(qrx), diag(ncol(WX)))

ind <- qrx$pivot
ind[ind] <- 1:length(ind)
Ri <- Ri[ind, ]
Vb <- B %*% Ri
Vb <- Vb %*% t(Vb)

edf <- sum(rowSums(Vb[tmp_ind, tmp_ind] * t(XVX)))
saveRDS(edf, file = "results/cognition_model/edf.rds")
v <- rowSums((Xp %*% Vb[tmp_ind, tmp_ind]) * Xp)

nmc <- 10000
betas <- mgcv::rmvn(nmc, as.numeric(beta_spline), as.matrix(Vb[tmp_ind, tmp_ind]))
fits <- Xp %*% t(betas)
pred <- c(Xp %*% beta_spline)
simDev <- fits - pred
# Find critical value for each of the three smooth functions
absDev <- map(levels(grid$itemgroup), function(igr) {
  apply(simDev[grid$itemgroup == igr, ], 2L, function(x) abs(x / sqrt(v[grid$itemgroup == igr])))
})
masd <- map(absDev, function(x) apply(x, 2L, max))
crit <- map_dbl(masd, function(x) quantile(x, prob = .95, type = 8))

plot_dat <- grid %>% 
  mutate(
    fit = pred,
    se = sqrt(as.numeric(v)),
    itemgroup = fct_recode(itemgroup, "Episodic memory" = "cvlt",
                           "Working memory" = "digitspan",
                           "Executive function" = "stroop"),
    crit = crit[as.integer(itemgroup)]
  )

p1 <- ggplot(plot_dat, aes(x = age, y = fit)) + 
  geom_line() +
  geom_ribbon(aes(ymin = fit + qnorm(.025) * se,
                  ymax = fit + qnorm(.975) * se), alpha = .3) + 
  geom_ribbon(aes(ymin = fit - crit * se,
                  ymax = fit + crit * se), alpha = .3) + 
  facet_wrap(vars(itemgroup), scales = "free_y") + 
  xlab("Age") +
  theme(axis.title.y = element_blank())

ggsave("results/cognition_model/figures/cognition_smooth_plot.pdf",
       plot = p1,
       width = 16, height = 6, units = "cm")

plot_dat %>% 
  group_by(itemgroup) %>% 
  filter(fit == max(fit))

# Find posterior distributions for age at max
posteriors <- map_dfr(levels(grid$itemgroup), function(x){
  ff <- fits[grid$itemgroup == x, ]
  tibble(
    domain = x,
    max_age = grid$age[grid$itemgroup == x][apply(ff, 2, which.max)]
  ) %>% 
    mutate(iteration = row_number())
}) %>% 
  mutate(
    itemgroup = fct_recode(domain, "Episodic memory" = "cvlt",
                           "Working memory" = "digitspan",
                           "Executive function" = "stroop")
  )

posteriors %>% group_by(itemgroup) %>% summarise(mean(max_age))

p1 <- ggplot(posteriors, aes(x = max_age, group = itemgroup, color = itemgroup, 
                             fill = itemgroup)) + 
  geom_density(alpha = .5) + 
  ggthemes::scale_color_colorblind() + 
  ggthemes::scale_fill_colorblind() + 
  labs(color = NULL, fill = NULL) + 
  theme(
    axis.title.y = element_blank(),
    legend.position = c(.65, .65)
  ) + 
  xlab("Age at maximum")

# Take an illustrative sample of curves from posterior distribution
examples <- map_dfr(levels(grid$itemgroup), function(x){
  ff <- fits[grid$itemgroup == x, 1:100]
  colnames(ff) <- paste0("fit", 1:ncol(ff))
  tibble(
    domain = x,
    age = grid$age[grid$itemgroup == x]
  ) %>% 
    bind_cols(as_tibble(ff)) %>% 
    pivot_longer(cols = starts_with("fit"))
})

p2 <- ggplot(examples, aes(x = age, y = value, group = interaction(domain, name),
                           color = domain)) + 
  geom_line(alpha = .2) +
  ggthemes::scale_color_colorblind() +
  theme(
    axis.title.y = element_blank(),
    legend.position = "none"
  ) + 
  xlab("Age")

p2 + p1 + plot_layout(nrow = 1)

ggsave("results/cognition_model/figures/cognition_posteriors.pdf",
       width = 16, height = 6, units = "cm")

posteriors %>% 
  arrange(iteration, max_age) %>% 
  group_by(iteration) %>% 
  summarise(
    order = paste(itemgroup, collapse = " - "), .groups = "drop"
  ) %>% 
  group_by(order) %>% 
  count() %>% 
  ungroup() %>% 
  mutate(prob = n / sum(n))

# Find estimates and standard errors of regression coefficients and factor loadings
test_vec <- levels(dat$test)
test_ref <- setdiff(test_vec, levels(mod$mdat$test))

lambda_est <- rep(1, length(test_vec))
lambda_est[!test_vec %in% test_ref] <- mod$opt$par[mod$lambda_inds]
lambda_se <- rep(0, length(test_vec))
lambda_se[!test_vec %in% test_ref] <- sqrt(diag(S)[length(mod$beta_inds) + seq_along(mod$lambda_inds)])
names(lambda_se) <- names(lambda_est) <- test_vec

beta_est <- mod$opt$par[mod$beta_inds]
beta_se <- sqrt(diag(S)[seq_along(mod$beta_inds)])
names(beta_se) <- names(beta_est) <- str_remove(colnames(mod$lmod$X), "^Xtest")

cbind(round(beta_est, 2), round(beta_se, 2))
cbind(round(lambda_est, 2), round(lambda_se, 2))

# Stroop effects must be transformed back to seconds
transf_df <- dat %>% 
  filter(str_detect(itemgroup, "stroop")) %>% 
  group_by(itemgroup_retest) %>% 
  summarise(location = mean(value), scale = sd(value))

saveRDS(transf_df, file = "results/cognition_model/transf_df.rds")

beta_est[str_detect(names(beta_est), "^teststroop.*[1-2]+$")] <- 
  -(beta_est[str_detect(names(beta_est), "^teststroop.*[1-2]+$")] * transf_df$scale[[1]] + transf_df$location[[1]])

beta_se[str_detect(names(beta_se), "^teststroop.*[1-2]+$")] <- 
  (beta_se[str_detect(names(beta_se), "^teststroop.*[1-2]+$")] * transf_df$scale[[1]])

beta_est[str_detect(names(beta_est), "^teststroop.*[3-4]+$")] <- 
  -(beta_est[str_detect(names(beta_est), "^teststroop.*[3-4]+$")] * transf_df$scale[[2]] + transf_df$location[[2]])

beta_se[str_detect(names(beta_se), "^teststroop.*[3-4]+$")] <- 
  (beta_se[str_detect(names(beta_se), "^teststroop.*[3-4]+$")] * transf_df$scale[[2]])

lambda_est[str_detect(names(lambda_est), "^stroop.*[1-2]+$")] <- lambda_est[str_detect(names(lambda_est), "^stroop.*[1-2]+$")] *
  transf_df$scale[[1]]

lambda_se[str_detect(names(lambda_se), "^stroop.*[1-2]+$")] <- lambda_se[str_detect(names(lambda_se), "^stroop.*[1-2]+$")] *
  transf_df$scale[[1]]

lambda_est[str_detect(names(lambda_est), "^stroop.*[3-4]+$")] <- lambda_est[str_detect(names(lambda_est), "^stroop.*[3-4]+$")] *
  transf_df$scale[[2]]

lambda_se[str_detect(names(lambda_se), "^stroop.*[3-4]+$")] <- lambda_se[str_detect(names(lambda_se), "^stroop.*[3-4]+$")] *
  transf_df$scale[[2]]

round(cbind(beta_est, beta_se), 2)
round(cbind(lambda_est, lambda_se), 2)

saveRDS(list(beta_true = beta_est, beta_se = beta_se,
             lambda_true = lambda_est, lambda_se = lambda_se),
        "results/cognition_model/params.rds")

beta_retest <- cbind(beta_est, beta_se)[str_detect(names(beta_est), "itemgroup_retest"), ]
beta_retest <- cbind(beta_retest, OR = exp(beta_retest[, 1]))
beta_retest[3, ] <- beta_retest[3, ] * transf_df$scale[[1]]
beta_retest[4, ] <- beta_retest[4, ] * transf_df$scale[[2]]

beta_retest

# Variance components
varcomp_df <- map_dfr(seq_along(mod$lmod$reTrms$cnms), function(i){
  level <- names(mod$lmod$reTrms$cnms)[[i]]
  effect <- mod$lmod$reTrms$cnms[[i]]
  ind <- seq(from = mod$lmod$reTrms$Gp[[i]] + 1L, length.out = length(effect))
  Psi <- as.matrix((t(Lambdat[ind, ind]) %*% Lambdat[ind, ind]) * mod$final_mod$phi[[1]])
  tibble(
    level = level,
    effect = list(effect),
    Psi = list(Psi),
    Cor = list(cov2cor(!!Psi))
  )
})

varcomp_df %>% 
  slice(1:3) %>% 
  unnest(cols = c(effect, Psi, Cor))

varcomp_df %>% 
  slice(4) %>% 
  unnest(cols = c(effect, Psi, Cor))

varcomp_df %>% 
  slice(5:7) %>% 
  unnest(cols = c(effect, Psi)) %>% 
  mutate(lambda = 1 / Psi)

# Dispersion parameters
# Stroop condition 1 and 2
sqrt(mod$final_mod$phi[[1]]) * transf_df$scale[[1]]
# Stroop conditions 3 and 4
sqrt(mod$final_mod$phi[[1]] / mod$opt$par[mod$weights_inds]) * transf_df$scale[[2]]

# Proportion of zeros
1 - length(Zt@x) / prod(dim(Zt))

