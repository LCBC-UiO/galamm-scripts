# Usage: Rscript code/ses_simulation/run_simulation.R "iter" "level_2_theta" "0/1"
iter <- as.integer(commandArgs(trailingOnly = TRUE)[[1]])
level_2_theta <- as.numeric(commandArgs(trailingOnly = TRUE)[[2]])
only_opt <- as.integer(commandArgs(trailingOnly = TRUE)[[3]])

library(tidyverse)
library(galamm)
library(lme4)
library(memoise)
library(Matrix)
library(mgcv)
library(gamm4)
library(mvtnorm)
library(patchwork)

# Keep structure of data but simulate new responses
dat <- readRDS("data/cognition_model/dat.rds")
mod <- readRDS("results/cognition_model/cognition_model.rds")

# Update matrices with final parameters
Lambdat <- mod$lmod$reTrms$Lambdat
Lambdat@x <- mod$opt$par[mod$theta_inds][mod$lmod$reTrms$Lind]
Zt <- mod$lmod$reTrms$Zt
Zt@x <- Zt@x * c(1, mod$opt$par[mod$lambda_inds])[mod$lambda_mapping_Zt + 2L]
X <- mod$lmod$X
X <- X * c(1, mod$opt$par[mod$lambda_inds])[mod$lambda_mapping_X + 2L]

span_stroop <- seq(from = mod$lmod$reTrms$Gp[[5]] + 1L, to = mod$lmod$reTrms$Gp[[6]])
span_digitspan <- seq(from = mod$lmod$reTrms$Gp[[6]] + 1L, to = mod$lmod$reTrms$Gp[[7]])
span_cvlt <- seq(from = mod$lmod$reTrms$Gp[[7]] + 1L, to = mod$lmod$reTrms$Gp[[8]])

# Level 2 random effects
level_2_var <- c(level_2_theta, level_2_theta, mod$opt$par[mod$theta_inds[[3]]])^2 * mod$final_mod$phi[[1]]
names(level_2_var) <- map_chr(mod$lmod$reTrms$cnms[1:3], ~ str_extract(.x, "[:alpha:]+$"))
level_2_ranef <- dat %>%
  distinct(itemgroup, id, timepoint) %>%
  mutate(zeta_2 = rnorm(nrow(.), sd = sqrt(level_2_var[as.character(itemgroup)])))

# Level 3 random effects
Psi3 <- as.matrix(t(Lambdat[11116:11118, 11116:11118]) %*% Lambdat[11116:11118, 11116:11118]) *
  mod$final_mod$phi[[1]]

zeta_3 <- rmvnorm(length(unique(dat$id)), sigma = Psi3)
colnames(zeta_3) <- map_chr(mod$lmod$reTrms$cnms[[4]], ~ str_remove(.x, "itemgroup"))

fixed_part <- as.numeric(X %*% mod$opt$par[mod$beta_inds])

penalized_stroop <- as.numeric(t(Lambdat[span_stroop, span_stroop] %*% Zt[span_stroop, ]) %*% mod$final_mod$u[span_stroop])
penalized_digitspan <- as.numeric(t(Lambdat[span_digitspan, span_digitspan] %*% Zt[span_digitspan, ]) %*% mod$final_mod$u[span_digitspan])
penalized_cvlt <- as.numeric(t(Lambdat[span_cvlt, span_cvlt] %*% Zt[span_cvlt, ]) %*% mod$final_mod$u[span_cvlt])
penalized_part <- penalized_stroop + penalized_digitspan + penalized_cvlt

loading <- c(1, mod$opt$par[mod$lambda_inds])[mod$factor_mapping + 2L]

simdat <- dat %>%
  distinct(id) %>%
  bind_cols(as_tibble(zeta_3)) %>%
  pivot_longer(cols = c("cvlt", "digitspan", "stroop"), values_to = "zeta_3",
               names_to = "itemgroup") %>%
  mutate(itemgroup = factor(itemgroup, levels = levels(dat$itemgroup))) %>%
  inner_join(dat %>% mutate(rn = row_number()), by = c("id", "itemgroup")) %>%
  inner_join(level_2_ranef, by = c("id", "itemgroup", "timepoint")) %>% 
  arrange(rn) %>% 
  mutate(
    linpred = fixed_part + penalized_part + loading * (zeta_2 + zeta_3),
    value_z_orig = value_z
  ) %>% 
  nest_by(itemgroup_retest, .keep = TRUE) %>% 
  pmap_dfr(function(itemgroup_retest, data){
    if(itemgroup_retest %in% c("stroop_12", "stroop_34")){
      
      if(itemgroup_retest == "stroop_12"){
        residuals <- rnorm(nrow(data), sd = sqrt(mod$final_mod$phi[[1]]))
      } else {
        residuals <- rnorm(nrow(data), sd = sqrt(mod$final_mod$phi[[1]] / mod$opt$par[mod$weights_inds]))
      }
      data %>% 
        mutate(
          value_z = linpred + residuals#,
          #value_z = (value_z - mean(value_z)) / sd(value_z)
        )
    } else {
      data %>% 
        mutate(
          value_z = rbinom(nrow(.), size = 16L, prob = plogis(linpred))
        )
    }
  }) %>% 
  arrange(rn)


dat <- simdat
rm(simdat)

test_ref <- c("cvlt_a1", "digitspan_bwd", "stroop_time_1")
sm <- smoothCon(s(age_z, by = itemgroup, bs = "cr", k = 10), 
                data = dat, absorb.cons = TRUE)
re <- lapply(sm, smooth2random, 2, "")
names(re) <- names(sm) <- levels(dat$itemgroup)

mdat <- as.list(dat)
for(igr in levels(mdat$itemgroup)){
  eval(parse(text = paste0("mdat$itemgroup_", igr, "<- as.numeric(dat$itemgroup == '", igr, "')")))
}
# Add intercept for each test that has not been standardized to zero
mdat$Xtest <- model.matrix(~ 0 + test + retest:itemgroup_retest, data = dat)

mdat$Xf_cvlt <- re[["cvlt"]]$Xf
mdat$Xr_cvlt <- re[["cvlt"]]$rand$Xr
mdat$pseudoGroups_cvlt <- as.integer(as.factor(rep(re[["cvlt"]]$rind, length = nrow(dat))))
mdat$Xf_digitspan <- re[["digitspan"]]$Xf
mdat$Xr_digitspan <- re[["digitspan"]]$rand$Xr
mdat$pseudoGroups_digitspan <- as.integer(as.factor(rep(re[["digitspan"]]$rind, length = nrow(dat))))
mdat$Xf_stroop <- re[["stroop"]]$Xf
mdat$Xr_stroop <- re[["stroop"]]$rand$Xr
mdat$pseudoGroups_stroop <- as.integer(as.factor(rep(re[["stroop"]]$rind, length = nrow(dat))))

lmod <- glFormula(cbind(value, trials - value) ~ 0 + Xtest + Xf_cvlt + Xf_digitspan + 
                    Xf_stroop + (1 | pseudoGroups_cvlt) + 
                    (1 | pseudoGroups_digitspan) + 
                    (1 | pseudoGroups_stroop) + 
                    (0 + itemgroup | id) + (0 + itemgroup_cvlt | id:timepoint) + 
                    (0 + itemgroup_digitspan | id:timepoint) + 
                    (0 + itemgroup_stroop | id:timepoint), 
                  family = "binomial", data = mdat)

lmod$reTrms$Ztlist$`1 | pseudoGroups_stroop` <- as(t(as.matrix(mdat$Xr_stroop))[], class(lmod$reTrms$Zt))
lmod$reTrms$Ztlist$`1 | pseudoGroups_digitspan` <- as(t(as.matrix(mdat$Xr_digitspan))[], class(lmod$reTrms$Zt))
lmod$reTrms$Ztlist$`1 | pseudoGroups_cvlt` <- as(t(as.matrix(mdat$Xr_cvlt))[], class(lmod$reTrms$Zt))

lmod$reTrms$Zt <- do.call(rbind, lmod$reTrms$Ztlist)

mdat$test <- as.character(mdat$test)
mdat$test[mdat$test %in% c(test_ref)] <- "reference"
mdat$test <- fct_drop(factor(mdat$test, levels = c("reference", levels(dat$test))))
factor_mapping <- as.integer(mdat$test) - 2L

delta <- unique(diff(lmod$reTrms$Zt@p))
lambda_mapping_Zt <- rep(factor_mapping, each = delta)

lambda_mapping_X <- rep(-1L, length(lmod$X))
lambda_mapping_X[seq(from = (min(str_which(colnames(lmod$X), "^Xf_")) - 1L) * nrow(lmod$X) + 1L,
                     to = length(lmod$X))] <- rep(factor_mapping, 3)

# Apply weights for last two Stroop conditions, as they differ quite a lot from the first two
weights_mapping <- if_else(mdat$test %in% c("stroop_time_3", "stroop_time_4"), 0L, -1L)

stopifnot(length(lambda_mapping_Zt) == length(lmod$reTrms$Zt@x))
stopifnot(length(lambda_mapping_X) == length(lmod$X))

theta_inds <- seq_along(lmod$reTrms$theta)
beta_inds <- seq(from = max(theta_inds) + 1L, length.out = ncol(lmod$X))
lambda_inds <- seq(from = max(beta_inds) + 1L,
                   length.out = sum(unique(lambda_mapping_Zt) != -1))
weights_inds <- seq(from = max(lambda_inds) + 1L, length.out = sum(unique(weights_mapping != -1L)))

lbound <- c(lmod$reTrms$lower, rep(-Inf, length(beta_inds)),
            rep(-Inf, length(lambda_inds)), rep(.1, length(weights_inds)))

mlwrapper <- function(par, gradient, hessian) {
  marginal_likelihood(
    y = mdat$value_z,
    trials = mdat$trials,
    X = lmod$X,
    Zt = lmod$reTrms$Zt,
    Lambdat = lmod$reTrms$Lambdat,
    beta = par[beta_inds],
    theta = par[theta_inds],
    theta_mapping = lmod$reTrms$Lind - 1L,
    lambda = par[lambda_inds],
    lambda_mapping_Zt = lambda_mapping_Zt,
    lambda_mapping_X = lambda_mapping_X,
    weights = par[weights_inds],
    weights_mapping = weights_mapping,
    family = c("gaussian", "binomial"),
    family_mapping = if_else(mdat$itemgroup == "stroop", 0L, 1L),
    maxit_conditional_modes = 20L,
    gradient = gradient, 
    hessian = hessian
  )
}

mlmem <- memoise(mlwrapper)  
fn <- function(par, gradient = FALSE, hessian = FALSE) {
  mlmem(par = par, gradient = gradient, hessian = hessian)$logLik  
}
gr <- function(par, gradient = TRUE, hessian = FALSE) {
  mlmem(par = par, gradient = gradient, hessian = hessian)$gradient
}

# Model without weights, otherwise identical
par_init <- c(lmod$reTrms$theta, rep(0, length(beta_inds)), rep(1, length(lambda_inds)), rep(1, length(weights_inds)))

stopifnot(is.finite(fn(par_init)))
stopifnot(length(gr(par_init)) == length(par_init))

opt <- optim(par_init, fn, gr, gradient = TRUE,
             method = "L-BFGS-B", lower = lbound,
             control = list(fnscale = -1, trace = 3, REPORT = 1,
                            maxit = 5000, factr = 1e8))


if(only_opt == 0L){
  final_mod <- mlwrapper(opt$par, gradient = TRUE, hessian = TRUE)
  
  saveRDS(list(
    final_mod = final_mod, opt = opt,
    lmod = lmod, re = re, sm = sm, mdat = mdat,
    beta_inds = beta_inds, theta_inds = theta_inds,
    lambda_inds = lambda_inds, weights_inds = weights_inds,
    lambda_mapping_Zt = lambda_mapping_Zt
  ), paste0("results/cognition_model/simulation/model_fits/simulation_", iter, "_", level_2_theta, ".rds"))
} else {
  saveRDS(
    opt, 
    paste0("results/cognition_model/simulation/model_fits/simulation_opt_", iter, "_", level_2_theta, ".rds")
  )
}


