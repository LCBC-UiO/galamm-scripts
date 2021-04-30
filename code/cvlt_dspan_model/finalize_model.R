library(tidyverse)
library(numDeriv)
library(latex2exp)
options(bitmapType = "cairo")
theme_set(theme_bw())
library(gamm4)

opt <- readRDS("results/cvlt_dspan_model/opt.rds")
source("code/cvlt_dspan_model/model_setup.R")

lambda_est <- c(a1 = 1, opt$par[paste0("a", 2:5)],
                bwd = 1, fwd = opt$par[["fwd"]])
fixed_lambda <- which(lambda_est == 1)
dat$weight <- lambda_est[dat$item]

final_gamm <- gamm4(
  formula = formula, random = random,
  family = binomial(), data = dat, verbose = 2L)

beta_hat <- coef(final_gamm$gam)
cov_lambda <- solve(opt$hessian)
cov_beta_naive <- vcov(final_gamm$gam)

gammstart <- list(
  theta = getME(final_gamm$mer, "theta"),
  fixef = getME(final_gamm$mer, "fixef"))

beta_deriv <- jacobian(function(x){
  lambda <- c(a1 = 1, x[paste0("a", 2:5)],
              bwd = 1, fwd = x[["fwd"]])
  dat$weight <- lambda[dat$item]
  mod <- gamm4(formula = formula, random = random,
               family = binomial(), data = dat,
               control = glmerControl(calc.derivs = FALSE))
  
  coef(mod$gam)
}, x = lambda_est[-fixed_lambda], method = "Richardson")

cov_beta <- cov_beta_naive + beta_deriv %*% cov_lambda %*% t(beta_deriv)

res <- list(
  final_gamm = final_gamm,
  beta_deriv = beta_deriv,
  cov_beta_naive = cov_beta_naive,
  cov_beta = cov_beta,
  cov_lambda = cov_lambda,
  lambda_est = lambda_est
)

saveRDS(res, "results/cvlt_dspan_model/res.rds")
