library(gamm4)
library(dplyr)
library(numDeriv)

opt <- readRDS("results/cvlt_model/opt.rds")
source("code/cvlt_model/model_setup.R")

lambda_est <- c(a1 = 1, opt$par)
dat$weight <- lambda_est[dat$item]

final_gamm <- gamm4(
    formula = formula, random = random,
    family = binomial(), data = dat)

beta_hat <- coef(final_gamm$gam)
cov_lambda <- solve(opt$hessian)
cov_beta_naive <- vcov(final_gamm$gam)

gammstart <- list(
    theta = getME(final_gamm$mer, "theta"),
    fixef = getME(final_gamm$mer, "fixef"))

beta_deriv <- jacobian(function(x){
    lambda <- c(a1 = 1, x)
    dat$weight <- lambda[dat$item]
    mod <- gamm4(formula = formula, random = random,
                 family = binomial(), data = dat,
                 start = gammstart,
                 control = glmerControl(calc.derivs = FALSE))

    coef(mod$gam)
    }, x = lambda_est[-1], method = "simple")

cov_beta <- cov_beta_naive + beta_deriv %*% cov_lambda %*% t(beta_deriv)

res <- list(
    final_gamm = final_gamm,
    beta_deriv = beta_deriv,
    cov_beta_naive = cov_beta_naive,
    cov_beta = cov_beta,
    cov_lambda = cov_lambda,
    lambda_est = lambda_est
    )

saveRDS(res, "results/cvlt_model/res.rds")
