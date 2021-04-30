lambda_init <- c(1, runif(4, min = 1, max = 2))
names(lambda_init) <- paste0("a", 1:5)

formula <- as.formula(
  cbind(successes, failures) ~ 0 + a2 + a3 + a4 + a5 + retest +
    s(age_z, by = weight, k = 15, bs = "cr")
)
random <- ~ (0 + weight | id / timepoint)

run_simulation <- function(n, iter, sim_params){
  dat <- sample_data(n, sim_params)
  
  opt <- optim(lambda_init[-1], function(x){
    ll <- c(a1 = 1, x)
    dat$weight <- ll[dat$item]
    if(!exists("gammstart")) gammstart <- NULL
    
    mod <- gamm4(
      formula = formula, random = random, family = binomial(), 
      start = gammstart, data = dat, 
      control = glmerControl(calc.derivs = FALSE)
    )
    gammstart <<- list(theta = getME(mod$mer, "theta"),
                       fixef = getME(mod$mer, "fixef"))
    -as.numeric(logLik(mod$mer))
  }, method = "L-BFGS-B", control = list(trace = 3),
  hessian = TRUE)
  
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
    dd <- dat
    dd$weight <- lambda[dd$item]
    mod <- gamm4(formula = formula, random = random,
                 family = binomial(), data = dd,
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
}
