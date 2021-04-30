library(dplyr)
library(gamm4)

source("code/cvlt_model/model_setup.R")

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

message("Finished optimizing")

saveRDS(opt, file = "results/cvlt_model/opt.rds")
