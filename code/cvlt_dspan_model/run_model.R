library(tidyverse)
library(gamm4)

source("code/cvlt_dspan_model/model_setup.R")

lambda_init <- c(1, c(1.5, 2, 2,2.5), 1, .5)
names(lambda_init) <- c(paste0("a", 1:5), "bwd", "fwd")
fix_inds <- c(1, 6)
fix_nms <- names(lambda_init)[fix_inds]

opt <- optim(lambda_init[-(fix_inds)], function(x){
  ll <- x
  ll[fix_nms] <- 1
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

saveRDS(opt, file = "results/cvlt_dspan_model/opt.rds")
