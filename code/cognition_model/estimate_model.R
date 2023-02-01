library(tidyverse)
library(galamm)
library(lme4)
library(memoise)
library(Matrix)
library(mgcv)
library(gamm4)

dat <- readRDS("data/cognition_model/dat.rds")

# Reference levels for identifiability of factor loadings. 
# The corresponding loadings are set to 1.
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

final_mod <- mlwrapper(opt$par, gradient = TRUE, hessian = TRUE)

saveRDS(list(
  final_mod = final_mod, opt = opt,
  lmod = lmod, re = re, sm = sm, mdat = mdat,
  beta_inds = beta_inds, theta_inds = theta_inds,
  lambda_inds = lambda_inds, weights_inds = weights_inds,
  factor_mapping = factor_mapping,
  lambda_mapping_Zt = lambda_mapping_Zt,
  lambda_mapping_X = lambda_mapping_X,
), "results/cognition_model/cognition_model.rds")
