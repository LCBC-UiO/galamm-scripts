library(galamm)
library(memoise)
library(tidyverse)
library(gamm4)
library(furrr)
dat <- readRDS("data/ses_model/dat.rds")

base_formula <- "value_z ~ 0 + itemedu_father + itemedu_mother + itemedu_self +
  itemincome_father + itemincome_mother + itemincome_self + 
  itemgrouphippocampus:(siteousAvanto + siteousPrisma + siteousSkyra + icv_z + sexmale) +
  Xf + (1 | pseudoGroups) + (0 + itemgrouphippocampus | id)"


for(case in letters[7:1]){
  if(file.exists(paste0("results/ses_model/case_", case, ".rds"))){
    next
  }
  sm <- smoothCon(s(age_z, by = itemgrouphippocampus, k = 15, bs = "cr"), data = dat)[[1]]
  re <- smooth2random(sm, "", 2)
  
  mdat <- as.list(dat)
  mdat$Xf <- re$Xf
  mdat$Xr <- re$rand$Xr
  mdat$pseudoGroups <- rep(1:ncol(re$rand$Xr), length = nrow(dat))
  mdat$itemgroupses <- as.integer(mdat$itemgroup %in% c("edu", "income"))
  
  form <- paste(
    base_formula,
    if(case == "g") "(0 + itemgroupses | id)" else "(1 | id)", sep = " + ")
  
  lmod <- lFormula(as.formula(form), data = mdat)
  
  lmod$reTrms$Ztlist$`1 | pseudoGroups` <- as(t(as.matrix(mdat$Xr))[], class(lmod$reTrms$Zt))
  lmod$reTrms$Zt <- do.call(rbind, lmod$reTrms$Ztlist)
  delta <- diff(lmod$reTrms$Zt@p)
  
  itemvec <- mdat$item
  if(case %in% c("a", "b")){
    factor_mapping <- as.integer(itemvec) - 2L
  } else if(case %in% c("c", "d")){
    factor_mapping <- as.integer(
      factor(
        str_remove(itemvec, "_father|_mother"),
        levels = c("edu", "edu_self", "income", "income_self", "hippocampus")
      )) - 2L
  } else if(case %in% c("e", "f", "g")){
    factor_mapping <- as.integer(
      factor(
        str_remove(itemvec, "_father|_mother|_self"),
        levels = c("edu", "income", "hippocampus")
      )
    ) - 2L
    
    if(case == "g"){
      factor_mapping[itemvec == "hippocampus"] <- -1L
    }
  }
  
  lambda_mapping_Zt <- unlist(map(seq_along(delta), ~ c(factor_mapping[[.x]], rep(-1L, delta[[.x]] - 1))))
  
  if(case %in% c("a", "c", "e")){
    mm <- max(factor_mapping)
    lambda_mapping_Zt <- map(lambda_mapping_Zt, ~ if(.x == mm) c(mm, mm + 1L) else .x)
    age_vec <- mdat$age_z[mdat$itemgroup == "hippocampus"]
    lambda_mapping_Zt_covs <- vector(mode = "list", length = length(lambda_mapping_Zt))
    j <- 1
    for(i in seq_along(lambda_mapping_Zt)){
      if(length(lambda_mapping_Zt[[i]]) == 1) {
        lambda_mapping_Zt_covs[[i]] <- 1
      } else {
        lambda_mapping_Zt_covs[[i]] <- c(1, age_vec[[j]])
        j <- j + 1
      }
    }
    walk2(lambda_mapping_Zt, lambda_mapping_Zt_covs, function(x, y){
      stopifnot(length(x) == length(y))
    })
  } else {
    lambda_mapping_Zt_covs <- integer()
  }
 
  theta_inds <- seq(from = 1, to = length(lmod$reTrms$theta))
  beta_inds <- seq(from = max(theta_inds) + 1L, length.out = ncol(lmod$X))
  lambda_inds <- seq(from = max(beta_inds) + 1L, length.out = max(unlist(lambda_mapping_Zt)) + 1L)
  weights_mapping <- if_else(mdat$itemgroup == "hippocampus", 1L, if_else(mdat$itemgroup == "edu", 0L, -1L))
  weights_inds <- seq(from = max(lambda_inds) + 1L, length.out = length(unique(weights_mapping[weights_mapping != -1L])))
  
  lbound <- c(lmod$reTrms$lower, rep(-Inf, length(beta_inds)), 
              rep(-Inf, length(lambda_inds)), rep(0.1, length(weights_inds)))
  
  mlwrapper <- function(par, gradient, hessian){
    marginal_likelihood(
      y = mdat$value_z,
      X = lmod$X,
      Zt = lmod$reTrms$Zt,
      Lambdat = lmod$reTrms$Lambdat,
      beta = par[beta_inds],
      theta = par[theta_inds],
      theta_mapping = lmod$reTrms$Lind - 1L,
      lambda = par[lambda_inds],
      lambda_mapping_Zt = lambda_mapping_Zt,
      lambda_mapping_Zt_covs = lambda_mapping_Zt_covs,
      weights = par[weights_inds],
      weights_mapping = weights_mapping,
      gradient = gradient,
      hessian = hessian
    )
  }
  
  par_init <- tryCatch({
    if(case == "f"){
      tmp <- readRDS(file.path("results/ses_model", "case_g.rds"))$opt$par
      c(tmp[1:17], 0, tmp[18:19])
    } else if(case == "e"){
      tmp <- readRDS(file.path("results/ses_model", "case_f.rds"))$opt$par
      c(tmp[1:18], 0, tmp[19:20])
    } else if(case == "d"){
      tmp <- readRDS(file.path("results/ses_model", "case_f.rds"))$opt$par
      c(tmp[1:16], 1, tmp[17], tmp[17], tmp[18], tmp[19:20])
    } else if(case == "c"){
      tmp <- readRDS(file.path("results/ses_model", "case_d.rds"))$opt$par
      c(tmp[1:20], 0, tmp[21:22])
    } else if(case == "b"){
      tmp <- readRDS(file.path("results/ses_model", "case_d.rds"))$opt$par
      c(tmp[1:16], tmp[17], tmp[17], tmp[18], tmp[18], tmp[19], tmp[20], tmp[21:22])
    } else if(case == "a"){
      tmp <- readRDS(file.path("results/ses_model", "case_b.rds"))$opt$par
      c(tmp[1:22], 0, tmp[23:24])
    } else {
      stop()
    }
  }, 
  error = function(e){
    c(lmod$reTrms$theta, rep(0, ncol(lmod$X)), runif(length(lambda_inds)), c(1, 2))
  })
  
  mlmem <- memoise(mlwrapper)
  fn <- function(par, gradient = FALSE, hessian = FALSE){
    mlmem(par, gradient, hessian)$logLik
  }
  gr <- function(par, gradient = TRUE, hessian = FALSE){
    mlmem(par, gradient, hessian)$gradient
  }
  
  stopifnot(length(gr(par_init)) == length(par_init))
  
  opt <- optim(
    par = par_init, fn = fn, gr = gr, method = "L-BFGS-B", lower = lbound, 
    control = list(fnscale = -1, maxit = 5000, trace = 3, REPORT = 10))
  
  final_model <- mlwrapper(opt$par, TRUE, TRUE)
  S <- tryCatch({
    -solve(final_model$hessian)
  }, error = function(e) NULL)
  
  ret <- list(opt = opt, final_model = final_model,
              beta_inds = beta_inds, theta_inds = theta_inds,
              lambda_inds = lambda_inds, weights_inds = weights_inds,
              reTrms = lmod$reTrms, X = lmod$X,
              lambda_mapping_Zt = lambda_mapping_Zt,
              lambda_mapping_Zt_covs = lambda_mapping_Zt_covs,
              sm = sm, re = re, mdat = mdat, weights_mapping = weights_mapping,
              S = S)
  
  Lambdat <- lmod$reTrms$Lambdat
  Lambdat@x <- opt$par[theta_inds][lmod$reTrms$Lind]
  
  # Fitting parametrization, modulo Cholesky
  spline_inds <- (lmod$reTrms$Gp[[3]] + 1L) : lmod$reTrms$Gp[[4]]
  b <- as.numeric(Lambdat[spline_inds, spline_inds] %*% final_model$u[spline_inds])
  
  # Back to original parametrization
  beta_spline <- re$trans.U %*% (
      re$trans.D * c(b, opt$par[beta_inds[str_which(colnames(lmod$X), "Xf[:digit:]")]]))
  
  Zt <- lmod$reTrms$Zt
  if(length(lambda_mapping_Zt_covs) == 0){
    Zt@x <- c(1, opt$par[lambda_inds])[lambda_mapping_Zt + 2L]
  } else {
    Zt@x <- map2_dbl(lambda_mapping_Zt, lambda_mapping_Zt_covs, function(l, x){
      if(identical(l, -1L)){
        1
      } else {
        sum(opt$par[lambda_inds[l + 1L]] * x)
      }
    })
  }
  
  Xfp <- cbind(as(lmod$X, "dgCMatrix")[, -str_which(colnames(lmod$X), "Xf")], re$rand$Xr, re$Xf)
  B <- Matrix(0, ncol(Xfp), ncol(Xfp))
  diag(B) <- 1
  tmp_ind <- seq(from = ncol(Xfp) - length(beta_spline) + 1L, to = ncol(Xfp))
  B[tmp_ind, tmp_ind] <- t(re$trans.D * t(re$trans.U))
  
  V <- Matrix::Diagonal(
    n = length(mdat$value), x = final_model$phi / final_model$V) + 
    crossprod(Lambdat[-spline_inds, -spline_inds] %*% Zt[-spline_inds, ]) * 
    final_model$phi
  
  R <- Matrix::chol(V, pivot = TRUE); piv <- attr(R, "pivot")
  
  WX <- as(solve(t(R), Xfp), "matrix")
  Sp <- matrix(0, ncol(Xfp), ncol(Xfp))
  diag(Sp)[c(rep(FALSE, min(tmp_ind) - 1L), (re$pen.ind == 1))] <- 1 / opt$par[theta_inds[[3]]]
  
  qrx <- qr(rbind(WX, Sp), LAPACK = TRUE)
  Ri <- backsolve(qr.R(qrx), diag(ncol(WX)))
  
  ind <- qrx$pivot;ind[ind] <- 1:length(ind)## qrx$pivot
  Ri <- Ri[ind,] ## unpivoted square root of cov matrix in fitting parameterization Ri Ri' = cov
  
  Vb <- B%*%Ri; Vb <- Vb%*%t(Vb)
  
  ret$spline_coefs <- list(
    beta_spline = beta_spline, Vb = Vb
  )
  
  saveRDS(ret, paste0("results/ses_model/case_", case, ".rds"))
}

# Fit case f model with factor loading for income set to zero, to compute likelihood
# ratio test
sm <- smoothCon(s(age_z, by = itemgrouphippocampus, k = 15, bs = "cr"), data = dat)[[1]]
re <- smooth2random(sm, "", 2)

mdat <- as.list(dat)
mdat$Xf <- re$Xf
mdat$Xr <- re$rand$Xr
mdat$pseudoGroups <- rep(1:ncol(re$rand$Xr), length = nrow(dat))
mdat$itemgroupses <- as.integer(mdat$itemgroup %in% c("edu", "income"))

form <- paste(base_formula, "(1 | id)", sep = " + ")

lmod <- lFormula(as.formula(form), data = mdat)
lmod$reTrms$Ztlist$`1 | pseudoGroups` <- as(t(as.matrix(mdat$Xr))[], class(lmod$reTrms$Zt))
lmod$reTrms$Zt <- do.call(rbind, lmod$reTrms$Ztlist)
delta <- diff(lmod$reTrms$Zt@p)

itemvec <- mdat$item

factor_mapping <- as.integer(factor(
      str_remove(itemvec, "_father|_mother|_self"),
      levels = c("edu", "income", "hippocampus")
    )) - 2L

lambda_mapping_Zt <- unlist(map(seq_along(delta), ~ c(factor_mapping[[.x]], rep(-1L, delta[[.x]] - 1))))
lambda_mapping_Zt_covs <- integer()

theta_inds <- seq(from = 1, to = length(lmod$reTrms$theta))
beta_inds <- seq(from = max(theta_inds) + 1L, length.out = ncol(lmod$X))
lambda_inds <- seq(from = max(beta_inds) + 1L, length.out = 1L)
weights_mapping <- if_else(mdat$itemgroup == "hippocampus", 1L, if_else(mdat$itemgroup == "edu", 0L, -1L))
weights_inds <- seq(from = max(lambda_inds) + 1L, length.out = length(unique(weights_mapping[weights_mapping != -1L])))

lbound <- c(lmod$reTrms$lower, rep(-Inf, length(beta_inds)), 
            rep(-Inf, length(lambda_inds)), rep(0.1, length(weights_inds)))

# Run model, fixing income loading to 0
mlwrapper <- function(par, gradient, hessian){
  marginal_likelihood(
    y = mdat$value_z,
    X = lmod$X,
    Zt = lmod$reTrms$Zt,
    Lambdat = lmod$reTrms$Lambdat,
    beta = par[beta_inds],
    theta = par[theta_inds],
    theta_mapping = lmod$reTrms$Lind - 1L,
    lambda = c(0, par[lambda_inds]),
    lambda_mapping_Zt = lambda_mapping_Zt,
    lambda_mapping_Zt_covs = lambda_mapping_Zt_covs,
    weights = par[weights_inds],
    weights_mapping = weights_mapping,
    gradient = gradient,
    hessian = hessian
  )
}

# Initialize with case f fit, except loading for income
par_init <- readRDS(file.path("results/ses_model", "case_f.rds"))$opt$par[-17]

mlmem <- memoise(mlwrapper)
fn <- function(par, gradient = FALSE, hessian = FALSE){
  mlmem(par, gradient, hessian)$logLik
}
gr <- function(par, gradient = TRUE, hessian = FALSE){
  mlmem(par, gradient, hessian)$gradient[-17] # Remove the 17th element of the gradient
}

stopifnot(length(gr(par_init)) == length(par_init))

opt <- optim(
  par = par_init, fn = fn, gr = gr, method = "L-BFGS-B", lower = lbound, 
  control = list(fnscale = -1, maxit = 5000, trace = 3, REPORT = 10))

saveRDS(opt, file = "results/ses_model/case_f_no_income_loading.rds")