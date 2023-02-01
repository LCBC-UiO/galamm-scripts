# Usage: Rscript code/ses_simulation/run_simulation.R "iter"
iter <- as.integer(commandArgs(trailingOnly = TRUE)[[1]])
library(tidyverse)
library(memoise)
library(gamm4)
library(galamm)

sim_params <- readRDS("data/ses_model/sim_params.rds")

outfile <- file.path(
  "results/ses_model/simulation/model_fits", paste0("simulation_", iter, ".rds"))

sims <- map(sim_params$lambda_interaction, 
  function(lambda_interaction){
  
  dat <- tibble(
    id = seq_len(sim_params$n),
    age_at_baseline = sim_params$age_at_baseline_dist(runif(sim_params$n)),
    item = map(seq_len(sim_params$n), function(i){
      v1 <- c("edu_father", "edu_mother", "edu_self", 
              "income_father", "income_mother", "income_self")
      v2 <- rep("hippocampus", sample(seq_along(sim_params$timepoints_dist), 1, 
                                      prob = sim_params$timepoints_dist))
        c(v1, v2)}),
      eta1 = rnorm(sim_params$n, sd = sim_params$psi1),
      eta2 = rnorm(sim_params$n, sd = sim_params$psi2)
    ) %>% 
    unnest(cols = item) %>% 
    filter(
      (age_at_baseline < 20 & !str_detect(item, "self")) |
        (age_at_baseline >= 20 & str_detect(item, "self")) |
        item == "hippocampus"
    ) %>% 
    mutate(time_interval = sim_params$interval_dist(runif(nrow(.)))) %>% 
    group_by(id, item) %>% 
    mutate(
      age = age_at_baseline + cumsum(c(0, time_interval[-1]))
    ) %>% 
    ungroup() %>% 
    filter(between(age, sim_params$age_min, sim_params$age_max)) %>% 
    select(-age_at_baseline, -time_interval) %>% 
    mutate(
      residual = case_when(
        str_detect(item, "edu") ~ 
          rnorm(nrow(.), sd = sim_params$sigma_edu),
        str_detect(item, "income") ~ 
          rnorm(nrow(.), sd = sim_params$sigma_income),
        item == "hippocampus" ~
          rnorm(nrow(.), sd = sim_params$sigma_hippocampus)
      ),
      age_z = (age - sim_params$age_mean) / sim_params$age_sd,
      smooth_age = (item == "hippocampus") * sim_params$smooth_age(age_z),
      intercept = if_else(item == "hippocampus", 0, sim_params$beta[paste0("item", item)]),
      loading = sim_params$lambda[item] * eta1 + 
        (item == "hippocampus") * (!!lambda_interaction * age_z * eta1 + eta2),
      value_z = smooth_age + intercept + loading + residual,
      itemgroup = str_extract(item, "edu|income|hippocampus")
    ) %>% 
      bind_cols(
        as_tibble(model.matrix(~ 0 + item + itemgroup, data = .))
      )
  
  form <- value_z ~ 0 + itemedu_father + itemedu_mother + itemedu_self +
  itemincome_father + itemincome_mother + itemincome_self +
  Xf + (1 | pseudoGroups) + (0 + itemgrouphippocampus | id) + (1 | id)

  cases <- map(c("e", "f"), function(case){
    sm <- smoothCon(s(age_z, by = itemgrouphippocampus, k = 15, bs = "cr"), data = dat)[[1]]
    re <- smooth2random(sm, "", 2)
    
    mdat <- as.list(dat)
    mdat$Xf <- re$Xf
    mdat$Xr <- re$rand$Xr
    mdat$pseudoGroups <- rep(1:ncol(re$rand$Xr), length = nrow(dat))
    mdat$itemgroupses <- as.integer(mdat$itemgroup %in% c("edu", "income"))
    
    lmod <- lFormula(form, data = mdat)
    
    lmod$reTrms$Ztlist$`1 | pseudoGroups` <- as(t(as.matrix(mdat$Xr))[], class(lmod$reTrms$Zt))
    lmod$reTrms$Zt <- do.call(rbind, lmod$reTrms$Ztlist)
    delta <- diff(lmod$reTrms$Zt@p)
   
    factor_mapping <- as.integer(
      factor(str_remove(mdat$item, "_father|_mother|_self"),
        levels = c("edu", "income", "hippocampus"))) - 2L
    
    lambda_mapping_Zt <- unlist(map(seq_along(delta), ~ c(factor_mapping[[.x]], rep(-1L, delta[[.x]] - 1))))
    
    if(case == "e"){
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
  
    tmp <- readRDS(paste0("results/ses_model/case_", case, ".rds"))
    
    par_init <- c(tmp$opt$par[tmp$theta_inds],
                  tmp$opt$par[tmp$beta_inds[which(colnames(tmp$X) %in% colnames(lmod$X))]],
                  tmp$opt$par[tmp$lambda_inds],
                  tmp$opt$par[tmp$weights_inds])
    
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
                lambda_mapping_Zt = lambda_mapping_Zt,
                lambda_mapping_Zt_covs = lambda_mapping_Zt_covs,
                sm = sm, re = re, weights_mapping = weights_mapping,
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
    Sp <- matrix(0, ncol(Xfp),ncol(Xfp))
    diag(Sp)[c(rep(FALSE, min(tmp_ind) - 1L), (re$pen.ind == 1))] <- 1 / opt$par[theta_inds[[3]]]
    
    qrx <- qr(rbind(WX, Sp), LAPACK = TRUE)
    Ri <- backsolve(qr.R(qrx), diag(ncol(WX)))
    
    ind <- qrx$pivot;ind[ind] <- 1:length(ind)## qrx$pivot
    Ri <- Ri[ind,] ## unpivoted square root of cov matrix in fitting parameterization Ri Ri' = cov
    
    Vb <- B%*%Ri; Vb <- Vb%*%t(Vb)
    
    ret$spline_coefs <- list(
      beta_spline = beta_spline, Vb = Vb
    )
    ret
  })
  names(cases) <- c("e", "f")
  cases
  })
names(sims) <- paste0("lambda_", sim_params$lambda_interaction)
saveRDS(sims, outfile)
