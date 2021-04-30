library(mgcv)
library(tidyverse)
library(numDeriv)
options(bitmapType = "cairo")

source("code/ses_model/model_setup.R")
dat <- readRDS("data/ses_model/simulated_data.rds")

# First do a model comparison
models <- list(
    "itemgroups_equal" = readRDS("results/ses_model/model_fits/opt_itemgroups_equal.rds"),
    "itemgroups_equal_no" = readRDS("results/ses_model/model_fits/opt_itemgroups_equal_no_intrx.rds"),
    "itemgroups_equal_nono" = readRDS("results/ses_model/model_fits/opt_itemgroups_equal_nono_intrx.rds"),
    "parents_equal" = readRDS("results/ses_model/model_fits/opt_parents_equal.rds"),
    "parents_equal_no" = readRDS("results/ses_model/model_fits/opt_parents_equal_no_intrx.rds"),
    "free_loadings" = readRDS("results/ses_model/model_fits/opt_free_loadings.rds"),
    "free_loadings_no" = readRDS("results/ses_model/model_fits/opt_free_loadings_no_intrx.rds")
)

# Compute fit for free loadings model to get the number of degrees of freedom
# for the parameters fitted by mgcv::gam
lambda_est <- c(edu_father = 1, models[["free_loadings"]]$par)
ref_mod <- modcall(formula, random, dat, weights, lambda_est)
ref_df <- length(ref_mod$lme$fixDF$X)

# Comparison with the full model
comp_df <- imap_dfr(models, function(mod, name){
    tibble(
        name = name,
        df = ref_df + length(mod$par),
        log_likelihood = -mod$value,
        aic = -2 * log_likelihood + 2 * df
    )
}) %>% 
    arrange(desc(df)) %>% 
    mutate(
        name = recode(
            name,
            free_loadings = "(a): Free loadings",
            free_loadings_no = "(b): (a) and no interaction, $\\lambda_{8}=0$",
            parents_equal = "(c): Parents equal, $\\lambda_{1}=\\lambda_{2}$ and $\\lambda_{4}=\\lambda_{5}$",
            parents_equal_no = "(d): (c) and no interaction, $\\lambda_{8}=0$",
            itemgroups_equal = "(e): Item groups equal, $\\lambda_{1}=\\lambda_{2}=\\lambda_{3}$ and $\\lambda_{4}=\\lambda_{5}=\\lambda_{6}$",
            itemgroups_equal_no = "(f): (e) and no interaction, $\\lambda_{8}=0$",
            itemgroups_equal_nono = "(g): (f) and no main effect, $\\lambda_{7}=\\lambda_{8}=0$"
            ),
        aic = aic - max(aic)
    )
comp_df

# Write model comparison to table
tab <- c(
    "\\begin{tabular}{llll}",
    "\\toprule",
    "Model & DF & Log-likelihood & AIC \\\\",
    "\\midrule",
    pmap_chr(comp_df, function(name, df, log_likelihood, aic){
        sprintf("%s & %d & %.f & %.2f \\\\", name, df, log_likelihood, aic)
    }),
    "\\bottomrule",
    "\\end{tabular}"
)

con <- file("results/ses_model/tables/ses_aic_table.txt")
writeLines(tab, con = con)
close(con)




opt <- models[["Group loadings equal no interaction"]]

loadings <- "itemgroups_equal"
inds_list <- indsfun(loadings)
fixed_inds <- inds_list$fixed_inds
equal_inds <- inds_list$equal_inds
lambda_est <- numeric(length(lambda_names))
names(lambda_est) <- lambda_names
zero_inds <- 8L
lambda_est[fixed_inds] <- 1
lambda_est[names(opt$par)] <- opt$par
lambda_est <- constrain_loadings(lambda_est, equal_inds, zero_inds)

final_gamm <- modcall(formula, random, dat, weights, lambda_est)

beta_hat <- coef(final_gamm$gam)
cov_lambda <- solve(opt$hessian)
cov_beta_naive <- vcov(final_gamm$gam)

beta_deriv <- jacobian(function(x){
    ll <- lambda_est
    ll[names(x)] <- x
    ll <- constrain_loadings(ll, equal_inds, zero_inds)
    
    mod <- modcall(formula, random, dat, weights, ll)
    coef(mod$gam)
    }, x = opt$par, method = "simple")

cov_beta <- cov_beta_naive + beta_deriv %*% cov_lambda %*% t(beta_deriv)

res <- list(
    final_gamm = final_gamm,
    beta_deriv = beta_deriv,
    cov_beta_naive = cov_beta_naive,
    cov_beta = cov_beta,
    cov_lambda = cov_lambda,
    lambda_est = lambda_est
    )

saveRDS(res, "results/ses_model/res.rds")
