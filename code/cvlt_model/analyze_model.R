library(tidyverse)
library(latex2exp)
options(bitmapType = "cairo")
theme_set(theme_bw())
library(gamm4)

mod <- readRDS("results/cvlt_model/res.rds")
dat <- readRDS("data/cvlt_model/simulated_data.rds")

dat %>% 
    distinct(id, timepoint) %>% 
    summarise(n_distinct(id), n())

cnm <- c(str_replace(colnames(mod$cov_lambda), "a", "lambda"), 
         colnames(mod$cov_beta))
nparams <- nrow(mod$cov_lambda) + nrow(mod$cov_beta)
covmat <- matrix(nrow = nparams, ncol = nparams,
                 dimnames = list(cnm, cnm))
lambda_inds <- str_which(colnames(covmat), "lambda")
beta_inds <- setdiff(1:ncol(covmat), lambda_inds)
covmat[lambda_inds, lambda_inds] <- mod$cov_lambda
covmat[beta_inds, beta_inds] <- mod$cov_beta
covmat[lambda_inds, beta_inds] <- mod$cov_lambda %*% t(mod$beta_deriv)
covmat[beta_inds, lambda_inds] <- t(covmat[lambda_inds, beta_inds])

ggplot(filter(dat, item %in% c("a1", "a3", "a5")), 
       aes(x = age, y = successes, group = id)) +
    geom_point(size = .3) +
    geom_line(alpha = .2) +
    facet_grid(
        cols = vars(item), labeller = 
            as_labeller(function(x) paste("Trial", str_extract(x, "[1-5]")))) +
    xlab("Age") + ylab("Words recalled") +
    theme(
        panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(),
        strip.background = element_blank()
    )

ggsave("results/cvlt_model/figures/cvlt_sample_plot.pdf", 
       width = 16, height = 6, units = "cm")


vars <- c(paste0("a", 2:5), "retest")
bhat <- coef(mod$final_gamm$gam)[vars]
bhat_se_naive <- sqrt(diag(mod$cov_beta_naive[vars, vars]))
bhat_se <- sqrt(diag(covmat[vars, vars]))

tab <- c("\\begin{tabular}{llll}", "\\toprule",
         "Parameter & Estimate & SE & SE-naive  \\\\",
         "\\midrule", "\\multicolumn{4}{l}{\\textit{Item effects}}\\\\")

for(i in 1:4){
    tab <- c(tab, 
             sprintf("\\hspace{3mm} $\\beta_{t%d}$ & %.2f & %.3f & %.3f \\\\", 
                     i+1L, bhat[[i]], bhat_se[[i]], bhat_se_naive[[i]]))
}

tab <- c(tab, 
         "\\midrule",
         "\\multicolumn{4}{l}{\\textit{Retest effect}}\\\\",
         sprintf("\\hspace{3mm} $\\beta_{r}$ & %.2f & %.3f & %.3f \\\\", 
                 bhat[[5]], bhat_se[[5]], bhat_se_naive[[5]]),
         "\\bottomrule",
         "\\end{tabular}")


con <- file("results/cvlt_model/tables/cvlt_parametric_terms.txt")
writeLines(tab, con)
close(con)

# Variance components
getME(mod$final_gamm$mer, "theta")

# Item characteristic curves
lambda_se <- c(0, sqrt(diag(mod$cov_lambda)))
names(lambda_hat) <- str_replace(names(lambda_hat), "a", "lambda")
coefs <- c(lambda_hat, c(a1 = 0, bhat[paste0("a", 2:5)]))
S <- matrix(0, nrow = length(coefs), ncol = length(coefs),
            dimnames = list(names(coefs), names(coefs)))
cnms <- intersect(colnames(S), colnames(covmat))
S[cnms, cnms] <- covmat[cnms, cnms]
grid <- crossing(item = paste0("a", 1:5), eta = seq(from = -3, to = 3, by = .01))
X <- model.matrix(~ 0 + item + item:eta, data = grid)
colnames(X) <- c(paste0("a", 1:5), paste0("lambda", 1:5))
X <- X[, names(coefs)]

df <- grid %>% 
    mutate(
        linpred = as.numeric(X %*% coefs),
        linpred_se = as.numeric(sqrt(diag(X %*% S %*% t(X)))),
        prob = plogis(linpred),
        prob_lower = plogis(linpred + qnorm(.025) * linpred_se),
        prob_upper = plogis(linpred + qnorm(.975) * linpred_se),
        item = factor(item)
    )

getformat <- function(i){
    sprintf("Trial %d: $\\hat{\\lambda}_{%d} = %.2f\\;(%.2f)", 
            i, i, lambda_hat[[i]], lambda_se[[i]])
}

ggplot(df, aes(x = eta, y = prob, group = item,
               ymin = prob_lower, ymax = prob_upper)) +
    geom_line(aes(color = item)) + 
    geom_ribbon(aes(fill = item), alpha = .3) +
    xlab("Episodic memory") +
    ylab("Prob. correct") +
    ggthemes::scale_color_colorblind(
        labels = unname(TeX(map_chr(1:5, getformat)))
    ) +
    ggthemes::scale_fill_colorblind() +
    guides(fill = FALSE) +
    labs(color = NULL, fill = NULL) +
    theme(
        panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(),
        legend.text = element_text(size = 6),
        legend.title = element_text(size = 6),
        legend.key.size = unit(0.1, 'cm'),
        legend.position = c(.77, .26), 
        legend.text.align = 0
    )

ggsave("results/cvlt_model/figures/cvlt_irc.pdf",
       width = 8, height = 6, units = "cm")


# Smooth function
grid <- tibble(
    a2 = 0, a3 = 0, a4 = 0, a5 = 0, retest = 0,
    age_z = seq(from = min(dat$age_z), to = max(dat$age_z), 
                length.out = 1000),
    weight = 1
)

Xp <- predict(mod$final_gamm$gam, newdata = grid, type = "lpmatrix")
cnms <- str_subset(colnames(Xp), "s\\(age\\_z\\)")
Xp <- Xp[, cnms]
S <- covmat[cnms, cnms]
bhat <- coef(mod$final_gamm$gam)[cnms]
v <- diag(Xp %*% S %*% t(Xp))

nmc <- 100000
betas <- rmvn(nmc, bhat, S)
fits <- Xp %*% t(betas)
pred <- c(Xp %*% bhat)
simDev <- fits - pred
absDev <- apply(simDev, 2L, function(x) abs(x / sqrt(v)))
masd <- apply(absDev, 2L, max)
crit <- quantile(masd, prob = .95, type = 8)
cat(crit, file = "results/cvlt_model/critical_value.txt")

df <- grid %>% 
    select(age_z) %>% 
    mutate(
        age = age_z * sd(dat$age) + mean(dat$age),
        f = as.numeric(Xp %*% bhat),
        se = sqrt(v),
        f = f, # Remove item 1 intercept
        f_low_pt = f + qnorm(.025) * se,
        f_high_pt = f + qnorm(.975) * se,
        f_low_sim = f - crit[[1]] * se,
        f_high_sim = f + crit[[1]] * se
    )

ggplot(df, aes(x = age, y = f)) +
    geom_line() +
    geom_ribbon(alpha = .3, aes(ymin = f_low_pt, ymax = f_high_pt)) +
    geom_ribbon(alpha = .3, aes(ymin = f_low_sim, ymax = f_high_sim)) +
    xlab("Age") +
    ylab(TeX("Episodic memory")) +
    theme_bw() +
    theme(
        panel.border = element_blank(),
        axis.line = element_line(),
        panel.grid = element_blank()
    )

ggsave(filename = "results/cvlt_model/figures/univariate_cvlt_smooth.pdf", 
       width = 8, height = 6, units = "cm")


