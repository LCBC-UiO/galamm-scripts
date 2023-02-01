library(galamm)
library(tidyverse)
theme_set(theme_bw())
theme_update(panel.grid = element_blank())
library(gamm4)

dat <- readRDS("data/ses_model/dat.rds")

# Go on analyzing the model with lowest AIC
cases <- sort(
  str_extract(list.files("results/ses_model/", pattern = "case_[a-z]+\\.rds"), 
              "(?<=case_)[:alpha:]"))
models <- map(cases, function(case){
  readRDS(paste0("results/ses_model/case_", case, ".rds"))
})
names(models) <- cases

case <- "f"
mod <- models$f

# Compute p-value for income loading
mod0 <- readRDS("results/ses_model/case_f_no_income_loading.rds")
1 - pchisq(2 * (mod$opt$value - mod0$value), 1)

beta_spline <- as.numeric(mod$spline_coefs$beta_spline)
Vb <- as.matrix(mod$spline_coefs$Vb)

# Proportion of structural zeroes:
1 - length(mod$reTrms$Zt@x) / prod(dim(mod$reTrms$Zt))
1 - length(mod$reTrms$Lambdat@x) / prod(dim(mod$reTrms$Lambdat))

# Comparison with the full model
comp_df <- imap_dfr(models, function(mod, name){
  tibble(
    name = name,
    df = length(mod$opt$par) + 1L,
    log_likelihood = mod$final_model$logLik,
    aic = -2 * log_likelihood + 2 * df
  )
}) %>% 
  arrange(desc(df)) %>% 
  mutate(
    name = recode(
      name,
      a = "(a): Free loadings",
      b = "(b): (a) and no interaction, $\\lambda_{8}=0$",
      c = "(c): Parents equal, $\\lambda_{1}=\\lambda_{2}$ and $\\lambda_{4}=\\lambda_{5}$",
      d = "(d): (c) and no interaction, $\\lambda_{8}=0$",
      e = "(e): Item groups equal, $\\lambda_{1}=\\lambda_{2}=\\lambda_{3}$ and $\\lambda_{4}=\\lambda_{5}=\\lambda_{6}$",
      f = "(f): (e) and no interaction, $\\lambda_{8}=0$",
      g = "(g): (f) and no main effect, $\\lambda_{7}=\\lambda_{8}=0$"
    ),
    aic = aic - max(aic)
  )

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

comp_df %>% 
  filter(df <= 21) %>% 
  summarise(likelihood_ratio = 2 * (last(log_likelihood) - first(log_likelihood)))

dat %>% 
  filter(str_detect(item, "hippocampus")) %>% 
  ggplot(aes(x = age, y = value, group = id)) +
  geom_line(alpha = .3, size = .2) +
  geom_point(size = .05, aes(color = site)) +
  theme(
    
    legend.position = c(.9, .87),
    legend.text = element_text(size = 6),
    legend.title = element_text(size = 6),
    legend.key.size = unit(0.1, 'cm'),
    panel.border = element_blank(),
    axis.line = element_line()
  ) +
  xlab("Age") +
  ylab("Hippocampal volume") +
  labs(color = "Scanner") +
  ggthemes::scale_color_colorblind()

ggsave("results/ses_model/figures/hippocampus_plot.pdf",
       width = 8, height = 6, units = "cm")

# Create table of coefficients
transf_df <- dat %>% 
  group_by(itemgroup) %>% 
  summarise(
    location = mean(value),
    scale = sd(value),
    .groups = "drop"
  )

scale_beta <- function(x, y, loc = TRUE){
  if(str_detect(y, "income")){
    as.numeric(
      (loc == TRUE) * transf_df[transf_df$itemgroup == "income", "location"] +
        transf_df[transf_df$itemgroup == "income", "scale"] * x)
  } else if (str_detect(y, "edu")){
    as.numeric(
      (loc == TRUE) * transf_df[transf_df$itemgroup == "edu", "location"] +
        transf_df[transf_df$itemgroup == "edu", "scale"] * x)
  } else
    if (str_detect(y, "hippocampus") & !str_detect(y, "icv")) {
      as.numeric(transf_df[transf_df$itemgroup == "hippocampus", "scale"] * x)
    } else if (str_detect(y, "hippocampus") & str_detect(y, "icv")){
      as.numeric(transf_df[transf_df$itemgroup == "hippocampus", "scale"] * x / 
                   sd(dat$icv))
    } else {
      as.numeric(x)
    }
}

item_translation <- c(
  itemedu_father = "Father's education, $\\beta_{s1}$",
  itemedu_mother = "Mother's education, $\\beta_{s2}$",
  itemedu_self = "Education, $\\beta_{s3}$",
  itemincome_father = "Father's income, $\\beta_{s4}$",
  itemincome_mother = "Mother's income, $\\beta_{s5}$",
  itemincome_self = "Income, $\\beta_{s6}$",    
  `itemgrouphippocampus:siteousAvanto` = "Scanner ousAvanto, $\\beta_{h1}$",
  `itemgrouphippocampus:siteousPrisma` = "Scanner ousPrisma, $\\beta_{h2}$",
  `itemgrouphippocampus:siteousSkyra` = "Scanner ousSkyra, $\\beta_{h3}$",
  `itemgrouphippocampus:icv_z` = "Total intracranial volume, $\\beta_{h4}$",
  `itemgrouphippocampus:sexmale` = "Sex=Male, $\\beta_{h5}$"
)

units <- c(
  itemincome_mother = "$\\log(\\text{NOK})$",
  itemincome_father = "$\\log(\\text{NOK})$",
  itemincome_self = "$\\log(\\text{NOK})$",
  itemedu_mother = "$\\log(\\text{years})$",
  itemedu_father = "$\\log(\\text{years})$",
  itemedu_self = "$\\log(\\text{years})$",
  `itemgrouphippocampus:siteousAvanto` = "mm$^{3}$",
  `hippocampus:siteousAvanto` = "mm$^{3}$",
  `itemgrouphippocampus:siteousPrisma` = "mm$^{3}$",
  `itemgrouphippocampus:siteousSkyra` = "mm$^{3}$",
  `itemgrouphippocampus:icv_z` = "mm$^{3}$/mm$^{3}$",
  `itemgrouphippocampus:sexmale` = "mm$^{3}$"
)

beta_df <- tibble(
  item = colnames(mod$X),
  estimate = mod$opt$par[mod$beta_inds],
  se = sqrt(diag(mod$S)[seq_along(mod$beta_inds)])
) %>% 
  filter(!str_detect(item, "^Xf")) %>% 
  mutate(
    estimate = map2_dbl(estimate, item, scale_beta),
    se = map2_dbl(se, item, scale_beta, loc = FALSE),
    units = units[item],
    item_tr = factor(item_translation[item], levels = item_translation),
    tab_string = case_when(
      item %in% c("mother_income", "father_income",
                  "mother_edu", "father_edu") ~ 
        sprintf(
          "\\hspace{3mm} %s & %.3G & %.2G & %s \\\\",
          item_tr, estimate, se, units),
      item == "hippocampus:icv_z" ~ sprintf(
        "\\hspace{3mm} %s & %.1e & %.1e & %s \\\\",
        item_tr, estimate, se, units),
      TRUE ~ sprintf(
        "\\hspace{3mm} %s & %.3G & %.3G & %s \\\\",
        item_tr, estimate, se, units)
    )
  )

lambda_translation <- c(
  edu_father = "Education, $\\lambda_{1}=\\lambda_{2}=\\lambda_{3}$",
  income_father = "Income, $\\lambda_{4}=\\lambda_{5}=\\lambda_{6}$",
  hippocampus = "Hippocampus, $\\lambda_{7}$"
)

lambda_units <- c(
  edu_father = "$\\log(\\text{years})$",
  edu_mother = "$\\log(\\text{years})$",
  edu_self = "$\\log(\\text{years})$",
  income_father = "$\\log(\\text{NOK})$",
  income_mother = "$\\log(\\text{NOK})$",
  income_self = "$\\log(\\text{NOK})$",
  hippocampus = "mm$^{3}$",
  interaction = "mm$^{3}$/year"
)

lambda_df <- tibble(
  item = names(lambda_translation),
  estimate = c(1, mod$opt$par[mod$lambda_inds]),
  se = c(NA_real_, sqrt(diag(mod$S)[length(mod$beta_inds) + seq_along(mod$lambda_inds)]))
) %>% 
  mutate(
    across(c(estimate, se), ~ case_when(
      str_detect(item, "income") ~ . * 
        transf_df$scale[transf_df$itemgroup == "income"],
      str_detect(item, "edu") ~ . * 
        transf_df$scale[transf_df$itemgroup == "edu"],
      item == "hippocampus" ~ . * 
        transf_df$scale[transf_df$itemgroup == "hippocampus"],
      item == "interaction" ~ . * 
        transf_df$scale[transf_df$itemgroup == "hippocampus"] /
        sd(dat$age),
      TRUE ~ NA_real_
    )),
    units = lambda_units[item],
    item = lambda_translation[item],
    tab_string = sprintf(
      "\\hspace{3mm} %s & %.3g & %.3g &  %s \\\\",
      item, estimate, se, units)
  )

vc <- list()

# Standard deviation of SES latent variable
vc$psi1 <- mod$opt$par[mod$theta_inds][[1]] * sqrt(mod$final_model$phi)
# Standard deviation of hippocampus random effect
vc$psi2 <- mod$opt$par[mod$theta_inds][[2]] * 
  transf_df$scale[transf_df$itemgroup == "hippocampus"] * 
  sqrt(mod$final_model$phi)

vc$sigma_edu <- sqrt(mod$final_model$phi / mod$opt$par[mod$weights_inds][[1]]) * 
  transf_df$scale[transf_df$itemgroup == "edu"]

vc$sigma_income <- sqrt(mod$final_model$phi) * 
  transf_df$scale[transf_df$itemgroup == "income"]

vc$sigma_hippocampus <- sqrt(mod$final_model$phi / mod$opt$par[mod$weights_inds][[2]]) * 
  transf_df$scale[transf_df$itemgroup == "hippocampus"]

vc_df <- tibble(
  item = c("Socioeconomic status, $\\surd\\psi_{1}$", 
           "Hippocampus, $\\surd\\psi_{2}$", 
           "Income residual, $\\sigma_{1}$", 
           "Education residual, $\\sigma_{2}$", 
           "Hippocampus residual, $\\sigma_{3}$"),
  estimate = with(vc, {c(psi1, psi2, sigma_income, sigma_edu, sigma_hippocampus)}),
  units = c("-", "mm$^{3}$", "$\\log(\\text{NOK})$", "$\\log(\\text{years})$", "mm$^{3}$")
) %>% 
  mutate(
    tab_string = sprintf("\\hspace{3mm} %s & %.3g & - &  %s \\\\",
                         item, estimate, units)
  )

tab <- c(
  "\\begin{tabular}{lllll}",
  "\\toprule",
  "Parameter & Estimate & SE & Units \\\\",
  "\\midrule",
  "\\multicolumn{4}{l}{\\textit{Effects on hippocampal volume}}\\\\",
  beta_df$tab_string[7:11],
  "\\midrule",
  "\\multicolumn{4}{l}{\\textit{Factor loadings}}\\\\",
  lambda_df$tab_string,
  "\\midrule",
  "\\multicolumn{4}{l}{\\textit{Variance components}}\\\\",
  vc_df$tab_string,
  "\\bottomrule",
  "\\end{tabular}"
)

con <- file("results/ses_model/tables/ses_param_table.txt")
writeLines(tab, con = con)
close(con)

supptab <- c(
  "\\begin{tabular}{lllll}",
  "\\toprule",
  "Parameter & Estimate & SE & Units \\\\",
  "\\midrule",
  "\\multicolumn{4}{l}{\\textit{Item intercepts}}\\\\",
  beta_df$tab_string[1:6],
  "\\bottomrule",
  "\\end{tabular}"
)

con <- file("results/ses_model/tables/ses_param_supp.txt")
writeLines(supptab, con = con)
close(con)

# Now compute the smooth term estimates
# This requires combining the betas with the lambdas
grid <- crossing(
  sdevs = c(-2, 0, 2),
  age_z = seq(from = min(dat$age_z), to = max(dat$age_z), 
              length.out = 100),
  item = "hippocampus"
) %>% 
  mutate(
    eta = sdevs * vc$psi1,
    itemedu_father = 0,         
    itemedu_mother = 0, 
    itemedu_self = 0,
    itemincome_father = 0, 
    itemincome_mother = 0, 
    itemincome_self = 0,
    siteousAvanto = 0, 
    siteousPrisma = 0,
    siteousSkyra = 0,
    icv_z = 0, sexmale = 0,
    lambda = 1, itemgrouphippocampus = 1
  )

Xp <- PredictMat(mod$sm, data = grid)
Xp <- cbind(Xp, grid$eta)
covmat <- rbind(cbind(Vb[12:26, 12:26], 0),
                cbind(matrix(0, ncol = 15), mod$S[length(mod$beta_inds) + seq_along(mod$lambda_inds)[[2]], 
                                                  length(mod$beta_inds) + seq_along(mod$lambda_inds)[[2]]]))

bhat <- c(beta_spline, mod$opt$par[mod$lambda_inds[[2]]])
grid$pred <- as.numeric(Xp %*% bhat)
grid$pred_se <- as.numeric(sqrt(rowSums((Xp %*% covmat) * Xp)))

df <- grid %>% 
  select(-starts_with("item"), -starts_with("site"), -icv_z) %>% 
  mutate(
    pred = as.numeric(pred) * 
      transf_df$scale[transf_df$itemgroup == "hippocampus"] +
      transf_df$location[transf_df$itemgroup == "hippocampus"],
    pred_se = as.numeric(pred_se) * 
      transf_df$scale[transf_df$itemgroup == "hippocampus"],
    age = age_z * sd(dat$age) + mean(dat$age),
    sdevs = factor(sdevs),
    sdevs = fct_relabel(
      sdevs, 
      ~ if_else(.x == 0, "Mean", sprintf("Mean %+d SD", as.numeric(.x)))),
    pred_se = if_else(sdevs == "Mean", 0, pred_se)
  )

ggplot(df, aes(x = age_z, y = pred, group = factor(sdevs))) + 
  geom_line() + 
  geom_ribbon(aes(ymin = pred + qnorm(.025) * pred_se,
                  ymax = pred + qnorm(.975) * pred_se, fill = factor(sdevs)),
              alpha = .3) +
  geom_line(aes(color = factor(sdevs))) +
  ggthemes::scale_color_colorblind() +
  ggthemes::scale_fill_colorblind() +
  labs(color = "SES", fill = "SES") +
  theme(
    panel.grid = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(),
    legend.position = c(.2, .2),
    legend.text = element_text(size = 6),
    legend.title = element_text(size = 6),
    legend.key.size = unit(0.2, 'cm')
  ) +
  ylab("Hippocampal volume") +
  xlab("Age")

ggsave("results/ses_model/figures/ses_smooth.pdf",
       width = 8, height = 6, units = "cm")

