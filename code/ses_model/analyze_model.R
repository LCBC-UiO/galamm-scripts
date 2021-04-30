library(tidyverse)
library(latex2exp)
options(bitmapType = "cairo")
theme_set(theme_bw())
library(mgcv)

mod <- readRDS("results/ses_model/res.rds")
dat <- readRDS("data/ses_model/simulated_data.rds")

dat %>% 
    filter(str_detect(item, "hippocampus")) %>% 
    ggplot(aes(x = age, y = value, group = id)) +
    geom_line(alpha = .3, size = .2) +
    geom_point(size = .05, aes(color = site)) +
    theme(
        panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(),
        strip.background = element_blank(),
        legend.position = c(.9, .87),
        legend.text = element_text(size = 6),
        legend.title = element_text(size = 6),
        legend.key.size = unit(0.1, 'cm')
    ) +
    xlab("Age") +
    ylab("Hippocampal volume") +
    labs(color = "Scanner") +
    ggthemes::scale_color_colorblind()

ggsave("results/ses_model/figures/hippocampus_plot.pdf",
       width = 8, height = 6, units = "cm")


# Counts for manuscripts
dat %>% 
    ungroup() %>% 
    distinct(id, age, site) %>% 
    group_by(id) %>% 
    mutate(ntp = n()) %>% 
    ungroup() %>% 
    summarise(
        timepoints = n(),
        subjects = n_distinct(id),
        min(age), max(age),
        min(ntp), max(ntp)
    )

dat %>% 
    group_by(id) %>% 
    filter(any(itemgroup %in% c("edu", "income"))) %>% 
    ungroup() %>% 
    filter(item != "hippocampus") %>% 
    distinct(id, item) %>% 
    group_by(id) %>% 
    mutate(items = n()) %>% 
    ungroup() %>% 
    summarise(n_distinct(id), mean(items))

dat %>% 
    filter(itemgroup != "hippocampus") %>% 
    group_by(id) %>% 
    summarise(
        has_edu = any(str_detect(item, "edu")),
        has_income = any(str_detect(item, "income")),
        .groups = "drop"
    ) %>% 
    summarise(
        sum(has_edu), sum(has_income)
    )


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
    item = names(coef(mod$final_gamm$gam)),
    estimate = coef(mod$final_gamm$gam),
    se = sqrt(diag(mod$cov_beta)),
    se_naive = sqrt(diag(mod$cov_beta_naive))
) %>% 
    filter(!str_detect(item, "s\\(age\\_"),
           !str_detect(item, "t2\\(age\\_")) %>% 
    mutate(
        estimate = map2_dbl(estimate, item, scale_beta),
        se = map2_dbl(se, item, scale_beta, loc = FALSE),
        se_naive = map2_dbl(se_naive, item, scale_beta, loc = FALSE),
        units = units[item],
        item_tr = factor(item_translation[item], levels = item_translation),
        tab_string = case_when(
            item %in% c("mother_income", "father_income",
                        "mother_edu", "father_edu") ~ 
                sprintf(
                    "\\hspace{3mm} %s & %.3G & %.2G & %.2G & %s \\\\",
                    item_tr, estimate, se, se_naive, units),
            item == "hippocampus:icv_z" ~ sprintf(
                "\\hspace{3mm} %s & %.1e & %.1e & %.1e & %s \\\\",
                item_tr, estimate, se, se_naive, units),
            TRUE ~ sprintf(
                "\\hspace{3mm} %s & %.3G & %.3G & %.3G & %s \\\\",
                item_tr, estimate, se, se_naive, units)
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

psi1 <- as.numeric(VarCorr(mod$final_gamm$lme)["weight", "StdDev"])

lambda_df <- tibble(
    item = names(lambda_translation),
    estimate = mod$lambda_est[item],
    se = c(NA_real_, sqrt(diag(mod$cov_lambda)))
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
            "\\hspace{3mm} %s & %.3g & %.3g & - & %s \\\\",
            item, estimate, se, units),
        tab_string = str_replace(tab_string, "NA", "-")
    )

psi2 <- as.numeric(VarCorr(mod$final_gamm$lme)["itemgrouphippocampus", "StdDev"]) * 
    transf_df$scale[transf_df$itemgroup == "hippocampus"]

residual_multipliers <- coef(mod$final_gamm$lme$modelStruct$varStruct, 
                             unconstrained = FALSE)
residual_multipliers[setdiff(c("edu", "income", "hippocampus"), 
                             names(residual_multipliers))] <- 1

sigma_edu <- sigma(mod$final_gamm$lme) * 
    transf_df$scale[transf_df$itemgroup == "edu"] *
    residual_multipliers[["edu"]]

sigma_income <- sigma(mod$final_gamm$lme) * 
    transf_df$scale[transf_df$itemgroup == "income"] *
    residual_multipliers[["income"]]

sigma_hippocampus <- sigma(mod$final_gamm$lme) * 
    transf_df$scale[transf_df$itemgroup == "hippocampus"] *
    residual_multipliers[["hippocampus"]]

vc_df <- tibble(
    item = c("Socioeconomic status, $\\surd\\psi_{1}^{(2)}$", 
             "Hippocampus, $\\surd\\psi_{2}^{(2)}$", 
             "Income residual, $\\sigma_{1}$", 
             "Education residual, $\\sigma_{2}$", 
             "Hippocampus residual, $\\sigma_{3}$"),
    estimate = c(psi1, psi2, sigma_income, sigma_edu, sigma_hippocampus),
    units = c("-", "mm$^{3}$", "$\\log(\\text{NOK})$", "$\\log(\\text{years})$", "mm$^{3}$")
) %>% 
    mutate(
        tab_string = sprintf("\\hspace{3mm} %s & %.3g & - & - & %s \\\\",
                             item, estimate, units)
    )

tab <- c(
    "\\begin{tabular}{lllll}",
    "\\toprule",
    "Parameter & Estimate & SE & SE-naive & Units \\\\",
    "\\midrule",
    #"\\multicolumn{5}{l}{\\textit{Item intercepts}}\\\\",
    # beta_df$tab_string[1:6],
    #"\\midrule",
    "\\multicolumn{5}{l}{\\textit{Effects on hippocampal volume}}\\\\",
    beta_df$tab_string[7:11],
    "\\midrule",
    "\\multicolumn{5}{l}{\\textit{Factor loadings}}\\\\",
    lambda_df$tab_string,
    "\\midrule",
    "\\multicolumn{5}{l}{\\textit{Variance components}}\\\\",
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
    "Parameter & Estimate & SE & SE-naive & Units \\\\",
    "\\midrule",
    "\\multicolumn{5}{l}{\\textit{Item intercepts}}\\\\",
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
        eta = sdevs * psi1,
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
        weight = 1, itemgrouphippocampus = 1
    )

inds <- str_detect(names(coef(mod$final_gamm$gam)), "s\\(age")

bhat <- c(lambda7 = mod$lambda_est[["hippocampus"]],
          lambda8 = mod$lambda_est[["interaction"]],
          coef(mod$final_gamm$gam)[inds])

Slambda <- matrix(c(mod$cov_lambda[["hippocampus", "hippocampus"]], 0, 0, 0),
                  nrow = 2)
Sbeta <- mod$cov_beta[inds, inds]
Slambda_beta <- Slambda %*% t(cbind(mod$beta_deriv[inds, 2], 0))
Sbeta_lambda <- t(Slambda_beta)

S <- rbind(
    cbind(Slambda, Slambda_beta),
    cbind(Sbeta_lambda, Sbeta))

Xp <- cbind(
    model.matrix(~ 0 + eta + eta:age_z, data = grid),
    predict(mod$final_gamm$gam, newdata = grid, type = "lpmatrix")[, inds]
)

fit <- Xp %*% bhat
se.fit <- sqrt(rowSums((Xp %*% S) * Xp))

df <- grid %>% 
    select(-starts_with("item")) %>% 
    mutate(
        fit = as.numeric(fit) * 
            transf_df$scale[transf_df$itemgroup == "hippocampus"] +
            transf_df$location[transf_df$itemgroup == "hippocampus"],
        se = as.numeric(se.fit) * 
            transf_df$scale[transf_df$itemgroup == "hippocampus"],
        age = age_z * sd(dat$age) + mean(dat$age),
        sdevs = factor(sdevs),
        sdevs = fct_relabel(
            sdevs, 
            ~ if_else(.x == 0, "Mean", sprintf("Mean %+d SD", as.numeric(.x))))
    ) %>% 
    group_by(sdevs, age) %>% 
    summarise(
        fit = sum(fit),
        se = sqrt(sum(se^2)),
        .groups = "drop"
    ) %>% 
    mutate(se = if_else(sdevs == "Mean", 0, se))

ggplot(df, aes(x = age, y = fit, group = sdevs)) + 
    geom_ribbon(aes(ymin = fit + qnorm(.025) * se,
                    ymax = fit + qnorm(.975) * se, fill = sdevs),
                alpha = .3) +
    geom_line(aes(color = sdevs)) +
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


dat %>% 
    group_by(id) %>% 
    summarise(age = mean(age), .groups = "drop") %>% 
    mutate(ses = ranef(mod$final_gamm$lme)$id$weight) %>% 
    ggplot(aes(x = age, y = ses)) + 
    geom_point()


# Effect of increasing by two years more than mean
mean_edu <- dat %>% 
    filter(itemgroup == "edu") %>% 
    summarise(mean(exp(value))) %>% 
    pull()

# Effect of 1 SD increase in SES on education
exp(2.81 + .169*.67) - exp(2.81)

# Effect of 1 SD increase in SES on income
exp(13.1 + .265*.67) - exp(13.1)

# Effect of 1 SD increase in SES on hippocampus
58.2 * .67

# Find annual rate of decline for hippocampal volume
df %>% 
    filter(sdevs == "Mean") %>% 
    mutate(
        dt = age - lag(age),
        df = fit - lag(fit),
        rate = df / dt
    ) %>% View()
    ggplot(aes(x = age, y = rate)) + 
    geom_line()
