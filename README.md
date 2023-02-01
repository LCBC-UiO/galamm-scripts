# Supplementary Code

This directory contains supplementary code to the paper ["Longitudinal modeling 
of age-dependent latent traits with generalized additive latent and mixed models"](https://arxiv.org/abs/2105.02488). The actual data cannot be shared due to privacy, so we have instead created simulated datasets resembling the actual data, for which all analyses and simulations can be run. The results will obviously differ from those presented in the paper.

All paths are relative to the
top directory. If `galamm-scripts.Rproj` is opened in RStudio, the top
directory will automatically be the working directory, and all scripts should
run with no modifications. The scripts require R version 4.0.0 or higher; due
to changes in the `approxfun()` function the code will not run with older versions
of R.

The `galamm` package is under rapid development, and hence the code is unlikely to work with the latest version of this package. Instead, use the version corresponding to [the latest commit at the time of writing the paper](https://github.com/LCBC-UiO/galamm/tree/0431f9826776340dd2016e169657e7655f07cc1d). You might first have to install the `remotes` package, and if on Windows, you might need to install Rtools.

```
remotes::install_github("LCBC-UiO/galamm", "0431f9826776340dd2016e169657e7655f07cc1d")
```

There are three subdirectories, `code/`, `data/`, and `results/`. Each of
these again contain two subdirectories. The `code/`
directories are described in the table below. 

| Path | Content | 
| ---- | ----        |
| `code/cognition_model` | Latent response model, Section 5. |
| `code/ses_model` | Latent covariates model, Section 6. |


