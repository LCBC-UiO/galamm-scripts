# Supplementary Code

This directory contains supplementary code to the paper "Longitudinal modeling 
of age-dependent latent traits with generalized additive latent and mixed models".

All paths are relative to the
top directory. If `latent-variable-gamm.Rproj` is opened in RStudio, the top
directory will automatically be the working directory, and all scripts should
run with no modifications. The scripts require R version 4.0.0 or higher; due
to changes in the `approxfun()` function the code will not run with older versions
of R.

There are three subdirectories, `code/`, `data/`, and `results/`. Each of
these again contain six subdirectories. The `code/`
directories are described in the table below. 

| Path | Content | 
| ---- | ----        |
| `code/cvlt_model` | Latent response model for episodic memory, Sec. 3.1. |
| `code/cvlt_simulation` | Simulations with latent response model, Sec. 1 of Supplement B. |
| `code/cvlt_dspan_model` | Factor-by-curve interaction model, Sec. 3.2. | 
| `code/cvlt_dspan_simulation` | Simulations with factor-by-curve model, Sec. 2 of Supplement B. |
| `code/ses_model` | Latent covariates model for socioeconomic status and hippocampal volume, Sec. 3.3. |
| `code/ses_simulation` | Simulations with latent covariates model, Sec. 3 of Supplement B. |

The `data/` directories have the following contents. These need not be directly accessed, as the scripts
in the `code/` directories contain paths to the right data.

| Path | Content | 
| ---- | ----        |
| `data/cvlt_model` | Simulated dataset allowing code for latent response model for episodic memory Sec. 3.1. to be run. |
| `data/cvlt_simulation` | Simulation parameters for simulations with latent response model, Sec. 1 of Supplement B.  |
| `data/cvlt_dspan_model` | Simulated dataset allowing code for factor-by-curve interaction model in Sec. 3.2 to be run. | 
| `data/cvlt_dspan_simulation` | Simulation parameters for simulations with factor-by-curve model, Sec. 2 in Supplement B. |
| `data/ses_model` | Simulated dataset allowing code for latent covariates model for socioeconomic status and hippocampal volume Sec. 3.3 to be run. |
| `data/ses_simulation` | Simulation parameters for simulations with latent covariates model, Sec. 3 of Supplement B. |


In the code directories related to the data analyses (`code/cvlt_model`, `code/cvlt_dspan_model`, `code/ses_model`), the file `run_model.R` runs the profile likelihood algorithm which finds the optimal factor loadings, and hence should be run first. Next, the file `finalize_model.R` computes asymptotic covariance matrices of the final fit, and creates a list with all necessary information for analyzing the models. Finally, `analyze_model.R` creates all the figures and tables presented in the manuscript, based on the results saved by `finalize_model.R`.

In the code directories related to simulation studies (`code/cvlt_simulation`, `code/cvlt_dspan_simulation`, `code/ses_simulation`), the files named `run_simulation.R` should be run first, to perform the simulation experiment. Next, the files `analyze_simulation.R` should be run to produce the figures and tables presented in the manuscript. In `code/cvlt_simulation` there are additional files starting with `run_` and `analyze_`, which produce additional simulation results related to variance components and basis dimensions. 
