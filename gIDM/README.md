# Repository Overview

This repository contains code, processed data, and model outputs for the manuscript:

> *Geometric overdispersion facilitates the integration of ecological data*

The primary workflow is contained in `gIDM_code.Rmd`, which includes data manipulation, simulation studies, and model fitting for the American robin and brook trout analyses.

---

# Repository Structure

## `gIDM_code.Rmd`

Main R Markdown workflow for:

- data preparation,
- simulation studies,
- model fitting,
- posterior analysis,
- model assessment,
- and visualization.

### Sections within `gIDM_code.Rmd`

#### Set up

Loads packages, sources functions, and sets working directory. 

#### Data preparation for brook trout analysis

Code for preparing brook trout covariates and processed data products.

> Raw brook trout data are not publicly distributed because of privacy restrictions.  
> Code is provided for review and reproducibility purposes only.

#### Fit brook trout model

Fits the proposed geometric integrated distribution model to brook trout data.

> This section can be run directly if the processed `data` object has already been created.

#### Assess model fit (RPS)

Computes ranked probability scores (RPS) and other posterior predictive summaries.

#### Plotting

Generates figures for:
- covariate associations,
- spatial intensity surfaces,
- and posterior summaries.

> This section can be run directly if `data` and `fit_list` are available.

#### Simulation study in parallel (geostatistical model)

Conducts the simulation study described in the manuscript using parallel computation.

#### Data preparation for American robin analysis

Prepares data for the integrated distance sampling analysis.

> Portions of this section are adapted from `IDS1_CaseStudy.R` archived at Zenodo:  
> https://doi.org/10.5281/zenodo.10666980

#### Fit models (Kéry et al. 2025 data)

Fits integrated distribution models to the American robin data from Kéry et al. (2025).

#### Visualize/map American robin density in Oregon

Produces spatial visualizations of estimated robin density across Oregon.

#### Compare posterior covariate effects

Compares posterior distributions of regression coefficients across models and data sources.

---

## `algorithms/`

Contains hard-coded MCMC algorithms using Pólya--Gamma augmentation as well as Stan programs used for Hamiltonian Monte Carlo (HMC) fitting.

These files are called throughout `gIDM_code.Rmd`.

---

## `Outputs/`

Storage directory for processed outputs, fitted objects, and simulation study results referenced throughout `gIDM_code.Rmd`.

Simulation study outputs are included for convenience and reproducibility.

Processed data objects required for fitting the integrated distribution models are also stored here.

---

## `Kery_2024/`

Contains data and code associated with the integrated distance sampling model described in:

Kéry et al. (2025). *Ecology*.

The files are included here for convenience but are also publicly archived on Zenodo:

https://doi.org/10.5281/zenodo.10666980

---

# Data Availability

Raw brook trout data are not publicly archived because of privacy restrictions.

Processed and summarized data objects necessary to fit the models presented in the manuscript are provided in `Outputs/`.

Additional raw data may be available from the authors upon reasonable request.