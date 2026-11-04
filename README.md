# icSurvML

Methods for Interval-Censored Survival Analysis, Imputation, Conditional
Survival Estimation, and Progressive Illness–Death Models.

The `icSurvML` package provides a unified framework for the analysis of
interval-censored time-to-event data. It includes nonparametric estimation of
survival functions, comparison of survival curves, imputation strategies for
interval-censored event times, conditional survival estimation, and
landmark-based estimators of transition probabilities in progressive
illness–death models.

The package combines classical statistical methodology with modern
machine-learning approaches, allowing researchers to analyse
interval-censored data using flexible and reproducible workflows.

## Main Features

### Survival Analysis under Interval Censoring

- Turnbull estimator for interval-censored survival data (`get_intervals()`, `survival_ic()`);
- Bootstrap confidence intervals;
- Smoothed survival estimation, with survival probabilities redistributed linearly within each Turnbull interval;
- Comparison of survival curves between groups (`test_ic()`).

### Imputation of Event Times

Event times can be imputed with `impute_ic()` and then analysed with the
Kaplan–Meier estimator (`fit_ic()`). The package implements the following approaches:

- Classical imputation: left-point, midpoint, and right-point;
- Parametric Accelerated Failure Time (AFT) models: log-normal, Weibull, exponential, and Gaussian;
- Machine-learning methods:
  - Extreme Gradient Boosting (XGBoost);
  - Random Survival Forests (RSF);
  - Survival Support Vector Machines (SVM);
  - Generalized Linear Model Boosting (GLMBoost);
  - Generalized Additive Model Boosting (GAMBoost).

All model-based procedures ensure that imputed event times remain within the
observed intervals \((L_i,R_i]\), through the Scaled Linear Redistribution
Method (SLRM), which preserves the variability of the model-based predictions.

### Conditional Survival Estimation

- `cond_surv_ic()`: survival conditional on a continuous covariate, using a
  kernel-weighted conditional Turnbull estimator (`method = "Turnbull"`) or a
  midpoint-based Beran estimator (`method = "Beran"`);
- `cond_aZb_surv_ic()`: survival of a terminal event conditional on an
  interval-censored intermediate event falling within a time window \((a,b]\),
  using a Turnbull-based reconstruction estimator (`method = "Turnbull"`) or a
  midpoint-based benchmark (`method = "Beran"`).

### Progressive Illness–Death Models

For multi-state settings with interval-censored intermediate events, the
package implements landmark estimators of transition probabilities
(`tp_ic()`):

- Interpolated Aalen–Johansen (IAJ);
- Interval Landmark Method (ILM);
- Landmark Turnbull estimator (LMTB).

These methods extend landmarking ideas to progressive illness–death models
where the transition from the healthy state to the illness state is
interval-censored.

## Data

The package includes the `breastCancer` dataset, a semi-synthetic dataset based
on the German Breast Cancer Study Group data, with interval-censored recurrence
times.

## Installation

The release version of the package is available on CRAN:

```r
install.packages("icSurvML")
```

The development version can be installed from GitHub:

```r
# install.packages("remotes")
remotes::install_github("gsoutinho/icSurvML")
```

## Authors and Maintainer

- Gustavo Soutinho
- Luís Meira-Machado (lmachado@math.uminho.pt)

Maintainer: Gustavo Soutinho (gustavosoutinho@sapo.pt)