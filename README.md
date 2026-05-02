# LongFPLS <img src="https://img.shields.io/badge/R-%3E=3.5.0-1f425f.svg" alt="R (>= 3.5.0)" align="right" height="20"/>

**LongFPLS** provides tools for fitting **longitudinal functional partial least squares regression** models for **longitudinal scalar-on-function data**. The package is designed for settings where a scalar response and a densely observed functional predictor are repeatedly measured over time for each subject.

The package implements a supervised dimension-reduction framework that decomposes longitudinal functional predictors into subject-specific functional intercept, subject-specific functional slope, and visit-specific deviation components, and then constructs partial least squares directions that are directly associated with the scalar response.

The package includes simulation-based data generation, pilot decomposition of longitudinal functional predictors, residualized supervised partial least squares construction, REML-based mixed-model fitting, smoothing of reconstructed coefficient functions, and BIC-based tuning of LFPLS components.

---

## 🚀 Key Features

- **Longitudinal scalar-on-function regression:** fits regression models where the response is scalar and the predictor is a repeatedly observed function.

- **Longitudinal functional decomposition:** decomposes the functional predictor into
  - a smooth population mean surface,
  - subject-specific functional intercept component,
  - subject-specific functional slope component,
  - visit-specific functional deviation component.

- **Supervised dimension reduction:** constructs functional partial least squares directions that maximize association with the scalar response, rather than merely explaining variation in the functional predictor.

- **Mixed-model-aware estimation:** incorporates within-subject dependence through a scalar random-intercept mixed model and REML-based variance-component estimation.

- **Residualized LFPLS construction:** adjusts for scalar covariates and time effects before extracting supervised functional directions.

- **Component-specific coefficient functions:** estimates separate coefficient functions for the subject-level and visit-level functional components.

- **BIC-based tuning:** provides automatic tuning over the number of LFPLS components and smoothing parameters.

- **Built-in simulation tools:** generates longitudinal scalar-on-function datasets under different functional-effect mechanisms.

- **User manual included:** the package is accompanied by the manual file `LongFPLS_1.0.0.pdf`.

---

## 📦 Installation

You can install the development version of **LongFPLS** from GitHub:

```r
install.packages("remotes")
remotes::install_github("UfukBeyaztas/LongFPLS")

Then load the package:

library(LongFPLS)

📘 Main Functions

The package contains four main user-facing functions:

Function	Description
simulate_data()	Simulates longitudinal scalar-on-function data under different functional-effect scenarios.
prepare_lfpls()	Prepares the data object required for LFPLS fitting.
lfpls_iterate()	Fits LFPLS for fixed tuning parameters.
tune_lfpls()	Selects LFPLS tuning parameters by BIC-based grid search.


📄 Manual

A user manual is included with the package:

LongFPLS_1.0.0.pdf

The manual contains detailed descriptions of the exported functions, their arguments, returned objects, and example workflows.
