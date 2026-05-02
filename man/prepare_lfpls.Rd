\name{prepare_lfpls}
\alias{prepare_lfpls}

\title{
Prepare Data for Longitudinal Functional PLS
}

\description{
Constructs the prepared data object required for longitudinal functional
partial least squares (LFPLS). The function estimates the pilot
longitudinal functional decomposition, builds the scalar and time-effect
design matrices, computes numerical integration weights, fits an initial
random-intercept mixed model, and stores all quantities required by
\code{\link{lfpls_iterate}} and \code{\link{tune_lfpls}}.
}

\usage{
prepare_lfpls(y, subject, times, X, tps,
              scalar_cov = NULL,
              presmooth_X = TRUE,
              df_spline = 10,
              df_time_eta = 5,
              df_time_Z = 4,
              lme_maxiter = 80)
}

\arguments{
  \item{y}{
A numeric vector of scalar responses. Its length must be equal to the
number of rows of \code{X}.
}

  \item{subject}{
A vector of subject identifiers. Observations with the same identifier are
treated as repeated visits from the same subject.
}

  \item{times}{
A numeric vector of visit times. The function internally standardizes these
times before estimating the longitudinal functional decomposition and
constructing the time-effect basis.
}

  \item{X}{
A numeric matrix containing the observed functional predictor. Rows
correspond to subject-visit observations and columns correspond to grid
points over the functional domain.
}

  \item{tps}{
A numeric vector giving the grid points of the functional domain. Its length
must be equal to \code{ncol(X)}. The vector is used for optional
presmoothing and for constructing numerical integration weights.
}

  \item{scalar_cov}{
Optional scalar covariates. This may be \code{NULL}, a numeric vector, a
numeric matrix, or a data frame. If a data frame is supplied, it is expanded
using \code{model.matrix(~ . - 1, data = scalar_cov)}. The default is
\code{NULL}.
}

  \item{presmooth_X}{
Logical. If \code{TRUE}, each observed functional predictor curve is
presmoothed over \code{tps} using \code{smooth.spline}. The default is
\code{TRUE}.
}

  \item{df_spline}{
Degrees of freedom used in \code{smooth.spline} when
\code{presmooth_X = TRUE}. The default is \code{10}.
}

  \item{df_time_eta}{
Degrees of freedom for the spline basis used to estimate the smooth
mean surface \eqn{\eta(s,t)} in the pilot longitudinal functional
decomposition. The default is \code{5}.
}

  \item{df_time_Z}{
Degrees of freedom for the spline basis used to model the smooth fixed
effect of visit time in the scalar response model. The default is \code{4}.
}

  \item{lme_maxiter}{
Maximum number of iterations passed to the \code{lme} optimizer when
fitting the initial random-intercept mixed model. The default is \code{80}.
}
}

\details{
The function prepares the data required to fit LFPLS for longitudinal
scalar-on-function regression. The method assumes that, for subject
\eqn{i} and visit \eqn{j}, the observed functional predictor can be
represented through the longitudinal decomposition
\deqn{
X_{ij}(s) =
\eta(s,t_{ij}) + B_{i0}(s) + t_{ij}B_{i1}(s) + U_{ij}(s),
}
where \eqn{\eta(s,t)} is a smooth population mean surface,
\eqn{B_{i0}(s)} is a subject-specific functional intercept process,
\eqn{B_{i1}(s)} is a subject-specific functional slope process, and
\eqn{U_{ij}(s)} is a visit-specific functional deviation process.

The preparation step proceeds as follows. First, the response, subject
identifiers, visit times, and functional predictor matrix are coerced to
the required formats, and the function checks that \code{nrow(X)} equals
\code{length(y)}. If \code{presmooth_X = TRUE}, each curve is smoothed over
the functional grid \code{tps}.

Second, visit times are standardized and a smooth mean surface
\eqn{\eta(s,t)} is estimated using a spline basis in the visit-time
direction. The centered curves are then decomposed subject by subject into
empirical estimates of the functional intercept \eqn{B_{i0}(s)}, functional
slope \eqn{B_{i1}(s)}, and visit-specific deviation \eqn{U_{ij}(s)}.

Third, numerical integration weights are computed over \code{tps} using the
trapezoidal rule. These weights are used to construct the weighted
functional blocks required by the LFPLS score construction:
\eqn{B_{i0}}, \eqn{\tilde t_{ij}B_{i1}}, and \eqn{U_{ij}}. In the returned
object these are stored as \code{B0_wtd}, \code{B1t_wtd}, and
\code{U_wtd}. The weighted blocks approximate the functional inner products
used in the scalar-on-function model.

Fourth, the scalar fixed-effect design is constructed. A spline basis in
standardized visit time is included to represent the smooth population-level
time effect. Optional scalar covariates supplied through \code{scalar_cov}
are appended to this time-basis design. The resulting matrix is stored as
\code{Z_full}.

Finally, an initial random-intercept mixed model is fitted by REML:
\deqn{
Y_{ij} =
\gamma^\top z(t_{ij}) + \theta^\top Z_{ij}
+ b_i + \epsilon_{ij},
}
where \eqn{b_i} is a subject-specific scalar random intercept. The resulting
initial estimates of the residual and random-intercept standard deviations
are stored as \code{sigma_e0} and \code{sigma_b0}. These estimates are used
to initialize the whitening step in \code{\link{lfpls_iterate}}.

The returned object is not itself a final LFPLS fit. It is a structured
input object for \code{\link{lfpls_iterate}}, which fits LFPLS for fixed
tuning parameters, and \code{\link{tune_lfpls}}, which searches over tuning
parameter grids.
}

\value{
A list with the following components:

  \item{y}{
Numeric vector of scalar responses.
}

  \item{subject}{
Vector of subject identifiers.
}

  \item{subject_factor}{
Factor version of \code{subject}, used in mixed-model fitting.
}

  \item{n}{
Total number of subject-visit observations.
}

  \item{D}{
Number of grid points in the functional predictor.
}

  \item{tstd}{
Standardized visit times.
}

  \item{w}{
Numerical integration weights over the functional grid, computed using the
trapezoidal rule.
}

  \item{Z_full}{
Full scalar fixed-effect design matrix, including the spline basis for
standardized visit time and any additional scalar covariates.
}

  \item{z_names}{
Column names of \code{Z_full}.
}

  \item{df0}{
Data frame used to fit the initial reduced mixed model.
}

  \item{formula0}{
Formula object for the initial reduced mixed model.
}

  \item{fit0}{
Initial \code{lme} fit of the reduced random-intercept mixed model.
}

  \item{sigma_e0}{
Initial REML estimate of the residual standard deviation.
}

  \item{sigma_b0}{
Initial REML estimate of the scalar random-intercept standard deviation.
}

  \item{B0_wtd}{
Weighted functional block corresponding to the subject-specific functional
intercept component \eqn{B_{i0}(s)}.
}

  \item{B1t_wtd}{
Weighted functional block corresponding to the standardized-time-weighted
functional slope component \eqn{\tilde t_{ij}B_{i1}(s)}.
}

  \item{U_wtd}{
Weighted functional block corresponding to the visit-specific deviation
component \eqn{U_{ij}(s)}.
}

  \item{penalty}{
Discrete second-difference penalty matrix used later to smooth the
reconstructed coefficient functions.
}

  \item{decomposition}{
A list containing the pilot functional decomposition, including the
presmoothed/used functional predictor matrix, standardized times, estimated
mean surface, centered curves, subject-specific functional intercepts,
subject-specific functional slopes, and visit-specific deviations.
}
}

\note{
The object returned by \code{prepare_lfpls} is designed to be passed
directly to \code{\link{lfpls_iterate}} or \code{\link{tune_lfpls}}. The
functional predictor matrix \code{X} should be observed on a common dense
grid \code{tps}. If \code{presmooth_X = TRUE}, \code{tps} must be supplied
and must have length equal to \code{ncol(X)}.
}

\seealso{
\code{\link{lfpls_iterate}},
\code{\link{tune_lfpls}},
\code{\link{simulate_data}}
}

\examples{
\donttest{
## Simulate a small longitudinal scalar-on-function dataset
dat <- simulate_data(
  case = "1",
  sigma_x = 0.05,
  sigma_y = 0.35,
  sigma_b = 0.50
)

X <- dat$X
y <- dat$y
Z <- dat$Z
subject <- dat$subject
times <- dat$times
tps <- dat$tps

## Prepare the LFPLS input object
prep <- prepare_lfpls(
  y = y,
  subject = subject,
  times = times,
  X = X,
  tps = tps,
  scalar_cov = Z,
  presmooth_X = TRUE,
  df_spline = 10,
  df_time_eta = 5,
  df_time_Z = 4,
  lme_maxiter = 80
)

names(prep)
prep$n
prep$D

## Initial variance-component estimates
prep$sigma_e0
prep$sigma_b0

## Weighted functional blocks used by LFPLS
dim(prep$B0_wtd)
dim(prep$B1t_wtd)
dim(prep$U_wtd)

## Fit LFPLS using the prepared object
fit <- lfpls_iterate(
  prep = prep,
  n_comp_b0 = 1,
  n_comp_b1 = 1,
  n_comp_u = 1,
  lambda_b = 1e-2,
  max_iter = 2
)

names(fit)
}
}
