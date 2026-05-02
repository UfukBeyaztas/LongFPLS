\name{lfpls_iterate}
\alias{lfpls_iterate}

\title{
Fit Longitudinal Functional Partial Least Squares
}

\description{
Fits the longitudinal functional partial least squares (LFPLS) model for
a fixed choice of tuning parameters. The function takes an object produced
by \code{\link{prepare_lfpls}}, extracts supervised component-specific
PLS directions, refits the score mixed model by REML, and reconstructs
the component-specific coefficient functions.
}

\usage{
lfpls_iterate(prep,
              n_comp_b0 = 3,
              n_comp_b1 = 3,
              n_comp_u = 3,
              ortho = TRUE,
              ridge = 1e-10,
              lambda_b = 1e-2,
              max_iter = 8,
              tol = 1e-4,
              lme_maxiter = 80)
}

\arguments{
  \item{prep}{
An object returned by \code{\link{prepare_lfpls}}. It contains the scalar
response, subject identifiers, standardized visit times, scalar covariate
design matrix, numerical integration weights, initial variance-component
estimates, and the estimated longitudinal functional decomposition.
}

  \item{n_comp_b0}{
A positive integer giving the number of LFPLS components extracted from
the between-subject functional intercept component \eqn{B_{i0}(s)}.
The default is \code{3}.
}

  \item{n_comp_b1}{
A positive integer giving the number of LFPLS components extracted from
the between-subject functional slope component \eqn{t_{ij}B_{i1}(s)}.
The default is \code{3}.
}

  \item{n_comp_u}{
A positive integer giving the number of LFPLS components extracted from
the visit-specific deviation component \eqn{U_{ij}(s)}. The default is
\code{3}.
}

  \item{ortho}{
Logical. If \code{TRUE}, the Krylov-type PLS basis vectors are
orthonormalized before constructing the score variables. The default is
\code{TRUE}.
}

  \item{ridge}{
A small nonnegative ridge constant used to stabilize matrix inversions
during residualization against scalar covariates and previously extracted
score components. The default is \code{1e-10}.
}

  \item{lambda_b}{
Nonnegative smoothing parameter used in the discrete Sobolev-type
regularization of the reconstructed coefficient functions. Larger values
produce smoother coefficient estimates. The default is \code{1e-2}.
}

  \item{max_iter}{
Maximum number of outer LFPLS iterations. Each iteration updates the
whitened PLS directions, the REML mixed-model fit, the variance components,
and the reconstructed coefficient functions. The default is \code{8}.
}

  \item{tol}{
Convergence tolerance for the weighted squared change in the reconstructed
coefficient functions. The default is \code{1e-4}.
}

  \item{lme_maxiter}{
Maximum number of iterations passed to the \code{lme} optimizer
when fitting the score mixed model. The default is \code{80}.
}
}

\details{
The function implements the fixed-tuning LFPLS fitting step for
longitudinal scalar-on-function regression. The underlying working model
is based on the decomposition
\deqn{
X_{ij}(s) = \eta(s,t_{ij}) + B_{i0}(s) + t_{ij}B_{i1}(s) + U_{ij}(s),
}
where \eqn{\eta(s,t)} is a smooth mean surface, \eqn{B_{i0}(s)} and
\eqn{B_{i1}(s)} are subject-specific functional intercept and slope
processes, and \eqn{U_{ij}(s)} is a visit-specific deviation process.

After preprocessing by \code{\link{prepare_lfpls}}, LFPLS fits the
component-wise mixed model
\deqn{
Y_{ij} =
\gamma^\top z(t_{ij}) + \theta^\top Z_{ij}
+ \langle B_{i0},\beta_0\rangle
+ \tilde t_{ij}\langle B_{i1},\beta_1\rangle
+ \langle U_{ij},\beta_U\rangle
+ b_i + \epsilon_{ij}.
}
The inferential target is the triple of component-specific coefficient
functions \eqn{(\beta_0,\beta_1,\beta_U)}.

At each outer iteration, the function first forms the stacked functional
design based on the weighted empirical components
\eqn{B_{i0}}, \eqn{\tilde t_{ij}B_{i1}}, and \eqn{U_{ij}}. It then whitens
the response, scalar covariates, and functional blocks using the current
random-intercept covariance estimates. For subject \eqn{i}, the working
covariance has the form
\deqn{
V_i = \sigma_e^2 I_{J_i} + \sigma_b^2 {\bf 1}_{J_i}{\bf 1}_{J_i}^\top.
}
This whitening step accounts for within-subject correlation when extracting
the supervised functional directions.

For the between-subject blocks, residualized alternative partial least
squares (RAPLS) directions are constructed after projecting both the
response and functional blocks away from the scalar covariate design. The
directions are generated through a Krylov-type recursion based on the
residualized empirical covariance and cross-covariance operators. In the
discretized implementation, this is equivalent to repeatedly applying
cross-products of the residualized functional design to the residualized
response.

The visit-specific block \eqn{U_{ij}} is treated sequentially. Its PLS
directions are extracted after residualizing against both the scalar
covariates and the already extracted between-subject score variables. This
prevents the within-subject component from re-capturing predictive signal
already explained by the subject-level intercept and slope components.

Given the extracted score variables, the function fits a linear mixed model
with a scalar random intercept by REML using \code{lme}. The fixed
effects associated with the LFPLS score variables are then used to
reconstruct the coefficient functions:
\deqn{
\hat\beta_0(s) =
\sum_{\ell=1}^{p_0}\hat\delta_{0\ell}\psi^{(0)}_\ell(s),
\quad
\hat\beta_1(s) =
\sum_{\ell=1}^{p_1}\hat\delta_{1\ell}\psi^{(1)}_\ell(s),
\quad
\hat\beta_U(s) =
\sum_{\ell=1}^{p_U}\hat\delta_{U\ell}\psi^{(U)}_\ell(s).
}
The raw reconstructed coefficient vectors are smoothed using a discrete
second-difference penalty, corresponding to a Sobolev-type roughness
regularization. A diminishing-step averaging update is used to stabilize
successive coefficient-function estimates.

This function is intended for fitting LFPLS with fixed values of
\code{n_comp_b0}, \code{n_comp_b1}, \code{n_comp_u}, and \code{lambda_b}.
For automatic tuning over grids of these parameters, use
\code{\link{tune_lfpls}}.
}

\value{
A list containing the fitted LFPLS model and associated estimates:

  \item{fit}{
The fitted \code{lme} object from the final REML score mixed model.
}

  \item{b0}{
Numeric vector containing the estimated coefficient function
\eqn{\hat\beta_0(s)} for the between-subject functional intercept component,
evaluated on the functional grid used in \code{\link{prepare_lfpls}}.
}

  \item{b1}{
Numeric vector containing the estimated coefficient function
\eqn{\hat\beta_1(s)} for the between-subject functional slope component,
evaluated on the functional grid used in \code{\link{prepare_lfpls}}.
}

  \item{bU}{
Numeric vector containing the estimated coefficient function
\eqn{\hat\beta_U(s)} for the visit-specific deviation component, evaluated
on the functional grid used in \code{\link{prepare_lfpls}}.
}

  \item{n_comp_b0}{
Number of LFPLS components used for the \eqn{B_{i0}} block.
}

  \item{n_comp_b1}{
Number of LFPLS components used for the \eqn{t_{ij}B_{i1}} block.
}

  \item{n_comp_u}{
Number of LFPLS components used for the \eqn{U_{ij}} block.
}

  \item{lambda_b}{
Smoothing parameter used for the coefficient-function regularization.
}

  \item{iter}{
Number of outer LFPLS iterations completed.
}

  \item{w}{
Numerical integration weights over the functional grid.
}

  \item{sigma_e}{
Final REML estimate of the residual standard deviation.
}

  \item{sigma_b}{
Final REML estimate of the random-intercept standard deviation.
}

If no score mixed model can be fitted successfully, the returned object may
be \code{NULL}.
}

\note{
The input object \code{prep} must be created by \code{\link{prepare_lfpls}}.
The estimated coefficient vectors \code{b0}, \code{b1}, and \code{bU} are
returned on the same functional grid used to construct the input functional
predictor matrix. The function uses a scalar random-intercept mixed model
for the response-side dependence structure.
}

\seealso{
\code{\link{prepare_lfpls}},
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

## Fit LFPLS with fixed tuning parameters
fit <- lfpls_iterate(
  prep = prep,
  n_comp_b0 = 1,
  n_comp_b1 = 1,
  n_comp_u = 1,
  lambda_b = 1e-2,
  max_iter = 2
)

names(fit)
fit$iter
fit$sigma_e
fit$sigma_b

## Estimated component-specific coefficient functions
plot(dat$tps, fit$b0, type = "l",
     xlab = "Functional domain",
     ylab = "Coefficient",
     main = "Estimated beta_0")
}
}
