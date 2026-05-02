\name{tune_lfpls}
\alias{tune_lfpls}

\title{
Tune Longitudinal Functional Partial Least Squares
}

\description{
Selects the tuning parameters of longitudinal functional partial least
squares (LFPLS) by grid search. The function fits LFPLS over candidate
numbers of component-specific PLS directions and smoothing parameters,
then returns the model with the smallest Bayesian information criterion
(BIC).
}

\usage{
tune_lfpls(prep,
           n_comp_b0_grid = 1:5,
           n_comp_b1_grid = 1:5,
           n_comp_u_grid = 1:5,
           lambda_b_grid = c(1e-3, 1e-2),
           max_iter = 6,
           tol = 1e-4,
           ortho = TRUE,
           ridge = 1e-10,
           lme_maxiter = 80)
}

\arguments{
  \item{prep}{
An object returned by \code{\link{prepare_lfpls}}. It contains the prepared
longitudinal functional decomposition, scalar design matrix, integration
weights, initial variance-component estimates, and all quantities required
to fit LFPLS.
}

  \item{n_comp_b0_grid}{
Integer vector giving the candidate numbers of LFPLS components extracted
from the between-subject functional intercept component \eqn{B_{i0}(s)}.
The default is \code{1:5}.
}

  \item{n_comp_b1_grid}{
Integer vector giving the candidate numbers of LFPLS components extracted
from the between-subject functional slope component
\eqn{\tilde t_{ij}B_{i1}(s)}. The default is \code{1:5}.
}

  \item{n_comp_u_grid}{
Integer vector giving the candidate numbers of LFPLS components extracted
from the visit-specific deviation component \eqn{U_{ij}(s)}. The default is
\code{1:5}.
}

  \item{lambda_b_grid}{
Numeric vector of candidate nonnegative smoothing parameters for the
component-specific coefficient functions. The default is
\code{c(1e-3, 1e-2)}.
}

  \item{max_iter}{
Maximum number of outer LFPLS iterations used in each candidate fit. The
default is \code{6}.
}

  \item{tol}{
Convergence tolerance passed to \code{\link{lfpls_iterate}}. The default is
\code{1e-4}.
}

  \item{ortho}{
Logical. If \code{TRUE}, the Krylov-type PLS basis vectors are
orthonormalized when constructing the LFPLS directions. The default is
\code{TRUE}.
}

  \item{ridge}{
Small nonnegative ridge constant used to stabilize residualization and
matrix inversions inside \code{\link{lfpls_iterate}}. The default is
\code{1e-10}.
}

  \item{lme_maxiter}{
Maximum number of iterations passed to the \code{lme} optimizer for
each candidate score mixed model. The default is \code{80}.
}
}

\details{
The function performs automatic tuning for the LFPLS estimator. For each
candidate tuple
\deqn{
(p_0,p_1,p_U,\lambda)
\in
\mathcal P_0 \times \mathcal P_1 \times \mathcal P_U \times \Lambda,
}
where \eqn{\mathcal P_0}, \eqn{\mathcal P_1}, and \eqn{\mathcal P_U} are
the supplied component-number grids and \eqn{\Lambda} is the supplied
smoothing-parameter grid, \code{tune_lfpls} calls
\code{\link{lfpls_iterate}}.

For a fixed candidate tuple, \code{\link{lfpls_iterate}} fits the LFPLS
working model
\deqn{
Y_{ij} =
\gamma^\top z(t_{ij}) + \theta^\top Z_{ij}
+ \langle B_{i0},\beta_0\rangle
+ \tilde t_{ij}\langle B_{i1},\beta_1\rangle
+ \langle U_{ij},\beta_U\rangle
+ b_i+\epsilon_{ij},
}
where \eqn{B_{i0}}, \eqn{B_{i1}}, and \eqn{U_{ij}} are obtained from the
pilot longitudinal functional decomposition stored in \code{prep}. The
parameters \eqn{p_0}, \eqn{p_1}, and \eqn{p_U} determine the dimensions of
the LFPLS Krylov subspaces used for the between-subject intercept,
between-subject slope, and visit-specific deviation blocks, respectively.
The parameter \eqn{\lambda} controls the Sobolev-type smoothing of the
reconstructed coefficient functions.

Each candidate fit produces a REML score mixed model. The BIC of this
mixed model is computed using \code{BIC(fit$fit)}. Candidate fits that fail,
return \code{NULL}, or do not produce a finite BIC are skipped. The selected
model is the successful candidate with the smallest finite BIC:
\deqn{
(\widehat p_0,\widehat p_1,\widehat p_U,\widehat\lambda)
=
\arg\min_{(p_0,p_1,p_U,\lambda)}
\operatorname{BIC}\{ \widehat{\mathcal M}(p_0,p_1,p_U,\lambda) \}.
}
The selected BIC value is stored in the returned object as
\code{best_bic}.

This function is useful when the analyst does not want to manually choose
the number of LFPLS components for each longitudinal functional block or
the smoothing strength for the reconstructed coefficient functions. For
large grids or high-resolution functional predictors, the search may be
computationally intensive because \code{\link{lfpls_iterate}} is refitted
for every candidate combination.
}

\value{
A fitted LFPLS object corresponding to the candidate tuning parameters with
the smallest finite BIC. The returned object is the list produced by
\code{\link{lfpls_iterate}}, with an additional component:

  \item{fit}{
The fitted \code{lme} object from the selected REML score mixed model.
}

  \item{b0}{
Estimated coefficient function \eqn{\hat\beta_0(s)} for the
between-subject functional intercept component.
}

  \item{b1}{
Estimated coefficient function \eqn{\hat\beta_1(s)} for the
between-subject functional slope component.
}

  \item{bU}{
Estimated coefficient function \eqn{\hat\beta_U(s)} for the visit-specific
functional deviation component.
}

  \item{n_comp_b0}{
Selected number of LFPLS components for the \eqn{B_{i0}} block.
}

  \item{n_comp_b1}{
Selected number of LFPLS components for the \eqn{\tilde t_{ij}B_{i1}} block.
}

  \item{n_comp_u}{
Selected number of LFPLS components for the \eqn{U_{ij}} block.
}

  \item{lambda_b}{
Selected smoothing parameter.
}

  \item{iter}{
Number of outer LFPLS iterations completed for the selected fit.
}

  \item{w}{
Numerical integration weights over the functional grid.
}

  \item{sigma_e}{
Final REML estimate of the residual standard deviation for the selected
fit.
}

  \item{sigma_b}{
Final REML estimate of the scalar random-intercept standard deviation for
the selected fit.
}

  \item{best_bic}{
Smallest finite BIC value found over the tuning grid.
}

If all candidate fits fail, the function prints the last error message,
when available, and stops with an error.
}

\note{
The object supplied to \code{prep} must first be created by
\code{\link{prepare_lfpls}}. The total number of candidate fits equals
\code{length(n_comp_b0_grid) * length(n_comp_b1_grid) *
length(n_comp_u_grid) * length(lambda_b_grid)}. For computationally
expensive applications, start with a small grid and then refine the grid
near the best candidate.
}

\seealso{
\code{\link{prepare_lfpls}},
\code{\link{lfpls_iterate}},
\code{\link{simulate_data}}
}

\examples{
\donttest{

## Simulate a small longitudinal scalar-on-function dataset
dat <- simulate_data(
      case = "2",
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

## Tune LFPLS over a small grid
fit_tuned <- tune_lfpls(
      prep = prep,
      n_comp_b0_grid = 1:5,
      n_comp_b1_grid = 1:5,
      n_comp_u_grid  = 1:5,
      lambda_b_grid  = c(1e-2, 5e-2, 1e-1),
      max_iter = 6,
      tol = 1e-4,
      ortho = TRUE,
      ridge = 1e-10,
      lme_maxiter = 80
    )

## Selected tuning parameters
fit_tuned$n_comp_b0
fit_tuned$n_comp_b1
fit_tuned$n_comp_u
fit_tuned$lambda_b
fit_tuned$best_bic

## Estimated component-specific coefficient functions
plot(dat$tps, fit_tuned$b0, type = "l",
     xlab = "Functional domain",
     ylab = "Coefficient",
     main = "Tuned LFPLS: beta_0")
}
}

\keyword{models}
\keyword{regression}
\concept{longitudinal functional data}
\concept{functional partial least squares}
\concept{model selection}
