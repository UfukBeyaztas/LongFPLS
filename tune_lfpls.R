tune_lfpls <- function(prep,
                       n_comp_b0_grid = 1:5,
                       n_comp_b1_grid = 1:5,
                       n_comp_u_grid = 1:5,
                       lambda_b_grid = c(1e-3, 1e-2),
                       max_iter = 6,
                       tol = 1e-4,
                       ortho = TRUE,
                       ridge = 1e-10,
                       lme_maxiter = 80) {
  
  best <- NULL
  best_bic <- Inf
  last_err <- NULL
  
  for (lam in lambda_b_grid) {
    for (p0 in n_comp_b0_grid) {
      for (p1 in n_comp_b1_grid) {
        for (pU in n_comp_u_grid) {
          fit <- tryCatch(
            lfpls_iterate(
              prep = prep,
              n_comp_b0 = p0,
              n_comp_b1 = p1,
              n_comp_u = pU,
              ortho = ortho,
              ridge = ridge,
              lambda_b = lam,
              max_iter = max_iter,
              tol = tol,
              lme_maxiter = lme_maxiter
            ),
            error = function(e) { last_err <<- e; NULL }
          )
          if (is.null(fit) || is.null(fit$fit)) next
          
          bic <- tryCatch(BIC(fit$fit), error = function(e) Inf)
          if (is.finite(bic) && bic < best_bic) {
            best_bic <- bic
            best <- fit
            best$best_bic <- best_bic
          }
        }
      }
    }
  }
  
  if (is.null(best)) {
    cat("tune_lfpls: all fits failed.\n")
    if (!is.null(last_err)) cat("Last error:\n", conditionMessage(last_err), "\n")
    stop("tune_lfpls: all fits failed.")
  }
  best
}