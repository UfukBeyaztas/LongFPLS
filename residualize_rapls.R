residualize_rapls <- function(Xw, yw, Zw, n_comp, ortho = TRUE, ridge = 1e-10) {
  Xw <- as.matrix(Xw)
  yw <- as.numeric(yw)
  Zw <- as.matrix(Zw)
  d <- ncol(Xw)
  if (d == 0) return(matrix(0, d, 0))
  
  q <- ncol(Zw)
  if (q > 0) {
    ZtZ <- crossprod(Zw) + ridge * diag(q)
    ZtZ_inv <- solve(ZtZ)
    
    B_X <- ZtZ_inv %*% crossprod(Zw, Xw) 
    B_y <- as.numeric(ZtZ_inv %*% crossprod(Zw, yw))
    
    X_res <- Xw - Zw %*% B_X
    y_res <- yw - as.numeric(Zw %*% B_y)
  } else {
    X_res <- Xw
    y_res <- yw
  }
  
  n_comp <- min(n_comp, ncol(X_res))
  if (n_comp <= 0) return(matrix(0, d, 0))
  
  K <- krylov_pls_basis(X_res, y_res, n_comp)
  if (ncol(K) == 0) return(matrix(0, d, 0))
  if (ortho) K <- orthonormalize(K)
  
  keep <- colSums(K^2) > 1e-14
  K[, keep, drop = FALSE]
}