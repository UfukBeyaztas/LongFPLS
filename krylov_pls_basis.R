krylov_pls_basis <- function(X, y, n_comp) {
  X <- as.matrix(X)
  y <- as.numeric(y)
  n <- length(y)
  d <- ncol(X)
  n_comp <- min(n_comp, d)
  if (n_comp <= 0 || d == 0) return(matrix(0, d, 0))
  
  K <- matrix(0, d, n_comp)
  K[, 1] <- as.numeric(crossprod(X, y)) / n
  if (n_comp > 1) {
    for (i in 2:n_comp) {
      K[, i] <- as.numeric(crossprod(X, X %*% K[, i - 1])) / n
    }
  }
  K
}