robust_inverse <- function(A, ridge = 1e-10) {
  A <- as.matrix(A)
  A <- A + ridge * diag(ncol(A))
  ginv(A)
}
