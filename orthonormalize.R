orthonormalize <- function(basis) {
  basis <- as.matrix(basis)
  n <- nrow(basis)
  p <- ncol(basis)
  if (p == 0) return(basis)
  
  ortho <- matrix(0, n, p)
  nrm1 <- sqrt(sum(basis[, 1]^2))
  if (nrm1 < 1e-14) return(matrix(0, n, 0))
  ortho[, 1] <- basis[, 1] / nrm1
  
  if (p > 1) {
    for (i in 2:p) {
      v <- basis[, i]
      for (j in 1:(i - 1))
        v <- v - sum(v * ortho[, j]) * ortho[, j]
      nrm <- sqrt(sum(v^2))
      ortho[, i] <- if (nrm < 1e-14) 0 else v / nrm
    }
  }
  
  keep <- colSums(ortho^2) > 1e-14
  ortho[, keep, drop = FALSE]
}