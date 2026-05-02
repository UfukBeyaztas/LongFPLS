shrink_coefficients <- function(b, penalty, lambda) {
  b <- as.numeric(b)
  if (lambda <= 0) return(b)
  D <- length(b)
  A <- diag(D) + lambda * penalty
  as.numeric(solve(A, b))
}
