second_diff_penalty <- function(D) {
  if (D <= 2) return(diag(D))
  D2 <- matrix(0, D - 2, D)
  for (i in 1:(D - 2))
    D2[i, i:(i + 2)] <- c(1, -2, 1)
  crossprod(D2)
}
