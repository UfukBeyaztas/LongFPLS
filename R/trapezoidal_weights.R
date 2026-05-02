trapezoidal_weights <- function(t_grid) {
  t_grid <- as.numeric(t_grid)
  D <- length(t_grid)
  if (D < 2) return(rep(1, D))
  dt <- diff(t_grid)
  w <- numeric(D)
  w[1] <- dt[1] / 2
  w[D] <- dt[D - 1] / 2
  if (D > 2)
    w[2:(D - 1)] <- (dt[1:(D - 2)] + dt[2:(D - 1)]) / 2
  w
}
