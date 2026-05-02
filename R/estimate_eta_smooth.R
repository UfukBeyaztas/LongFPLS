estimate_eta_smooth <- function(X, tstd, df_time = 5) {
  X <- as.matrix(X)
  Zt <- cbind(1, ns(tstd, df = df_time))
  qrZ <- qr(Zt)
  theta <- qr.coef(qrZ, X)
  eta_hat <- Zt %*% theta 
  list(eta_mat = eta_hat, Zt = Zt)
}
