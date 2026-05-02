whiten_by_subject <- function(y, Z, X, subject, sigma_b, sigma_e) {
  subject <- as.vector(subject)
  u_id <- unique(subject)
  
  yw <- numeric(length(y))
  Zw <- matrix(0, nrow(Z), ncol(Z))
  Xw <- matrix(0, nrow(X), ncol(X))
  
  for (id in u_id) {
    idx <- which(subject == id)
    yw[idx]   <- as.numeric(whiten_block_ri(y[idx], sigma_b, sigma_e))
    Zw[idx, ] <- whiten_block_ri(Z[idx, , drop = FALSE], sigma_b, sigma_e)
    Xw[idx, ] <- whiten_block_ri(X[idx, , drop = FALSE], sigma_b, sigma_e)
  }
  list(yw = yw, Zw = Zw, Xw = Xw)
}