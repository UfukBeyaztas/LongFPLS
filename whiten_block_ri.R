whiten_block_ri <- function(block, sigma_b, sigma_e, eps = 1e-12) {
  block <- as.matrix(block)
  ni <- nrow(block)
  if (ni == 0) return(block)
  
  se2 <- max(sigma_e^2, eps)
  sb2 <- max(sigma_b^2, 0)
  
  u <- rep(1, ni) / sqrt(ni)
  proj <- as.numeric(u %*% block)
  block_par <- u %*% matrix(proj, 1)
  block_perp <- block - block_par
  
  lam_par <- se2 + ni * sb2
  lam_perp <- se2
  
  block_perp / sqrt(lam_perp) + block_par / sqrt(lam_par)
}