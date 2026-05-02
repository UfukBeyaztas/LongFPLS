simulate_data <- function(case = c("1", "2"),
                          nI = 100,
                          nV = 4,
                          tps = seq(0.5, 99.5, by = 1),
                          sigma_x = 0.05,
                          sigma_y = 0.35,
                          sigma_b = 0.50,
                          q_scalar = 5,
                          theta_true = NULL) {
  
  case <- match.arg(case)
  
  if (is.null(theta_true))
    theta_true <- c(0.6, -0.4, 0.3, 0.2, 0.1)
  
  generate_random_times <- function(j) {
    tt <- 0
    if (j > 1) {
      for (jj in 2:j) tt <- c(tt, sum(tt) + runif(1))
    }
    tt - mean(tt)
  }
  
  w_orth <- function(B, w) {
    B <- as.matrix(B)
    Wsqrt <- diag(sqrt(w), nrow = length(w), ncol = length(w))
    Z <- Wsqrt %*% B
    Q <- qr.Q(qr(Z))
    Bout <- solve(Wsqrt, Q)
    keep <- colSums(Bout^2) > 1e-12
    Bout[, keep, drop = FALSE]
  }
  
  normalize_curve <- function(v, w) {
    v <- as.numeric(v)
    nrm <- sqrt(sum(w * v^2))
    if (nrm < 1e-12) return(v)
    v / nrm
  }
  
  toeplitz_cov <- function(sd_vec, rho) {
    p <- length(sd_vec)
    R <- outer(1:p, 1:p, function(i, j) rho^abs(i - j))
    diag(sd_vec, p, p) %*% R %*% diag(sd_vec, p, p)
  }
  
  D <- length(tps)
  w <- trapezoidal_weights(tps)
  s <- (tps - min(tps)) / (max(tps) - min(tps))
  
  dict_raw <- cbind(
    bs(s, df = 5, degree = 3, intercept = TRUE),
    sin(2 * pi * s),
    cos(2 * pi * s),
    sin(4 * pi * s),
    cos(4 * pi * s),
    exp(-70  * (s - 0.25)^2),
    exp(-100 * (s - 0.55)^2),
    exp(-140 * (s - 0.80)^2)
  )
  
  Psi <- w_orth(dict_raw, w)
  K <- ncol(Psi)
  
  ncl <- c(1:nI, sample(1:nI, (nV - 1) * nI, replace = TRUE))
  Jvec <- as.numeric(sort(table(ncl), decreasing = TRUE))
  n <- sum(Jvec)
  
  subject <- rep(1:nI, Jvec)
  
  times <- unlist(lapply(Jvec, FUN = generate_random_times))
  times <- times / sqrt(var(times))
  
  subject_index <- match(subject, 1:nI)
  
  sd_B0 <- c(2.50, 2.00, 1.60, 1.20, 0.45, 0.35, 0.28, 0.22, 0.18, 0.14, 0.10, 0.08)[1:K]
  sd_B1 <- c(1.80, 1.40, 1.10, 0.85, 0.30, 0.24, 0.18, 0.14, 0.12, 0.10, 0.08, 0.06)[1:K]
  sd_U  <- c(1.70, 1.35, 1.05, 0.80, 0.55, 0.42, 0.30, 0.22, 0.18, 0.14, 0.10, 0.08)[1:K]
  
  Sigma_B0 <- toeplitz_cov(sd_B0, rho = 0.45)
  Sigma_B1 <- toeplitz_cov(sd_B1, rho = 0.35)
  Sigma_U  <- toeplitz_cov(sd_U,  rho = 0.30)
  
  A0 <- mvrnorm(n = nI, mu = rep(0, K), Sigma = Sigma_B0)
  A1 <- mvrnorm(n = nI, mu = rep(0, K), Sigma = Sigma_B1)
  C  <- mvrnorm(n = n,  mu = rep(0, K), Sigma = Sigma_U)
  
  eta_base <- outer(
    times, s,
    FUN = function(tt, ss) {
      0.30 * tt * cos(2 * pi * ss) +
        0.12 * (tt^2 - mean(times^2)) * (ss - 0.5)
    }
  )
  
  B0 <- A0[subject_index, , drop = FALSE] %*% t(Psi)
  B1 <- A1[subject_index, , drop = FALSE] %*% t(Psi)
  time_mat <- matrix(times, nrow = n, ncol = D)
  
  if (q_scalar > 0) {
    Z <- matrix(rnorm(n * q_scalar), n, q_scalar)
    colnames(Z) <- paste0("Z", 1:q_scalar)
    Z_effect <- as.numeric(Z %*% theta_true)
  } else {
    Z <- matrix(0, n, 0)
    Z_effect <- rep(0, n)
  }
  
  if (case == "1") {
    B_comp <- 0.10 * (B0 + B1 * time_mat)
    U_comp <- 1.00 * (C %*% t(Psi))
    eta_mat <- 0.03 * eta_base
    
    W_latent <- eta_mat + B_comp + U_comp
    X <- W_latent + matrix(rnorm(n * D, mean = 0, sd = sigma_x), nrow = n, ncol = D)
    
    beta_true <- 1.10 * dnorm(s, mean = 0.30, sd = 0.18) -
      0.85 * dnorm(s, mean = 0.72, sd = 0.16)
    
    beta_true <- beta_true / sqrt(sum(w * beta_true^2))
    
    mu_true <- as.numeric(
      rep(rnorm(nI, 0, sigma_b), Jvec) +
        Z_effect +   
        1.50 * W_latent %*% (w * beta_true)
    )
    
    truth <- list(
      beta_true = beta_true,
      beta_B_true = NULL,
      beta_U_true = NULL
    )
  }
  
  if (case == "2") {
    B_comp <- 1.20 * (B0 + B1 * time_mat)
    U_comp <- 0.80 * (C %*% t(Psi))
    eta_mat <- 0.10 * eta_base
    
    W_latent <- eta_mat + B_comp + U_comp
    X <- W_latent + matrix(rnorm(n * D, mean = 0, sd = sigma_x), nrow = n, ncol = D)
    
    beta_B_true <- 1.00 * dnorm(s, mean = 0.25, sd = 0.20) +
      0.75 * dnorm(s, mean = 0.78, sd = 0.18) -
      0.55 * dnorm(s, mean = 0.55, sd = 0.24)
    
    beta_U_true <- -0.80 * dnorm(s, mean = 0.38, sd = 0.16) +
      0.95 * dnorm(s, mean = 0.82, sd = 0.15)
    
    beta_B_true <- beta_B_true / sqrt(sum(w * beta_B_true^2))
    beta_U_true <- beta_U_true / sqrt(sum(w * beta_U_true^2))
    
    mu_true <- as.numeric(
      rep(rnorm(nI, 0, sigma_b), Jvec) +
        Z_effect +   
        1.60 * B_comp %*% (w * beta_B_true) +
        1.10 * U_comp %*% (w * beta_U_true)
    )
    
    truth <- list(
      beta_true = NULL,
      beta_B_true = beta_B_true,
      beta_U_true = beta_U_true
    )
  }
  
  y <- mu_true + rnorm(n, mean = 0, sd = sigma_y)
  
  list(
    case = case,
    nI = nI,
    n = n,
    Jvec = Jvec,
    subject = subject,
    times = times,
    tps = tps,
    X = X,
    y = y,
    mu_true = mu_true,
    
    Z = Z,
    theta_true = theta_true,
    
    W_latent = W_latent,
    B_comp = B_comp,
    U_comp = U_comp,
    eta_mat = eta_mat,
    Psi = Psi,
    truth = truth
  )
}
