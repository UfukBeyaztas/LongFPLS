prepare_lfpls <- function(y, subject, times, X, tps,
                          scalar_cov = NULL,
                          presmooth_X = TRUE,
                          df_spline = 10,
                          df_time_eta = 5,
                          df_time_Z = 4,
                          lme_maxiter = 80) {
  
  y <- as.numeric(y)
  subject <- as.vector(subject)
  times <- as.numeric(times)
  X <- as.matrix(X)
  
  n <- length(y)
  D <- ncol(X)
  if (nrow(X) != n) stop("nrow(X) must equal length(y).")
  
  dec <- decompose_fx(
    X = X,
    subject = subject,
    times = times,
    presmooth = presmooth_X,
    tps = tps,
    df_spline = df_spline,
    df_time_eta = df_time_eta
  )
  
  tstd <- dec$tstd
  w <- trapezoidal_weights(tps)
  
  Ws <- diag(w, length(w))
  
  u_id <- dec$u_id
  subject_index <- match(subject, u_id)
  
  B0_wtd <- dec$B0_subj[subject_index, , drop = FALSE] %*% Ws
  B1_wtd <- dec$B1_subj[subject_index, , drop = FALSE] %*% Ws
  B1t_wtd <- B1_wtd * matrix(tstd, nrow = n, ncol = D)
  U_wtd <- dec$U %*% Ws
  
  Z_time <- cbind(1, ns(tstd, df = df_time_Z))
  colnames(Z_time) <- paste0("Z", seq_len(ncol(Z_time)))  # Z1, Z2, ...
  
  Z_scalar <- as_scalar_matrix(scalar_cov, n)
  if (ncol(Z_scalar) > 0) {
    colnames(Z_scalar) <- paste0("S_", make.names(colnames(Z_scalar), unique = TRUE))
  }
  
  Z_full <- cbind(Z_time, Z_scalar)
  z_names <- colnames(Z_full)
  
  df0 <- data.frame(y = y, subject = factor(subject))
  for (nm in z_names) df0[[nm]] <- Z_full[, nm]
  
  formula0 <- as.formula(paste0("y ~ -1 + ", paste(z_names, collapse = " + ")))
  
  fit0 <- lme(
    fixed = formula0,
    random = ~1 | subject,
    data = df0,
    method = "REML",
    control = lmeControl(msMaxIter = lme_maxiter, opt = "optim")
  )
  
  vc0 <- VarCorr(fit0)
  sigma_e0 <- as.numeric(vc0["Residual", "StdDev"])
  sigma_b0 <- as.numeric(vc0[1, "StdDev"])
  if (!is.finite(sigma_b0)) sigma_b0 <- 0
  
  P <- second_diff_penalty(D)
  
  list(
    y = y,
    subject = subject,
    subject_factor = factor(subject),
    n = n,
    D = D,
    tstd = tstd,
    w = w,
    
    Z_full = Z_full,
    z_names = z_names,
    df0 = df0,
    formula0 = formula0,
    fit0 = fit0,
    
    sigma_e0 = sigma_e0,
    sigma_b0 = sigma_b0,
    
    B0_wtd = B0_wtd,
    B1t_wtd = B1t_wtd,
    U_wtd = U_wtd,
    
    penalty = P,
    decomposition = dec
  )
}
