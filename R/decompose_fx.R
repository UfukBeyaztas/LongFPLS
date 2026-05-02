decompose_fx <- function(X, subject, times,
                         presmooth = TRUE,
                         tps = NULL,
                         df_spline = 10,
                         df_time_eta = 5) {
  X <- as.matrix(X)
  subject <- as.vector(subject)
  times <- as.numeric(times)
  n <- nrow(X)
  D <- ncol(X)
  
  if (presmooth) {
    if (is.null(tps)) stop("presmooth = TRUE requires tps.")
    X <- t(apply(X, 1, function(x) smooth.spline(tps, x, df = df_spline)$y))
  }
  
  tstd <- (times - mean(times)) / sqrt(var(times))
  eta <- estimate_eta_smooth(X, tstd, df_time = df_time_eta)
  X_tilde <- X - eta$eta_mat
  
  u_id <- unique(subject)
  I <- length(u_id)
  
  B0_subj <- matrix(0, I, D)
  B1_subj <- matrix(0, I, D)
  U_mat <- matrix(0, n, D)
  
  for (ii in seq_len(I)) {
    id <- u_id[ii]
    idx <- which(subject == id)
    Ji <- length(idx)
    
    Ti <- cbind(1, tstd[idx])
    Xt_i <- X_tilde[idx, , drop = FALSE]
    
    if (Ji >= 2) {
      Ginv <- robust_inverse(crossprod(Ti), ridge = 1e-10)
      Bi <- Ginv %*% crossprod(Ti, Xt_i)
    } else {
      Bi <- rbind(Xt_i[1, ], rep(0, D))
    }
    
    B0_subj[ii, ] <- Bi[1, ]
    B1_subj[ii, ] <- Bi[2, ]
    
    B_i <- Ti %*% Bi
    U_mat[idx, ] <- Xt_i - B_i
  }
  
  list(
    X_used = X,
    tstd = tstd,
    eta_mat = eta$eta_mat,
    X_tilde = X_tilde,
    u_id = u_id,
    B0_subj = B0_subj,
    B1_subj = B1_subj,
    U = U_mat
  )
}
