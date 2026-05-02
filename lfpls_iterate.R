lfpls_iterate <- function(prep,
                          n_comp_b0 = 3,
                          n_comp_b1 = 3,
                          n_comp_u = 3,
                          ortho = TRUE,
                          ridge = 1e-10,
                          lambda_b = 1e-2,
                          max_iter = 8,
                          tol = 1e-4,
                          lme_maxiter = 80) {
  
  y <- prep$y
  subject <- prep$subject
  n <- prep$n
  D <- prep$D
  w <- prep$w
  
  Z_full  <- prep$Z_full
  z_names <- prep$z_names
  df0 <- prep$df0
  
  B0_wtd <- prep$B0_wtd
  B1t_wtd <- prep$B1t_wtd
  U_wtd <- prep$U_wtd
  P <- prep$penalty
  
  sigma_e_cur <- prep$sigma_e0
  sigma_b_cur <- prep$sigma_b0
  
  b0_old <- rep(0, D); b1_old <- rep(0, D); bU_old <- rep(0, D)
  best <- NULL
  
  for (it in 1:max_iter) {
    
    X_stack <- cbind(B0_wtd, B1t_wtd, U_wtd)
    
    wb <- whiten_by_subject(
      y = y,
      Z = Z_full,
      X = X_stack,
      subject = subject,
      sigma_b = sigma_b_cur,
      sigma_e = sigma_e_cur
    )
    
    yw <- wb$yw
    Zw <- wb$Zw
    Xw_all <- wb$Xw
    
    B0_w  <- Xw_all[, 1:D, drop = FALSE]
    B1t_w <- Xw_all[, (D + 1):(2 * D), drop = FALSE]
    U_w   <- Xw_all[, (2 * D + 1):(3 * D), drop = FALSE]
    
    K0 <- residualize_rapls(B0_w,  yw, Zw, n_comp = n_comp_b0, ortho = ortho, ridge = ridge)
    K1 <- residualize_rapls(B1t_w, yw, Zw, n_comp = n_comp_b1, ortho = ortho, ridge = ridge)
    
    compsB_w <- compsB <- NULL
    
    if (ncol(K0) > 0) {
      comps0_w <- B0_w %*% K0
      comps0   <- B0_wtd %*% K0
      colnames(comps0_w) <- paste0("PLS0w_", seq_len(ncol(comps0_w)))
      colnames(comps0)   <- paste0("PLS0_",  seq_len(ncol(comps0)))
    } else { comps0_w <- comps0 <- NULL }
    
    if (ncol(K1) > 0) {
      comps1_w <- B1t_w %*% K1
      comps1   <- B1t_wtd %*% K1
      colnames(comps1_w) <- paste0("PLS1w_", seq_len(ncol(comps1_w)))
      colnames(comps1)   <- paste0("PLS1_",  seq_len(ncol(comps1)))
    } else { comps1_w <- comps1 <- NULL }
    
    if (!is.null(comps0_w) || !is.null(comps1_w)) {
      compsB_w <- cbind(comps0_w, comps1_w)
      compsB   <- cbind(comps0,   comps1)
    }
    
    Zw_ext <- Zw
    if (!is.null(compsB_w)) Zw_ext <- cbind(Zw, compsB_w)
    
    KU <- residualize_rapls(U_w, yw, Zw_ext, n_comp = n_comp_u, ortho = ortho, ridge = ridge)
    compsU <- if (ncol(KU) > 0) U_wtd %*% KU else NULL
    if (!is.null(compsU)) colnames(compsU) <- paste0("PLSU_", seq_len(ncol(compsU)))
    
    if (is.null(compsB) && is.null(compsU)) break
    
    df1 <- df0
    comp_names <- character(0)
    
    if (!is.null(compsB)) {
      for (nm in colnames(compsB)) df1[[nm]] <- compsB[, nm]
      comp_names <- c(comp_names, colnames(compsB))
    }
    if (!is.null(compsU)) {
      for (nm in colnames(compsU)) df1[[nm]] <- compsU[, nm]
      comp_names <- c(comp_names, colnames(compsU))
    }
    
    rhs <- paste(z_names, collapse = " + ")
    if (length(comp_names) > 0) rhs <- paste(rhs, paste(comp_names, collapse = " + "), sep = " + ")
    formula1 <- as.formula(paste0("y ~ -1 + ", rhs))
    
    fit1 <- tryCatch(
      nlme::lme(
        fixed = formula1,
        random = ~1 | subject,
        data = df1,
        method = "REML",
        control = lmeControl(msMaxIter = lme_maxiter, opt = "optim")
      ),
      error = function(e) NULL
    )
    if (is.null(fit1)) break
    
    vc <- VarCorr(fit1)
    sigma_e_cur <- as.numeric(vc["Residual", "StdDev"])
    sigma_b_cur <- as.numeric(vc[1, "StdDev"])
    if (!is.finite(sigma_b_cur)) sigma_b_cur <- 0
    
    coef_all <- fixed.effects(fit1)
    
    keep0 <- grep("^PLS0_", names(coef_all), value = TRUE)
    keep1 <- grep("^PLS1_", names(coef_all), value = TRUE)
    keepU <- grep("^PLSU_", names(coef_all), value = TRUE)
    
    b0_raw <- rep(0, D); b1_raw <- rep(0, D); bU_raw <- rep(0, D)
    
    if (length(keep0) > 0 && ncol(K0) > 0) {
      b0_raw <- as.numeric(K0[, seq_along(keep0), drop = FALSE] %*% as.numeric(coef_all[keep0]))
    }
    if (length(keep1) > 0 && ncol(K1) > 0) {
      b1_raw <- as.numeric(K1[, seq_along(keep1), drop = FALSE] %*% as.numeric(coef_all[keep1]))
    }
    if (length(keepU) > 0 && ncol(KU) > 0) {
      bU_raw <- as.numeric(KU[, seq_along(keepU), drop = FALSE] %*% as.numeric(coef_all[keepU]))
    }
    
    b0_tmp <- shrink_coefficients(b0_raw, P, lambda_b)
    b1_tmp <- shrink_coefficients(b1_raw, P, lambda_b)
    bU_tmp <- shrink_coefficients(bU_raw, P, lambda_b)
    
    a_m <- 1 / it
    b0 <- (1 - a_m) * b0_old + a_m * b0_tmp
    b1 <- (1 - a_m) * b1_old + a_m * b1_tmp
    bU <- (1 - a_m) * bU_old + a_m * bU_tmp
    
    err <- sum(w * (b0 - b0_old)^2) + sum(w * (b1 - b1_old)^2) + sum(w * (bU - bU_old)^2)
    
    best <- list(
      fit = fit1,
      b0 = b0, b1 = b1, bU = bU,
      n_comp_b0 = n_comp_b0,
      n_comp_b1 = n_comp_b1,
      n_comp_u  = n_comp_u,
      lambda_b = lambda_b,
      iter = it,
      w = w,
      sigma_e = sigma_e_cur,
      sigma_b = sigma_b_cur
    )
    
    if (is.finite(err) && err < tol) break
    b0_old <- b0; b1_old <- b1; bU_old <- bU
  }
  
  best
}