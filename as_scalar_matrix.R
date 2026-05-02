as_scalar_matrix <- function(scalar_cov, n) {
  if (is.null(scalar_cov)) return(matrix(0, n, 0))
  
  if (is.data.frame(scalar_cov)) {
    mm <- model.matrix(~ . - 1, data = scalar_cov)
    if (nrow(mm) != n) stop("scalar_cov has wrong number of rows.")
    return(mm)
  }
  
  if (is.vector(scalar_cov)) {
    scalar_cov <- matrix(scalar_cov, ncol = 1)
  }
  
  scalar_cov <- as.matrix(scalar_cov)
  if (nrow(scalar_cov) != n) stop("scalar_cov has wrong number of rows.")
  
  if (is.null(colnames(scalar_cov))) {
    colnames(scalar_cov) <- paste0("S", seq_len(ncol(scalar_cov)))
  }
  scalar_cov
}