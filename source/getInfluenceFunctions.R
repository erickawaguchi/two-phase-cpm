getBetaInfluenceFunctions <- function(fit, getC1 = FALSE) {
  # Quick checks
  if(!(all(class(fit) %in% c("orm", "rms")))) stop("Fitted object must be of class 'orm' from the rms package")
  if(is.null(fit$mscore)) stop("orm must be fit with 'mscore = TRUE'")
  if(strsplit(as.character(packageVersion("rms")), "[.]")[[1]][1] < 8) stop("Requires an orm fitted object from rms version 8.0.0 or higher")

  Na <- fit$non.slopes # Number of intercepts
  Np <- length(fit$coefficients) - Na # Number of slope parameters
  
  a  <- Matrix::bandSparse(Na, k = c(0, 1), diagonals = fit$info.matrix$a, symmetric = TRUE)
  ab <- fit$info.matrix$ab                        
  b  <- fit$info.matrix$b
  
  
  #B  <- fit$mscore
  B1 <- fit$mscore[, 1:Na]
  B2 <- fit$mscore[, (Na + 1):(Na + Np)]
  
  a_inv_ab <- Matrix::solve(a, ab)
  
  #b - t(ab) %*% solve(a, ab) is the Schur complement
  #C2 <- t(Matrix::solve(b - Matrix::t(ab) %*% Matrix::solve(a, ab), Matrix::t(B2) - t(ab) %*% Matrix::solve(a, Matrix::t(B1))))
  
  C2 <- solve(b - t(ab) %*% a_inv_ab, t(B2 - B1 %*% a_inv_ab))
  
  if(getC1) {
    C1 <- -t(Matrix::solve(a, t(B1) - ab %*% C2))
  } else {
    C1 <- NULL
  }
  return(list(C1 = C1, C2 = -t(C2)))
}