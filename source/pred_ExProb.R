pred_ExProb <- function(f, v, d, alpha = 0.05, method = "discrete") {
  X <- as.matrix(d)
  # Manual (most of this code is derived from orm.s)
  # See https://rdrr.io/cran/rms/src/R/orm.s
  
  # Get necessary components from fitted model "f"
  ns   <- f$non.slopes                    # number of intercepts
  ints <- f$coef[1:ns]                    # intercept vector (cutpoints)
  bet  <- f$coef[-(1:ns)]                 # slopes
  vals <- f$yunique; if(!length(vals)) vals <- names(f$freq)
  vals <- as.numeric(vals)                  # ordered outcome values (numeric)
  inverse <- eval(f$famfunctions[2])
  cumprob <- eval(f$famfunctions[1])
  yname <- f$yname
  info <- f$info.matrix
  
  if(!(all(names(bet) %in% colnames(X)))) stop("Column names in d should match with coefficient names in model") # Verify all coefficient names are correct and included in X

  lp <- c(X %*% bet); 
  #lp - (lp0 - f$coefficients[f$interceptRef])
  y <- NULL
  
  # From here on out is from ExProb
  prob <- cumprob(sapply(c(1e+30, ints), "+", lp)); prob
  
  dim(prob) <- c(length(lp), length(vals))
  if (!length(y)) {
    colnames(prob) <- paste("Prob(Y>=", vals, ")", sep = "")
    y <- vals
    result <- structure(list(y = vals, prob = prob, yname = yname), 
                        class = "ExProb")
  }
  
  # Still need to do this part here...(to get conf intervals)
  conf.int = 1 - alpha
  
  index <- sapply(y, FUN = function(x) {
    if (x <= min(vals)) 
      result <- 1
    else if (x >= max(vals)) 
      result <- length(vals)
    else which(x <= vals)[1] - 1
  })
  
  idx <- which(names(c(ints, bet)) %in% colnames(X))
  dlb.dtheta <- as.matrix(cbind(1, X))
  lb.se <- sapply(1:length(y), function(i) Matrix::diag(dlb.dtheta %*% v[c(index[i], idx), c(index[i], idx)] %*% t(dlb.dtheta)))
  lb.se <- matrix(sqrt(lb.se), ncol = length(y))
  m.alpha <- c(ints, bet)[index]
  lb <- matrix(sapply(m.alpha, "+", lp), ncol = length(y))
  ci.ub <- matrix(sapply(1:length(y), FUN = function(i) {
    cumprob(lb[, i] + qnorm((1 + conf.int)/2) * lb.se[, 
                                                      i])
  }), ncol = length(y))
  ci.lb <- matrix(sapply(1:length(y), FUN = function(i) {
    cumprob(lb[, i] - qnorm((1 + conf.int)/2) * lb.se[, 
                                                      i])
  }), ncol = length(y))
  ci.ub[, which(y <= min(vals))] <- ci.lb[, which(y <= 
                                                    min(vals))] <- 1
  ci.ub[, which(y >= max(vals))] <- ci.lb[, which(y >= 
                                                    max(vals))] <- 0
  if (length(y) > 1) 
    colnames(ci.ub) <- colnames(ci.lb) <- colnames(result$prob)

    return(
    list(
      y = y,
      p = prob,
      lower = ci.lb,
      upper = ci.ub
    )
  )
}