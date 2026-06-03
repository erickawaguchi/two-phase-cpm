pred_quantile <- function(f, v, d, tau = 0.5, alpha = 0.05, method = "discrete") {
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
  
  if(!(all(names(bet) %in% colnames(X)))) stop("Column names in d should match with coefficient names in model") # Verify all coefficient names are correct and included in X
  
  lp <- X %*% bet 
  lb <- matrix(sapply(ints, `+`, lp), ncol = ns)
  
  if (method == "interpolated") {
    m_cdf <- cbind(1 - cumprob(lb), 1)
    
    # Get interpolated fitted values
    cp <- cbind(cumprob(lb), 0)
    
    # build y-grid for interpolation, row by row
    m_yvals <- matrix(NA, nrow(lb), ns + 2)
    for (j in seq_len(nrow(lb))) {
      # normalized weights between first and last cutpoint
      ws <- c(0, (cp[j, -(ns+1)] - cp[j, 1]) / (cp[j, ns] - cp[j, 1]), 1)
      m_yvals[j, ] <- (1 - ws) * c(vals[1], vals) + ws * c(vals, vals[ns + 1])
    }
    
    # map tau onto the same probability scale and interpolate y
    target <- cumprob(inverse(1 - tau))
    
    z <- vapply(seq_len(nrow(lb)), function(i)
      approx(x = c(1, cp[i, ]), y = m_yvals[i, ], xout = target, rule = 2)$y,
      numeric(1)) 
  }
  
  if (method == "discrete") {
    m.cdf <- cbind(1 - cumprob(lb), 1)
    id <- apply(m.cdf, 1, FUN = function(x) {
      min(which(x >= tau))[1]
    })
    z <- vals[id]
  }
  
  # Verify standard errors
  lb.se        <- matrix(NA, ncol = ns, nrow = nrow(X))
  idx          <- which(names(c(ints, bet)) %in% colnames(X))
  dlb.dtheta   <- as.matrix(cbind(1, X))

  # Doing this the original way since info.inverse will be calculated via survey package
  for(i in 1:ns){
    v.i <- v[c(i, idx), c(i, idx)]
    lb.se[, i] <- sqrt(diag(dlb.dtheta %*% v.i %*% t(dlb.dtheta)))
    # Compute (i, idx) portion of info inverse, multiplied by t(dlb.dtheta)
    #v.i        <- infoMxop(f$info.matrix, i=c(i, idx), B=t(dlb.dtheta))
    #lb.se[, i] <- sqrt(Matrix::diag(dlb.dtheta %*% v.i))
  }
  
  w <- qnorm(1 - alpha / 2)
  
  ci.ub <- matrix(sapply(1:ns, FUN=function(i) {1 - cumprob(lb[, i] - w * lb.se[, i])}), ncol = ns)
  ci.lb <- matrix(sapply(1:ns, FUN=function(i) {1 - cumprob(lb[, i] + w * lb.se[, i])}), ncol = ns)
  
  if (method == "interpolated") {
  z.ub <- sapply(1:nrow(lb),
                 function(i) approx(c(1, 1 - ci.lb[i,], 0), m_yvals[i,],
                                    xout = cumprob(inverse(1 - tau)),
                                    rule = 2)$y)
  z.lb <- sapply(1:nrow(lb),
                 function(i) approx(c(1, 1 - ci.ub[i,], 0), m_yvals[i,],
                                    xout = cumprob(inverse(1 - tau)),
                                    rule = 2)$y)
  }
  
  if (method == "discrete") {
    id <- apply(cbind(ci.lb, 1), 1, FUN = function(x) {
      min(which(x >= tau))[1]
    })
    z.ub <- vals[id]
    id <- apply(cbind(ci.ub, 1), 1, FUN = function(x) {
      min(which(x >= tau))[1]
    })
    z.lb <- vals[id]
  }
  
  return(
    list(
      p = z,
      lower = z.lb,
      upper = z.ub
    )
  )
}
