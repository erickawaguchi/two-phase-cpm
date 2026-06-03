# Simulation for Supplemental Table A1
rm(list = ls())
library(rms)
library(dplyr)
library(survey)
library(foreach)
library(readxl)
library(sandwich)

data_dir <- "/Users/ekawaguc/GitHub/svycpm/"

# Source files:
source(paste0(data_dir, "source/getInfluenceFunctions.R"))
source(paste0(data_dir, "source/pred_quantile.R"))
source(paste0(data_dir, "source/pred_ExProb.R"))

simList <- read_xlsx(paste0(data_dir, "simulations/simulationList.xlsx"))

simIdx <- 1
simVals <- simList[simIdx, ]

B     <- 1000 # Number of Monte Carlo replicates
N     <- 3000 # Phase 1 size
n     <- simVals$n # Phase 2 size
rho_x <- simVals$r_x
rho_y <- simVals$r_y

s_true <- 1
s_x    <- sqrt(s_true^2 * (1 - rho_x^2) / rho_x^2)
s_y    <- sqrt(s_true^2 * (1 - rho_y^2) / rho_y^2)

p_true  <- 0.35          
Se      <- simVals$Se
Sp      <- simVals$Sp

b0     <- c(0.5, -0.5)

# Generate "population"
set.seed(10000)
NN <- 1000000
y  <-  0.1 + sum(c(-1, 1) * b0) + rnorm(n = NN)
y2 <- y + rnorm(NN, mean = 0, sd = s_y)
cor(y, y2)

# Simulate error-prone latent y (and error-prone transformed y)
yt  <- qchisq(pnorm(y), df = 5)
yt2 <- qchisq(pnorm(y2), df = 5)
cor(yt, yt2)

qY1 <- median(yt); qY1
qY2 <- quantile(yt, probs = 0.8); qY2
ep3     <- mean(yt >= 3); ep3
ep6     <- mean(yt >= 6); ep6

set.seed(simVals$seed) + 200
B <- 1000
out <- foreach(b = 1:B, .combine = 'rbind') %do% {
  cat(paste(b, "out of", B, "\n"))
  x1 <- rnorm(N, mean = 0, sd = s_true)
  x2 <- rbinom(N, size = 1, prob = p_true)
  X  <- cbind(x1, x2)
  
  # Generate latent y
  y  <-  0.1 + X %*% b0 + rnorm(n = N)
  
  # Simulate error-prone latent y
  y2 <- y + rnorm(N, mean = 0, sd = s_y)
  
  # Simulate error-prone latent y (and error-prone transformed y)
  yt <- qchisq(pnorm(y), df = 5)
  yt2 <- qchisq(pnorm(y2), df = 5)
  
  # Simulate "error prone phase 1 covariate data"
  x1.p1 <- x1 + rnorm(N, mean = 0, sd = s_x)
  x2.p1 <- ifelse(x2 == 1, rbinom(N, 1, Se), rbinom(N, 1, 1 - Sp))
  
  d <- data.frame(
    id = 1:N,
    y = y,
    y2 = y2,
    yt.true = yt,
    yt.p1 = yt2,
    x1.true = x1,
    x2.true = x2,
    x1.p1 = x1.p1,
    x2.p1 = x2.p1,
    x1.p2 = x1,
    x2.p2 = x2
  )
  
  # Get phase 2 subcohort as a stratified sample based on Y and X2
  qbreaks <- quantile(d$yt.p1, probs = c(0, 0.25, 0.75, 1), na.rm = TRUE)
  
  
  d <- d %>% mutate(
    yt.p1.cat = cut(yt.p1, breaks = qbreaks,
                    include.lowest = TRUE,
                    labels = c("0", "1", "2")),
    strata = as.numeric(interaction(x2.p1, yt.p1.cat, drop = TRUE))
  )
  
  strata.names <- unique(d$strata)
  strata.n     <- n / length(strata.names)
  
  a <- foreach(i = 1:length(strata.names), .combine = "rbind") %do% {
    data.frame(id = d %>% filter(strata == strata.names[i[]]) %>% 
                 slice_sample(n = strata.n) %>% pull(id), 
               strata_num = i,
               in.p2 = 1)
  }
  
  
  d <- left_join(d, a) %>%
    mutate(
      in.p2 = ifelse(is.na(in.p2), 0, 1),
      x1.p2 = ifelse(in.p2 == 0, NA, x1),
      x2.p2 = ifelse(in.p2 == 0, NA, x2),
      yt.p2 = ifelse(in.p2 == 0, NA, yt),
    )
  
  d %>% group_by(strata_num) %>%
    summarise(n = n())
  
  # Create deciles and percentiles for Phase 1 y
  d <- d %>% mutate(
    yt10 = as.numeric(cut(yt.p1, breaks = quantile(yt.p1, probs = seq(0, 1, 0.1), na.rm = TRUE),
                          include.lowest = TRUE,
                          labels = FALSE)),
    yt20 = as.numeric(cut(yt.p1, breaks = quantile(yt.p1, probs = seq(0, 1, 0.05), na.rm = TRUE),
                          include.lowest = TRUE,
                          labels = FALSE)),
    yt5 = as.numeric(cut(yt.p1, breaks = quantile(yt.p1, probs = seq(0, 1, 0.2), na.rm = TRUE),
                          include.lowest = TRUE,
                          labels = FALSE))
  )
  
  # Cumulative probability models
  
  # True model
  fit.true <- orm(yt.true ~ x1 + x2, data = d,
                  family = "probit")
  
  # Fit on only phase 1
  fit.p1 <- orm(yt.p1 ~ x1.p1 + x2.p1, data = d, mscore= TRUE,
                family = "probit")
  
  fit.p11 <- orm(yt5 ~ x1.p1 + x2.p1, data = d, mscore= TRUE,
                 family = "probit")
  
  fit.p12 <- orm(yt10 ~ x1.p1 + x2.p1, data = d, mscore= TRUE,
                 family = "probit")
  
  fit.p13 <- orm(yt20 ~ x1.p1 + x2.p1, data = d, mscore= TRUE,
                 family = "probit")
  
  
  d1 <- d2 <- d3 <- d4 <- d
  
  # 1) Calibrate on beta only
  d1$naiveInfl.x1 <- getBetaInfluenceFunctions(fit.p1)$C2[, 1]
  d1$naiveInfl.x2 <- getBetaInfluenceFunctions(fit.p1)$C2[, 2]
  
  designO <- twophase(id=list(~1,~1), strata = list(NULL, ~strata), subset = ~in.p2 ==1,
                      data = d1, method = "approx")
  Ninflcal.orm <- survey::calibrate(designO,formula=~naiveInfl.x1 + naiveInfl.x2,phase=2,calfun="raking") 
  sdat <- model.frame(Ninflcal.orm)
  sdat$.survey.prob.weights <- (1 / Ninflcal.orm$prob) / mean(1 / Ninflcal.orm$prob)
  
  fit1 <- orm(yt.true ~ x1.p2 + x2.p2, 
              data = sdat,
              mscore = TRUE,
              weights = .survey.prob.weights,
              family = "probit")
  
  inflFun <- getBetaInfluenceFunctions(fit1, getC1 = TRUE)
  var1  <- twophasevar(cbind(inflFun$C1, inflFun$C2), Ninflcal.orm)
  
  # 2) Calibrate on alpha AND beta (deciles)
  tmp <- getBetaInfluenceFunctions(fit.p11, getC1 = TRUE)
  tmp1 <- as.matrix(cbind(tmp$C1, tmp$C2))
  colnames(tmp1) <- c(paste0("naiveInfl.a", 1:dim(tmp$C1)[2]), paste0("naiveInfl.x", 1:dim(tmp$C2)[2]))
  d2 <- cbind(d2, tmp1)
  
  designO <- twophase(id=list(~1,~1), strata = list(NULL, ~strata), subset = ~in.p2 ==1,
                      data = d2, method = "approx")
  
  fml <- paste0("~", paste(colnames(tmp1), collapse = "+"))
  Ninflcal.orm <- survey::calibrate(designO,formula=as.formula(fml),phase=2,calfun="raking") 
  sdat <- model.frame(Ninflcal.orm)
  sdat$.survey.prob.weights <- (1 / Ninflcal.orm$prob) / mean(1 / Ninflcal.orm$prob)
  
  fit2 <- orm(yt.true ~ x1.p2 + x2.p2, 
              data = sdat,
              mscore = TRUE,
              weights = .survey.prob.weights,
              family = "probit")
  
  inflFun <- getBetaInfluenceFunctions(fit2, getC1 = TRUE)
  var2  <- twophasevar(cbind(inflFun$C1, inflFun$C2), Ninflcal.orm)
  
  # 3) Calibrate on alpha AND beta (deciles)
  tmp <- getBetaInfluenceFunctions(fit.p12, getC1 = TRUE)
  tmp1 <- as.matrix(cbind(tmp$C1, tmp$C2))
  colnames(tmp1) <- c(paste0("naiveInfl.a", 1:dim(tmp$C1)[2]), paste0("naiveInfl.x", 1:dim(tmp$C2)[2]))
  d3 <- cbind(d3, tmp1)
  
  designO <- twophase(id=list(~1,~1), strata = list(NULL, ~strata), subset = ~in.p2 ==1,
                      data = d3, method = "approx")
  
  fml <- paste0("~", paste(colnames(tmp1), collapse = "+"))
  Ninflcal.orm <- survey::calibrate(designO,formula=as.formula(fml),phase=2,calfun="raking") 
  sdat <- model.frame(Ninflcal.orm)
  sdat$.survey.prob.weights <- (1 / Ninflcal.orm$prob) / mean(1 / Ninflcal.orm$prob)
  
  fit3 <- orm(yt.true ~ x1.p2 + x2.p2, 
              data = sdat,
              mscore = TRUE,
              weights = .survey.prob.weights,
              family = "probit")
  
  inflFun <- getBetaInfluenceFunctions(fit3, getC1 = TRUE)
  var3  <- twophasevar(cbind(inflFun$C1, inflFun$C2), Ninflcal.orm)
  
  # 3) Calibrate on alpha AND beta (deciles)
  tmp <- getBetaInfluenceFunctions(fit.p13, getC1 = TRUE)
  tmp1 <- as.matrix(cbind(tmp$C1, tmp$C2))
  colnames(tmp1) <- c(paste0("naiveInfl.a", 1:dim(tmp$C1)[2]), paste0("naiveInfl.x", 1:dim(tmp$C2)[2]))
  d4 <- cbind(d4, tmp1)
  
  designO <- twophase(id=list(~1,~1), strata = list(NULL, ~strata), subset = ~in.p2 ==1,
                      data = d4, method = "approx")
  
  fml <- paste0("~", paste(colnames(tmp1), collapse = "+"))
  Ninflcal.orm <- survey::calibrate(designO,formula=as.formula(fml),phase=2,calfun="raking") 
  sdat <- model.frame(Ninflcal.orm)
  sdat$.survey.prob.weights <- (1 / Ninflcal.orm$prob) / mean(1 / Ninflcal.orm$prob)
  
  fit4 <- orm(yt.true ~ x1.p2 + x2.p2, 
              data = sdat,
              mscore = TRUE,
              weights = .survey.prob.weights,
              family = "probit")
  
  inflFun <- getBetaInfluenceFunctions(fit4, getC1 = TRUE)
  var4  <- twophasevar(cbind(inflFun$C1, inflFun$C2), Ninflcal.orm)
  
  # Get predicted quantiles (0.5 and 0.8) and predicted EPs when setting x1 = -1 and x2 = 1
  d0 <- data.frame(x1.p2 = -1, x2.p2 = 1)
  
  fit1.q1  <- pred_quantile(fit1, var1, d0, tau = 0.5, method = "discrete")
  fit1.q2  <- pred_quantile(fit1, var1, d0, tau = 0.8, method = "discrete")
  fit1.ep  <- pred_ExProb(fit1, var1, d0)
  
  fit2.q1  <- pred_quantile(fit2, var2, d0, tau = 0.5, method = "discrete")
  fit2.q2  <- pred_quantile(fit2, var2, d0, tau = 0.8, method = "discrete")
  fit2.ep  <- pred_ExProb(fit2, var2, d0)
  
  fit3.q1  <- pred_quantile(fit3, var3, d0, tau = 0.5, method = "discrete")
  fit3.q2  <- pred_quantile(fit3, var3, d0, tau = 0.8, method = "discrete")
  fit3.ep  <- pred_ExProb(fit3, var3, d0)
  
  fit4.q1  <- pred_quantile(fit4, var4, d0, tau = 0.5, method = "discrete")
  fit4.q2  <- pred_quantile(fit4, var4, d0, tau = 0.8, method = "discrete")
  fit4.ep  <- pred_ExProb(fit4, var4, d0)
  
  qu  <- Quantile(fit.true)                       
  true.q1  <- Predict(fit.true, x1 = -1, x2 = 1, fun = function(lp) qu(0.5, lp = lp, method = "discrete"))
  true.q2  <- Predict(fit.true, x1 = -1, x2 = 1, fun = function(lp) qu(0.8, lp = lp, method = "discrete"))
  
  dd <- ExProb(fit.true)
  lp0 <- predict(fit.true, newdata = data.frame(x1 = -1, x2 = 1))
  w <- dd(lp0, X = data.frame(x1 = -1, x2 = 1), conf.int = 0.95)
  
  c(
    # To confirm with original simulation
    tail(fit.true$coefficients, 2),
    tail(fit1$coefficients, 2),
    tail(fit2$coefficients, 2),
    tail(fit3$coefficients, 2),
    tail(fit4$coefficients, 2),
    # New results
    true.q1$yhat,
    true.q1$lower,
    true.q1$upper,
    fit1.q1$p,
    fit1.q1$lower,
    fit1.q1$upper,
    fit2.q1$p,
    fit2.q1$lower,
    fit2.q1$upper,
    fit3.q1$p,
    fit3.q1$lower,
    fit3.q1$upper,
    fit4.q1$p,
    fit4.q1$lower,
    fit4.q1$upper,
    true.q2$yhat,
    true.q2$lower,
    true.q2$upper,
    fit1.q2$p,
    fit1.q2$lower,
    fit1.q2$upper,
    fit2.q2$p,
    fit2.q2$lower,
    fit2.q2$upper,
    fit3.q2$p,
    fit3.q2$lower,
    fit3.q2$upper,
    fit4.q2$p,
    fit4.q2$lower,
    fit4.q2$upper,  
    w$prob[min(which(w$y >= 3))],
    attr(w, "limits")$lower[min(which(w$y >= 3))],
    attr(w, "limits")$upper[min(which(w$y >= 3))],
    fit1.ep$p[min(which(fit1.ep$y >= 3))],
    fit1.ep$lower[min(which(fit1.ep$y >= 3))],
    fit1.ep$upper[min(which(fit1.ep$y >= 3))],
    fit2.ep$p[min(which(fit2.ep$y >= 3))],
    fit2.ep$lower[min(which(fit2.ep$y >= 3))],
    fit2.ep$upper[min(which(fit2.ep$y >= 3))],
    fit3.ep$p[min(which(fit3.ep$y >= 3))],
    fit3.ep$lower[min(which(fit3.ep$y >= 3))],
    fit3.ep$upper[min(which(fit3.ep$y >= 3))],
    fit4.ep$p[min(which(fit4.ep$y >= 3))],
    fit4.ep$lower[min(which(fit4.ep$y >= 3))],
    fit4.ep$upper[min(which(fit4.ep$y >= 3))],
    w$prob[min(which(w$y >= 3))],
    attr(w, "limits")$lower[min(which(w$y >= 6))],
    attr(w, "limits")$upper[min(which(w$y >= 6))],
    fit1.ep$p[min(which(fit1.ep$y >= 6))],
    fit1.ep$lower[min(which(fit1.ep$y >= 6))],
    fit1.ep$upper[min(which(fit1.ep$y >= 6))],
    fit2.ep$p[min(which(fit2.ep$y >= 6))],
    fit2.ep$lower[min(which(fit2.ep$y >= 6))],
    fit2.ep$upper[min(which(fit2.ep$y >= 6))],
    fit3.ep$p[min(which(fit3.ep$y >= 6))],
    fit3.ep$lower[min(which(fit3.ep$y >= 6))],
    fit3.ep$upper[min(which(fit3.ep$y >= 6))],
    fit4.ep$p[min(which(fit4.ep$y >= 6))],
    fit4.ep$lower[min(which(fit4.ep$y >= 6))],
    fit4.ep$upper[min(which(fit4.ep$y >= 6))],
    qY1, qY2, ep3, ep6,
    tail(diag(var1), 2),
    tail(diag(var2), 2),
    tail(diag(var3), 2),
    tail(diag(var4), 2)
    )
}

# Save rds for later use
#saveRDS(out, file = "simulations/results/A1.rds")
#out <- readRDS("")

# Using results from above, create Table A1.
dat1 <- data.frame(
  y = c(out[, 1:10]),
  method = rep(rep(c("True", "GR (0)", "GR (5)", "GR (10)", "GR (20)"), each = 2000), times = 1),
  type = rep(rep(c("beta1", "beta2"), each = 1000), times = 5)
)

dat1 %>% group_by(type, method) %>% summarise(
  x = round(mean(y), 3),
  sd = round(sd(y), 3))


# Get coverage probabilities

#beta1
mean(0.5 >= out[, 3] - qnorm(0.975) * sqrt(out[, 75]) & 0.5 <= out[, 3] + qnorm(0.975) * sqrt(out[, 75]))
mean(0.5 >= out[, 5] - qnorm(0.975) * sqrt(out[, 77]) & 0.5 <= out[, 5] + qnorm(0.975) * sqrt(out[, 77]))
mean(0.5 >= out[, 7] - qnorm(0.975) * sqrt(out[, 79]) & 0.5 <= out[, 7] + qnorm(0.975) * sqrt(out[, 79]))
mean(0.5 >= out[, 9] - qnorm(0.975) * sqrt(out[, 81]) & 0.5 <= out[, 9] + qnorm(0.975) * sqrt(out[, 81]))

#beta2
mean(-0.5 >= out[, 4] - qnorm(0.975) * sqrt(out[, 76]) & -0.5 <= out[, 4] + qnorm(0.975) * sqrt(out[, 76]))
mean(-0.5 >= out[, 6] - qnorm(0.975) * sqrt(out[, 78]) & -0.5 <= out[, 6] + qnorm(0.975) * sqrt(out[, 78]))
mean(-0.5 >= out[, 8] - qnorm(0.975) * sqrt(out[, 80]) & -0.5 <= out[, 8] + qnorm(0.975) * sqrt(out[, 80]))
mean(-0.5 >= out[, 10] - qnorm(0.975) * sqrt(out[, 82]) & -0.5 <= out[, 10] + qnorm(0.975) * sqrt(out[, 82]))

# Median
mean(out[, 71] >= out[, 15] & out[, 71] <= out[, 16])
mean(out[, 71] >= out[, 18] & out[, 71] <= out[, 19])
mean(out[, 71] >= out[, 21] & out[, 71] <= out[, 22])
mean(out[, 71] >= out[, 24] & out[, 71] <= out[, 25])

# 80th percentile
mean(out[, 72] >= out[, 30] & out[, 72] <= out[, 31])
mean(out[, 72] >= out[, 33] & out[, 72] <= out[, 34])
mean(out[, 72] >= out[, 36] & out[, 72] <= out[, 37])
mean(out[, 72] >= out[, 39] & out[, 72] <= out[, 40])

# Pr(Y > 3)
mean(out[, 73] >= out[, 45] & out[, 73] <= out[, 46])
mean(out[, 73] >= out[, 48] & out[, 73] <= out[, 49])
mean(out[, 73] >= out[, 51] & out[, 73] <= out[, 52])
mean(out[, 73] >= out[, 54] & out[, 73] <= out[, 55])


# Pr(Y > 6)
mean(out[, 74] >= out[, 60] & out[, 74] <= out[, 61])
mean(out[, 74] >= out[, 63] & out[, 74] <= out[, 64])
mean(out[, 74] >= out[, 66] & out[, 74] <= out[, 67])
mean(out[, 74] >= out[, 69] & out[, 74] <= out[, 70])

