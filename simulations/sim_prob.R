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

set.seed(simVals$seed)
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
  
  # Fit following models
  
  # Cumulative probability models
  
  # True model
  fit.m1 <- orm(yt.true ~ x1 + x2, data = d,
                  family = "probit")
  
  # Fit on only phase 1
  fit.p1 <- orm(yt.p1 ~ x1.p1 + x2.p1, data = d, mscore= TRUE,
                family = "probit")
  
  d$naiveInfl.x1 <- getBetaInfluenceFunctions(fit.p1)$C2[, 1]
  d$naiveInfl.x2 <- getBetaInfluenceFunctions(fit.p1)$C2[, 2]
  
  # Fit two-phase design
  designO <- twophase(id=list(~1,~1), strata = list(NULL, ~strata), subset = ~in.p2 ==1,
                      data = d, method = "approx")
  
  # IPW
  sdat <- model.frame(designO)
  sdat$.survey.prob.weights <- (1 / designO$prob) / mean(1 / designO$prob)
  
  fit.ipw <- orm(yt.true ~ x1.p2 + x2.p2, 
                 data = sdat,
                 mscore = TRUE,
                 weights = .survey.prob.weights,
                 family = "probit")
  inflFun  <- getBetaInfluenceFunctions(fit.ipw, getC1 = TRUE)
  ipw.var <- twophasevar(cbind(inflFun$C1, inflFun$C2), designO)
  
  # Generalized raking
  Ninflcal.orm <- survey::calibrate(designO,formula=~naiveInfl.x1 + naiveInfl.x2,phase=2,calfun="raking") 
  sdat <- model.frame(Ninflcal.orm)
  sdat$.survey.prob.weights <- (1 / Ninflcal.orm$prob) / mean(1 / Ninflcal.orm$prob)
  
  fit.gr <- orm(yt.true ~ x1.p2 + x2.p2, 
                data = sdat,
                mscore = TRUE,
                weights = .survey.prob.weights,
                family = "probit")
  
  inflFun <- getBetaInfluenceFunctions(fit.gr, getC1 = TRUE)
  gr.var  <- twophasevar(cbind(inflFun$C1, inflFun$C2), Ninflcal.orm)
  
  
  # Get predicted quantiles and exceedence probs based on X1 = -1 and X2 = 1
  d <- data.frame(x1.p2 = -1, x2.p2 = 1)
  
  ipw.q1  <- pred_quantile(fit.ipw, ipw.var, d, tau = 0.5, method = "discrete")
  ipw.q2  <- pred_quantile(fit.ipw, ipw.var, d, tau = 0.8, method = "discrete")
  ipw.ep  <- pred_ExProb(fit.ipw, ipw.var, d)
  
  gr.q1  <- pred_quantile(fit.gr, gr.var, d, tau = 0.5, method = "discrete")
  gr.q2  <- pred_quantile(fit.gr, gr.var, d, tau = 0.8, method = "discrete")
  gr.ep  <- pred_ExProb(fit.gr, gr.var, d)
  
  c(
    # To confirm with original simulation
    tail(fit.gr$coefficients, 2),
    tail(fit.ipw$coefficients, 2),
    tail(sqrt(diag(gr.var)), 2),
    tail(sqrt(diag(ipw.var)), 2),
    # New results
    gr.q1$p,
    gr.q1$lower,
    gr.q1$upper,
    ipw.q1$p,
    ipw.q1$lower,
    ipw.q1$upper,
    gr.q2$p,
    gr.q2$lower,
    gr.q2$upper,
    ipw.q2$p,
    ipw.q2$lower,
    ipw.q2$upper,
    gr.ep$p[min(which(gr.ep$y >= 3))],
    gr.ep$lower[min(which(gr.ep$y >= 3))],
    gr.ep$upper[min(which(gr.ep$y >= 3))],
    ipw.ep$p[min(which(ipw.ep$y >= 3))],
    ipw.ep$lower[min(which(ipw.ep$y >= 3))],
    ipw.ep$upper[min(which(ipw.ep$y >= 3))],  
    gr.ep$p[min(which(gr.ep$y >= 6))],
    gr.ep$lower[min(which(gr.ep$y >= 6))],
    gr.ep$upper[min(which(gr.ep$y >= 6))],
    ipw.ep$p[min(which(ipw.ep$y >= 6))],
    ipw.ep$lower[min(which(ipw.ep$y >= 6))],
    ipw.ep$upper[min(which(ipw.ep$y >= 6))],
    qY1,
    qY2,
    ep3,
    ep6)
}
saveRDS(out, file = paste0(data_dir, "simulations/results/sim_p.rds"))


