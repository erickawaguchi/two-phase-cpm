rm(list = ls())
library(rms)
library(dplyr)
library(survey)
library(foreach)
library(readxl)
library(sandwich)


# Update directory here
data_dir <- "/Users/ekawaguc/GitHub/svycpm/"

# Source files:
source(paste0(data_dir, "source/getInfluenceFunctions.R"))
simList <- read_xlsx(paste0(data_dir, "simulations/simulationList.xlsx"))


simIdx <- 1
simVals <- simList[simIdx, ]

# simIdx = 1: Base case simulation
# 2, 3: Varying n scenario (Figure 6 in the manuscript)
# 4, 5: Varying Phase 1 / Phase 2 correlation (Table 3 in the manuscript)

set.seed(simVals$seed)

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

out <- foreach(b = 1:B, .combine = 'rbind') %do% {
  
  cat(paste(b, "out of", B, "\n"))
  
  # Generate covariates, (true) X1 and X2
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
  if(rho_x != 0) {
    x1.p1 <- x1 + rnorm(N, mean = 0, sd = s_x)
  } else {
    x1.p1 <- rnorm(N, mean = 0, sd = s_true)
  }
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
  
  # Fit models for comparison
  
  # Cumulative probability models
  
  # Model 1 (TRUE)
  fit.m1 <- orm(yt.true ~ x1 + x2, data = d,
                  family = "probit")
  
  # Model 2 (Phase 1 only)
  fit.m2 <- orm(yt.p1 ~ x1.p1 + x2.p1, data = d, mscore= TRUE,
                family = "probit")
  
  # Model 3 (Phase 2 only)
  fit.m3 <- orm(yt.true ~ x1.p2 + x2.p2, data = d, 
                subset = (in.p2 == 1),
                family = "probit")
  
  d$naiveInfl.x1 <- getBetaInfluenceFunctions(fit.m2)$C2[, 1]
  d$naiveInfl.x2 <- getBetaInfluenceFunctions(fit.m2)$C2[, 2]
  
  # Fit two-phase design
  designO <- twophase(id=list(~1,~1), strata = list(NULL, ~strata), subset = ~in.p2 ==1,
                      data = d, method = "approx")
  
  # Fit IPW model (using SRS weights)
  sdat <- model.frame(designO)
  sdat$.survey.prob.weights <- (1 / designO$prob) / mean(1 / designO$prob)
  
  
  # IPW model
  fit.m4 <- orm(yt.true ~ x1.p2 + x2.p2, 
                 data = sdat,
                 mscore = TRUE,
                 weights = .survey.prob.weights,
                 family = "probit")
  
  inflFun  <- getBetaInfluenceFunctions(fit.m4, getC1 = FALSE)
  m4.var   <- twophasevar(inflFun$C2, designO)
  
  
  # Proposed GR model
  Ninflcal.orm <- survey::calibrate(designO,formula=~naiveInfl.x1 + naiveInfl.x2,phase=2,calfun="raking") 
  sdat <- model.frame(Ninflcal.orm)
  sdat$.survey.prob.weights <- (1 / Ninflcal.orm$prob) / mean(1 / Ninflcal.orm$prob)
  
  fit.m5 <- orm(yt.true ~ x1.p2 + x2.p2, 
                data = sdat,
                mscore = TRUE,
                weights = .survey.prob.weights,
                family = "probit")
  
  inflFun <- getBetaInfluenceFunctions(fit.m5, getC1 = FALSE)
  m5.var  <- twophasevar(inflFun$C2, Ninflcal.orm)

  # Design-based linear regression, proper and improper transform
  
  # Improper
  fit.lm1 <- glm(log(yt.p1) ~ x1.p1 + x2.p1, data = d)
  
  IF  <- estfun(fit.lm1) %*% t(bread(fit.lm1))
  d$if1 <- IF[, 2]
  d$if2 <- IF[, 3]
  
  designO <- twophase(id=list(~1,~1), strata = list(NULL, ~strata), subset = ~in.p2 ==1,
                      data = d, method = "approx")
  Ninflcal.lm1 <-survey::calibrate(designO,formula=~if1 + if2,phase=2,calfun="raking") 
  
  fit.m6 <- svyglm(log(yt.true) ~ x1.p2 + x2.p2, design = Ninflcal.lm1)
  
  
  # Proper
  fit.lm2 <- glm(y2 ~ x1.p1 + x2.p1, data = d)
  IF  <- estfun(fit.lm2) %*% t(bread(fit.lm2))
  d$if1 <- IF[, 2]
  d$if2 <- IF[, 3]
  designO <- twophase(id=list(~1,~1), strata = list(NULL, ~strata), subset = ~in.p2 ==1,
                      data = d, method = "approx")
  Ninflcal.lm2 <-survey::calibrate(designO,formula=~if1 + if2,phase=2,calfun="raking") 
  
  fit.m7 <- svyglm(y ~ x1.p2 + x2.p2, design = Ninflcal.lm2)

  # Model 8 (true gold standard)
  fit.m8 <- glm(y ~ x1.true + x2.true, data = d)
  
  tmp <- c(table(x2, x2.p1))
  
  # Store all results:
  c(
    # Betas
    tail(fit.m1$coefficients, 2),
    tail(fit.m2$coefficients, 2),
    tail(fit.m3$coefficients, 2),
    tail(fit.m4$coefficients, 2),
    tail(fit.m5$coefficients, 2),
    fit.m6$coefficients[-1],
    fit.m7$coefficients[-1],
    fit.m8$coefficients[-1],
    # SE
    sqrt(diag(infoMxop(fit.m1$info.matrix, i = "x"))),
    sqrt(diag(infoMxop(fit.m2$info.matrix, i = "x"))),
    sqrt(diag(infoMxop(fit.m3$info.matrix, i = "x"))),
    sqrt(diag(m4.var)),
    sqrt(diag(m5.var)),
    SE(fit.m6)[-1],
    SE(fit.m7)[-1],
    SE(fit.m8)[-1],
    # Information about data
    cor(x1, x1.p1),
    cor(x2, x2.p1),
    tmp[4] / sum(tmp[c(2, 4)]),
    tmp[1] / sum(tmp[c(1, 3)]),
    cor(y, y2),
    cor(yt, yt2)
  )
}

saveRDS(out, file = paste0(data_dir, "simulations/results/sim_", simVals$model, ".rds"))
