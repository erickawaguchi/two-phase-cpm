rm(list = ls())
library(rms)
library(dplyr)
library(foreach)
library(readxl)
library(ggplot2)
library(ggpubr)

# Update directory here
data_dir <- "/Users/ekawaguc/GitHub/svycpm/"

idx = 1
# idx = 1: Base case simulation
# 2, 3: Varying n scenario (Figure 6 in the manuscript)
# 4, 5: Varying Phase 1 / Phase 2 correlation (Table 3 in the manuscript)

out <- readRDS(paste0(data_dir, "simulations/results/sim_", idx, ".rds"))

apply(out, 2, mean)

# beta estimates
d1 <- data.frame(
  y = c(out[, c(1:16)]),
  y2 = c(out[, c(1:16)]) - rep(rep(c(0.5, -0.5), each = 1000), times = 8),
  method = rep(c("True CPM", "P1 Only", "P2 Only",
                 "IPW CPM",
                 "GR CPM",
                 "GR LR (log)", 
                 "GR LR (C)",
                 "True LR"
                 ),
               each = 2000
  ),
  x = rep(rep(c("hat(beta)[1]", "hat(beta)[2]"), each = 1000), times = 8)
)


d1$method <- factor(d1$method, levels = c("True CPM", "True LR",
                                          "P1 Only", "P2 Only",
                                          "IPW CPM",
                                          "GR CPM",
                                          "GR LR (log)", 
                                          "GR LR (C)"))

p1 <- ggplot(data = d1, aes(y = y2, x = method)) +
  geom_boxplot() +
  facet_grid(~x, labeller = label_parsed) + 
  #ylab(expression(beta~"-"~hat(beta))) +
  ylab("Bias") +
  ylim(c(-0.4, 0.4)) +
  geom_abline(slope = 0, linetype = "dashed", col = "darkgrey", linewidth = 1) +
  xlab("")+
  theme_bw() +
  theme(
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 12),
    strip.text = element_text(size = 14)
  )

p1

# SE estimates
d2 <- data.frame(
  y = c(out[, c(17, 18, 23, 24, 25, 26, 29, 30, 31, 32)]),
  method = rep(c("True CPM",
                 "IPW CPM",
                 "GR CPM",
                 "GR LR (C)",
                 "True LR"),
               each = 2000
  ),
  x = rep(rep(c("hat(beta)[1]", "hat(beta)[2]"), each = 1000), times = 5)
)


d2$method <- factor(d2$method, levels = c("True CPM",
                                          "True LR",
                                          "IPW CPM",
                                          "GR CPM",
                                          "GR LR (C)"))

p2 <- ggplot(data = d2, aes(y = y, x = method)) +
  geom_boxplot() +
  facet_grid(~x, labeller = label_parsed) + 
  ylab(expression(SE)) +
  xlab("") + 
  #abline(h = 0.5) +
  theme_bw() +
  theme(
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    #axis.title.x = element_blank(),
    axis.title.y = element_text(size = 12),
    strip.text = element_text(size = 14)
  )

# 1800 x 900
ggarrange(p1, p2, nrow = 2, labels = c("A)", "B)"))

#########################################
#########################################
# Get figures for comparing beta values across scenarios
data_dir <- "/Users/ekawaguc/Desktop/CPM/Survey\ Sampling/"

simList <- read_xlsx(paste0(data_dir, "simulations/simulationList.xlsx"))

res <- foreach(j = c(1, 2, 3), .combine = 'rbind') %do% {
  simVals <- simList[j, ]
  out     <- readRDS(paste0(data_dir, "sims/results/sim_", j, ".rds"))
  data.frame(
    y = c(out[, 7:10]),  
    y2 = c(out[, 7:10]) - rep(rep(c(0.5, -0.5), each = 1000), times = 2),
    
    x = rep(rep(c("hat(beta)[1]", "hat(beta)[2]"), each = 1000), times = 2),
    method = rep(c("IPW", "GR"), each = 2000),
    n = factor(rep(simVals$n, 4000)))
}

res$n <- factor(res$n, levels = c(300, 600, 900))
res %>% group_by(method, n, x) %>% summarise(sd = sd(y)) %>% filter(method == "GR")

# Figure 6
ggplot(res %>% filter(method == "GR"), aes(y = y2, x = n)) + 
  geom_boxplot() +
  facet_wrap(~x, labeller = label_parsed)  + 
  ylab("Bias") +
  scale_x_discrete(
    labels = c(
      "300" = "n = 300",
      "600" = "n = 600",
      "900" = "n = 900")
  ) +
  #xlab("Phase 2 Sample Size (n)") +
  #abline(h = 0.5) +
  theme_bw() +
  theme(
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 12),
    strip.text = element_text(size = 14)
  )
# 800 x 600


res <- foreach(j = c(5, 4, 1), .combine = 'rbind') %do% {
  simVals <- simList[j, ]
  out     <- readRDS(paste0(data_dir, "sims/results/sim_", j, ".rds"))
  data.frame(
    y = c(out[, 7:10]),  
    x = rep(rep(c("beta1", "beta2"), each = 1000), times = 2),
    method = rep(c("IPW", "GR"), each = 2000),
    n = factor(rep(simVals$model, 4000)),
    m = paste0("rho_x = ", simVals$r_x, "; Se = ", simVals$Se, 
               "; Sp = ", simVals$Sp))
}

# Table 3
res %>% group_by(method, x, m) %>% summarise(sd = round(sd(y), 3))

#########################################
#########################################
# Table 2

# Load rds file from sim_prob.R
#out <- readRDS(".../sim_p.rds")
out <- as.data.frame(out)

#beta1 

# GR
mean(out[, 1] - 0.5)
sd(out[, 1])
mean((out[, 1] - qnorm(0.975) * out[, 5]) < 0.5 & 
      0.5 < (out[, 1] + qnorm(0.975) * out[, 5]))

# IPW
mean(out[, 3] - 0.5)
sd(out[, 3])
mean((out[, 3] - qnorm(0.975) * out[, 7]) < 0.5 & 
       0.5 < (out[, 3] + qnorm(0.975) * out[, 7]))

# GR
mean(out[, 2] - -0.5)
sd(out[, 2])
mean((out[, 2] - qnorm(0.975) * out[, 6]) < -0.5 & 
       -0.5 < (out[, 2] + qnorm(0.975) * out[, 6]))

# IPW
mean(out[, 4] - -0.5)
sd(out[, 4])
mean((out[, 4] - qnorm(0.975) * out[, 8]) < -0.5 & 
       -0.5 < (out[, 4] + qnorm(0.975) * out[, 8]))


q5 <- out[, 33]
q8 <- out[, 34]
ep3 <- out[, 35]
ep6 <- out[, 36]

# Median:

#GR
mean(out[, 9] - q5)
sd(out[, 9])
mean(out[, 10] < q5 & out[, 11] > q5)

#IPW
mean(out[, 12] - q5)
sd(out[, 12])
mean(out[, 13] < q5 & out[, 14] > q5)

# 80th percentile:

#GR
mean(out[, 15] - q8)
sd(out[, 15])
mean(out[, 16] < q8 & out[, 17] > q8)

#IPW
mean(out[, 18] - q8)
sd(out[, 18])
mean(out[, 19] < q8 & out[, 20] > q8)


# Pr(Y > 3)

#GR
mean(out[, 21] - ep3)
sd(out[, 21])
mean(out[, 22] < ep3 & out[, 23] > ep3)

#IPW
mean(out[, 24] - ep3)
sd(out[, 24])
mean(out[, 25] < ep3 & out[, 26] > ep3)

# Pr(Y > 6)

#GR
mean(out[, 27] - ep6)
sd(out[, 27])
mean(out[, 28] < ep6 & out[, 29] > ep6)

#IPW
mean(out[, 30] - ep6)
sd(out[, 30])
mean(out[, 31] < ep6 & out[, 32] > ep6)
