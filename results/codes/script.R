library(dplyr)

# simulate data -----------------------------------------------------------


data <- simBAMMIT(V = 4,
                  Q = 3,
                  Bv = c(12,10, 4, 2),
                  mu = 100,
                  lambda = c(8,10, 12),
                  sb = 1,
                  sB = 1,
                  sy = 1)
# N = B1*B2*...*Bv
#saveRDS(data, file = "Simulated data/V4_Q3_N960_L81012.rds")

# read the data

V2_Q1_N120_L10 <- readRDS("~/Documents/GitHub/bammit/Simulated data/V2_Q1_N120_L10.rds")
V3_Q1_N480_L10 <- readRDS("~/Documents/GitHub/bammit/Simulated data/V3_Q1_N480_L10.rds")
V4_Q1_N960_L10 <- readRDS("~/Documents/GitHub/bammit/Simulated data/V4_Q1_N960_L10.rds")
V2_Q2_N120_L1012 <- readRDS("~/Documents/GitHub/bammit/Simulated data/V2_Q2_N120_L1012.rds")
V3_Q2_N480_L1012 <- readRDS("~/Documents/GitHub/bammit/Simulated data/V3_Q2_N480_L1012.rds")
V4_Q2_N960_L1012 <- readRDS("~/Documents/GitHub/bammit/Simulated data/V4_Q2_N960_L1012.rds")
V2_Q3_N120_L81012 <- readRDS("~/Documents/GitHub/bammit/Simulated data/V2_Q3_N120_L81012.rds")
V3_Q3_N480_L81012 <- readRDS("~/Documents/GitHub/bammit/Simulated data/V3_Q3_N480_L81012.rds")
V4_Q3_N960_L81012 <- readRDS("~/Documents/GitHub/bammit/Simulated data/V4_Q3_N960_L81012.rds")

# jags code ---------------------------------------------------------------

data <- V2_Q1_N120_L10
model <- bammitJags(data = data,
                    Q = 2,
                    mmu = 100,
                    smu = 10,
                    a = 0.1,
                    b = 0.1,
                    nthin = 1,
                    nburnin = 2000)

saveRDS(model, file = "Running models/V2_Q2_N960_L10.rds")

plot(data$int, model$BUGSoutput$mean$int)


# plots --------------------------------------------------------------
library(reshape2)
library(ggplot2)
# compare V's simulated Q = 1 and run model Q = 3

V2_Q3_N120_L10 <- readRDS("~/Documents/GitHub/bammit/Running models/V2_Q3_N120_L10.rds")
V3_Q3_N480_L10 <- readRDS("~/Documents/GitHub/bammit/Running models/V3_Q3_N480_L10.rds")
V4_Q3_N960_L10 <- readRDS("~/Documents/GitHub/bammit/Running models/V4_Q3_N960_L10.rds")

# compare bilinear
VQ3int <- list(V2 = V2_Q3_N120_L10$BUGSoutput$mean$int,
            V3 = V3_Q3_N480_L10$BUGSoutput$mean$int,
            V4 = V4_Q3_N960_L10$BUGSoutput$mean$int) |>
  melt() |>
  select(value, L1)

quantVQ3int <- list(
  V2 = data.frame(q1 = apply(V2_Q3_N120_L10$BUGSoutput$sims.list$int, 2, function(x) quantile(x, 0.05)),
                  q2 = apply(V2_Q3_N120_L10$BUGSoutput$sims.list$int, 2, function(x) quantile(x, 0.95))),
  V3 = data.frame(q1 = apply(V3_Q3_N480_L10$BUGSoutput$sims.list$int, 2, function(x) quantile(x, 0.05)),
                  q2 = apply(V3_Q3_N480_L10$BUGSoutput$sims.list$int, 2, function(x) quantile(x, 0.95))),
  V4 = data.frame(q1 = apply(V4_Q3_N960_L10$BUGSoutput$sims.list$int, 2, function(x) quantile(x, 0.05)),
                  q2 = apply(V4_Q3_N960_L10$BUGSoutput$sims.list$int, 2, function(x) quantile(x, 0.95)))
) |> plyr::ldply()


# quantVQ3int <- list(
#   V2 = data.frame(q1 = apply(V2_Q3_N120_L10$BUGSoutput$sims.list$int, 2, function(x) 1.96*sd(x)),
#                   q2 = apply(V2_Q3_N120_L10$BUGSoutput$sims.list$int, 2, function(x) 1.96*sd(x))),
#   V3 = data.frame(q1 = apply(V3_Q3_N480_L10$BUGSoutput$sims.list$int, 2, function(x) 1.96*sd(x)),
#                   q2 = apply(V3_Q3_N480_L10$BUGSoutput$sims.list$int, 2, function(x) 1.96*sd(x))),
#   V4 = data.frame(q1 = apply(V4_Q3_N960_L10$BUGSoutput$sims.list$int, 2, function(x) 1.96*sd(x)),
#                   q2 = apply(V4_Q3_N960_L10$BUGSoutput$sims.list$int, 2, function(x) 1.96*sd(x)))
# ) |> plyr::ldply()

VQ3int <- cbind(VQ3int, quantVQ3int$q1, quantVQ3int$q2)
colnames(VQ3int) <- c("value", "L1", "q1", "q2")

VQ3Dataint <- list(V2 = V2_Q1_N120_L10$int,
                      V3 = V3_Q1_N480_L10$int,
                      V4 = V4_Q1_N960_L10$int) |>
  melt() |>
  select(value, L1)

dfVQ3int <- data.frame(v = VQ3int$L1, est = VQ3int$value, true = VQ3Dataint$value,
                        q1 = VQ3int$q1, q2 = VQ3int$q2)

# LVQ3int <- list(VQ3int, VQ3Dataint)
# names(LVQ3int) <- c("Estimated", "True")
# dfVQ3int <- plyr::ldply(LVQ3int)
# colnames(dfVQ3int) <- c("data","value", "v")


lab <- function(x) c("V = 2", "V = 3", "V = 4")
dfVQ3int |>
  ggplot(aes(x = true, y = est)) +
  geom_abline(size = 0.3, col = "gray") +
  geom_linerange(aes(ymin =  q1, ymax = q2), alpha = 0.2, size = 0.1) +
  geom_point(colour =  "steelblue", size = 0.2) +
  facet_wrap(~v,
             labeller = as_labeller(lab)) +
  labs(x = "True", y = "Estimated") +
  theme_bw(base_size = 14) +
  theme(strip.background = element_blank(),
        panel.spacing.x = unit(5,"mm"),
        legend.position = c(0.75, 0.04))


# considering just V = 4

V4_Q3_N960_L10 <- readRDS("~/Documents/GitHub/bammit/Running models/V4_Q3_N960_L10.rds")
V4_Q1_N960_L10 <- readRDS("~/Documents/GitHub/bammit/Simulated data/V4_Q1_N960_L10.rds")

dfEstV4 <- list(b1 = V4_Q3_N960_L10$BUGSoutput$mean$b1,
                b2 = V4_Q3_N960_L10$BUGSoutput$mean$b2,
                b3 = V4_Q3_N960_L10$BUGSoutput$mean$b3,
                b4 = V4_Q3_N960_L10$BUGSoutput$mean$b4) |>
  melt() |> select(L1, value)


dfDataV4 <- list(b1 = V4_Q1_N960_L10$bv[[1]],
                 b2 = V4_Q1_N960_L10$bv[[2]],
                 b3 = V4_Q1_N960_L10$bv[[3]],
                 b4 = V4_Q1_N960_L10$bv[[4]]) |>
  melt() |> select(L1, value)

dfMainEffV4 <- cbind(dfEstV4, dfDataV4)[,c(1,2,4)]
colnames(dfMainEffV4) <- c("var", "est", "true")

quantV4mainEff <- list(
  b1 = data.frame(q1 = apply(V4_Q3_N960_L10$BUGSoutput$sims.list$b1, 2, function(x) quantile(x, 0.05)),
                  q2 = apply(V4_Q3_N960_L10$BUGSoutput$sims.list$b1, 2, function(x) quantile(x, 0.95))),
  b2 = data.frame(q1 = apply(V4_Q3_N960_L10$BUGSoutput$sims.list$b2, 2, function(x) quantile(x, 0.05)),
                  q2 = apply(V4_Q3_N960_L10$BUGSoutput$sims.list$b2, 2, function(x) quantile(x, 0.95))),
  b3 = data.frame(q1 = apply(V4_Q3_N960_L10$BUGSoutput$sims.list$b3, 2, function(x) quantile(x, 0.05)),
                  q2 = apply(V4_Q3_N960_L10$BUGSoutput$sims.list$b3, 2, function(x) quantile(x, 0.95))),
  b4 = data.frame(q1 = apply(V4_Q3_N960_L10$BUGSoutput$sims.list$b4, 2, function(x) quantile(x, 0.05)),
                  q2 = apply(V4_Q3_N960_L10$BUGSoutput$sims.list$b4, 2, function(x) quantile(x, 0.95)))
) |> plyr::ldply()

dfMainEffV4 <- cbind(dfMainEffV4, quantV4mainEff$q1, quantV4mainEff$q2 )
colnames(dfMainEffV4) <- c("var", "est", "true", "q1", "q2")


dfMainEffV4$facet <- factor(dfMainEffV4$var,
                            labels = c(expression(bold(b)^{(1)}),
                                       expression(bold(b)^{(2)}),
                                       expression(bold(b)^{(3)}),
                                       expression(bold(b)^{(4)})))

dfMainEffV4 |>
  ggplot(aes(x = true, y = est)) +
  geom_abline(size = 0.3, col = "gray") +
  geom_linerange(aes(ymin =  q1, ymax = q2), alpha = 0.5, size = 0.2) +
  geom_point(colour =  "steelblue", size = 1) +
  facet_wrap(~facet, ncol = 4, labeller = label_parsed) +
  labs(x = "True", y = "Estimated") +
  theme_bw(base_size = 14) +
  theme(strip.background = element_blank(),
        panel.spacing.x = unit(5,"mm"),
        legend.position = c(0.75, 0.04))


# plot V and Q

V2_Q1_N120_L10R <- V2_Q1_N120_L10
V3_Q1_N480_L10R <- V3_Q1_N480_L10
V4_Q1_N960_L10R <- V4_Q1_N960_L10
V2_Q2_N120_L1012R  <- V2_Q2_N120_L1012
V3_Q2_N480_L1012R <- V3_Q2_N480_L1012
V4_Q2_N960_L1012R <- V4_Q2_N960_L1012
V2_Q3_N120_L81012R <- V2_Q3_N120_L81012
V3_Q3_N480_L81012R <- V3_Q3_N480_L81012
V4_Q3_N960_L81012R <- V4_Q3_N960_L81012

# load models
V2_Q1_N120_L10 <- readRDS("~/Documents/GitHub/bammit/Running models/V2_Q1_N120_L10.rds")
V2_Q2_N120_L10 <- readRDS("~/Documents/GitHub/bammit/Running models/V2_Q2_N960_L10.rds")
V2_Q3_N120_L10 <- readRDS("~/Documents/GitHub/bammit/Running models/V2_Q3_N120_L10.rds")
V3_Q1_N480_L10 <- readRDS("~/Documents/GitHub/bammit/Running models/V3_Q1_N960_L10.rds")
V3_Q2_N480_L10 <- readRDS("~/Documents/GitHub/bammit/Running models/V3_Q2_N960_L10.rds")
V3_Q3_N480_L10 <- readRDS("~/Documents/GitHub/bammit/Running models/V3_Q3_N480_L10.rds")
V4_Q1_N960_L10 <- readRDS("~/Documents/GitHub/bammit/Running models/V4_Q1_N960_L10.rds")
V4_Q2_N960_L10 <- readRDS("~/Documents/GitHub/bammit/Running models/V4_Q2_N960_L10.rds")
V4_Q3_N960_L10 <- readRDS("~/Documents/GitHub/bammit/Running models/V4_Q3_N960_L10.rds")



# data.frame(V = rep("V2", 120),
#            Q = rep("Q1", 120),
#            int = V2_Q1_N120_L10$BUGSoutput$mean$int,
#            q1 = apply(V2_Q1_N120_L10$BUGSoutput$sims.list$int, 2, function(x) quantile(x, 0.05)),
#            q2 = apply(V2_Q1_N120_L10$BUGSoutput$sims.list$int, 2, function(x) quantile(x, 0.95)))





bV2Q1 <- data.frame(V = rep("V2", 120),
           Q = rep("Q1", 120),
           int = V2_Q1_N120_L10$BUGSoutput$mean$int,
           intR = V2_Q1_N120_L10R$int)

bV2Q2 <- data.frame(V = rep("V2", 120),
                    Q = rep("Q2", 120),
                    int = V2_Q2_N120_L10$BUGSoutput$mean$int,
                    intR = V2_Q1_N120_L10R$int)

bV2Q3 <- data.frame(V = rep("V2", 120),
                    Q = rep("Q3", 120),
                    int = V2_Q3_N120_L10$BUGSoutput$mean$int,
                    intR = V2_Q1_N120_L10R$int)


bV3Q1 <- data.frame(V = rep("V3", 480),
                    Q = rep("Q1", 480),
                    int = V3_Q1_N480_L10$BUGSoutput$mean$int,
                    intR = V3_Q1_N480_L10R$int)

bV3Q2 <- data.frame(V = rep("V3", 480),
                    Q = rep("Q2", 480),
                    int = V3_Q2_N480_L10$BUGSoutput$mean$int,
                    intR = V3_Q1_N480_L10R$int)

bV3Q3 <- data.frame(V = rep("V3", 480),
                    Q = rep("Q3", 480),
                    int = V3_Q3_N480_L10$BUGSoutput$mean$int,
                    intR = V3_Q1_N480_L10R$int)

bV4Q1 <- data.frame(V = rep("V4", 960),
                    Q = rep("Q1", 960),
                    int = V4_Q1_N960_L10$BUGSoutput$mean$int,
                    intR = V4_Q1_N960_L10R$int)

bV4Q2 <- data.frame(V = rep("V4", 960),
                    Q = rep("Q2", 960),
                    int = V4_Q2_N960_L10$BUGSoutput$mean$int,
                    intR = V4_Q1_N960_L10R$int)

bV4Q3 <- data.frame(V = rep("V4", 960),
                    Q = rep("Q3", 960),
                    int = V4_Q3_N960_L10$BUGSoutput$mean$int,
                    intR = V4_Q1_N960_L10R$int)


dfb <- rbind(bV2Q1, bV2Q2, bV2Q3, bV3Q1, bV3Q2, bV3Q3, bV4Q1, bV4Q2, bV4Q3)
lab <- function(x) c("Q = 1", "Q = 2", "Q = 3")
dfb$V <- factor(dfb$V,labels = c("2", "3", "4"))

dfb |>
  ggplot(aes(x = int, col = V, fill = V)) +
  geom_histogram(aes(x = intR, y=..density..), colour="black", fill="white") +
  geom_density(alpha = 0.2) +
  #geom_density(aes(x = intR, col = "gray"), alpha = 0.2) +
  facet_wrap(~Q, labeller = as_labeller(lab)) +
  labs(x = "Interaction", y = "Density") +
  theme_bw(base_size = 14) +
  theme(strip.background = element_blank(),
        panel.spacing.x = unit(5,"mm"))



dfb |>
  ggplot(aes(x = int, col = V, fill = V)) +
  geom_boxplot(aes(x = intR, col = V, fill = V), colour="black", fill="white") +
  geom_boxplot(alpha = 0.2) +
  #geom_density(aes(x = intR, col = "gray"), alpha = 0.2) +
  facet_wrap(~Q, labeller = as_labeller(lab)) +
  labs(x = "Interaction", y = "Density") +
  theme_bw(base_size = 14) +
  theme(strip.background = element_blank(),
        panel.spacing.x = unit(5,"mm"))




# sum of the interactions -------------------------------------------------



bVQ <- data.frame(V = c("V2", "V2", "V2", "V3", "V3", "V3", "V4", "V4", "V4"),
                    Q = c("Q1", "Q2", "Q3", "Q1", "Q2", "Q3", "Q1", "Q2", "Q3"),
                    int = c(sum(V2_Q1_N120_L10$BUGSoutput$mean$int),
                             sum(V2_Q2_N120_L10$BUGSoutput$mean$int),
                             sum(V2_Q3_N120_L10$BUGSoutput$mean$int),
                             sum(V3_Q1_N480_L10$BUGSoutput$mean$int),
                             sum(V3_Q2_N480_L10$BUGSoutput$mean$int),
                             sum(V3_Q3_N480_L10$BUGSoutput$mean$int),
                             sum(V4_Q1_N960_L10$BUGSoutput$mean$int),
                             sum(V4_Q2_N960_L10$BUGSoutput$mean$int),
                             sum(V4_Q3_N960_L10$BUGSoutput$mean$int)),
                    intR = c(sum(V2_Q1_N120_L10R$int),
                              sum(V2_Q1_N120_L10R$int),
                              sum(V2_Q1_N120_L10R$int),
                              sum(V3_Q1_N480_L10R$int),
                              sum(V3_Q1_N480_L10R$int),
                              sum(V3_Q1_N480_L10R$int),
                              sum(V4_Q1_N960_L10R$int),
                              sum(V4_Q1_N960_L10R$int),
                              sum(V4_Q1_N960_L10R$int)))


p1 <- bVQ |>
  ggplot(aes(x = int, y = intR, col = V, shape = Q)) +
  geom_abline(size = 0.3, col = "gray") +
  geom_point(size = 4) +
  labs(x = "True", y = "Estimated") +
  theme_bw(base_size = 14)

p2 <- bVQ |>
ggplot(aes(x = int, y = intR, col = V)) +
  geom_abline(size = 0.3, col = "gray") +
  geom_point(size = 4) +
  labs(x = "True", y = "Estimated") +
  theme_bw(base_size = 14) +
  facet_wrap(~Q)


gridExtra::grid.arrange(p1,p2)

