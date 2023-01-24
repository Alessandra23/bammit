# Simulation study of AR-BAMMIT model

library(R2jags)
library(ggplot2)
library(reshape2)
library(patchwork)
library(bammit)

set.seed(02)
data <- simArBAMMIT(
  V = 3,
  Q = 1,
  Bv = c(12, 10, 4),
  mu = 10,
  lambda = c(12),
  sb = 1,
  sB = 1,
  seta = 1,
  somega = 1,
  sy = 1
)

ARbammitModel <- ARbammitJags(
  data = data, Q = 1, mmu = 10, smu = 10, nthin = 2, nburnin = 100, stheta = 1, somega = 1, a = 0.001, b = 0.001
)



# Plots for the interaction
dataInt <- data$int
modelInt <- ARbammitModel$BUGSoutput$mean$int
ciInt <- data.frame(
  q1 = apply(ARbammitModel$BUGSoutput$sims.list$int, 2, function(x) quantile(x, 0.05)),
  q2 = apply(ARbammitModel$BUGSoutput$sims.list$int, 2, function(x) quantile(x, 0.95))
)


dfint <- data.frame(est = modelInt, true = dataInt, q1 = ciInt$q1, q2 = ciInt$q2)

pBlin <- dfint |>
  ggplot(aes(x = true, y = est)) +
  geom_abline(size = 0.3, col = "gray") +
  geom_linerange(aes(ymin = q1, ymax = q2), alpha = 0.2, linewidth = 0.1) +
  geom_point(colour = "steelblue", size = 0.2) +
  labs(x = "True", y = "Estimated") +
  theme_bw(base_size = 14)


# plot posterior of precicion and teh true precision
dfSd <- as.data.frame(ARbammitModel$BUGSoutput$sims.list$sy)
dfSdT <- data.frame(true = 1 / (data$sy^2))
pSd <- ggplot(dfSd, aes(x = V1)) +
  stat_density(colour = "steelblue", fill = "steelblue", alpha = 0.2, geom = "area") +
  labs(x = expression(1 / sigma["y"]^"2"), y = "Density") +
  geom_vline(data = dfSdT, aes(xintercept = true), colour = "firebrick") +
  theme_bw(base_size = 14)
pSd

pBlin + pSd


# plots for the main effects

g <- data.frame(true = data$bv[[1]], est = ARbammitModel$BUGSoutput$median$b1)
e <- data.frame(true = data$bv[[2]], est = ARbammitModel$BUGSoutput$median$b2)
t <- data.frame(true = data$bv[[3]], est = ARbammitModel$BUGSoutput$median$b3)

quantmainEff <- list(
  g = data.frame(q1 = apply(ARbammitModel$BUGSoutput$sims.list$b1, 2, function(x) quantile(x, 0.05)),
                 q2 = apply(ARbammitModel$BUGSoutput$sims.list$b1, 2, function(x) quantile(x, 0.95))),
  e = data.frame(q1 = apply(ARbammitModel$BUGSoutput$sims.list$b2, 2, function(x) quantile(x, 0.05)),
                 q2 = apply(ARbammitModel$BUGSoutput$sims.list$b2, 2, function(x) quantile(x, 0.95))),
  t = data.frame(q1 = apply(ARbammitModel$BUGSoutput$sims.list$b3, 2, function(x) quantile(x, 0.05)),
                 q2 = apply(ARbammitModel$BUGSoutput$sims.list$b3, 2, function(x) quantile(x, 0.95)))
) |> plyr::ldply()

mainEff <- reshape2::melt(list(g = g,e = e,t = t), id=c("true","est"))
mainEff <- cbind(mainEff, quantmainEff)

mainEff$facet <- factor(mainEff$L1,
                        labels = c(expression(bold(b)^{(1)}),
                                   expression(bold(b)^{(2)}),
                                   expression(bold(b)^{(3)})))

mainEff |> ggplot(aes(x = true, y = est)) +
  geom_abline(size = 0.3, col = "gray") +
  geom_linerange(aes(ymin =  q1, ymax = q2), alpha = 0.5, size = 0.4) +
  geom_point(colour =  "steelblue", size = 0.4) +
  facet_wrap(~facet, ncol = 4, labeller = label_parsed, scales = 'free') +
  labs(x = "True", y = "Estimated") +
  theme_bw(base_size = 14) +
  theme(strip.background = element_blank(),
        panel.spacing.x = unit(5,"mm"),
        legend.position = c(0.75, 0.04))


# comparison with bammit


# ------ V = 2

set.seed(02)
data <- simArBAMMIT(
  V = 2,
  Q = 1,
  Bv = c(12, 10),
  mu = 10,
  lambda = c(12),
  sb = 1,
  sB = 1,
  seta = 1,
  somega = 1,
  sy = 1
)

ARbammitModel <- ARbammitJags(
  data = data, Q = 1, mmu = 10, smu = 10, nthin = 2, nburnin = 100, stheta = 1, somega = 1, a = 0.001, b = 0.001
)

plot(ARbammitModel)

bammitModel <- bammitJags(data = data,
                          Q = 1,
                          mmu = 10,
                          smu = 10,
                          a = 0.001,
                          b = 0.001,
                          nthin = 2,
                          nburnin = 100)
plot(bammitARdata)

qplot(data$y, ARbammitModel$BUGSoutput$median$mu) +
  geom_abline()+
  theme_light() +
  labs(x = expression(hat(y)), y = 'y')

qplot(data$y, bammitModel$BUGSoutput$median$mu) +
  geom_abline()+
  theme_light() +
  labs(x = expression(hat(y)), y = 'y')

round(caret::RMSE(data$y, ARbammitModel$BUGSoutput$median$mu), 4)
round(caret::RMSE(data$y, bammitModel$BUGSoutput$median$mu), 4)




# ------ V = 3

set.seed(02)
data <- simArBAMMIT(
  V = 3,
  Q = 1,
  Bv = c(12, 10, 4),
  mu = 10,
  lambda = c(12),
  sb = 1,
  sB = 1,
  seta = 1,
  somega = 1,
  sy = 1
)

ARbammitModel <- ARbammitJags(
  data = data, Q = 1, mmu = 10, smu = 10, nthin = 2, nburnin = 100, stheta = 1, somega = 1, a = 0.001, b = 0.001
)

plot(ARbammitModel)

bammitModel <- bammitJags(data = data,
                          Q = 1,
                          mmu = 10,
                          smu = 10,
                          a = 0.001,
                          b = 0.001,
                          nthin = 2,
                          nburnin = 100)
plot(bammitARdata)

qplot(data$y, ARbammitModel$BUGSoutput$median$mu) +
  geom_abline()+
  theme_light() +
  labs(x = expression(hat(y)), y = 'y')

qplot(data$y, bammitModel$BUGSoutput$median$mu) +
  geom_abline()+
  theme_light() +
  labs(x = expression(hat(y)), y = 'y')

round(caret::RMSE(data$y, ARbammitModel$BUGSoutput$median$mu), 4)
round(caret::RMSE(data$y, bammitModel$BUGSoutput$median$mu), 4)


# ------ V = 4
set.seed(02)
data <- simArBAMMIT(
  V = 4,
  Q = 1,
  Bv = c(12, 10, 4, 2),
  mu = 10,
  lambda = c(12),
  sb = 1,
  sB = 1,
  seta = 1,
  somega = 1,
  sy = 1
)

ARbammitModel <- ARbammitJags(
  data = data, Q = 1, mmu = 10, smu = 10, nthin = 2, nburnin = 100, stheta = 1, somega = 1, a = 0.001, b = 0.001
)



bammitModel <- bammitJags(data = data,
                            Q = 1,
                            mmu = 10,
                            smu = 10,
                            a = 0.001,
                            b = 0.001,
                            nthin = 2,
                            nburnin = 100)
plot(bammitARdata)

qplot(data$y, ARbammitModel$BUGSoutput$median$mu) +
  geom_abline()+
  theme_light() +
  labs(y = expression(hat(y)), x = 'y')+

qplot(data$y, bammitModel$BUGSoutput$median$mu) +
  geom_abline()+
  theme_light() +
  labs(y = expression(hat(y)), x = 'y')

round(caret::RMSE(data$y, ARbammitModel$BUGSoutput$median$mu), 4)
round(caret::RMSE(data$y, bammitModel$BUGSoutput$median$mu), 4)
