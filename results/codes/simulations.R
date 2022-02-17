# Simulate data from ammi equation
library(R2jags)
library(tidyverse)

data <- simulateDataAmmi(I = 25, J = 6, mu = 100, sg = 10, se = 10, sy = 2, lambda = 10)
data <- simulateDataAmmi(I = 6, J = 4, mu = 100, sg = 10, se = 10, sy = 2, lambda = c(10,12))
data <- simulateDataAmmi(I = 12, J = 6, mu = 100, sg = 10, se = 10, sy = 0.5, lambda = c(10,12,25))

josseModel <- josseJags(data = data,  mmu = 90, smu = 10,
                        sg = 10, se = 10, slambda = 1,
                        a = 0.1, b = 0.1, nthin = 1, nburnin = 1000)
crossaModel <- crossaJags(data = data, mmu = 90, smu = 10,
                          mug = 0, sg = 10, mue = 0, se = 10,
                          mulambda = 10, slambda = 1, a = 0.1, b = 0.1, stheta = 1,
                          nthin = 1, nburnin = 1000)

crossaModel$BUGSoutput$sims.list$gamma[,,2]
qqplot(data$blin, crossaModel$BUGSoutput$sims.list$blin)
abline(0,1)


qplot(data$blin, crossaModel$BUGSoutput$mean$blin, ylab = "Estimated bilinear - Crossa", xlab = "True bilinear") + geom_abline() +
  geom_abline() + theme_bw() + geom_point(colour = "red", size = 0.6)

qplot(data$blin, josseModel$BUGSoutput$mean$blin, ylab = "Estimated bilinear - Josse", xlab = "True bilinear") + geom_abline() +
  geom_abline() + theme_bw() + geom_point(colour = "red", size = 0.6)


# Simulations BAMMIT

dataT <- simulateDataBammit(I = 6, J = 4, K = 2, mu = 100, sg = 10, se = 10, st = 10, sy = 2, lambda = 10)
bammitModel <- bammitJags(data = dataT, mmu = 90, smu = 10, mug = 0, sg = 10, mue = 0, se = 10, mut = 0, st = 10,
                          mulambda = 10, slambda = 1, stheta = 10, a = 0.1, b = 0.1, nthin = 1, nburnin = 100)

qplot(dataT$blin, bammitModel$BUGSoutput$mean$blin, ylab = "Estimated bilinear - BAMMIT", xlab = "True bilinear") + geom_abline() +
  geom_abline() + theme_bw() + geom_point(colour = "red", size = 0.6)

yhatB <- predictionBAMMIT(bammitModel, dataT)

bammbammitModel$BUGSoutput$mean$blin
dataT$blin
data <- dataT

qplot(yhatB, dataT$y, xlab = expression(hat(y)), ylab = "y") + geom_abline() +
  geom_abline() + theme_bw() + geom_point(colour = "red", size = 0.6)
