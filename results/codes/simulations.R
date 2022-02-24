# Simulate data from ammi equation
library(R2jags)
library(tidyverse)
library(reshape2)
library(RColorBrewer)
library(ggplot2)
library(gridExtra)


data <- simulateDataAmmi(I = 25, J = 6, mu = 100, sg = 10, se = 10, sy = 2, lambda = 10)
data <- simulateDataAmmi(I = 6, J = 4, mu = 100, sg = 10, se = 10, sy = 2, lambda = c(10, 12))
data <- simulateDataAmmi(I = 12, J = 6, mu = 100, sg = 10, se = 10, sy = 0.5, lambda = c(2, 12, 25))

set.seed(2022)
data <- simulateDataAmmi(I = 12, J = 6, mu = 100, sg = 10, se = 10, sy = 0.5, lambda = 10)
josseModel <- josseJags(
  data = data, mmu = 90, smu = 10,
  sg = 10, se = 10, slambda = 1,
  a = 0.1, b = 0.1, nthin = 1, nburnin = 1000
)
qplot(data$blin, josseModel$BUGSoutput$mean$blin, ylab = "Estimated bilinear - Josse", xlab = "True bilinear") + geom_abline() +
  geom_abline() + theme_bw() + geom_point(colour = "red", size = 0.6)

crossaModel <- crossaJags(
  data = data, mmu = 90, smu = 10,
  mug = 0, mue = 0, a = 0.1, b = 0.1, stheta = 10,
  nthin = 1, nburnin = 500
)


qqplot(data$blin, crossaModel$BUGSoutput$sims.list$blin)
abline(0, 1)


qplot(data$blin, crossaModel$BUGSoutput$mean$blin, ylab = "Estimated bilinear - Crossa", xlab = "True bilinear") + geom_abline() +
  geom_abline() + theme_bw() + geom_point(colour = "red", size = 0.6)




compSigmas <- sigmaComp(sigma = c(0.5, 2, 10), lambda = 10, slambda = 1, ncol = 3)
compSigmas$plot[[1]]
compSigmas$plot[[2]]
compSigmas$plot[[3]]
compSigmas$plot[[4]]
compSigmas$plot[[5]]
grid.arrange(
  compSigmas$plot[[1]],
  compSigmas$plot[[2]]
)


compSigmas2 <- sigmaComp(sigma = c(0.5, 2, 10), lambda = c(2,10, 25),slambda = 1,  ncol = 3)
compSigmas2$plot[[1]]
compSigmas2$plot[[2]]
compSigmas2$plot[[3]]
compSigmas2$plot[[4]]
compSigmas2$plot[[5]]

grid.arrange(
  compSigmas$plot[[1]],
  compSigmas2$plot[[1]]
)


grid.arrange(
  compSigmas$plot[[2]],
  compSigmas2$plot[[2]]
)

grid.arrange(
  compSigmas$plot[[3]],
  compSigmas2$plot[[3]]
)



## Visualize gamma and delta simulated

set.seed(2022)
data <- simulateDataAmmi(I = 25, J = 6, mu = 100, sg = 10, se = 10, sy = 0.5, lambda = c(2))




test <- vizSimAmmi(data)
test$boxplotGamma
test$boxplotDelta
test$boxplotGammaDelta
test$geBoxplot




# Simulations BAMMIT

data <- simulateDataBammit(I = 6, J = 4, K = 4, mu = 100, sg = 10, se = 10, st = 10, sy = 2,
                           lambda = c(2,12, 25))
data <- simulateDataBammit(I = 6, J = 4, K = 2, mu = 100, sg = 10, se = 10, st = 10, sy = 2,
                           lambda = 12)

vizBammit <- vizSimBammit(data)
vizBammit$boxplotGamma
vizBammit$boxplotDelta
vizBammit$boxplotKappa
vizBammit$boxplotGammaDeltaKappa
vizBammit$getBoxplot

bammitModel <- bammitJags(
  data = data, mmu = 90, smu = 10, stheta = 10, a = 0.1, b = 0.1, nthin = 1, nburnin = 100
)

qplot(dataT$blin, bammitModel$BUGSoutput$mean$blin, ylab = "Estimated bilinear - BAMMIT", xlab = "True bilinear") + geom_abline() +
  geom_abline() + theme_bw() + geom_point(colour = "red", size = 0.6)

yhatB <- predictionBAMMIT(bammitModel, dataT)

bammbammitModel$BUGSoutput$mean$blin
dataT$blin
data <- dataT

qplot(yhatB, dataT$y, xlab = expression(hat(y)), ylab = "y") + geom_abline() +
  geom_abline() + theme_bw() + geom_point(colour = "red", size = 0.6)
