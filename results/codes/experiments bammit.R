# Simulation of bammit
library(tidyverse)
library(reshape2)
library(R2jags)
library(bammit)


set.seed(2022)

data <- simulateDataBammit(I = 12, J = 6, K = 2, mu = 100, sg = 1, se = 1, st = 1, sy = 1, stheta = 1,
                           lambda = c(12))

hist(data$y)

bammitModel <- bammitJags(
  data = data, mmu = 100, smu = 10,  nthin = 1, nburnin = 500, a = 0.1, b = 0.1,  stheta = 1
)

curve(dgamma(x, 0.1, 0.5), 0, 10)

# hist(bammitModel$BUGSoutput$sims.list$sy)
# abline(v = 1/(10^2), col = "red")

# if the sum  is zero
# nn <-  nrow(bammitModel$BUGSoutput$sims.list$gamma[,,1])
# summm <- vector()
# for (i in 1:nn) {
#   summm[i] <- (matrix(rep(1,6), nrow = 1))%*%matrix(bammitModel$BUGSoutput$sims.list$gamma[,,1][i,], ncol = 1)
# }
# head(summm)
#
#
# # if the sum of square is zero
# t(bammitModel$BUGSoutput$sims.list$gamma[,,1][1,])%*%bammitModel$BUGSoutput$sims.list$gamma[,,1][1,]
# summ <- vector()
# for (i in 1:nn) {
#   summ[i] <- t(bammitModel$BUGSoutput$sims.list$gamma[,,1][i,])%*%bammitModel$BUGSoutput$sims.list$gamma[,,1][i,]
# }
# head(summ)


pBlin <- qplot(data$blin,bammitModel$BUGSoutput$mean$blin, ylab = "Estimated", xlab = "True") +
  geom_abline() + theme_bw(base_size = 14) + geom_point(colour = "steelblue", size = 2)
pBlin

# plot posterior of precicion and teh true precision
dfSd <- as.data.frame(bammitModel$BUGSoutput$sims.list$sy)
dfSdT <- data.frame(true = 1/(data$sy^2))
pSd <- ggplot(dfSd, aes(x = V1)) +
  geom_density(colour = "steelblue", fill = "steelblue", alpha = 0.7) +
  labs(x = expression(1/sigma["y"]^"2"), y = "Frequency") +
  geom_vline(data = dfSdT, aes(xintercept = true), colour = "firebrick")+
  theme_bw(base_size = 14)
pSd

gridExtra::grid.arrange(pBlin, pSd, nrow = 1)


dfG <- as.data.frame(bammitModel$BUGSoutput$sims.list$g)
names(dfG) <- paste0("g", seq(1:ncol(dfG)))
dfG <- melt(dfG)
dfGT <- data.frame(variable = levels(dfG$variable), true = data$g)

pG <- ggplot(dfG, aes(x = value)) +
  geom_density(colour = "steelblue", fill = "steelblue", alpha = 0.4) +
  facet_wrap(~variable) +
  labs(x = "Genotype", y = "Frequency") +
  geom_vline(data = dfGT, aes(xintercept = true), colour = "firebrick")+
  theme_bw(base_size = 14) +
  theme(strip.background = element_rect(fill="white"),
        panel.spacing.x = unit(5,"mm"),
        legend.position = c(0.75, 0.04))
pG


qplot(data$g,bammitModel$BUGSoutput$mean$g, ylab = "Estimated", xlab = "True") +
  geom_abline() + theme_bw(base_size = 14) + geom_point(colour = "steelblue", size = 2)



dfE <- as.data.frame(bammitModel$BUGSoutput$sims.list$e)
names(dfE) <- paste0("e", seq(1:ncol(dfE)))
dfE <- melt(dfE)
dfET <- data.frame(variable = levels(dfE$variable), true = data$e)

pE <- ggplot(dfE, aes(x = value)) +
  geom_density(colour = "steelblue", fill = "steelblue", alpha = 0.4) +
  facet_wrap(~variable) +
  labs(x = "Environment", y = "Frequency") +
  geom_vline(data = dfET, aes(xintercept = true), colour = "firebrick")+
  theme_bw(base_size = 14)+
  theme(strip.background = element_rect(fill="white"),
        panel.spacing.x = unit(5,"mm"),
        legend.position = c(0.75, 0.04))
pE


dfT <- as.data.frame(bammitModel$BUGSoutput$sims.list$t)
names(dfT) <- paste0("t",  seq(1:ncol(dfT)))
dfT <- melt(dfT)
dfTT <- data.frame(variable = levels(dfT$variable), true = data$t)

pT <- ggplot(dfT, aes(x = value)) +
  geom_density(colour = "steelblue", fill = "steelblue", alpha = 0.4) +
  facet_wrap(~variable) +
  labs(x = "Year posterior", y = "Frequency") +
  geom_vline(data = dfTT, aes(xintercept = true), colour = "firebrick")+
  theme_bw(base_size = 14)+
  theme(strip.background = element_rect(fill="white"),
        panel.spacing.x = unit(5,"mm"),
        legend.position = c(0.75, 0.04))
pT


pG
pE
pT



## one data set

dfMainEff <- melt(list(dfG, dfE, dfT))
dfMainEff$L1 <- as.factor(dfMainEff$L1)
levels(dfMainEff$L1) <- c("Genotype", "Environment", "Time")


dfMainEff <- dfMainEff |> group_by(variable) |>
  mutate(mean = mean(value))

dfMainEffTrue <- melt(list(dfGT, dfET, dfTT))
dfMainEffTrue$L1 <- as.factor(dfMainEffTrue$L1)
levels(dfMainEffTrue$L1) <- c("Genotype", "Environment", "Time")
dfMainEffTrue$est <- unique(dfMainEff$mean)

dfMainEff$true = dfMainEffTrue$value[ match(dfMainEff$variable, dfMainEffTrue$variable) ]

ggplot(dfMainEffTrue, aes(x = est, y = value, group = variable)) +
  geom_point(colour = "steelblue", size = 2.5) +
  geom_abline() +
  facet_wrap(~L1) +
  theme_bw(base_size = 16) +
  theme(strip.background = element_rect(fill = "white"),
        panel.spacing.x = unit(5,"mm")) +
  labs(x = "Estimated", y = "True")






qplot(data$g,bammitModel$BUGSoutput$mean$g, ylab = "Estimated", xlab = "True") +
  geom_abline() + theme_bw(base_size = 14) + geom_point(colour = "steelblue", size = 2)



# compare with multiple Q

data <- simulateDataBammit(I = 6, J = 4, K = 4, mu = 100, sg = 10, se = 10, st = 10, sy = 1, stheta = 1,
                           lambda = c(10,12))

bammitModel <- bammitJags(
  data = data, mmu = 90, smu = 10,  nthin = 1, nburnin = 100, a = 4, b = 6,  stheta = 1
)








