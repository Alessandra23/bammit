# simulate data
library(bammit)
library(tidyverse)
library(reshape2)
library(R2jags)


set.seed(2022)
data <- simulateDataAmmi(I = 12, J = 10, mu = 100, sg = 10, se = 10, sy = 2,
                         lambda = c(12), stheta = 1)

crossaModel <- crossaJags(
  data = data, mmu = 90, smu = 10,
  mug = 0, mue = 0, a = 4, b = 10, stheta = 1,
  nthin = 1, nburnin = 100
)

qplot(data$blin, crossaModel$BUGSoutput$mean$blin, ylab = "Estimated", xlab = "True") + geom_abline() +
  geom_abline() + theme_bw() + geom_point(colour = "red", size = 0.6)


# compare gammas and deltas
# data$gamma %>% var()
# data$delta %>% var()


getGammaDelta <- function(I, J, lambda = c(12),  sy = 0.5, mu = 100, sg = 10, se = 10,
                          stheta = 0.1){

  data <-  simulateDataAmmi(I = I, J = J, mu = mu, sg = sg, se = se, sy = sy,
                            lambda = lambda, stheta = stheta)
  getGamma <- data$gamma
  getDelta <- data$delta
  getKappa <- data$kappa

  return(list(gamma = getGamma,
              delta = getDelta,
              kappa = getKappa))

}


getGammaDeltaDelta <- function(I, J, K, lambda = c(12),  sy = 0.5, mu = 100, sg = 10, se = 10, st = 10,
                          stheta = 0.1){

  data <-  simulateDataBammit(I = I, J = J, K = K, mu = mu, sg = sg, se = se, st = st, sy = sy, stheta = stheta,
                              lambda = lambda)
  getGamma <- data$gamma
  getDelta <- data$delta
  getKappa <- data$kappa

  return(list(gamma = getGamma,
              delta = getDelta,
              kappa = getKappa))

}

getGammaDeltaDelta(I = 12, J = 10, K = 4)

m <- 2000
I = 10
J = 6
K = 4

getGamma <- replicate(m, getGammaDeltaDelta(I = I, J = J, K = K, stheta = 0.1)$gamma)
getGammadf <- as.data.frame(getGamma)
(varGamma <- apply(getGammadf, 1, var))
#varGamma <- data.frame(gamma = factor(1:I), variances = varGamma)


getGammamelt <- melt(t(getGammadf))
gammaLabs <- paste0('\u03B3', 1:I)
getGammamelt %>% ggplot(aes(x = as.factor(Var2), y = value)) +
  geom_boxplot(fill = "steelblue") + #scale_fill_brewer(palette="RdBu", labels = c("Crossa", "Josse", "Real")) +
  #geom_point(varGamma, aes( y = variances)) +
  labs(x = " ", y = "value") +
  theme_minimal(base_size = 14) +
  scale_x_discrete(labels= gammaLabs) #+
  #geom_jitter(data=varGamma, aes(x=gamma, y=variances),
              #position=position_jitter(.2), color="blue", size=1.5, pch=20)


# qaundo aumento o tamanho de m, a

## teste


generateBlinT <- function(index, Q, stheta = 1){
  theta <- matrix(NA, nrow = index, ncol = Q)
  variable <- matrix(NA, nrow = index, ncol = Q)
  sqrtTheta <-  vector()
  for(q in 1:Q){
    for (i in 1:(index)) {
      theta[i,q] <- rnorm(1, 0, stheta)
    }
    #theta[index, q] <- -sum(theta[1:(index-1), q])
    m <- apply(theta, 2, mean)
    thetaN <- apply(theta, 1, function(x){x-m}) %>% as.matrix()

    if(Q>1){
      thetaN <- t(thetaN)
    }

    sqrtTheta[q] <- sqrt(1/sum(thetaN[,q]^2))

    for (i in 1:index) {
      variable[i,q] <- (thetaN[i,q])*sqrtTheta[q]
    }
  }
  return(list(theta, thetaN,  variable))
}
tt <-  generateBlinT(6,1,1)
sum(tt[[2]])
sum(tt[[3]]^2)


for(q in 1:Q){
  for(i in 1:I){
    thetaG[i,q] ~ dnorm(0,stheta)
  }
  mG[q] = sum(thetaG[1:I,q])/I
  for(i in 1:I){
    thetaGNew[i,q] = thetaG[i,q] - mG[q]
  }
  sqrtThetaG[q] = sqrt(1/(sum(thetaGNew[1:I,q]^2 + 0.000001)))
  for(i in 1:I){
    gamma[i,q] = thetaGNew[i,q]*sqrtThetaG[q]
  }
}


# if the sum is zero
(matrix(rep(1,6), nrow = 1))%*%matrix(bammitModel$BUGSoutput$sims.list$gamma[,,1][1,], ncol = 1)

bammitModel$BUGSoutput$sims.matrix
sum(bammitModel$BUGSoutput$mean$gamma^2)
sum(bammitModel$BUGSoutput$mean$thetaGNew^2)

# if the sum of square is zero
t(bammitModel$BUGSoutput$sims.list$gamma[,,1][1,])%*%bammitModel$BUGSoutput$sims.list$gamma[,,1][1,]
nn <-  nrow(bammitModel$BUGSoutput$sims.list$gamma[,,1])
summ <- vector()
for (i in 1:nn) {
  summ[i] <- t(bammitModel$BUGSoutput$sims.list$gamma[,,1][i,])%*%bammitModel$BUGSoutput$sims.list$gamma[,,1][i,]
}
summ


sqrtThetaG[q] = sqrt(1/(sum(thetaGNew[1:I,q]^2 + 0.000001)))

generateBlin <- function(index, Q, stheta){
  theta <- matrix(NA, nrow = index, ncol = Q)
  variable <- matrix(NA, nrow = index, ncol = Q)
  for(q in 1:Q){
    for (i in 1:(index-1)) {
      theta[i,q] <- rnorm(1, 0, stheta)
    }
    theta[index, q] <- -sum(theta[1:(index-1), q])
    sqrtTheta <- sqrt(1/sum(theta[,q]^2))
    for (i in 1:index) {
      variable[i,q] <- (theta[i,q])*sqrtTheta
    }
  }
  return(list(theta, variable))
}

generateBlin(6,1,1)

generateBlinT(6, 1, 1)
I = 6
m = 500
stheta = 1
testVar <- function(m, I, stheta){
  getReplicate <- replicate(m, generateBlin(I, 1, stheta))
  getVariable <- getReplicate[2,]
  getTheta <- getReplicate[1,]

  getThetadf <- as.data.frame(getTheta)
  colnames(getThetadf) <- paste0("V", 1:m)
  (varTheta <- apply(getThetadf, 1, var))
  mean(varTheta)

  mean(unlist(lapply(getTheta, var)))
  mean(unlist(lapply(getVariable, var)))

  getVariabledf <- as.data.frame(getVariable)
  colnames(getVariabledf) <- paste0("V", 1:m)
  (varVariable <- apply(getVariabledf, 1, var))
  mean(varVariable)

  lbs <- c("\u03B3", "\u03B8")
  custLab <- function (x){
    lbs
  }

  vardf <- data.frame(variable = varVariable, theta = varTheta)
  vardfmelt <- melt(vardf)
  varPlot <- vardfmelt %>% ggplot(aes(y = value)) +
    geom_boxplot() + #scale_fill_brewer(palette="RdBu", labels = c("Crossa", "Josse", "Real")) +
    #geom_point(varGamma, aes( y = variances)) +
    labs(x = " ", y = "value") +
    theme_bw(base_size = 14) +
    facet_wrap(~variable,scales = 'free', labeller = as_labeller(custLab))+
    theme(legend.position="bottom",
          axis.ticks.x = element_blank(),
          axis.text.x = element_blank())+
    theme(strip.background =element_rect(fill="white"))
  varPlot
  # par(mfrow = c(1,2))
  # boxplot(varTheta[-10])
  # boxplot(varVariable[-10])

  # plot gamma
  getVariablemelt <- melt(t(getVariabledf))
  varLabs <- paste0('\u03B3', 1:I) # gamma
  #varLabs <- paste0('\u03B4', 1:J) # delta
  variablePlot <- getVariablemelt %>% ggplot(aes(x = as.factor(Var2), y = value)) +
    geom_boxplot() + #scale_fill_brewer(palette="RdBu", labels = c("Crossa", "Josse", "Real")) +
    #geom_point(varGamma, aes( y = variances)) +
    labs(x = " ", y = "value") +
    theme_minimal(base_size = 14) +
    scale_x_discrete(labels= varLabs)
  variablePlot
  # plot theta
  getThetamelt <- melt(t(getThetadf))
  varLabs <- paste0('\u03B8', 1:I) # theta
  thetaPlot <- getThetamelt %>% ggplot(aes(x = as.factor(Var2), y = value)) +
    geom_boxplot() + #scale_fill_brewer(palette="RdBu", labels = c("Crossa", "Josse", "Real")) +
    #geom_point(varGamma, aes( y = variances)) +
    labs(x = " ", y = "value") +
    theme_minimal(base_size = 14) +
    scale_x_discrete(labels= varLabs)

  return(list(varPlot = varPlot,
              variablePlot = variablePlot,
              thetaPlot = thetaPlot,
              vardf = vardf,
              getVariablemelt = getVariablemelt))
}

I = 6
m = 500
test1 <- testVar(m = m, I = I, stheta = 0.1)
test2 <- testVar(m = m, I = I, stheta = 10)
test3 <- testVar(m = m, I = I, stheta = 100)
varcomp <- data.frame(V1 = test1$vardf$variable[-I],
                      V2 = test2$vardf$variable[-I],
                      V3 = test3$vardf$variable[-I])
varcompmelt <- melt(varcomp)
lbs <- paste0('\u03c3', " = ", c(0.1, 10, 100))
custLab <- function (x){
  lbs
}
varcompmelt %>% ggplot(aes(y = value)) +
  geom_boxplot() + #scale_fill_brewer(palette="RdBu", labels = c("Crossa", "Josse", "Real")) +
  #geom_point(varGamma, aes( y = variances)) +
  labs(x = "Distribution of variances", y = "value") +
  theme_bw(base_size = 14) +
  facet_wrap(~variable, labeller = as_labeller(custLab))+
  theme(legend.position="bottom",
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank())+
  theme(strip.background =element_rect(fill="white"))


variableComp <- list(V1 = test1$getVariablemelt,
                     V2 = test3$getVariablemelt,
                     V3 = test3$getVariablemelt)

variableCompmelt <- melt(variableComp, c("Var2", "value"))

varLabs <- paste0('\u03B3', 1:I) # gamma

variableCompmelt %>% ggplot(aes(x = as.factor(Var2), y = value)) +
  geom_boxplot() + #scale_fill_brewer(palette="RdBu", labels = c("Crossa", "Josse", "Real")) +
  #geom_point(varGamma, aes( y = variances)) +
  labs(x = " ", y = "value") +
  theme_bw(base_size = 14) +
  facet_wrap(~L1, labeller = as_labeller(custLab))+
  theme(legend.position="bottom",
        axis.ticks.x = element_blank())+
        #axis.text.x = element_blank())+
  theme(strip.background =element_rect(fill="white"))+
  scale_x_discrete(labels= varLabs)




var(test$vardf$variable)
test$variablePlot
test$vardf
test$varPlot
test$thetaPlot



generateBlin <- function(index, Q, stheta){
  theta <- matrix(NA, nrow = index, ncol = Q)
  variable <- matrix(NA, nrow = index, ncol = Q)
  for(q in 1:Q){
    for (i in 1:(index-1)) {
      theta[i,q] <- rnorm(1, 0, stheta)
    }
    theta[index, q] <- -sum(theta[1:(index-1), q])
    sqrtTheta <- sqrt(1/sum(theta[,q]^2))
    for (i in 1:index) {
      variable[i,q] <- (theta[i,q])*sqrtTheta
    }
  }
  return(variable)
}

I = 6
m = 500
stheta = 1

getReplicate <- replicate(m, generateBlin(I, 1, stheta))
getVariable <- getReplicate

getVariabledf <- as.data.frame(getVariable)
colnames(getVariabledf) <- paste0("V", 1:m)
(varVariable <- apply(getVariabledf, 1, var))
mean(varVariable)

getVariablemelt <- melt(t(getVariabledf))
varLabs <- paste0('\u03B3', 1:I) # gamma
#varLabs <- paste0('\u03B4', 1:J) # delta
variablePlot <- getVariablemelt %>% ggplot(aes(x = as.factor(Var2), y = value)) +
  geom_boxplot() + #scale_fill_brewer(palette="RdBu", labels = c("Crossa", "Josse", "Real")) +
  #geom_point(varGamma, aes( y = variances)) +
  labs(x = " ", y = "value") +
  theme_minimal(base_size = 14) +
  scale_x_discrete(labels= varLabs)
variablePlot



generateBlinT <- function(index, Q, stheta){
  theta <- matrix(NA, nrow = index, ncol = Q)
  variable <- matrix(NA, nrow = index, ncol = Q)
  for(q in 1:Q){
    for (i in 1:(index)) {
      theta[i,q] <- rnorm(1, 0, stheta)
    }
    #theta[index, q] <- -sum(theta[1:(index-1), q])
    m <- apply(theta, 2, mean)
    thetaN <- apply(theta, 1, function(x){x-m}) %>% as.matrix()
    sqrtTheta <- sqrt(1/sum(thetaN[,q]^2))

    for (i in 1:index) {
      variable[i,q] <- (thetaN[i,q])*sqrtTheta
    }
  }
  return(list(thetaN, variable))
}

aa = generateBlinT(6,1,1)
sum(aa[[1]])
sum(aa[[2]])
sum(aa[[2]]^2)


I = 12
m = 500
stheta = 1

getReplicate <- replicate(m, generateBlinT(I, 1, stheta))
getVariable <- getReplicate[2,]

getReplicate[1,]

getVariabledf <- as.data.frame(getVariable)
colnames(getVariabledf) <- paste0("V", 1:m)
(varVariable <- apply(getVariabledf, 1, var))
mean(varVariable)

getVariablemelt <- melt(t(getVariabledf))
varLabs <- paste0('\u03B3', 1:I) # gamma
#varLabs <- paste0('\u03B4', 1:J) # delta
variablePlot <- getVariablemelt %>% ggplot(aes(x = as.factor(Var2), y = value)) +
  geom_boxplot() + #scale_fill_brewer(palette="RdBu", labels = c("Crossa", "Josse", "Real")) +
  #geom_point(varGamma, aes( y = variances)) +
  labs(x = " ", y = "value") +
  theme_minimal(base_size = 14) +
  scale_x_discrete(labels= varLabs)
variablePlot


getTheta <- getReplicate[1,]
getThetadf <- as.data.frame(getTheta)
colnames(getThetadf) <- paste0("V", 1:m)
(varTheta <- apply(getThetadf, 1, var))
mean(varTheta)

lbs <- c("\u03B3", "\u03B8")
custLab <- function (x){
  lbs
}

vardf <- data.frame(variable = varVariable, theta = varTheta)
vardfmelt <- melt(vardf)
varPlot <- vardfmelt %>% ggplot(aes(y = value)) +
  geom_boxplot() + #scale_fill_brewer(palette="RdBu", labels = c("Crossa", "Josse", "Real")) +
  #geom_point(varGamma, aes( y = variances)) +
  labs(x = " ", y = "value") +
  theme_bw(base_size = 14) +
  facet_wrap(~variable,scales = 'free', labeller = as_labeller(custLab))+
  theme(legend.position="bottom",
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank())+
  theme(strip.background =element_rect(fill="white"))
varPlot


abc <- c('x1', 'x2', 'x3')
values <- '#c3c3c3'
setNames(abc,rep(values, length(abc)))

names(abc) <- c("name1", "name2", "name3")
abc




