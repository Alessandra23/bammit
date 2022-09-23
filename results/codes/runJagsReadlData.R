load("Real data/ireland.RData")
library(dplyr)
library(R2jags)
library(bayesplot)
library(tidyr)
library(ggplot2)
library(ggridges)
library(MCMCvis)


bammitModelRealData <- readRDS("~/Documents/GitHub/bammit/Real data/bammitModelRealData.rds")
MCMCtrace(bammitModelRealData,
          params = c('g[1]', 'g[2]', 'g[3]', 'g[4]', 'g[5]', 'g[6]'),
          ISB = FALSE,
          exact = TRUE,
          pdf = FALSE)

MCMCtrace(bammitModelRealData,
          params = c('blin[1]', 'blin[2]', 'blin[3]', 'blin[4]', 'blin[5]', 'blin[6]'),
          ISB = FALSE,
          exact = TRUE,
          pdf = FALSE)

MCMCtrace(bammitModelRealData)

MCMCtrace(bammitModelRealData,
          params = c('gamma[1]', 'blin[2]', 'blin[3]', 'blin[4]', 'blin[5]', 'blin[6]'),
          ISB = FALSE,
          exact = TRUE,
          pdf = FALSE)



test <- ireland |>
  group_by(Year, Genotype, Environment) |>
  summarise(mean = mean(Yield))
test |> filter(Genotype == "g3")  |>
  ggplot(aes(x = Year, y = mean)) +
  geom_point() +
  geom_line(aes(x = Year, y = mean))
ireland %>% filter(Environment == "e1")
ireland %>% filter(Genotype == "g3")

data <- ireland #|> filter(!(Year == "2019"))
valid <- ireland #|> filter((Year == "2019"))

levels(data$Genotype) <-(1:length(levels(data$Genotype)))
levels(data$Environment) <- 1:length(levels(ireland$Environment))

data <- data[1:217,]


bammitModelRealData <- bammitJagsRealData(
  data = data, mmu = 90, smu = 10, stheta = 1/100, a = 0.1, b = 0.1,
  nthin = 1, nburnin = 10, Q = 1
)


# Check convergence

bammitModelRealData$BUGSoutput$summary
traceplot(bammitModelRealData, parameters = c("g[1]"))
MCMCtrace(bammitModelRealData,
          params = c('g[1]', 'g[2]', 'g[3]', 'g[4]', 'g[5]', 'g[6]'),
          ISB = FALSE,
          exact = TRUE,
          pdf = FALSE)


data <- ireland |> dplyr::filter((Year == "2010"))
#save(bammitModelRealData, file = "Real data/bammitRealResults.RData")
yhatRealData <- predictionBAMMITRealData(bammitModelRealData, data)
RMSE(data$Mean, yhatRealData)
RMSE(data$Mean, yhat)

qplot(data$Mean, yhat, ylab = expression(hat(y)), xlab = "y") +
  geom_abline() + theme_bw(base_size = 14) + geom_point(colour = "steelblue", size = 1)

bammitModelRealData$BUGSoutput$sims.list$g
bammitModelRealData$BUGSoutput$sims.array

color_scheme_set("blue")
listPost <- bammitModelRealData$BUGSoutput


fitcmcm <- as.mcmc(bammitModelRealData)
fitmat <- as.matrix(fitcmcm)
fitdf <- as.data.frame(fitmat)


fitg <- fitdf[, grep(x = colnames(fitdf), pattern = "g[", fixed = TRUE)]
fite <- fitdf[, grep(x = colnames(fitdf), pattern = "e[", fixed = TRUE)]
fitt <- fitdf[, grep(x = colnames(fitdf), pattern = "t[", fixed = TRUE)]
fitblin <- fitdf[, grep(x = colnames(fitdf), pattern = "blin[", fixed = TRUE)]
fitgamma <- fitdf[, grep(x = colnames(fitdf), pattern = "gamma[", fixed = TRUE)]
#get all
#gPostCoeff <- select(fitdf,paste0("g[", seq(1, ncol(fitg), by = 1), "]"))
gPostCoeff <- select(fitdf,paste0("g[", seq(1, 10, by = 1), "]"))
names(gPostCoeff) <- paste0("g", 1:ncol(gPostCoeff))
gPostCoeffLong <- gather(gPostCoeff)


g1 <- ggplot(data = gPostCoeffLong, aes(x = value, y = key)) +
  stat_density_ridges(quantile_lines = TRUE,
                      quantiles = c(0.025, 0.5, 0.975),
                      alpha = 0.7,
                      colour = "steelblue",
                      fill = "steelblue") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  theme_bw() + xlab("Posterior estimate") + ylab("")
g1

gPostCoeffLongSum <- summarize(group_by(gPostCoeffLong, key),
                                        median_coef = median(value),
                                        lower_coef = quantile(value, probs = c(0.025)),
                                        upper_coef = quantile(value, probs = c(0.975)))

g2 <- ggplot(data = gPostCoeffLongSum, aes(x = median_coef, y = key)) +
  geom_errorbarh(aes(xmin = lower_coef, xmax = upper_coef), height = 0.1) +
  geom_point() +
  geom_vline(xintercept = 0, linetype = "dashed") +
  xlab("Posterior estimate") + ylab("") +
  theme_bw()
g2
plotsG <- gridExtra::grid.arrange(g1,g2, nrow = 1)

# e

ePostCoeff <- select(fitdf,paste0("e[", seq(1, 10, by = 1), "]"))
names(ePostCoeff) <- paste0("e", 1:ncol(ePostCoeff))
ePostCoeffLong <- gather(ePostCoeff)


e1 <- ggplot(data = ePostCoeffLong, aes(x = value, y = key)) +
  stat_density_ridges(quantile_lines = TRUE,
                      quantiles = c(0.025, 0.5, 0.975),
                      alpha = 0.7,
                      #alpha = 0.2,
                      colour = "steelblue",
                      fill = "steelblue") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  theme_bw() + xlab("Posterior estimate") + ylab("")
e1

ePostCoeffLongSum <- summarize(group_by(ePostCoeffLong, key),
                               median_coef = median(value),
                               lower_coef = quantile(value, probs = c(0.025)),
                               upper_coef = quantile(value, probs = c(0.975)))

e2 <- ggplot(data = ePostCoeffLongSum, aes(x = median_coef, y = key)) +
  geom_errorbarh(aes(xmin = lower_coef, xmax = upper_coef), height = 0.1) +
  geom_point() +
  geom_vline(xintercept = 0, linetype = "dashed") +
  xlab("Posterior estimate") + ylab("") +
  theme_bw()

plotsE <- gridExtra::grid.arrange(e1,e2, nrow = 1)

## t

tPostCoeff <- select(fitdf,paste0("t[", seq(1, 9, by = 1), "]"))
names(gPostCoeff) <- paste0("t", 1:ncol(tPostCoeff))
tPostCoeffLong <- gather(tPostCoeff)


t1 <- ggplot(data = tPostCoeffLong, aes(x = value, y = key)) +
  stat_density_ridges(quantile_lines = TRUE,
                      quantiles = c(0.025, 0.5, 0.975),
                      alpha = 0.7,
                      colour = "steelblue",
                      fill = "steelblue") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  theme_bw() + xlab("Posterior estimate") + ylab("")
t1

tPostCoeffLongSum <- summarize(group_by(tPostCoeffLong, key),
                               median_coef = median(value),
                               lower_coef = quantile(value, probs = c(0.025)),
                               upper_coef = quantile(value, probs = c(0.975)))

t2 <- ggplot(data = tPostCoeffLongSum, aes(x = median_coef, y = key)) +
  geom_errorbarh(aes(xmin = lower_coef, xmax = upper_coef), height = 0.1) +
  geom_point() +
  geom_vline(xintercept = 0, linetype = "dashed") +
  xlab("Posterior estimate") + ylab("") +
  theme_bw()

plotsT <- gridExtra::grid.arrange(t1,t2,nrow = 1)


# blin

blinPostCoeff <- select(fitdf,paste0("blin[", seq(1, 10, by = 1), "]"))
names(blinPostCoeff) <- paste0("blin", 1:ncol(blinPostCoeff))
blinPostCoeffLong <- gather(blinPostCoeff)


blin1 <- ggplot(data = blinPostCoeffLong, aes(x = value, y = key)) +
  stat_density_ridges(quantile_lines = TRUE,
                      quantiles = c(0.025, 0.5, 0.975),
                      alpha = 0.7,
                      colour = "steelblue",
                      fill = "steelblue") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  theme_bw() + xlab("Posterior estimate") + ylab("")


blinPostCoeffLongSum <- summarize(group_by(blinPostCoeffLong, key),
                               median_coef = median(value),
                               lower_coef = quantile(value, probs = c(0.025)),
                               upper_coef = quantile(value, probs = c(0.975)))

blin2 <- ggplot(data = blinPostCoeffLongSum, aes(x = median_coef, y = key)) +
  geom_errorbarh(aes(xmin = lower_coef, xmax = upper_coef), height = 0.1) +
  geom_point() +
  geom_vline(xintercept = 0, linetype = "dashed") +
  xlab("Posterior estimate") + ylab("") +
  theme_bw()

plotsBlin <- gridExtra::grid.arrange(blin1, blin2,nrow = 1)



## Gammas

#get all
#gPostCoeff <- select(fitdf,paste0("g[", seq(1, ncol(fitg), by = 1), "]"))
gPostCoeff <- select(fitdf,paste0("gamma[", seq(1, 10, by = 1),",", 1, "]"))
names(gPostCoeff) <- paste0("gamma", 1:ncol(gPostCoeff))
gPostCoeffLong <- gather(gPostCoeff)


g1 <- ggplot(data = gPostCoeffLong, aes(x = value, y = key)) +
  stat_density_ridges(quantile_lines = TRUE,
                      quantiles = c(0.025, 0.5, 0.975),
                      alpha = 0.7,
                      colour = "steelblue",
                      fill = "steelblue") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  theme_bw() + xlab("Posterior estimate") + ylab("")
g1

gPostCoeffLongSum <- summarize(group_by(gPostCoeffLong, key),
                               median_coef = median(value),
                               lower_coef = quantile(value, probs = c(0.025)),
                               upper_coef = quantile(value, probs = c(0.975)))

g2 <- ggplot(data = gPostCoeffLongSum, aes(x = median_coef, y = key)) +
  geom_errorbarh(aes(xmin = lower_coef, xmax = upper_coef), height = 0.1) +
  geom_point() +
  geom_vline(xintercept = 0, linetype = "dashed") +
  xlab("Posterior estimate") + ylab("") +
  theme_bw()
g2
plotsG <- gridExtra::grid.arrange(g1,g2, nrow = 1)


## deltas

#get all
#gPostCoeff <- select(fitdf,paste0("g[", seq(1, ncol(fitg), by = 1), "]"))
gPostCoeff <- select(fitdf,paste0("delta[", seq(1, 10, by = 1),",", 1, "]"))
names(gPostCoeff) <- paste0("delta", 1:ncol(gPostCoeff))
gPostCoeffLong <- gather(gPostCoeff)


g1 <- ggplot(data = gPostCoeffLong, aes(x = value, y = key)) +
  stat_density_ridges(quantile_lines = TRUE,
                      quantiles = c(0.025, 0.5, 0.975),
                      alpha = 0.7,
                      colour = "steelblue",
                      fill = "steelblue") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  theme_bw() + xlab("Posterior estimate") + ylab("")
g1

gPostCoeffLongSum <- summarize(group_by(gPostCoeffLong, key),
                               median_coef = median(value),
                               lower_coef = quantile(value, probs = c(0.025)),
                               upper_coef = quantile(value, probs = c(0.975)))

g2 <- ggplot(data = gPostCoeffLongSum, aes(x = median_coef, y = key)) +
  geom_errorbarh(aes(xmin = lower_coef, xmax = upper_coef), height = 0.1) +
  geom_point() +
  geom_vline(xintercept = 0, linetype = "dashed") +
  xlab("Posterior estimate") + ylab("") +
  theme_bw()
g2
plotsG <- gridExtra::grid.arrange(g1,g2, nrow = 1)

## rhos

#get all
#gPostCoeff <- select(fitdf,paste0("g[", seq(1, ncol(fitg), by = 1), "]"))
gPostCoeff <- select(fitdf,paste0("rho[", seq(1, 10, by = 1),",", 1, "]"))
names(gPostCoeff) <- paste0("rho", 1:ncol(gPostCoeff))
gPostCoeffLong <- gather(gPostCoeff)


g1 <- ggplot(data = gPostCoeffLong, aes(x = value, y = key)) +
  stat_density_ridges(quantile_lines = TRUE,
                      quantiles = c(0.025, 0.5, 0.975),
                      alpha = 0.7,
                      colour = "steelblue",
                      fill = "steelblue") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  theme_bw() + xlab("Posterior estimate") + ylab("")
g1

gPostCoeffLongSum <- summarize(group_by(gPostCoeffLong, key),
                               median_coef = median(value),
                               lower_coef = quantile(value, probs = c(0.025)),
                               upper_coef = quantile(value, probs = c(0.975)))

g2 <- ggplot(data = gPostCoeffLongSum, aes(x = median_coef, y = key)) +
  geom_errorbarh(aes(xmin = lower_coef, xmax = upper_coef), height = 0.1) +
  geom_point() +
  geom_vline(xintercept = 0, linetype = "dashed") +
  xlab("Posterior estimate") + ylab("") +
  theme_bw()
g2
plotsG <- gridExtra::grid.arrange(g1,g2, nrow = 1)


## lambdas

#get all
#gPostCoeff <- select(fitdf,paste0("g[", seq(1, ncol(fitg), by = 1), "]"))
gPostCoeff <- matrix(fitdf$lambda, ncol = 1)
lambda <- fitdf$lambda
names(gPostCoeff) <- rep("lambda", nrow(gPostCoeff))
gPostCoeffLong <- data.frame(key = seq_along(1:nrow(gPostCoeff)), value = gPostCoeff)


ggplot(gPostCoeffLong, aes(x = key, y = value)) +
  geom_density_ridges()


names(gPostCoeffLong$key)

g1 <- ggplot(data = gPostCoeffLong, aes(x = value, y = key)) +
  stat_density_ridges(quantile_lines = TRUE,
                      quantiles = c(0.025, 0.5, 0.975),
                      alpha = 0.7,
                      colour = "steelblue",
                      fill = "steelblue") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  theme_bw() + xlab("Posterior estimate") + ylab("")
g1

gPostCoeffLongSum <- gPostCoeffLong |> summarize(median_coef = median(value),
                               lower_coef = quantile(value, probs = c(0.025)),
                               upper_coef = quantile(value, probs = c(0.975)))

g2 <- ggplot(data = gPostCoeffLongSum, aes(x = median_coef, y = key)) +
  geom_errorbarh(aes(xmin = lower_coef, xmax = upper_coef), height = 0.1) +
  geom_point() +
  geom_vline(xintercept = 0, linetype = "dashed") +
  xlab("Posterior estimate") + ylab("") +
  theme_bw()
g2
plotsG <- gridExtra::grid.arrange(g1,g2, nrow = 1)

mcmc_areas(
  bammitModelRealData,
  pars = c("lambda"),
  prob = 0.5, # 80% intervals
  prob_outer = 0.99, # 99%
  point_est = "mean"
)


gPostCoeffLong |> ggplot(aes(x = value)) +
  geom_density(alpha = 0.7,
               colour = "steelblue",
               fill = "steelblue")+
  theme_bw() +
  # geom_vline(xintercept = 0, linetype = "dashed") +
  # geom_vline(xintercept = gPostCoeffLongSum$median_coef, colour = "steelblue") +
  # geom_vline(xintercept = gPostCoeffLongSum$lower_coef, colour = "steelblue") +
  # geom_vline(xintercept = gPostCoeffLongSum$upper_coef, colour = "steelblue") +
  theme_bw() + xlab("Posterior estimate") + ylab("")






