load("Real data/ireland.RData")

test <- ireland %>%
  group_by(Year, Genotype) %>%
  summarise(mean = mean(Mean))
test %>%  filter(Genotype == "g3") %>%
  ggplot(aes(x = Year, y = mean)) +
  geom_point() +
  geom_line(aes(x = Year, y = mean))
ireland %>% filter(Environment == "e1")
ireland %>% filter(Genotype == "g3")

data <- ireland

levels(data$Genotype) <-(1:length(levels(data$Genotype)))
levels(data$Environment) <- 1:length(levels(ireland$Environment))


bammitModelRealData <- bammitJagsRealData(
  data = ireland, mmu = 90, smu = 10, stheta = 10, a = 0.1, b = 0.1,
  nthin = 1, nburnin = 100, Q = 2
)
#save(bammitModelRealData, file = "Real data/bammitRealResults.RData")
yhatRealData <- predictionBAMMITRealData(bammitModelRealData, data)
RMSE(data$Mean, yhatRealData)

bammitModelRealData$BUGSoutput$sims.list$g
bammitModelRealData$BUGSoutput$sims.array

color_scheme_set("blue")
listPost <- bammitModelRealData$BUGSoutput

spread_draws(bammitModelRealData$BUGSoutput$sims.array[blin])
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


library("ggridges")
g1 <- ggplot(data = gPostCoeffLong, aes(x = value, y = key)) +
  stat_density_ridges(quantile_lines = TRUE,
                      quantiles = c(0.025, 0.5, 0.975),
                      alpha = 0.7) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  theme_bw() + xlab("Posterior estimate") + ylab("")


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
plotsG <- gridExtra::grid.arrange(g1,g2, nrow = 1)

# e

ePostCoeff <- select(fitdf,paste0("e[", seq(1, 10, by = 1), "]"))
names(ePostCoeff) <- paste0("e", 1:ncol(ePostCoeff))
ePostCoeffLong <- gather(ePostCoeff)


library("ggridges")
e1 <- ggplot(data = ePostCoeffLong, aes(x = value, y = key)) +
  stat_density_ridges(quantile_lines = TRUE,
                      quantiles = c(0.025, 0.5, 0.975),
                      alpha = 0.7) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  theme_bw() + xlab("Posterior estimate") + ylab("")


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

tPostCoeff <- select(fitdf,paste0("t[", seq(1, 10, by = 1), "]"))
names(gPostCoeff) <- paste0("t", 1:ncol(tPostCoeff))
tPostCoeffLong <- gather(tPostCoeff)


library("ggridges")
t1 <- ggplot(data = tPostCoeffLong, aes(x = value, y = key)) +
  stat_density_ridges(quantile_lines = TRUE,
                      quantiles = c(0.025, 0.5, 0.975),
                      alpha = 0.7) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  theme_bw() + xlab("Posterior estimate") + ylab("")


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


library("ggridges")
blin1 <- ggplot(data = blinPostCoeffLong, aes(x = value, y = key)) +
  stat_density_ridges(quantile_lines = TRUE,
                      quantiles = c(0.025, 0.5, 0.975),
                      alpha = 0.7) +
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





