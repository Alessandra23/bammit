library(tidyverse)
library(dplyr)
library(ggplot2)
library(R2jags)
library(ggridges)

## Get a small data set of Ireland
load("Real data/ireland.RData")
years <- c("2017","2018", "2019")
irelandSmall <- ireland |> dplyr::filter(Year %in% years)
irelandSmall$Year <- factor(irelandSmall$Year)
irelandSmall$Genotype <- factor(irelandSmall$Genotype)
irelandSmall$Environment <- factor(irelandSmall$Environment)


test <-  irelandSmall |>
  group_by(Genotype, Environment) |>
  dplyr::summarise(mean=mean(Mean))

test$mean  = test$mean - mean(test$mean)

bp1 <- ggplot(test, aes(x = Genotype, y = mean)) +
  geom_boxplot() +theme_bw()
bp1

bp2 <- ggplot(test, aes(x = Environment, y = mean)) +
  geom_boxplot() +theme_bw()
bp2

bp3 <- ggplot(irelandSmall, aes(x = Genotype, y = Mean, colour=Year)) +
  geom_boxplot() + theme_bw() +
  scale_colour_manual(values=c("steelblue", "firebrick"))
bp3
# Test original code ------------------------------------------------------

model <- bammitJagsRealData(
  data = irelandSmall, mmu = 90, smu = 10, stheta = 1/100, a = 0.1, b = 0.1,
  nthin = 1, nburnin = 10, Q = 1
)

#saveRDS(model, file = "Real data/model")


# Plots -------------------------------------------------------------------

fitcmcm <- as.mcmc(model)
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
gPostCoeffLong <- gather(fitg)

gennames <- c("g[1]", 'g[20]', 'g[24]')
gPostCoeffLong <- gPostCoeffLong |> filter(key %in% gennames)


g1 <- ggplot(data = gPostCoeffLong, aes(x = value, y = key)) +
  stat_density_ridges(quantile_lines = TRUE,
                      quantiles = c(0.025, 0.5, 0.975),
                      alpha = 0.4,
                      colour = "steelblue",
                      fill = "steelblue") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  theme_bw(base_size = 16) + xlab("Genotypes") + ylab("")
g1



## env

ePostCoeff <- select(fitdf,paste0("e[", seq(1, 10, by = 1), "]"))
names(ePostCoeff) <- paste0("e", 1:ncol(ePostCoeff))
ePostCoeffLong <- gather(ePostCoeff)


e1 <- ggplot(data = ePostCoeffLong, aes(x = value, y = key)) +
  stat_density_ridges(quantile_lines = TRUE,
                      quantiles = c(0.025, 0.5, 0.975),
                      alpha = 0.4,
                      #alpha = 0.2,
                      colour = "steelblue",
                      fill = "steelblue") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  theme_bw(base_size = 16) + xlab("Environments") + ylab("")
e1


e2 <- ggplot(data = ePostCoeffLong, aes(x = value, y = key)) +
  stat_density_ridges(quantile_lines = TRUE,
                      quantiles = c(0.025, 0.5, 0.975),
                      alpha = 0.4,
                      #alpha = 0.2,
                      colour = "steelblue",
                      fill = "steelblue") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  theme_bw() + xlab("Posterior estimate") + ylab("")
e2


# time
tPostCoeff <- select(fitdf,paste0("t[", seq(1, 10, by = 1), "]"))
names(tPostCoeff) <- paste0("t", 1:ncol(tPostCoeff))
tPostCoeffLong <- gather(tPostCoeff)


t1 <- ggplot(data = tPostCoeffLong, aes(x = value, y = key)) +
  stat_density_ridges(quantile_lines = TRUE,
                      quantiles = c(0.025, 0.5, 0.975),
                      alpha = 0.7,
                      #alpha = 0.2,
                      colour = "steelblue",
                      fill = "steelblue") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  theme_bw() + xlab(" ") + ylab("")
t1

# blin
tPostCoeff <- select(fitdf,paste0("g[", seq(10, 20, by = 1), "]"))
names(tPostCoeff) <- paste0("Interaction", 1:ncol(tPostCoeff))
tPostCoeffLong <- gather(tPostCoeff)


t1 <- ggplot(data = tPostCoeffLong, aes(x = value, y = key)) +
  stat_density_ridges(quantile_lines = TRUE,
                      quantiles = c(0.025, 0.5, 0.975),
                      alpha = 0.4,
                      #alpha = 0.2,
                      colour = "steelblue",
                      fill = "steelblue") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  theme_bw() + xlab("Interaction") + ylab("")
t1


gridExtra::grid.arrange(g1, e1, nrow = 1)



####


## get the predictions

data <- list(y = irelandSmall$Mean, x = irelandSmall[, c("Genotype", "Environment", "Year")])
names(data$x) <- c("g", "e", "t")
yhat <- predictionBAMMIT(model = model, data = data)
length(yhat)

dfG <- as.data.frame(model$BUGSoutput$sims.list$g)
names(dfG) <- unique(irelandSmall$Genotype)
dfE <- as.data.frame(model$BUGSoutput$sims.list$e)
names(dfE) <-  unique(irelandSmall$Environment)
dfT <- as.data.frame(model$BUGSoutput$sims.list$t)
names(dfT) <-  unique(irelandSmall$Year)

g <- data$x$g
e <- data$x$e
t <- data$x$t


# create a df
dat <- data.frame(yhat = yhat, e = e, g = g, t = t)
intPal = rev(diverging_hcl(palette = "Blue-Red 3", n = 100))

# factor levels
dat$e <- factor(dat$e, levels = unique(dat$e))
dat$g <- factor(dat$g, levels = unique(dat$g))
dat$t <- factor(dat$t, levels = unique(dat$t))

# plot limits
name <- "\u0177"
intLims <- range(dat$yhat)
limitsInt <- range(labeling::rpretty(intLims[1], intLims[2]))

# plot
ggplot(data = dat, aes(x = reorder(e, -yhat), y = reorder(g, yhat), fill = yhat))  +
  geom_tile() +
  #geom_text(aes(label = round(yhat,2))) +# just added this line to confirm the values are correct
  scale_fill_gradientn(
    limits = limitsInt,
    colors = intPal,
    guide = guide_colorbar(
      frame.colour = "black",
      ticks.colour = "black"
    ),
    name = name
  ) +
  xlab("Environment") +
  ylab("Genotype") +
  facet_wrap(~t, scales = "free") +
  theme_classic(base_size = 12)
  #theme(strip.background = element_blank(), panel.spacing.x = unit(5,"mm"))


# no sorted

ggplot(data = dat, aes(x = e, y = g, fill = yhat))  +
  geom_tile() +
  #geom_text(aes(label = round(yhat,2))) +# just added this line to confirm the values are correct
  scale_fill_gradientn(
    limits = limitsInt,
    colors = intPal,
    guide = guide_colorbar(
      frame.colour = "black",
      ticks.colour = "black"
    ),
    name = name
  ) +
  xlab("Environment") +
  ylab("Genotype") +
  facet_wrap(~t, scales = "free") +
  theme_classic(base_size = 16)


dat |> filter(t == "2010") |> ggplot(aes(x = reorder(e, -yhat), y = reorder(g, yhat), fill = yhat))  +
  geom_tile() +
  #geom_text(aes(label = round(yhat,2))) +# just added this line to confirm the values are correct
  scale_fill_gradientn(
    limits = limitsInt,
    colors = intPal,
    guide = guide_colorbar(
      frame.colour = "black",
      ticks.colour = "black"
    ),
    name = name
  ) +
  xlab("Environment") +
  ylab("Genotype") +
  theme_classic(base_size = 16)



# network plot ------------------------------------------------------------






