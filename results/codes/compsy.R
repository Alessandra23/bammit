library("ggridges")
library(tidybayes)
library(bayesplot)
color_scheme_set("blue")

sy <- c(0.1, 0.5, 2, 10)

data1 <- simulateDataAmmi(I = 25, J = 6, mu = 100, sg = 10, se = 10, sy = sy[1],
                          lambda = 12, stheta1 = 1, stheta2 = 1)
data2 <- simulateDataAmmi(I = 25, J = 6, mu = 100, sg = 10, se = 10, sy = sy[2],
                          lambda = 12, stheta1 = 1, stheta2 = 1)
data3 <- simulateDataAmmi(I = 25, J = 6, mu = 100, sg = 10, se = 10, sy = sy[3],
                          lambda = 12, stheta1 = 1, stheta2 = 1)
data4 <- simulateDataAmmi(I = 25, J = 6, mu = 100, sg = 10, se = 10, sy = sy[4],
                          lambda = 12, stheta1 = 1, stheta2 = 1)

crossaModel1 <- crossaJags(
  data = data1, mmu = 90, smu = 10,
  mug = 0, mue = 0, a = 1, b = 10, stheta = 10,
  nthin = 1, nburnin = 100
)


crossaModel2 <- crossaJags(
  data = data2, mmu = 90, smu = 10,
  mug = 0, mue = 0, a = 1, b = 10, stheta = 10,
  nthin = 1, nburnin = 100
)

crossaModel3 <- crossaJags(
  data = data3, mmu = 90, smu = 10,
  mug = 0, mue = 0, a = 1, b = 10, stheta = 10,
  nthin = 1, nburnin = 500
)

crossaModel4 <- crossaJags(
  data = data4, mmu = 90, smu = 10,
  mug = 0, mue = 0, a = 1, b = 10, stheta = 10,
  nthin = 1, nburnin = 500
)


qplot(data4$blin, crossaModel3$BUGSoutput$mean$blin, ylab = "Estimated bilinear - Josse", xlab = "True bilinear") + geom_abline() +
  geom_abline() + theme_bw() + geom_point(colour = "red", size = 0.6)



fitcmcm <- as.mcmc(crossaModel4)
fitmat <- as.matrix(fitcmcm)
fitdf <- as.data.frame(fitmat)

blinPostCoeff <- select(fitdf,paste0("blin[", seq(1, 150, by = 1), "]"))
names(blinPostCoeff) <- paste0("blin", 1:ncol(blinPostCoeff))
blinPostCoeffLong <- gather(blinPostCoeff)

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

blinPostCoeffLongSum$key <- as.numeric(gsub('blin', '', blinPostCoeffLongSum$key))
blinPostCoeffLongSum <- blinPostCoeffLongSum[order(blinPostCoeffLongSum$key), ]
blinPostCoeffLongSum <- cbind(blin = data4$blin, blinPostCoeffLongSum, meann = crossaModel4$BUGSoutput$mean$blin)


ggplot(data = blinPostCoeffLongSum, aes(x = meann, y = blin)) +
    geom_errorbarh(aes(xmin = lower_coef, xmax = upper_coef),
                   height = 0.1) +
    geom_point(colour = "steelblue") +
    geom_abline()+
    #geom_vline(xintercept = 0, linetype = "dashed") +
    xlab("Posterior estimate") + ylab("True bilinear") +
    theme_bw()







