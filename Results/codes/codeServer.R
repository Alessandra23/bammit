## All codes on the server


# bammit real data --------------------------------------------------------


bammitJagsRealData <- function(data, Q = 2, mmu, smu, stheta, a, b,
                               nthin = 1, nburnin = 2000){


  # Q <- Q
  # Y <- data$Mean
  # K <- length(levels(data$Year))
  # I <- length(levels(data$Genotype))
  # J <- length(levels(data$Environment))
  # N <- length(data$Mean)
  #
  # genotype <- data[, "Genotype"]
  # environment <- data[, "Environment"]
  # time <- data[, "Year"]

  # levels(data$Genotype) <- (1:length(levels(data$Genotype)))
  # levels(data$Environment) <- 1:length(levels(data$Environment))

  # Q <- Q
  # Y <- data$Mean
  # K <- length(levels(data$Year))
  # I <- length(levels(data$Genotype))
  # J <- length(levels(data$Environment))
  # N <- length(data$Mean)

  Q <- Q
  Y <- data$Mean
  K <- length(levels(as.factor(as.numeric(data$Year))))
  I <- length(levels(as.factor(as.numeric(data$Genotype))))
  J <- length(levels(as.factor(as.numeric(data$Environment))))
  N <- length(data$Mean)

  genotype <- as.factor(as.numeric(data$Genotype))
  environment <- as.factor(as.numeric(data$Environment))
  time <- as.factor(as.numeric(data$Year))

  modelCode <- "
  model
  {
  # Likelihood
   for (n in 1:N) {
    Y[n] ~ dnorm(mu[n], sy)
    mu[n] = muall + g[genotype[n]] + e[environment[n]] + t[time[n]] + blin[n]
    blin[n] = sum(lambda[1:Q] * gamma[genotype[n],1:Q] * delta[environment[n],1:Q]*rho[time[n],1:Q])
   }

  # Priors
  # Prior on grand mean
   muall ~ dnorm(mmu, smu^-2)

  # Prior on genotype effect
  for(i in 1:I) {
  g[i] ~ dnorm(0, sg^-2) # Prior on genotype effect
  }

  # Prior on environment effect
  for(j in 1:J) {
  e[j] ~ dnorm(0, se^-2) # Prior on environment effect
  }

  # Prior on time effect
  for(k in 1:K) {
  t[k] ~ dnorm(0, st^-2) # Prior on time effect
  }


   # Priors on gamma
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

   # Priors on delta
   for(q in 1:Q){
    for(j in 1:J){
      thetaD[j,q] ~ dnorm(0,stheta)
    }
    mD[q] = sum(thetaD[1:J,q])/J
    for(j in 1:J){
    thetaDNew[j,q] = thetaD[j,q] - mD[q]
    }
    sqrtThetaD[q] = sqrt(1/(sum(thetaDNew[1:J,q]^2+ 0.000001)))
    for(j in 1:J){
      delta[j,q] = thetaDNew[j,q]*sqrtThetaD[q]
    }
  }

   # Priors on rho
   for(q in 1:Q){
    for(k in 1:K){
      thetaK[k,q] ~ dnorm(0,stheta)
    }
    mK[q] = sum(thetaK[1:K,q])/K
    for(k in 1:K){
    thetaKNew[k,q] = thetaK[k,q] - mK[q]
    }
    sqrtThetaK[q] = sqrt(1/(sum(thetaKNew[1:K,q]^2+ 0.000001)))
    for(k in 1:K){
      rho[k,q] = thetaKNew[k,q]*sqrtThetaK[q]
    }
  }

  # Prior on eigenvalues
  for(q in 1:Q) {
    lambda_raw[q] ~ dnorm(0, slambda^-2)T(0,)
  }
  lambda = sort(lambda_raw)

  #Priors
  sg ~ dt(0, 1, 1)T(0,)
  se ~ dt(0, 1, 1)T(0,)
  st ~ dt(0, 1, 1)T(0,)
  slambda ~ dt(0, 1, 1)T(0,)

  # Prior on residual standard deviation
   sy ~ dgamma(a, b) # inverse of sy
  }
  "

  # Set up the data
  modelData <- list(
    N = N,
    Y = Y,
    I = I,
    J = J,
    K = K,
    Q = Q,
    genotype = genotype,
    environment = environment,
    time = time,
    mmu = mmu,
    smu = smu,
    # mug = mug,
    # mue = mue,
    # mut = mut,
    #mulambda = mulambda,
    # sg = sg,
    # se = se,
    # st = st,
    # slambda = slambda,
    stheta = stheta,
    a = a,
    b = b
  )

  # Choose the parameters to watch
  modelParameters <- c(
    "g", "e", "t", "lambda", "gamma", "delta", "rho", "sy", "muall", "blin"
  )

  # Run the model
  modelRun <- jags(
    data = modelData,
    parameters.to.save = modelParameters,
    model.file = textConnection(modelCode),
    #progress.bar = "none",
    n.thin = nthin,
    n.burnin = nburnin,
    n.iter = 2000 * nthin + nburnin
  )

  return(modelRun)
}

library(tidyr)
library(dplyr)
library(ggplot2)
library(R2jags)
library(ggridges)

load("~/ireland.RData")
years <- c("2018")
irelandSmall <- ireland %>% dplyr::filter(Year %in% years) %>% as.data.frame()
irelandSmall %>% str()


length(levels(as.factor(as.numeric(irelandSmall$Genotype))))

model <- bammitJagsRealData(
  data = irelandSmall, mmu = 90, smu = 10, stheta = 1/100, a = 0.1, b = 0.1,
  nthin = 1, nburnin = 10, Q = 1
)



fitcmcm <- as.mcmc(model)
fitmat <- as.matrix(fitcmcm)
fitdf <- as.data.frame(fitmat)

#get all
#gPostCoeff <- select(fitdf,paste0("g[", seq(1, ncol(fitg), by = 1), "]"))
gPostCoeff <- select(fitdf,paste0("g[", seq(1, 15, by = 1), "]"))
names(gPostCoeff) <- paste0("g", 1:ncol(gPostCoeff))
gPostCoeffLong <- gather(gPostCoeff)

ggplot(data = gPostCoeffLong, aes(y = value, x = key))+
  geom_boxplot() + theme_bw() +
  labs(x = " ", y = "")


g1 <- ggplot(data = gPostCoeffLong, aes(x = value, y = key)) +
  stat_density_ridges(quantile_lines = TRUE,
                      quantiles = c(0.025, 0.5, 0.975),
                      alpha = 0.7,
                      colour = "steelblue",
                      fill = "steelblue") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  theme_bw() + xlab("Posterior estimate") + ylab("")
g1


ePostCoeff <- select(fitdf,paste0("e[", seq(1, 8, by = 1), "]"))
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

ggplot(data = ePostCoeffLong, aes(y = value, x = key))+
  geom_boxplot() + theme_bw() +
  labs(x = "Genotype (2015)", y = "")


# data set 2 --------------------------------------------------------------


load("~/data2.RData")
years <- c("2014")
data2small <- data2%>% dplyr::filter(YEAR %in% years) %>%
  rename(Genotype = GENOTYPE,
         Environment = LOCATION,
         Year = YEAR,
         Mean = MEAN) %>%
  as.data.frame()

data2small %>% str()

modelData2 <- bammitJagsRealData(
  data = data2small, mmu = 90, smu = 10, stheta = 1/100, a = 0.1, b = 0.1,
  nthin = 1, nburnin = 10, Q = 1
)


fitcmcm <- as.mcmc(modelData2)
fitmat <- as.matrix(fitcmcm)
fitdf <- as.data.frame(fitmat)

fitg <- fitdf[, grep(x = colnames(fitdf), pattern = "g[", fixed = TRUE)]
fite <- fitdf[, grep(x = colnames(fitdf), pattern = "e[", fixed = TRUE)]
fitt <- fitdf[, grep(x = colnames(fitdf), pattern = "t[", fixed = TRUE)]
fitblin <- fitdf[, grep(x = colnames(fitdf), pattern = "blin[", fixed = TRUE)]
fitgamma <- fitdf[, grep(x = colnames(fitdf), pattern = "gamma[", fixed = TRUE)]
#get all
#gPostCoeff <- select(fitdf,paste0("g[", seq(1, ncol(fitg), by = 1), "]"))
gPostCoeff <- select(fitdf,paste0("g[", seq(10, 20, by = 1), "]"))
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

ggplot(data = gPostCoeffLong, aes(y = value, x = key))+
  geom_boxplot() + theme_bw() +
  labs(x = " ", y = "")



ePostCoeff <- select(fitdf,paste0("e[", seq(1, 7, by = 1), "]"))
names(ePostCoeff) <- paste0("e", 1:ncol(ePostCoeff))
ePostCoeffLong <- gather(ePostCoeff)

ggplot(data = ePostCoeffLong, aes(y = value, x = key))+
  geom_boxplot() + theme_bw() +
  labs(x = " ", y = "")

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

levels(unique(data$Year))
as.factor(as.numeric(data$Year))



# ambarti -----------------------------------------------------------------


# plot ambarti function

library(colorspace)
library(ggplot2)
library(patchwork)
library(tidyverse)





# This function is set up to work with the .RData files in the
# zip files containing the results for ambarti on Estevao's GitHub.
# for example: Ireland_VCU_2010_AMBARTI.RData.zip

# To begin, load the data from GitHub.

# The arguments for the plotting function are:
# data: Results from GitHub. e.g: Ireland_VCU_2010_AMBARTI.RData
# fill: use either the bart componant or ambarti to fill the principal heat map
# mainEffPal: Colour palette to colour the main effects
# intPal: Colour palette to colour the principal heat map
# title: Title of the plot




# Plotting function -------------------------------------------------------


plot_ambarti <- function(data,
                         fill = c("bart", "ambarti"),
                         mainEffPal = rev(diverging_hcl(palette = "Blue-Red 3", n = 100)),
                         intPal = rev(diverging_hcl(palette = "Blue-Red 3", n = 100)),
                         reorder = FALSE) {

  # set up names
  features <- attributes(data$x)
  g_names <- colnames(features$contrasts$g)
  e_names <- colnames(features$contrasts$e)
  e <- rep(e_names, times = length(g_names))
  g <- rep(g_names, each = length(e_names))


  g_names <- colnames(dfG)
  e_names <- colnames(dfE)
  e <- rep(e_names, times = length(g_names))
  g <- rep(g_names, each = length(e_names))

  # get yhats
  yhatBart <- apply(data$y_hat_bart, 2, mean) # predict values for bart
  yhat <- y

  # create a df
  df.ambarti <- data.frame(yhat = yhat, e = e, g = g)
  df.ambarti$e <- factor(df.ambarti$e, levels = unique(df.ambarti$e))
  df.ambarti$g <- factor(df.ambarti$g, levels = unique(df.ambarti$g))



  # used for continuous bars in plot
  df.ambarti$e <- match(df.ambarti$e, sort(unique(df.ambarti$e)))
  df.ambarti$g <- match(df.ambarti$g, sort(unique(df.ambarti$g)))

  # set axis tick labels
  labels_e <- paste0("e", 1:length(e_names))
  labels_g <- paste0("g", 1:length(g_names))

  # get the main effect (beta hat) means
  mainEffects <- cbind(dfG, dfE)
  bHatMeans <- apply(mainEffects, 2, mean)# - data$y
  mainLims <- range(bHatMeans)
  mainLims <- range(labeling::rpretty(mainLims[1], mainLims[2]))

  # choose between bart and ambarti
  if (fill == "bart") {
    plotValue = TRUE
    fill <- yhatBart
    name <- "BART\nInteraction\\\nMain Effects"
    intLims <- range(yhatBart)
    limitsInt <- c(floor(intLims[1]), ceiling(intLims[2]))
  } else if (fill == "ambarti") {
    plotValue = FALSE
    fill <- yhat
    name <- "\u0177"
    intLims <- range(yhat)
    limitsInt <- c(floor(intLims[1]), ceiling(intLims[2]))
  }

  # get lengths for plotting main effect bars
  NoOfEnv <- length(unique(df.ambarti$e))
  NoOfGene <- length(unique(df.ambarti$g))
  NoOfBeta <- length(bHatMeans)

  # Plotting ----------------------------------------------------------------


  # main heatmap

  if(reorder){
    breakPattern_g <- as.numeric(str_extract_all(levels(reorder(g, fill)), "[0-9]+"))
    breakPattern_e <- as.numeric(str_extract_all(levels(reorder(e, fill)), "[0-9]+"))

    p_main <- ggplot(df.ambarti, aes(reorder(e, fill), reorder(g, fill))) +
      geom_tile(aes(fill = fill)) +
      scale_y_discrete(breaks = breakPattern_g, labels = levels(reorder(g, fill))) +
      scale_x_discrete(breaks =  breakPattern_e, labels = levels(reorder(e, fill)))
  } else{
    p_main <- ggplot(df.ambarti, aes(e, g)) +
      geom_tile(aes(fill = fill)) +
      scale_y_continuous(breaks = c(1:length(labels_g)), labels = labels_g) +
      scale_x_continuous(breaks = c(1:length(labels_e)), labels = labels_e)
  }

  p_main <- p_main +
    scale_fill_gradientn(
      colors = intPal,
      guide = guide_colorbar(
        frame.colour = "black",
        ticks.colour = "black"
      ),
      name = name
    ) +
    theme_classic() +
    labs(x = "Environment", y = "Genotype")

  if(reorder){
    data_e <- tibble(e = breakPattern_e, bHatMeans = bHatMeans[1:length(breakPattern_e)])
    data_e <- tibble(e = breakPattern_e, bHatMeans = bHatMeans[data_e$e])
    data_g <- tibble(g = breakPattern_g, bHatMeans = bHatMeans[(length(breakPattern_e) + 1):length(bHatMeans)])
    data_g <- tibble(g = breakPattern_g, bHatMeans = bHatMeans[(data_g$g)+length(breakPattern_e)])

    data_e$e <- factor(data_e$e, levels = unique(data_e$e))
    data_g$g <- factor(data_g$g, levels = unique(data_g$g))

  } else{
    data_e <- tibble(e = 1:NoOfEnv, bHatMeans = bHatMeans[1:NoOfEnv])
    data_g <- tibble(g = (NoOfEnv) + 1:NoOfGene, bHatMeans = bHatMeans[(NoOfEnv) + 1:NoOfGene])
  }

  # bottom bar
  p_bottom <- ggplot() +
    geom_tile(
      data = data_e,
      aes(x = e, y = 0, fill = bHatMeans)
    ) +
    scale_fill_gradientn(
      limits = mainLims,
      colors = mainEffPal,
      guide = guide_colorbar(
        frame.colour = "black", ticks.colour = "black"
      )
    ) +
    theme_void() +
    theme(legend.position = "none")

  # left bar
  p_left <- ggplot() +
    geom_tile(
      data = data_g,
      aes(x = 0, y = g, fill = bHatMeans)
    ) +
    scale_fill_gradientn(
      limits = mainLims,
      colors = mainEffPal,
      guide = guide_colorbar(
        frame.colour = "black",
        ticks.colour = "black"
      ),
      name = "Main\neffects"
    ) +
    theme_void() +
    theme(legend.position = "none")

  p_left <- p_left + theme(
    legend.key.size = unit(0.7, "cm"),
    legend.key.width = unit(0.8, "cm")
  )


  # add them all together
  p <- p_left + p_main + plot_spacer() + p_bottom +
    plot_layout(
      guides = "collect",
      heights = c(1, 0.07),
      widths = c(0.07, 1)
    )

  if(plotValue == FALSE){
    return(p_main)
  }else {
    return(p)
  }

}


# Example -------------------------------------------------------------------
plot_ambarti(ambarti, fill = "bart", reorder = T)

# using different colour palette
plot_ambarti(ambarti,
             fill = "ambarti",
             intPal = rev(sequential_hcl(palette = "Plasma", n = 100)), reorder = F)




# ammi josse --------------------------------------------------------------

library(R2jags)
library(reshape2)

generateBlin <- function(index, Q, stheta = 1){
  theta <- matrix(NA, nrow = index, ncol = Q)
  variable <- matrix(NA, nrow = index, ncol = Q)
  sqrtTheta <-  vector()
  for(q in 1:Q){
    for (i in 1:index) {
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

  return(variable)
}



simulateDataAmmi <- function(I, J, mu, sg, se, sy, lambda, stheta) {
  # number of obs
  N <- I * J
  Q <- length(lambda)

  # generate main effects

  g <- rnorm(I, 0, sg)
  e <- rnorm(J, 0, se)

  # generate lambda, gamma, delta and kappa
  gamma <- generateBlin(I, Q)
  delta <- generateBlin(J, Q)

  # gamma <- generate_gamma_delta(I, Q, stheta = stheta1)
  # delta <- generate_gamma_delta(J, Q, stheta = stheta2)

  x <- expand.grid(1:I, 1:J)
  names(x) <- c("g", "e")
  x$g <- as.factor(x$g)
  x$e <- as.factor(x$e)

  # generate bilinear term
  blin <- rep(0, N)
  for (k in 1:length(lambda)) {
    blin <- blin + lambda[k] * gamma[x[, "g"], k] * delta[x[, "e"], k]
  }

  mu_ij <- mu + g[x[, "g"]] + e[x[, "e"]] + blin

  # generate y
  y <- rnorm(N, mu_ij, sy)

  return(list(
    y = y,
    g = g,
    e = e,
    lambda = lambda,
    gamma = gamma,
    delta = delta,
    blin = blin,
    I = I,
    J = J,
    Q = Q,
    x = x,
    sy = sy
  ))
}


josseJags <- function(data, mmu, smu, sg, se, slambda, a, b, nthin = 1, nburnin = 2000){
  modelCode <- '
  model
  {
  # Likelihood
   for (k in 1:N) {
    Y[k] ~ dnorm(mu[k], sy^-2)
    mu[k] = muall + g[genotype[k]] + e[environment[k]] + blin[k]
    blin[k] = sum(lambda[1:Q] * gamma[genotype[k],1:Q] * delta[environment[k],1:Q])
   }

  # Priors
  # Prior on grand mean
   muall ~ dnorm(mmu, smu^-2)

  # Prior on genotype effect
  for(i in 1:I) {
  g[i] ~ dnorm(0, sg^-2) # Prior on genotype effect
  }

  # Prior on environment effect
  for(j in 1:J) {
  e[j] ~ dnorm(0, se^-2) # Prior on environment effect
  }

  # Priors on gamma
  for(q in 1:Q) {
  gamma[1, q] ~ dnorm(0, 1)T(0,) # First one is restriced to be positive
  for(i in 2:I) {
  gamma[i, q] ~ dnorm(0, 1) # Prior on genotype interactions
  }
  }

  # Priors on delta
  for(q in 1:Q) {
  for(j in 1:J) {
  delta[j, q] ~ dnorm(0, 1) # Prior on environment interactions
  }
  }

  # Prior on eigenvalues
  for(q in 1:Q) {
  lambda_raw[q] ~ dnorm(0, slambda^-2)T(0,)
  }
  lambda = sort(lambda_raw)


  # Prior on residual standard deviation
  sy ~ dgamma(a, b)
  }
  '
  # from data
  Y <- data$y
  I <- data$I
  J <- data$J
  Q <- data$Q
  N <- I*J
  genotype <- data$x[, "g"]
  environment <- data$x[, "e"]

  # Set up the data
  modelData <- list(
    N = N,
    Y = Y,
    I = I,
    J = J,
    Q = Q,
    genotype = genotype,
    environment = environment,
    mmu = mmu,
    smu = smu,
    sg = sg,
    se = se,
    slambda = slambda,
    a = a,
    b = b
  )

  # Choose the parameters to watch
  modelParameters <- c(
    "g", "e", "lambda", "gamma", "delta", "sy", "muall", "blin"
  )

  # Run the model
  modelRun <- jags(
    data = modelData,
    parameters.to.save = modelParameters,
    model.file = textConnection(modelCode),
    #progress.bar = "none",
    n.thin = nthin,
    n.burnin = nburnin,
    n.iter = 2000 * nthin + nburnin
  )

  return(modelRun)
}


data <- simulateDataAmmi(I = 12, J = 6, mu = 100, sg = 10, se = 10, sy = 1.5, lambda = 10, stheta = 1)

josseModel <- josseJags(
  data = data, mmu = 90, smu = 10,
  sg = 10, se = 10, slambda = 1,
  a = 0.1, b = 0.1, nthin = 1, nburnin = 1000
)

predictionAMMI <- function(model, data, p) {
  muhat <- model$BUGSoutput$mean$muall
  ghat <- apply(josseModel$BUGSoutput$sims.list$g, 2, function(x){
    quantile(x, p)
  })
  ehat <- apply(josseModel$BUGSoutput$sims.list$e, 2, function(x){
    quantile(x, p)
  })
  blinhat <- apply(josseModel$BUGSoutput$sims.list$blin, 2, function(x){
    quantile(x, p)
  })

  N <- length(data$y)

  yhat <- rep(muhat, N) + ghat[data$x$g] + ehat[data$x$e] +  blinhat

  return(yhat)
}


y05 <- predictionAMMI(model = josseModel,data, p = 0.05)
y50 <- predictionAMMI(josseModel,data, p = 0.5)
y95 <- predictionAMMI(josseModel,data, p = 0.95)

yhat <- y50



dfG <- as.data.frame(josseModel$BUGSoutput$sims.list$g)
names(dfG) <- unique(data$Genotype)
dfE <- as.data.frame(josseModel$BUGSoutput$sims.list$e)
names(dfE) <-  unique(data$Environment)

g_names <- colnames(dfG)
e_names <- colnames(dfE)
e <- rep(e_names, times = length(g_names))
g <- rep(g_names, each = length(e_names))

# create a df
dat <- data.frame(yhat = yhat, e = e, g = g)
intPal = rev(diverging_hcl(palette = "Blue-Red 3", n = 100))

# factor levels
dat$e <- factor(dat$e, levels = unique(dat$e))
dat$g <- factor(dat$g, levels = unique(dat$g))

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
  theme_classic(base_size = 16)




dfG <- as.data.frame(josseModel$BUGSoutput$sims.list$g)
names(dfG) <- paste0("g", seq(1:ncol(dfG)))
dfG <- melt(dfG)
dfGT <- data.frame(variable = levels(dfG$variable), true = data$g)

pG <- ggplot(dfG, aes(x = value)) +
  geom_density(colour = "steelblue", fill = "steelblue", alpha = 0.7) +
  facet_wrap(~variable) +
  labs(x = "Genotype", y = "Frequency") +
  geom_vline(data = dfGT, aes(xintercept = true), colour = "firebrick")+
  theme_bw(base_size = 14) +
  theme(strip.background = element_rect(fill="white"),
        panel.spacing.x = unit(5,"mm"),
        legend.position = c(0.75, 0.04))
pG



dfE <- as.data.frame(josseModel$BUGSoutput$sims.list$e)
names(dfE) <- paste0("e", seq(1:ncol(dfE)))
dfE <- melt(dfE)
dfET <- data.frame(variable = levels(dfE$variable), true = data$e)






g_names <- colnames(dfG)
e_names <- colnames(dfE)
e <- rep(e_names, times = length(g_names))
g <- rep(g_names, each = length(e_names))

yhat <- y95

# create a df
df.ambarti <- data.frame(yhat = yhat, e = e, g = g)
df.ambarti$e <- factor(df.ambarti$e, levels = unique(df.ambarti$e))
df.ambarti$g <- factor(df.ambarti$g, levels = unique(df.ambarti$g))



# used for continuous bars in plot
df.ambarti$e <- match(df.ambarti$e, sort(unique(df.ambarti$e)))
df.ambarti$g <- match(df.ambarti$g, sort(unique(df.ambarti$g)))

# set axis tick labels
labels_e <- paste0("e", 1:length(e_names))
labels_g <- paste0("g", 1:length(g_names))


plotValue = FALSE
fill <- yhat
name <- "\u0177"
intLims <- range(yhat)
limitsInt <- c(floor(intLims[1]), ceiling(intLims[2]))


p_main <- ggplot(df.ambarti, aes(e, g)) +
  geom_tile(aes(fill = fill)) +
  scale_y_continuous(breaks = c(1:length(labels_g)), labels = labels_g) +
  scale_x_continuous(breaks = c(1:length(labels_e)), labels = labels_e) +
  scale_fill_gradientn(
    limits = limitsInt,
    colors = intPal,
    guide = guide_colorbar(
      frame.colour = "black",
      ticks.colour = "black"
    ),
    name = name
  ) +
  theme_classic() +
  labs(x = "Environment", y = "Genotype")
p_main



# simulate data bammit ----------------------------------------------------



library(dplyr)
library(R2jags)
library(ggplot2)

# Functions to simulated data ---------------------------------------------


genBlin <- function(index, Q){
  theta <- matrix(NA, nrow = index, ncol = Q)
  variable <- matrix(NA, nrow = index, ncol = Q)
  sqrtTheta <-  vector()
  for(q in 1:Q){
    for (i in 1:index) {
      theta[i,q] <- rnorm(1)
    }
    m <- apply(theta, 2, mean)
    thetaN <- apply(theta, 1, function(x){x-m}) %>% as.matrix()

    if(Q >1){
      thetaN <- t(thetaN)
    }

    sqrtTheta[q] <- sqrt(1/sum(thetaN[,q]^2))

    for (i in 1:index) {
      variable[i,q] <- (thetaN[i,q])*sqrtTheta[q]
    }
  }

  return(variable)
}


simBAMMIT <- function(V = 3,
                      Q = 1,
                      Bv = c(3,2,2),
                      mu = 100,
                      lambda = c(10),
                      sb = 1,
                      sB = 10,
                      sy = 1){

  N <- Reduce("*", Bv)

  # generate main effects
  bv <- lapply(Bv, function(x) {
    bv <- rnorm(x, 0, sb)
    bv <- bv - mean(bv)
    return(bv)
  })
  #names(bv) <- paste0("b", 1:length(Bv))

  meff <- rowSums(rev(expand.grid(rev(bv))))

  # generate bilinear term
  Beta <- vector(mode = "list", length = V)
  for (i in 1:V) {
    Beta[[i]] <- genBlin(Bv[i], Q)
  }

  k <- as.list(rep(1, Q))
  for (j in 1:length(k)) {
    for (i in 1:length(Beta)) {
      k[[j]] <- kronecker(k[[j]], Beta[[i]][,j])
    }
    k[[j]] <- k[[j]]*lambda[j]
  }

  blin <- Reduce("+", k)

  # generate y
  m <- mu + meff + blin
  y <- rnorm(N, m, sy)

  return(list(mu = mu,
              bv = bv,
              blin = blin,
              y = y,
              Bv = Bv))

}

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
saveRDS(data, file = "Synthetic data/V4_Q3_N960_L81012.rds")


Bv <- list(c(10,5,2),
           c(25,5,4),
           c(50,10,2),
           c(100,10,5),
           c(50,20,10))

genDataSets <- function(Bv){
  for (i in 1:10) {
    data <- simBAMMIT(V = 3,
                      Q = 2,
                      Bv = Bv,
                      mu = 100,
                      lambda = c(8,10, 12),
                      sb = 1,
                      sB = 1,
                      sy = 1)
    N <- Reduce('*', Bv)
    fileName <- paste0("bammit/data/N", N, "Rep", i, ".rds")
    saveRDS(data, file = fileName)
  }
}
lapply(Bv, genDataSets)


# Run jags ----------------------------------------------------------------

bammitJags <- function(data, Q, mmu, smu, a,b, stheta = 1, nthin, nburnin){

  Y <- data$y
  V <- length(data$Bv)
  Bv <- data$Bv
  N <- length(data$y)

  if(V == 2){

    modelCode <- "
    model{
      # Likelihood
      for (n in 1:N) {
          Y[n] ~ dnorm(mu[n], sy)
          mu[n] = muall + b1[var1[n]] + b2[var2[n]] + blin[n]
          blin[n] = sum(lambda[1:Q] * beta1[var1[n],1:Q] * beta2[var2[n],1:Q])
      }

      # Priors

      muall ~ dnorm(mmu, smu^-2) # grand mean

      for(i in 1:B1) {
          b1[i] ~ dnorm(0, sb1^-2) # Prior on effect 1
      }

      for(i in 1:B2) {
          b2[i] ~ dnorm(0, sb2^-2) # Prior on effect 2
      }

      for(q in 1:Q){
          for(i in 1:B1){
            thetaBeta1[i,q] ~ dnorm(0,stheta)
          }
          mBeta1[q] = sum(thetaBeta1[1:B1,q])/B1

          for(i in 1:B1){
            thetaBeta1New[i,q] = thetaBeta1[i,q] - mBeta1[q]
          }

          sqrtThetaBeta1[q] = sqrt(1/(sum(thetaBeta1New[1:B1,q]^2 + 0.000001)))
          for(i in 1:B1){
            beta1[i,q] = thetaBeta1New[i,q]*sqrtThetaBeta1[q]
          }
      }


       for(q in 1:Q){
          for(i in 1:B2){
            thetaBeta2[i,q] ~ dnorm(0,stheta)
          }
          mBeta2[q] = sum(thetaBeta2[1:B2,q])/B2

          for(i in 1:B2){
            thetaBeta2New[i,q] = thetaBeta2[i,q] - mBeta2[q]
          }

          sqrtThetaBeta2[q] = sqrt(1/(sum(thetaBeta2New[1:B2,q]^2 + 0.000001)))
          for(i in 1:B2){
            beta2[i,q] = thetaBeta2New[i,q]*sqrtThetaBeta2[q]
          }
      }

      # Prior on eigenvalues
      for(q in 1:Q) {
        lambda_raw[q] ~ dnorm(0, slambda^-2)T(0,)
      }
      lambda = sort(lambda_raw)

      #Priors
      sb1 ~ dt(0, 1, 1)T(0,)
      sb2 ~ dt(0, 1, 1)T(0,)
      slambda ~ dt(0, 1, 1)T(0,)

      sy ~ dgamma(a, b) # Prior on residual standard deviation - inverse of sy

    }
    "

    B1 <- Bv[1]
    B2 <- Bv[2]

    x <- rev(expand.grid(lapply(rev(Bv), function(x) 1:x)))
    colnames(x) <- paste0('var', 1:2)
    var1 <- x$var1
    var2 <- x$var2

    # Set up the data
    modelData <- list(
      N = N,
      Y = Y,
      B1 = B1,
      B2 = B2,
      Q = Q,
      var1 = var1,
      var2 = var2,
      mmu = mmu,
      smu = smu,
      stheta = stheta,
      a = a,
      b = b
    )

    # Choose the parameters to watch
    modelParameters <- c(
      "b1", "b2",  "lambda", "beta1", "beta2", "sy", "muall", "blin"
    )


  }

  if(V == 3){

    modelCode <- "
    model{
      # Likelihood
      for (n in 1:N) {
          Y[n] ~ dnorm(mu[n], sy)
          mu[n] = muall + b1[var1[n]] + b2[var2[n]] + b3[var3[n]] + blin[n]
          blin[n] = sum(lambda[1:Q] * beta1[var1[n],1:Q] * beta2[var2[n],1:Q]* beta3[var3[n],1:Q])
      }

      # Priors

      muall ~ dnorm(mmu, smu^-2) # grand mean

      for(i in 1:B1) {
          b1[i] ~ dnorm(0, sb1^-2) # Prior on effect 1
      }

      for(i in 1:B2) {
          b2[i] ~ dnorm(0, sb2^-2) # Prior on effect 2
      }

      for(i in 1:B3) {
          b3[i] ~ dnorm(0, sb3^-2) # Prior on effect 3
      }

      for(q in 1:Q){
          for(i in 1:B1){
            thetaBeta1[i,q] ~ dnorm(0,stheta)
          }
          mBeta1[q] = sum(thetaBeta1[1:B1,q])/B1

          for(i in 1:B1){
            thetaBeta1New[i,q] = thetaBeta1[i,q] - mBeta1[q]
          }

          sqrtThetaBeta1[q] = sqrt(1/(sum(thetaBeta1New[1:B1,q]^2 + 0.000001)))
          for(i in 1:B1){
            beta1[i,q] = thetaBeta1New[i,q]*sqrtThetaBeta1[q]
          }
      }


       for(q in 1:Q){
          for(i in 1:B2){
            thetaBeta2[i,q] ~ dnorm(0,stheta)
          }
          mBeta2[q] = sum(thetaBeta2[1:B2,q])/B2

          for(i in 1:B2){
            thetaBeta2New[i,q] = thetaBeta2[i,q] - mBeta2[q]
          }

          sqrtThetaBeta2[q] = sqrt(1/(sum(thetaBeta2New[1:B2,q]^2 + 0.000001)))
          for(i in 1:B2){
            beta2[i,q] = thetaBeta2New[i,q]*sqrtThetaBeta2[q]
          }
       }

      for(q in 1:Q){
          for(i in 1:B3){
            thetaBeta3[i,q] ~ dnorm(0,stheta)
          }
          mBeta3[q] = sum(thetaBeta3[1:B3,q])/B3

          for(i in 1:B3){
            thetaBeta3New[i,q] = thetaBeta3[i,q] - mBeta3[q]
          }

          sqrtThetaBeta3[q] = sqrt(1/(sum(thetaBeta3New[1:B3,q]^2 + 0.000001)))
          for(i in 1:B3){
            beta3[i,q] = thetaBeta3New[i,q]*sqrtThetaBeta3[q]
          }
      }

      # Prior on eigenvalues
      for(q in 1:Q) {
        lambda_raw[q] ~ dnorm(0, slambda^-2)T(0,)
      }
      lambda = sort(lambda_raw)

      #Priors
      sb1 ~ dt(0, 1, 1)T(0,)
      sb2 ~ dt(0, 1, 1)T(0,)
      sb3 ~ dt(0, 1, 1)T(0,)
      slambda ~ dt(0, 1, 1)T(0,)

      sy ~ dgamma(a, b) # Prior on residual standard deviation - inverse of sy

    }
    "

    B1 <- Bv[1]
    B2 <- Bv[2]
    B3 <- Bv[3]

    x <- rev(expand.grid(lapply(rev(Bv), function(x) 1:x)))
    colnames(x) <- paste0('var', 1:3)
    var1 <- x$var1
    var2 <- x$var2
    var3 <- x$var3

    # Set up the data
    modelData <- list(
      N = N,
      Y = Y,
      B1 = B1,
      B2 = B2,
      B3 = B3,
      Q = Q,
      var1 = var1,
      var2 = var2,
      var3 = var3,
      mmu = mmu,
      smu = smu,
      stheta = stheta,
      a = a,
      b = b
    )

    # Choose the parameters to watch
    modelParameters <- c(
      "b1", "b2", "b3" ,"lambda", "beta1", "beta2", "beta3", "sy", "muall", "blin"
    )

  }

  if(V == 4){

    modelCode <- "
    model{
      # Likelihood
      for (n in 1:N) {
          Y[n] ~ dnorm(mu[n], sy)
          mu[n] = muall + b1[var1[n]] + b2[var2[n]] + b3[var3[n]] + b4[var4[n]]  + blin[n]
          blin[n] = sum(lambda[1:Q] * beta1[var1[n],1:Q] * beta2[var2[n],1:Q] * beta3[var3[n],1:Q] * beta4[var4[n],1:Q])
      }

      # Priors

      muall ~ dnorm(mmu, smu^-2) # grand mean

      for(i in 1:B1) {
          b1[i] ~ dnorm(0, sb1^-2) # Prior on effect 1
      }

      for(i in 1:B2) {
          b2[i] ~ dnorm(0, sb2^-2) # Prior on effect 2
      }

      for(i in 1:B3) {
          b3[i] ~ dnorm(0, sb3^-2) # Prior on effect 3
      }

      for(i in 1:B4) {
          b4[i] ~ dnorm(0, sb4^-2) # Prior on effect 4
      }

      for(q in 1:Q){
          for(i in 1:B1){
            thetaBeta1[i,q] ~ dnorm(0,stheta)
          }
          mBeta1[q] = sum(thetaBeta1[1:B1,q])/B1

          for(i in 1:B1){
            thetaBeta1New[i,q] = thetaBeta1[i,q] - mBeta1[q]
          }

          sqrtThetaBeta1[q] = sqrt(1/(sum(thetaBeta1New[1:B1,q]^2 + 0.000001)))
          for(i in 1:B1){
            beta1[i,q] = thetaBeta1New[i,q]*sqrtThetaBeta1[q]
          }
      }


       for(q in 1:Q){
          for(i in 1:B2){
            thetaBeta2[i,q] ~ dnorm(0,stheta)
          }
          mBeta2[q] = sum(thetaBeta2[1:B2,q])/B2

          for(i in 1:B2){
            thetaBeta2New[i,q] = thetaBeta2[i,q] - mBeta2[q]
          }

          sqrtThetaBeta2[q] = sqrt(1/(sum(thetaBeta2New[1:B2,q]^2 + 0.000001)))
          for(i in 1:B2){
            beta2[i,q] = thetaBeta2New[i,q]*sqrtThetaBeta2[q]
          }
       }

      for(q in 1:Q){
          for(i in 1:B3){
            thetaBeta3[i,q] ~ dnorm(0,stheta)
          }
          mBeta3[q] = sum(thetaBeta3[1:B3,q])/B3

          for(i in 1:B3){
            thetaBeta3New[i,q] = thetaBeta3[i,q] - mBeta3[q]
          }

          sqrtThetaBeta3[q] = sqrt(1/(sum(thetaBeta3New[1:B3,q]^2 + 0.000001)))
          for(i in 1:B3){
            beta3[i,q] = thetaBeta3New[i,q]*sqrtThetaBeta3[q]
          }
      }

      for(q in 1:Q){
          for(i in 1:B4){
            thetaBeta4[i,q] ~ dnorm(0,stheta)
          }
          mBeta4[q] = sum(thetaBeta4[1:B4,q])/B4

          for(i in 1:B4){
            thetaBeta4New[i,q] = thetaBeta4[i,q] - mBeta4[q]
          }

          sqrtThetaBeta4[q] = sqrt(1/(sum(thetaBeta4New[1:B4,q]^2 + 0.000001)))
          for(i in 1:B4){
            beta4[i,q] = thetaBeta4New[i,q]*sqrtThetaBeta4[q]
          }
      }

      # Prior on eigenvalues
      for(q in 1:Q) {
        lambda_raw[q] ~ dnorm(0, slambda^-2)T(0,)
      }
      lambda = sort(lambda_raw)

      #Priors
      sb1 ~ dt(0, 1, 1)T(0,)
      sb2 ~ dt(0, 1, 1)T(0,)
      sb3 ~ dt(0, 1, 1)T(0,)
      sb4 ~ dt(0, 1, 1)T(0,)
      slambda ~ dt(0, 1, 1)T(0,)

      sy ~ dgamma(a, b) # Prior on residual standard deviation - inverse of sy

    }
    "

    B1 <- Bv[1]
    B2 <- Bv[2]
    B3 <- Bv[3]
    B4 <- Bv[4]

    x <- rev(expand.grid(lapply(rev(Bv), function(x) 1:x)))
    colnames(x) <- paste0('var', 1:4)
    var1 <- x$var1
    var2 <- x$var2
    var3 <- x$var3
    var4 <- x$var4

    # Set up the data
    modelData <- list(
      N = N,
      Y = Y,
      B1 = B1,
      B2 = B2,
      B3 = B3,
      B4 = B4,
      Q = Q,
      var1 = var1,
      var2 = var2,
      var3 = var3,
      var4 = var4,
      mmu = mmu,
      smu = smu,
      stheta = stheta,
      a = a,
      b = b
    )

    # Choose the parameters to watch
    modelParameters <- c(
      "b1", "b2", "b3", "b4", "lambda", "beta1", "beta2", "beta3", "beta4", "sy", "muall", "blin"
    )

  }

  # Run the model
  modelRun <- jags(
    data = modelData,
    parameters.to.save = modelParameters,
    model.file = textConnection(modelCode),
    #progress.bar = "none",
    n.thin = nthin,
    n.burnin = nburnin,
    n.iter = 2000 * nthin + nburnin
  )

  return(modelRun)

}

for (n in c(100,500,1000,5000, 10000)) {
  for (i in 1:10) {
    fileName <-  paste0("~/pineR/bammit/data/N", n, "Rep", i, ".rds")
    data <- readRDS(fileName)
    model <- bammitJags(data = data,
                        Q = 1,
                        mmu = 100,
                        smu = 10,
                        a = 0.1,
                        b = 0.1,
                        nthin = 1,
                        nburnin = 2000)

    saveRDS(model, file = paste0("~/pineR/bammit/Running models/N", n, "Rep", i, ".rds"))
  }
}


# read the data sets ------------------------------------------------------

data <- list("100" = list(),
             "500" = list(),
             "1000" = list(),
             "5000" = list(),
             "10000" = list())
for (n in c(100,500,1000,5000, 10000)) {
  nn <- as.character(n)
  for (i in 1:10) {
    fileName <-  paste0("~/pineR/bammit/data/N", n, "Rep", i, ".rds")
    data[[nn]][[i]] <- readRDS(fileName)
  }
}

model <- list("100" = list(),
              "500" = list(),
              "1000" = list(),
              "5000" = list(),
              "10000" = list())

for (n in c(100,500,1000,5000, 10000)) {
  nn <- as.character(n)
  for (i in 1:10) {
    fileName <-  paste0("~/pineR/bammit/Running models/N", n, "Rep", i, ".rds")
    modelAux <- readRDS(fileName)
    model[[nn]][[i]] <- modelAux$BUGSoutput$mean
  }
}

RMSE <- function(true, predicted){
  sqrt(mean(((true - predicted))^2))
}


# get each parameter ------------------------------------------------------

#get the data
d100 <- data[['100']]
d100Blin <- list()
for(i in 1:10){
  d100Blin[[i]] <- d100[[i]]$blin
}

d100b1 <- list()
for(i in 1:10){
  d100b1[[i]] <- d100[[i]]$bv[[1]]
}

d100b2 <- list()
for(i in 1:10){
  d100b2[[i]] <- d100[[i]]$bv[[2]]
}

d100b3 <- list()
for(i in 1:10){
  d100b3[[i]] <- d100[[i]]$bv[[3]]
}

# get the models
m100 <- model$`100`
m100Blin <- list()
for(i in 1:10){
  m100Blin[[i]] <- m100[[i]]$blin
}

m100b1 <- list()
for(i in 1:10){
  m100b1[[i]] <- m100[[i]][[1]]
}

m100b2 <- list()
for(i in 1:10){
  m100b2[[i]] <- m100[[i]][[2]]
}

m100b3 <- list()
for(i in 1:10){
  m100b3[[i]] <- m100[[i]][[3]]
}


# get the rmse
rmseBlin <- c()
for (i in 1:10) {
  rmseBlin[i] <- RMSE(d100Blin[[i]], m100Blin[[i]])
}

rmseb1 <- c()
for (i in 1:10) {
  rmseb1[i] <- RMSE(d100b1[[i]], m100b1[[i]])
}

rmseb2 <- c()
for (i in 1:10) {
  rmseb2[i] <- RMSE(d100b2[[i]], m100b2[[i]])
}

rmseb3 <- c()
for (i in 1:10) {
  rmseb3[i] <- RMSE(d100b3[[i]], m100b3[[i]])
}

# df of rmse
rmsedf100 <- data.frame(b1 = rmseb1, b2 = rmseb2 , b3 = rmseb3, blin = rmseBlin)

## 500


#get the data
d500 <- data[['500']]
d500Blin <- list()
for(i in 1:10){
  d500Blin[[i]] <- d500[[i]]$blin
}

d500b1 <- list()
for(i in 1:10){
  d500b1[[i]] <- d500[[i]]$bv[[1]]
}

d500b2 <- list()
for(i in 1:10){
  d500b2[[i]] <- d500[[i]]$bv[[2]]
}

d500b3 <- list()
for(i in 1:10){
  d500b3[[i]] <- d500[[i]]$bv[[3]]
}

# get the models
m500 <- model$`500`
m500Blin <- list()
for(i in 1:10){
  m500Blin[[i]] <- m500[[i]]$blin
}

m500b1 <- list()
for(i in 1:10){
  m500b1[[i]] <- m500[[i]][[1]]
}

m500b2 <- list()
for(i in 1:10){
  m500b2[[i]] <- m500[[i]][[2]]
}

m500b3 <- list()
for(i in 1:10){
  m500b3[[i]] <- m500[[i]][[3]]
}


# get the rmse
rmseBlin500 <- c()
for (i in 1:10) {
  rmseBlin500[i] <- RMSE(d500Blin[[i]], m500Blin[[i]])
}

rmseb1500 <- c()
for (i in 1:10) {
  rmseb1500[i] <- RMSE(d500b1[[i]], m500b1[[i]])
}

rmseb2500 <- c()
for (i in 1:10) {
  rmseb2500[i] <- RMSE(d500b2[[i]], m500b2[[i]])
}

rmseb3500 <- c()
for (i in 1:10) {
  rmseb3500[i] <- RMSE(d500b3[[i]], m500b3[[i]])
}

# df of rmse
rmsedf500 <- data.frame(b1 = rmseb1500, b2 = rmseb2500 , b3 = rmseb3500, blin = rmseBlin500)

## 1000


#get the data
d1000 <- data[['1000']]
d1000Blin <- list()
for(i in 1:10){
  d1000Blin[[i]] <- d1000[[i]]$blin
}

d1000b1 <- list()
for(i in 1:10){
  d1000b1[[i]] <- d1000[[i]]$bv[[1]]
}

d1000b2 <- list()
for(i in 1:10){
  d1000b2[[i]] <- d1000[[i]]$bv[[2]]
}

d1000b3 <- list()
for(i in 1:10){
  d1000b3[[i]] <- d1000[[i]]$bv[[3]]
}

# get the models
m1000 <- model$`1000`
m1000Blin <- list()
for(i in 1:10){
  m1000Blin[[i]] <- m1000[[i]]$blin
}

m1000b1 <- list()
for(i in 1:10){
  m1000b1[[i]] <- m1000[[i]][[1]]
}

m1000b2 <- list()
for(i in 1:10){
  m1000b2[[i]] <- m1000[[i]][[2]]
}

m1000b3 <- list()
for(i in 1:10){
  m1000b3[[i]] <- m1000[[i]][[3]]
}


# get the rmse
rmseBlin1000 <- c()
for (i in 1:10) {
  rmseBlin1000[i] <- RMSE(d1000Blin[[i]], m1000Blin[[i]])
}

rmseb11000 <- c()
for (i in 1:10) {
  rmseb11000[i] <- RMSE(d1000b1[[i]], m1000b1[[i]])
}

rmseb21000 <- c()
for (i in 1:10) {
  rmseb21000[i] <- RMSE(d1000b2[[i]], m1000b2[[i]])
}

rmseb31000 <- c()
for (i in 1:10) {
  rmseb31000[i] <- RMSE(d1000b3[[i]], m1000b3[[i]])
}

# df of rmse
rmsedf1000 <- data.frame(b1 = rmseb11000, b2 = rmseb21000 , b3 = rmseb31000, blin = rmseBlin1000)

## 5000


#get the data
d5000 <- data[['5000']]
d5000Blin <- list()
for(i in 1:10){
  d5000Blin[[i]] <- d5000[[i]]$blin
}

d5000b1 <- list()
for(i in 1:10){
  d5000b1[[i]] <- d5000[[i]]$bv[[1]]
}

d5000b2 <- list()
for(i in 1:10){
  d5000b2[[i]] <- d5000[[i]]$bv[[2]]
}

d5000b3 <- list()
for(i in 1:10){
  d5000b3[[i]] <- d5000[[i]]$bv[[3]]
}

# get the models
m5000 <- model$`5000`
m5000Blin <- list()
for(i in 1:10){
  m5000Blin[[i]] <- m5000[[i]]$blin
}

m5000b1 <- list()
for(i in 1:10){
  m5000b1[[i]] <- m5000[[i]][[1]]
}

m5000b2 <- list()
for(i in 1:10){
  m5000b2[[i]] <- m5000[[i]][[2]]
}

m5000b3 <- list()
for(i in 1:10){
  m5000b3[[i]] <- m5000[[i]][[3]]
}


# get the rmse
rmseBlin5000 <- c()
for (i in 1:10) {
  rmseBlin5000[i] <- RMSE(d5000Blin[[i]], m5000Blin[[i]])
}

rmseb15000 <- c()
for (i in 1:10) {
  rmseb15000[i] <- RMSE(d5000b1[[i]], m5000b1[[i]])
}

rmseb25000 <- c()
for (i in 1:10) {
  rmseb25000[i] <- RMSE(d5000b2[[i]], m5000b2[[i]])
}

rmseb35000 <- c()
for (i in 1:10) {
  rmseb35000[i] <- RMSE(d5000b3[[i]], m5000b3[[i]])
}

# df of rmse
rmsedf5000 <- data.frame(b1 = rmseb15000, b2 = rmseb25000 , b3 = rmseb35000, blin = rmseBlin5000)


ll <- list(
  l1 = rmsedf100,
  l2 = rmsedf500,
  l3 = rmsedf1000,
  l4 = rmsedf5000)

df <- melt(ll)
df$code <- rep(paste0("v",1:40))

dff <- df %>% group_by(code) %>% group_split()

# ggplot() +
# lapply(dff, function(x) {
#   geom_line(data = x, aes(L1, value, group = 1)) +
#     facet_wrap(~.variable)}
# )


xlabs <- c("100", "500", "1000", "5000")

vnames <-list(
  'b1' = bquote(b[1]),
  'b2' = bquote(b[2]),
  'b3' = bquote(b[3]),
  'int' = "Interaction")
lab <- function(variable,value){
  return(vnames[value])
}
ggplot(bind_rows(dff, .id="df"), aes(L1, value, group = df)) +
  geom_line(colour =  "steelblue") +
  labs(x = "N", y = "RMSE") +
  facet_wrap(~variable, nrow = 1,
             labeller = lab) +
  theme_bw(base_size = 14) +
  theme(strip.background = element_blank(),
        panel.spacing.x = unit(5,"mm"),
        legend.position = c(0.75, 0.04)) +
  scale_x_discrete(labels = xlabs)


geom_point(colour =  "steelblue", size = 0.2) +
  facet_wrap(~v,
             labeller = as_labeller(lab)) +
  labs(x = "True", y = "Estimated") +
  theme_bw(base_size = 14) +
  theme(strip.background = element_blank(),
        panel.spacing.x = unit(5,"mm"),
        legend.position = c(0.75, 0.04))


abc <- data.frame(values = c(rmsedf100$b1 %>% median(),
                             rmsedf500$b1 %>% median(),
                             rmsedf1000$b1 %>% median(),
                             rmsedf5000$b1 %>% median()),
                  sizes = c(1,2,3,4))

abc %>% ggplot(aes(x = sizes, y = values)) +
  geom_line()




#### The other data

V2_Q1_N120_L10 <- readRDS("~/pineR/bammit/data/V2_Q1_N120_L10.rds") # model 1
V3_Q1_N480_L10 <- readRDS("~/pineR/bammit/data/V3_Q1_N480_L10.rds") # model 2
V4_Q1_N960_L10 <- readRDS("~/pineR/bammit/data/V4_Q1_N960_L10.rds") # model 3
V2_Q2_N120_L1012 <- readRDS("~/pineR/bammit/data/V2_Q2_N120_L1012.rds") # model 4
V3_Q2_N480_L1012 <- readRDS("~/pineR/bammit/data/V3_Q2_N480_L1012.rds") # model 5
V4_Q2_N960_L1012 <- readRDS("~/pineR/bammit/data/V4_Q2_N960_L1012.rds") # model 6
V2_Q3_N120_L81012 <- readRDS("~/pineR/bammit/data/V2_Q3_N120_L81012.rds") # model 7
V3_Q3_N480_L81012 <- readRDS("~/pineR/bammit/data/V3_Q3_N480_L81012.rds") # model 8
V4_Q3_N960_L81012 <- readRDS("~/pineR/bammit/data/V4_Q3_N960_L81012.rds") # model 9

data <- list(V2_Q1_N120_L10, V3_Q1_N480_L10, V4_Q1_N960_L10, V2_Q2_N120_L1012, V3_Q2_N480_L1012,
             V4_Q2_N960_L1012 ,V2_Q3_N120_L81012,V3_Q3_N480_L81012,V4_Q3_N960_L81012
)
names(data) <- paste0("model", 1:9)

lapply(1:length(data), function(i){
  modelAux <- bammitJags(data = data[[i]],
                         Q = 2,
                         mmu = 100,
                         smu = 10,
                         a = 0.1,
                         b = 0.1,
                         nthin = 1,
                         nburnin = 2000)
  model <- modelAux$BUGSoutput
  saveRDS(model, file = paste0("~/pineR/bammit/Running models/model", i, ".rds"))
})



# Call the data  ----------------------------------------------------------

model <- list()
for (i in 1:9) {
  model[[i]] <- readRDS(paste0("~/pineR/bammit/Running models/model", i, ".rds"))
}

dataInt <- list()
for (i in 1:9) {
  dataInt[[i]] <- data[[i]]$blin
}

predInt <- list()
for (i in 1:9) {
  predInt[[i]] <- model[[i]]$mean$blin
}
names(predInt) <- paste0("model", 1:9)


quantInt <- vector("list", 9)
for (i in 1:9) {
  aux <- model[[i]]$sims.list$blin
  quantInt[[i]]$lower <- apply(aux, 2, function(x) quantile(x, 0.05))
  quantInt[[i]]$upper <- apply(aux, 2, function(x) quantile(x, 0.95))
}

for (i in 1:9) {
  aux <- model[[i]]$sims.list$blin
  quantInt[[i]]$lower <- apply(aux, 2, function(x) mean(x) - 1.96*sd(x))
  quantInt[[i]]$upper <- apply(aux, 2, function(x) mean(x) + 1.96*sd(x))
}
#names(quantInt) <-  paste0("model", 1:9)
uppers <- lapply(quantInt, function(x) x[["upper"]]) %>% melt() %>% dplyr::select(value)
lowers <- lapply(quantInt, function(x) x[["lower"]]) %>% melt() %>% dplyr::select(value)
colnames(uppers) <- "upper"
colnames(lowers) <- "lower"



predInt <- reshape2::melt(predInt) %>% dplyr::select(value, L1)
colnames(predInt) <- c("int", "model")
dataInt <- reshape2::melt(dataInt) %>% dplyr::select(value)
colnames(dataInt) <- "true"
Vs <- c(rep('V2', 120), rep('V3', 480), rep('V4', 960),
        rep('V2', 120), rep('V3', 480), rep('V4', 960),
        rep('V2', 120), rep('V3', 480), rep('V4', 960))
Qs <- c(rep("Q1", 1560),
        rep("Q2", 1560),
        rep("Q3", 1560))
predInt <- cbind(predInt, dataInt, uppers, lowers, Vs, Qs)

str(predInt)

predInt %>% filter(Qs == "Q1")
q3 <- predInt[predInt$Qs == "Q1", ]

lab <- function(x) c("V = 2", "V = 3", "V = 4")
q3 %>%
  ggplot(aes(x = true, y = int)) +
  geom_abline(size = 0.3, col = "gray") +
  geom_linerange(aes(ymin =  lower, ymax = upper), alpha = 0.2, size = 0.1) +
  geom_point(colour =  "steelblue", size = 0.2) +
  facet_wrap(~Vs,
             labeller = as_labeller(lab)) +
  labs(x = "True", y = "Estimated") +
  theme_bw(base_size = 14) +
  theme(strip.background = element_blank(),
        panel.spacing.x = unit(5,"mm"),
        legend.position = c(0.75, 0.04))


predInt$Qs <- factor(predInt$Qs,labels = c("1", "2", "3"))

lab <- function(x) c("V = 2", "V = 3", "V = 4")
predInt %>%
  ggplot(aes(x = true,
             y = int,
             fill = Qs))+
  #colour = Qs)) +
  # Add a ribbon with the confidence band
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha=0.4) +
  facet_wrap(~Vs,labeller = as_labeller(lab)) +
  labs(x = "True", y = "Estimated", fill = "Q") +
  theme_bw(base_size = 14) +
  theme(strip.background = element_blank(),
        panel.spacing.x = unit(5,"mm"),
        legend.position = "bottom")




# real data bammit and metrics --------------------------------------------



load("~/pineR/bammit/data/train_ireland.RData")
train
library(R2jags)
data <- train


R2 <- function(pred, obs){
  yobshat <- mean(obs)
  R2 <- 1 - (sum((obs - pred)^2)/sum((obs - yobshat)^2))
  return(R2)
}

bammitJagsRealData <- function(data, Q = 1, mmu = 10, smu = 2, stheta = 1, a = 1, b = 0.1,
                               nthin = 2, nburnin = 2000){

  Q <- Q
  Y <- data$Yield
  K <- length(levels(unique(data$Year)))
  I <- length(levels(unique(data$Genotype)))
  J <- length(levels(unique(data$Environment)))
  B <- length(levels(unique(data$Bloc)))
  N <- length(data$Yield)

  genotype <- data$Genotype
  environment <- data$Environment
  time <- data$Year
  bloc <- data$Bloc

  modelCode <- "
  model
  {
  # Likelihood
   for (n in 1:N) {
    Y[n] ~ dnorm(mu[n], sy)
    mu[n] = muall + g[genotype[n]] + e[environment[n]] + t[time[n]] + bl[bloc[n]] + blin[n]
    blin[n] = sum(lambda[1:Q] * gamma[genotype[n],1:Q] * delta[environment[n],1:Q] * rho[time[n], 1:Q] * kappa[bloc[n], 1:Q])
   }

  # Priors
  # Prior on grand mean
   muall ~ dnorm(mmu, smu^-2)

  # Prior on genotype effect
  for(i in 1:I) {
  g[i] ~ dnorm(0, sg^-2) # Prior on genotype effect
  }

  # Prior on environment effect
  for(j in 1:J) {
  e[j] ~ dnorm(0, se^-2) # Prior on environment effect
  }

  # Prior on time effect
  for(k in 1:K) {
  t[k] ~ dnorm(0, st^-2) # Prior on time effect
  }

   # Prior on bloc effect
  for(l in 1:B) {
  bl[l] ~ dnorm(0, sb^-2) # Prior on time effect
  }

   # Priors on gamma
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

   # Priors on delta
   for(q in 1:Q){
    for(j in 1:J){
      thetaD[j,q] ~ dnorm(0,stheta)
    }
    mD[q] = sum(thetaD[1:J,q])/J
    for(j in 1:J){
    thetaDNew[j,q] = thetaD[j,q] - mD[q]
    }
    sqrtThetaD[q] = sqrt(1/(sum(thetaDNew[1:J,q]^2+ 0.000001)))
    for(j in 1:J){
      delta[j,q] = thetaDNew[j,q]*sqrtThetaD[q]
    }
  }

   # Priors on rho
   for(q in 1:Q){
    for(k in 1:K){
      thetaK[k,q] ~ dnorm(0,stheta)
    }
    mK[q] = sum(thetaK[1:K,q])/K
    for(k in 1:K){
    thetaKNew[k,q] = thetaK[k,q] - mK[q]
    }
    sqrtThetaK[q] = sqrt(1/(sum(thetaKNew[1:K,q]^2+ 0.000001)))
    for(k in 1:K){
      rho[k,q] = thetaKNew[k,q]*sqrtThetaK[q]
    }
   }


   # Priors on kappa
   for(q in 1:Q){
    for(l in 1:B){
      thetaB[l,q] ~ dnorm(0,stheta)
    }
    mB[q] = sum(thetaB[1:B,q])/B
    for(l in 1:B){
    thetaBNew[l,q] = thetaB[l,q] - mB[q]
    }
    sqrtThetaB[q] = sqrt(1/(sum(thetaBNew[1:B,q]^2+ 0.000001)))
    for(l in 1:B){
      kappa[l,q] = thetaBNew[l,q]*sqrtThetaB[q]
    }
  }

  # Prior on eigenvalues
  for(q in 1:Q) {
    lambda_raw[q] ~ dnorm(0, slambda^-2)T(0,)
  }
  lambda = sort(lambda_raw)

  #Priors
  sg ~ dt(0, 1, 1)T(0,)
  se ~ dt(0, 1, 1)T(0,)
  st ~ dt(0, 1, 1)T(0,)
  sb ~ dt(0, 1, 1)T(0,)
  slambda ~ dt(0, 1, 1)T(0,)

  # Prior on residual standard deviation
   sy ~ dgamma(a, b) # inverse of sy
  }
  "

  # Set up the data
  modelData <- list(
    N = N,
    Y = Y,
    I = I,
    J = J,
    K = K,
    B = B,
    Q = Q,
    genotype = genotype,
    environment = environment,
    time = time,
    bloc = bloc,
    mmu = mmu,
    smu = smu,
    # mug = mug,
    # mue = mue,
    # mut = mut,
    #mulambda = mulambda,
    # sg = sg,
    # se = se,
    # st = st,
    # slambda = slambda,
    stheta = stheta,
    a = a,
    b = b
  )

  # Choose the parameters to watch
  modelParameters <- c(
    "g", "e", "t", "bl", "lambda", "gamma", "delta", "rho", "kappa","sy", "muall", "blin"
  )

  # Run the model
  modelRun <- jags(
    data = modelData,
    parameters.to.save = modelParameters,
    model.file = textConnection(modelCode),
    #progress.bar = "none",
    n.thin = nthin,
    n.burnin = nburnin,
    n.iter = 2000 * nthin + nburnin
  )

  return(modelRun)
}



trainQ1 <- bammitJagsRealData(data = train, Q = 1)
saveRDS(trainQ1, file = "~/pineR/bammit/Running models/trainQ1.rds")

trainQ2 <- bammitJagsRealData(data = train, Q = 2)
saveRDS(trainQ2, file = "~/pineR/bammit/Running models/trainQ2.rds")

trainQ3 <- bammitJagsRealData(data = train, Q = 3)
saveRDS(trainQ3, file = "~/pineR/bammit/Running models/trainQ3.rds")

trainQ1$BUGSoutput$median$blin %>% range()

model = trainQ1$BUGSoutput
saveRDS(model, file = "~/pineR/bammit/Running models/modeltrainQ1.rds")

predictionBAMMITReal <- function(model, data) {


  genotype <- data$Genotype
  environment <- data$Environment
  time <- data$Year
  bloc <- data$Bloc

  muhat <- model$mean$muall
  g <- model$mean$g
  e <- model$mean$e
  t <- model$mean$t
  bl <- model$mean$bl
  inthat <- model$mean$blin

  N <- length(data$Yield)
  yhat <- rep(muhat, N) + g[genotype] + e[environment] + t[time] + bl[bloc]  +  inthat


  return(yhat)
}

load("~/pineR/bammit/data/test_ireland.RData")
test$Yield
test

predR <- predictionBAMMITReal(trainQ1$BUGSoutput, test)
round(caret::RMSE(predR, test$Yield),2)
round(caret::R2(predR, test$Yield), 2)




# plots -------------------------------------------------------------------


V3_Q1_N480_L10 <- readRDS("~/pineR/bammit/data/V3_Q1_N480_L10.rds") # model 2
V3_Q2_N480_L1012 <- readRDS("~/pineR/bammit/data/V3_Q2_N480_L1012.rds") # model 5
V3_Q3_N480_L81012 <- readRDS("~/pineR/bammit/data/V3_Q3_N480_L81012.rds") # model 8

for (i in 1:3) {
  modelAux <- bammitJags(data = V3_Q2_N480_L1012,
                         Q = i,
                         mmu = 100,
                         smu = 10,
                         a = 0.1,
                         b = 0.1,
                         nthin = 1,
                         nburnin = 2000)
  model <- modelAux$BUGSoutput
  saveRDS(model, file = paste0("~/pineR/bammit/Running models/model5_", i, ".rds"))
}

# model2 <- list()
# predInt2 <- list()
# quantInt2 <-  vector("list", 3)
# for (i in 1:3) {
#   model2[[i]]  <- readRDS(paste0("~/pineR/bammit/Running models/model2_", i, ".rds"))
#   predInt2[[i]] <- model2[[i]] $mean$blin
#   aux <- model2[[i]]$sims.list$blin
#   quantInt2[[i]]$lower <- apply(aux, 2, function(x) quantile(x, 0.05))
#   quantInt2[[i]]$upper <- apply(aux, 2, function(x) quantile(x, 0.95))
# }
# names(model2) <- c("Q1", "Q2", "Q3")
# names(predInt2) <- c("Q1", "Q2", "Q3")
# names(quantInt2) <- c("Q1", "Q2", "Q3")
#
# model5 <- list()
# predInt5 <- list()
# quantInt5 <-  vector("list", 3)
# for (i in 1:3) {
#   model5[[i]] <- readRDS(paste0("~/pineR/bammit/Running models/model5_", i, ".rds"))
#   predInt5[[i]] <- model5[[i]]$mean$blin
#   aux <- model5[[i]]$sims.list$blin
#   quantInt5[[i]]$lower <- apply(aux, 2, function(x) quantile(x, 0.05))
#   quantInt5[[i]]$upper <- apply(aux, 2, function(x) quantile(x, 0.95))
# }
# names(model5) <- c("Q1", "Q2", "Q3")
# names(predInt5) <- c("Q1", "Q2", "Q3")
# names(quantInt5) <- c("Q1", "Q2", "Q3")
#
# model8 <- list()
# predInt8 <- list()
# quantInt8 <-  vector("list", 3)
# for (i in 1:3) {
#   model8[[i]] <- readRDS(paste0("~/pineR/bammit/Running models/model8_", i, ".rds"))
#   predInt8[[i]] <- model8[[i]]$mean$blin
#   aux <- model8[[i]]$sims.list$blin
#   quantInt8[[i]]$lower <- apply(aux, 2, function(x) quantile(x, 0.05))
#   quantInt8[[i]]$upper <- apply(aux, 2, function(x) quantile(x, 0.95))
# }
# names(model8) <- c("Q1", "Q2", "Q3")
# names(predInt8) <- c("Q1", "Q2", "Q3")
# names(quantInt8) <- c("Q1", "Q2", "Q3")




model<- list()
for (i in 1:3) {
  model[[i]] <- readRDS(paste0("~/pineR/bammit/Running models/model8_", i, ".rds"))
}

predInt <- list()
for (i in 1:3) {
  predInt[[i]] <- model[[i]]$median$blin
}
names(predInt) <- paste0("model", 1:3)


quantInt <- vector("list", 3)
for (i in 1:3) {
  aux <- model[[i]]$sims.list$blin
  quantInt[[i]]$lower <- apply(aux, 2, function(x) quantile(x, 0.05))
  quantInt[[i]]$upper <- apply(aux, 2, function(x) quantile(x, 0.95))
}

uppers <- lapply(quantInt, function(x) x[["upper"]]) %>% melt() %>% dplyr::select(value)
lowers <- lapply(quantInt, function(x) x[["lower"]]) %>% melt() %>% dplyr::select(value)
colnames(uppers) <- "upper"
colnames(lowers) <- "lower"


predInt <- reshape2::melt(predInt) %>% dplyr::select(value, L1)
colnames(predInt) <- c("int", "model")
#dataInt <- reshape2::melt(dataInt) %>% dplyr::select(value)
colnames(dataInt) <- "true"
Qs <- c(rep("Q1", 480),
        rep("Q2", 480),
        rep("Q3", 480))
predInt <- cbind(predInt, V3_Q3_N480_L81012$blin, uppers, lowers,  Qs)
colnames(predInt) <- c("int", "model", "true", "upper", "lower", "Q")


lab <- function(x) c("Q = 1", "Q = 2", "Q = 3")
predInt %>%
  ggplot(aes(x = true,
             y = int))+
  geom_point( size = 0.2) +
  geom_abline(size = 0.3, col = "gray") +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha=0.4, fill =  "steelblue") +
  facet_wrap(~Q,labeller = as_labeller(lab)) +
  labs(x = "True", y = "Estimated") +
  theme_bw(base_size = 14) +
  theme(strip.background = element_blank(),
        panel.spacing.x = unit(5,"mm"),
        legend.position = "bottom")


predInt %>%
  ggplot(aes(x = true, y = int)) +
  geom_abline(size = 0.3, col = "gray") +
  geom_linerange(aes(ymin =  lower, ymax = upper), alpha = 0.2, size = 0.1) +
  geom_point(colour =  "steelblue", size = 0.2) +
  facet_wrap(~Q,
             labeller = as_labeller(lab)) +
  labs(x = "True", y = "Estimated") +
  theme_bw(base_size = 14) +
  theme(strip.background = element_blank(),
        panel.spacing.x = unit(5,"mm"),
        legend.position = c(0.75, 0.04))



lab <- function(x) c("Q = 1", "Q = 2", "Q = 3")
predInt %>%
  ggplot(aes(x = true,
             y = int))+
  #geom_point(colour =  "steelblue", size = 0.2) +
  #geom_abline(size = 0.3, col = "gray") +
  geom_smooth(aes(ymin = lower, ymax = upper),
              stat = "identity", alpha = 0.4, size = 0.1, fill =  "steelblue", colour =  "steelblue") +
  facet_wrap(~Q,labeller = as_labeller(lab)) +
  labs(x = "True", y = "Estimated") +
  theme_bw(base_size = 14) +
  theme(strip.background = element_blank(),
        panel.spacing.x = unit(5,"mm"),
        legend.position = "bottom")


#predInt$Q <- factor(predInt$Q,labels = c("1", "2", "3"))
ggplot(predInt,
       aes(x = true,
           y = int,
           col = Q)) +
  geom_smooth(method = "loess", aes(fill = Q), alpha = 0.4, size = 0.4, se = FALSE)+
  geom_line(aes(x = true, y = lower, col = Q)) +
  #geom_smooth(aes(ymin = lower, ymax = upper,fill = Q, colour = Q),stat = "identity", alpha = 0.1, size = 0.1) +
  labs(x = "True", y = "Estimated")  +
  theme_bw(base_size = 14)



ggplot(predInt,
       aes(
         y = int,
         col = Q)) +
  geom_smooth(method = "lm", aes(fill = Q), alpha = 0.4, size = 0.4, se = FALSE)+
  geom_line(aes(x = true, y = lower, col = Q)) +
  #geom_line(aes(x = true, y = upper, col = Q)) +
  #geom_smooth(aes(ymin = lower, ymax = upper,fill = Q, colour = Q),stat = "identity", alpha = 0.1, size = 0.1) +
  labs(x = "True", y = "Estimated")  +
  theme_bw(base_size = 14)


test <- predInt[predInt$Q == "Q1", ]
test2 <- predInt[predInt$Q == "Q2", ]
test3 <- predInt[predInt$Q == "Q3", ]
indices = order(test$int)
indices2 = order(test2$int)
indices3 = order(test3$int)

df <- rbind(test[indices,], test2[indices2,], test3[indices3,])

df$id <- c(1:nrow(test), 1:nrow(test2), 1:nrow(test3))


df %>%
  ggplot(aes(x = id,
             y = int))+
  geom_line(size = 0.8, col = "steelblue") +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha=0.4, fill =  "steelblue") +
  geom_smooth(method = "loess", aes(x = id, y = true), alpha = 0.4, size = 0.8, se = FALSE, col = "firebrick")+
  geom_point(aes(x = id, y = true), size = 0.2)+
  facet_wrap(~Q,labeller = as_labeller(lab)) +
  labs(x = "Ordered indexes of interaction", y = "Estimated") +
  theme_bw(base_size = 14) +
  theme(strip.background = element_blank(),
        panel.spacing.x = unit(5,"mm"),
        legend.position = "bottom")





# ar-bammit ---------------------------------------------------------------


arbammitJagsRealData <- function(data, Q = 1, mmu = 10, smu = 2, stheta = 1, a = 0.1, b = 0.1,
                                 nthin = 2, nburnin = 3000){


  Q <- Q
  Y <- data$Yield
  B1 <- length(levels(data$Genotype))
  B2 <- length(levels(data$Environment))
  B3 <- length(levels(data$Year))
  B4 <- length(levels(data$Bloc))
  N <- length(data$Yield)

  genotype <- data$Genotype
  environment <- data$Environment
  time <- data$Year
  bloc <- data$Bloc

  modelCode <- "
  model
  {
  # Likelihood
   for (n in 1:N) {
    Y[n] ~ dnorm(mu[n], sy)
    mu[n] = muall + g[genotype[n]] + e[environment[n]] + t[time[n]] + bl[bloc[n]] + blin[n]
    blin[n] = sum(lambda[1:Q] * beta1[genotype[n],1:Q] * beta2[environment[n],1:Q] * beta3[time[n], 1:Q] * beta4[bloc[n], 1:Q])
   }

  # Priors
  # Prior on grand mean
   muall ~ dnorm(mmu, smu^-2)

  # Prior on genotype effect
  for(i in 1:B1) {
  g[i] ~ dnorm(0, sg^-2) # Prior on genotype effect
  }

  # Prior on environment effect
  for(j in 1:B2) {
  e[j] ~ dnorm(0, se^-2) # Prior on environment effect
  }

  # Prior on time effect
  tt[1] ~ dnorm(0,1)
  alpha ~ dunif(-1,1)
  for (k in 2:B3) {
    muk[k] = a0 + alpha * tt[k-1]
    tt[k] ~ dnorm(muk[k], st)
  }

  mut = sum(tt)/B3

  for(k in 1:B3){
    t[k] = tt[k] - mut
  }

   # Prior on bloc effect
  for(l in 1:B4) {
  bl[l] ~ dnorm(0, sb^-2) # Prior on time effect
  }

   # Priors on beta1
  for(q in 1:Q){
    for(i in 1:B1){
      thetaG[i,q] ~ dnorm(0,stheta)
    }
    mG[q] = sum(thetaG[1:B1,q])/B1
    for(i in 1:B1){
    thetaGNew[i,q] = thetaG[i,q] - mG[q]
    }
    sqrtThetaG[q] = sqrt(1/(sum(thetaGNew[1:B1,q]^2 + 0.000001)))
    for(i in 1:B1){
      beta1[i,q] = thetaGNew[i,q]*sqrtThetaG[q]
    }
  }

   # Priors on beta2
   for(q in 1:Q){
    for(j in 1:B2){
      thetaD[j,q] ~ dnorm(0,stheta)
    }
    mD[q] = sum(thetaD[1:B2,q])/B2
    for(j in 1:B2){
    thetaDNew[j,q] = thetaD[j,q] - mD[q]
    }
    sqrtThetaD[q] = sqrt(1/(sum(thetaDNew[1:B2,q]^2+ 0.000001)))
    for(j in 1:B2){
      beta2[j,q] = thetaDNew[j,q]*sqrtThetaD[q]
    }
  }

   # Priors on beta3
   for(q in 1:Q){
    thetaB3[1,q] ~ dnorm(0,stheta)
   }

   phiB3 ~ dunif(-1,1)
   for(q in 1:Q){
    for(k in 2:B3){
      thetaB3[k,q] ~ dnorm(b0 + phiB3 * thetaB3[k-1,q], stheta)
    }
    mB3[q] = sum(thetaB3[1:B3,q])/B3
    for(k in 1:B3){
    thetaB3New[k,q] = thetaB3[k,q] - mB3[q]
    }
    sqrtThetaB3[q] = sqrt(1/(sum(thetaB3New[1:B3,q]^2+ 0.000001)))
    for(k in 1:B3){
      beta3[k,q] = thetaB3New[k,q]*sqrtThetaB3[q]
    }
  }


   # Priors on beta4
   for(q in 1:Q){
    for(l in 1:B4){
      thetaB4[l,q] ~ dnorm(0,stheta)
    }
    mB4[q] = sum(thetaB4[1:B4,q])/B4
    for(l in 1:B4){
    thetaB4New[l,q] = thetaB4[l,q] - mB4[q]
    }
    sqrtThetaB4[q] = sqrt(1/(sum(thetaB4New[1:B4,q]^2+ 0.000001)))
    for(l in 1:B4){
      beta4[l,q] = thetaB4New[l,q]*sqrtThetaB4[q]
    }
  }

  # Prior on eigenvalues
  for(q in 1:Q) {
    lambda_raw[q] ~ dnorm(0, slambda^-2)T(0,)
  }
  lambda = sort(lambda_raw)

  #Priors
  sg ~ dt(0, 1, 1)T(0,)
  se ~ dt(0, 1, 1)T(0,)
  st ~ dt(0, 1, 1)T(0,)
  sb ~ dt(0, 1, 1)T(0,)
  slambda ~ dt(0, 1, 1)T(0,)
  a0 ~ dnorm(0, 10^-2)
  b0 ~ dnorm(0, 10^-2)

  # Prior on residual standard deviation
   sy ~ dgamma(a, b) # inverse of sy
  }
  "

  # Set up the data
  modelData <- list(
    N = N,
    Y = Y,
    B1 = B1,
    B2 = B2,
    B3 = B3,
    B4 = B4,
    Q = Q,
    genotype = genotype,
    environment = environment,
    time = time,
    bloc = bloc,
    mmu = mmu,
    smu = smu,
    stheta = stheta,
    a = a,
    b = b
  )

  # Choose the parameters to watch
  modelParameters <- c(
    "g", "e", "t", "bl", "lambda", "beta1", "beta2", "beta3", "beta4","sy", "muall", "blin"
  )

  # Run the model
  modelRun <- jags(
    data = modelData,
    parameters.to.save = modelParameters,
    model.file = textConnection(modelCode),
    #progress.bar = "none",
    n.thin = nthin,
    n.burnin = nburnin,
    n.iter = 2000 * nthin + nburnin
  )

  return(modelRun)
}

predictionBAMMITReal <- function(model, data) {

  genotype <- data$Genotype
  environment <- data$Environment
  time <- data$Year
  bloc <- data$Bloc

  muhat <- model$mean$muall
  g <- model$mean$g
  e <- model$mean$e
  t <- model$mean$t
  bl <- model$mean$bl
  inthat <- model$mean$blin

  N <- length(data$Yield)
  yhat <- rep(muhat, N) + g[genotype] + e[environment] + t[time] + bl[bloc]  +  inthat

  return(yhat)
}

load("~/pineR/bammit/data/train_ireland.RData")
load("~/pineR/bammit/data/test_ireland.RData")
train
data <- train

data <- train %>% dplyr::filter(Year %in% c(2010,2019))

trainQ1ar <- arbammitJagsRealData(data = train, Q = 1)
plot(trainQ1ar)
saveRDS(trainQ1ar$BUGSoutput, file = "~/pineR/bammit/Running models/trainQ3ar.rds")

trainQ1ar$BUGSoutput$median$blin
predRar <- predictionBAMMITReal(trainQ1ar$BUGSoutput, test)
round(caret::RMSE(predRar, test$Yield),2)
round(caret::R2(predRar, test$Yield),2)


train[train$Year == "2012",] %>%
  dplyr::group_by(Genotype,Environment) %>%
  dplyr::summarize(mean = mean(Yield))

dat <- data %>%
  dplyr::group_by(Year) %>%
  dplyr::summarize(mean = mean(Yield))
plot(as.numeric(dat$Year), dat$mean, type = 'l')
abline(reg=lm(dat$mean~time(dat$mean)))
acf(dat$mean)

# Reading the train data --------------------------------------------------------

base_dir <- '~/pineR/bammit/data/'
all_files <- list.files(base_dir, pattern="*.rds", full.names = TRUE)
train <- lapply(all_files, readRDS)
names(train) <- sub(pattern = "(.*)\\..*$", replacement = "\\1", basename(all_files))

# Running models for simulated data ---------------------------------------


runModel <- function(data, Q, name, base_dir = '~/pineR/bammit/sim run/'){

  V <- length(data$Bv)
  N <- Reduce('*', data$Bv)

  modelAux <- bammitJags(data = data,
                         Q = Q,
                         mmu = 100,
                         smu = 10,
                         a = 0.1,
                         b = 0.1,
                         nthin = 2,
                         nburnin = 2000)

  model <- modelAux$BUGSoutput
  fileName <- paste0(base_dir, 'model', name,'QM' , Q, '.rds')
  saveRDS(model, file = fileName)

}

Q = c(1,2,3)
for (q in Q) {
  for (i in 1:length(train)) {
    runModel(train[[i]], Q = q, name = names(train)[i])
  }
}





