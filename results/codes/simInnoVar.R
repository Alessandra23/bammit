library(R2jags)
library(ggplot2)
bammitJagssimInnoVar <- function(data, Q = 1, mmu, smu, stheta, a, b,
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
  # soil <- data[, "Year"]

  # levels(data$genotype) <-(1:length(levels(data$genotype)))
  # levels(data$environment) <- 1:length(levels(data$environment))

  Q <- Q
  Y <- data$y
  K <- length(levels(data$soil))
  I <- length(levels(data$genotype))
  J <- length(levels(data$environment))
  N <- length(data$y)

  genotype <- data[, "genotype"]
  environment <- data[, "environment"]
  soil <- data[, "soil"]

  modelCode <- "
  model
  {
  # Likelihood
   for (n in 1:N) {
    Y[n] ~ dnorm(mu[n], sy)
    mu[n] = muall + g[genotype[n]] + e[environment[n]] + t[soil[n]] + blin[n]
    blin[n] = sum(lambda[1:Q] * gamma[genotype[n],1:Q] * delta[environment[n],1:Q]*rho[soil[n],1:Q])
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

  # Prior on soil effect
  for(k in 1:K) {
  t[k] ~ dnorm(0, st^-2) # Prior on soil effect
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
    soil = soil,
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

simInnoVar <- readRDS("~/Documents/GitHub/bammit/Real data/simInnoVar.RDS")
simInnoVar <- simInnoVar |> select(target, genotype, environment, DranaigeClass)
colnames(simInnoVar) <- c("y", "genotype", "environment", "soil")

model <- bammitJagssimInnoVar(data = simInnoVar,  mmu = 5, smu = 1, stheta = 1/100,
                              a = 0.1, b = 0.1, nthin = 1, nburnin = 100)



dfBlin <- as.data.frame(model$BUGSoutput$sims.list$blin)
names(dfBlin) <- paste0("term", seq(1:ncol(dfBlin)))
dfBlin <- reshape2::melt(dfBlin[,1:8])

p <- ggplot(dfBlin, aes(x = value)) +
  geom_density(colour = "steelblue", fill = "steelblue", alpha = 0.4) +
  facet_wrap(~variable) +
  labs(x = "Genotype", y = "Frequency") +
  theme_bw(base_size = 14) +
  theme(strip.background = element_rect(fill="white"),
        panel.spacing.x = unit(5,"mm"),
        legend.position = c(0.75, 0.04))
p
