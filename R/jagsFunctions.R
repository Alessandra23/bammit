#' Jags code Jossi's model
#' @description A Jags code to estimated the parameters of AMMI model followuing Josse's approach
#' @param data A list containing the simulated values of y, the values of I, J, and Q and a mtrix od g and e.
#' @param mmu Mean of the grand mean
#' @param smu sd of grand mean
#' @param sg sd of genotypes
#' @param se sd of environments
#' @param slambda sd of $\lambda$
#' @param a shape parameter of a Gamma distribuion
#' @param b scale parameter of a Gamma distribuion
#' @param nthin thinning rate
#' @param nburnin length of burn in
#' @return The output of Jags function
#' @export
#' @importFrom R2jags 'jags'
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

#' @export
crossaJags <- function(data, mmu, smu, mug, sg, mue, se, mulambda, slambda, stheta, a, b, nthin = 1, nburnin = 2000){

  # from data
  Y <- data$y
  I <- data$I
  J <- data$J
  Q <- data$Q
  N <- I*J
  genotype <- data$x[, "g"]
  environment <- data$x[, "e"]

  modelCode <- "
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
  g[i] ~ dnorm(mug, sg^-2) # Prior on genotype effect
  }

  # Prior on environment effect
  for(j in 1:J) {
  e[j] ~ dnorm(mue, se^-2) # Prior on environment effect
  }


  # Priors on gamma
  for(q in 1:Q){
    for(i in 1:(I-1)){
      theta[i,q] ~ dnorm(0,stheta)
    }
    theta[I,q] = -sum(theta[1:(I-1),q])
    thetaSum[q] = sqrt(sum(theta[1:I,q]^2)) + 0.00001
    for(i in 1:I){
      gamma[i,q] = theta[i,q]/thetaSum[q]
    }
  }

  # Priors on delta
   for(q in 1:Q){
    for(j in 1:(J-1)){
      aux[j,q] ~ dnorm(0,stheta)
    }
    aux[J,q] = -sum(aux[1:(J-1),q])
    auxSum[q] = sqrt(sum(aux[1:J,q]^2)) + 0.000001
    for(j in 1:J){
      delta[j,q] = aux[j,q]/auxSum[q]
    }
  }

  # Prior on eigenvalues
  for(q in 1:Q) {
    lambda_raw[q] ~ dnorm(mulambda, slambda^-2)T(0,)
  }
  lambda = sort(lambda_raw)

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
    Q = Q,
    genotype = genotype,
    environment = environment,
    mmu = mmu,
    smu = smu,
    mug = mug,
    mue = mue,
    mulambda = mulambda,
    sg = sg,
    se = se,
    slambda = slambda,
    stheta = stheta,
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


bammitJags <- function(data, mmu, smu, mug, sg, mue, se, mut, st, mulambda, slambda, stheta, a, b, nthin = 1, nburnin = 2000){

  # from data
  Y <- data$y
  I <- data$I
  J <- data$J
  K <- data$K
  Q <- data$Q
  N <- I*J*K
  genotype <- data$x[, "g"]
  environment <- data$x[, "e"]
  time <- data$x[, "t"]

  modelCode <- "
  model
  {
  # Likelihood
   for (n in 1:N) {
    Y[n] ~ dnorm(mu[n], sy^-2)
    mu[n] = muall + g[genotype[n]] + e[environment[n]] + t[time[n]] + blin[n]
    blin[n] = sum(lambda[1:Q] * gamma[genotype[n],1:Q] * delta[environment[n],1:Q]*kappa[time[n],1:Q])
   }

  # Priors
  # Prior on grand mean
   muall ~ dnorm(mmu, smu^-2)

  # Prior on genotype effect
  for(i in 1:I) {
  g[i] ~ dnorm(mug, sg^-2) # Prior on genotype effect
  }

  # Prior on environment effect
  for(j in 1:J) {
  e[j] ~ dnorm(mue, se^-2) # Prior on environment effect
  }

  # Prior on time effect
  for(k in 1:K) {
  t[k] ~ dnorm(mut, st^-2) # Prior on time effect
  }


  # Priors on gamma
  for(q in 1:Q){
    for(i in 1:(I-1)){
      theta[i,q] ~ dnorm(0,stheta)
    }
    theta[I,q] = -sum(theta[1:(I-1),q])
    thetaSum[q] = sqrt(sum(theta[1:I,q]^2)) + 0.00001
    for(i in 1:I){
      gamma[i,q] = theta[i,q]/thetaSum[q]
    }
  }

  # Priors on delta
   for(q in 1:Q){
    for(j in 1:(J-1)){
      aux[j,q] ~ dnorm(0,stheta)
    }
    aux[J,q] = -sum(aux[1:(J-1),q])
    auxSum[q] = sqrt(sum(aux[1:J,q]^2)) + 0.000001
    for(j in 1:J){
      delta[j,q] = aux[j,q]/auxSum[q]
    }
  }

   # Priors on kappa
   for(q in 1:Q){
    for(k in 1:(K-1)){
      thetaK[k,q] ~ dnorm(0,stheta)
    }
    thetaK[K,q] = -sum(thetaK[1:(K-1),q])
    thetaKSum[q] = sqrt(sum(thetaK[1:K,q]^2)) + 0.000001
    for(k in 1:K){
      kappa[k,q] = thetaK[k,q]/thetaKSum[q]
    }
  }

  # Prior on eigenvalues
  for(q in 1:Q) {
    lambda_raw[q] ~ dnorm(mulambda, slambda^-2)T(0,)
  }
  lambda = sort(lambda_raw)

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
    mug = mug,
    mue = mue,
    mut = mut,
    mulambda = mulambda,
    sg = sg,
    se = se,
    st = st,
    slambda = slambda,
    stheta = stheta,
    a = a,
    b = b
  )

  # Choose the parameters to watch
  modelParameters <- c(
    "g", "e", "t", "lambda", "gamma", "delta", "kappa", "sy", "muall", "blin"
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
