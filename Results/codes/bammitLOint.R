# simulate the model with just the lower-interaction

library(R2jags)
library(ggplot2)
library(patchwork)

rm(list = ls())

# Jags code ---------------------------------------------------------------

# Model 1 - Interaction just between beta1 and beta 2

model_jags_1 <- function(data, Q, mmu, smu, a,b, stheta = 1, nthin, nburnin){

  Y <- data$y
  V <- length(data$Bv)
  Bv <- data$Bv
  N <- length(data$y)


  modelCode <- "
    model{
      # Likelihood
      for (n in 1:N) {
          Y[n] ~ dnorm(mu[n], sy)
          mu[n] = muall + b1[var1[n]] + b2[var2[n]] + b3[var3[n]] + int[n]
          int[n] = sum(lambda[1:Q] * beta1[var1[n],1:Q] * beta2[var2[n],1:Q])
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
            thetaBeta1[i,q] ~ dnorm(0, 1/stheta^2)
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
            thetaBeta2[i,q] ~ dnorm(0, 1/stheta^2)
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
      sb3 ~ dt(0, 1, 1)T(0,)
      slambda ~ dt(0, 1, 1)T(0,)

      sy ~ dgamma(a, b) # Prior on residual standard deviation - inverse of sy

    }
    "

  B1 <- Bv[1]
  B2 <- Bv[2]
  B3 <- Bv[3]

  x <- expand.grid(lapply(Bv, function(x) 1:x))
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
    "b1", "b2", "b3", "lambda", "beta1", "beta2", "sy", "muall",
    "int", "mu"
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

# Model 2 - Two parts: interaction between beta1 and beta 2  + between beta1, beta 2 and beta3

model_jags_2 <- function(data, Q, mmu, smu, a,b, stheta = 1, nthin, nburnin){

  Y <- data$y
  V <- length(data$Bv)
  Bv <- data$Bv
  N <- length(data$y)


  modelCode <- "
    model{
      # Likelihood
      for (n in 1:N) {
          Y[n] ~ dnorm(mu[n], sy)
          mu[n] = muall + b1[var1[n]] + b2[var2[n]] + b3[var3[n]] + int_p1[n] + int_p2[n]
          int_p1[n] = sum(lambda_1[1:Q] * beta1[var1[n],1:Q] * beta2[var2[n],1:Q])
          int_p2[n] = sum(lambda_2[1:Q] * beta1[var1[n],1:Q] * beta2[var2[n],1:Q]* beta3[var3[n],1:Q])
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
            thetaBeta1[i,q] ~ dnorm(0, 1/stheta^2)
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
            thetaBeta2[i,q] ~ dnorm(0, 1/stheta^2)
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
        lambda_raw_1[q] ~ dnorm(0, slambda^-2)T(0,)
      }
      lambda_1 = sort(lambda_raw_1)


      for(q in 1:Q) {
        lambda_raw_2[q] ~ dnorm(0, slambda^-2)T(0,)
      }
      lambda_2 = sort(lambda_raw_2)

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

  x <- expand.grid(lapply(Bv, function(x) 1:x))
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
    "b1", "b2", "b3", "lambda_1", 'lambda_2', "beta1", "beta2", 'beta3', "sy", "muall",
    "mu", 'int_p1', 'int_p2'
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

# Model 3 - original model, just the full interaction

model_jags_3 <- function(data, Q, mmu, smu, a,b, stheta = 1, nthin, nburnin){

  Y <- data$y
  V <- length(data$Bv)
  Bv <- data$Bv
  N <- length(data$y)


  modelCode <- "
    model{
      # Likelihood
      for (n in 1:N) {
          Y[n] ~ dnorm(mu[n], sy)
          mu[n] = muall + b1[var1[n]] + b2[var2[n]] + b3[var3[n]] + int[n]
          int[n] = sum(lambda[1:Q] * beta1[var1[n],1:Q] * beta2[var2[n],1:Q]* beta3[var3[n],1:Q])
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
            thetaBeta1[i,q] ~ dnorm(0, 1/stheta^2)
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
            thetaBeta2[i,q] ~ dnorm(0, 1/stheta^2)
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

  x <- expand.grid(lapply(Bv, function(x) 1:x))
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
    "b1", "b2", "b3", "lambda", "beta1", "beta2", 'beta3', "sy", "muall",
    "int", "mu"
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

# Model 4 - ammi models

model_jags_4 <- function(data, Q, mmu, smu, a,b, stheta = 1, nthin, nburnin){

  Y <- data$y
  V <- length(data$Bv)
  Bv <- data$Bv
  N <- length(data$y)


  modelCode <- "
    model{
      # Likelihood
      for (n in 1:N) {
          Y[n] ~ dnorm(mu[n], sy)
          mu[n] = muall + b1[var1[n]] + b2[var2[n]] + int[n]
          int[n] = sum(lambda[1:Q] * beta1[var1[n],1:Q] * beta2[var2[n],1:Q])
      }

      # Priors

      muall ~ dnorm(mmu, smu^-2) # grand mean

      for(i in 1:B1) {
          b1[i] ~ dnorm(0, sb1^-2) # Prior on effect 1
      }

      for(i in 1:B2) {
          b2[i] ~ dnorm(0, sb2^-2) # Prior on effect 2
      }

      # for(i in 1:B3) {
      #     b3[i] ~ dnorm(0, sb3^-2) # Prior on effect 3
      # }


      for(q in 1:Q){
          for(i in 1:B1){
            thetaBeta1[i,q] ~ dnorm(0, 1/stheta^2)
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
            thetaBeta2[i,q] ~ dnorm(0, 1/stheta^2)
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
      #sb3 ~ dt(0, 1, 1)T(0,)
      slambda ~ dt(0, 1, 1)T(0,)

      sy ~ dgamma(a, b) # Prior on residual standard deviation - inverse of sy

    }
    "

  B1 <- Bv[1]
  B2 <- Bv[2]

  x <- expand.grid(lapply(Bv, function(x) 1:x))
  colnames(x) <- paste0('var', 1:2)
  var1 <- x$var1
  var2 <- x$var2
  #var3 <- x$var3

  # Set up the data
  modelData <- list(
    N = N,
    Y = Y,
    B1 = B1,
    B2 = B2,
    #B3 = B3,
    Q = Q,
    var1 = var1,
    var2 = var2,
    #var3 = var3,
    mmu = mmu,
    smu = smu,
    stheta = stheta,
    a = a,
    b = b
  )

  # Choose the parameters to watch
  modelParameters <- c(
    "b1", "b2",  "lambda", "beta1", "beta2", "sy", "muall",
    "int", "mu"
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

# Model 5 - applying a transformation to sum 1

model_jags_5 <- function(data, Q, mmu, smu, a,b, nthin, nburnin){

  Y <- data$y
  V <- length(data$Bv)
  Bv <- data$Bv
  N <- length(data$y)


  modelCode <- "
    model{
      # Likelihood
      for (n in 1:N) {
          Y[n] ~ dnorm(mu[n], sy)
          mu[n] = muall + b1[var1[n]] + b2[var2[n]] + b3[var3[n]] + int[n]
          int[n] = sum(lambda[1:Q] * (beta1[var1[n],1:Q]) * (beta2[var2[n],1:Q])*(beta3[var3[n],1:Q]))
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
            thetaBeta1[i,q] ~ dnorm(0, sbeta1^-2)
          }
          mBeta1[q] = sum(thetaBeta1[1:B1,q])
          for(i in 1:B1){
            beta1[i,q] = thetaBeta1[i,q]/(mBeta1[q] + 0.0000001)
          }
      }

       for(q in 1:Q){
          for(i in 1:B2){
            thetaBeta2[i,q] ~ dnorm(0, sbeta2^-2)
          }
          mBeta2[q] = sum(thetaBeta2[1:B2,q])
          for(i in 1:B2){
            beta2[i,q] = thetaBeta2[i,q]/(mBeta2[q] + 0.0000001)
          }
       }

       for(q in 1:Q){
          for(i in 1:B3){
            thetaBeta3[i,q] ~ dnorm(0, sbeta3^-2)
          }
          mBeta3[q] = sum(thetaBeta3[1:B3,q])
          for(i in 1:B3){
            beta3[i,q] = thetaBeta3[i,q]/(mBeta3[q] + 0.0000001)
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

      sbeta1 ~ dt(0, 1, 1)T(0,)
      sbeta2 ~ dt(0, 1, 1)T(0,)
      sbeta3 ~ dt(0, 1, 1)T(0,)

      slambda ~ dt(0, 1, 1)T(0,)

      sy ~ dgamma(a, b) # Prior on residual standard deviation - inverse of sy

    }
    "

  B1 <- Bv[1]
  B2 <- Bv[2]
  B3 <- Bv[3]

  x <- expand.grid(lapply(Bv, function(x) 1:x))
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
    a = a,
    b = b
  )

  # Choose the parameters to watch
  modelParameters <- c(
    "b1", "b2", "b3", "lambda", "beta1", "beta2", 'beta3', "sy", "muall",
    "int", "mu"
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

# Model 7 - Andrew's code

model_jags_7 <- function(data, Q, mmu, smu, a,b, nthin, nburnin, alpha_beta){

  Y <- data$y
  V <- length(data$Bv)
  Bv <- data$Bv
  N <- length(data$y)



  modelCode <- "
model{
  # Likelihood
  for (n in 1:N) {
    Y[n] ~ dnorm(mu[n], sy^-2)
    mu[n] = muall + b1[var1[n]] + b2[var2[n]] + b3[var3[n]] + int[n]
    int[n] = sum(lambda[1:Q] * beta1[var1[n],1:Q] * beta2[var2[n],1:Q] * beta3[var3[n],1:Q])
  }

  # Priors
  muall ~ dnorm(mmu, smu^-2) # grand mean

  # Main effects
  for(i in 1:B1) {
    b1[i] ~ dnorm(0, sb1^-2) # Prior on effect 1
  }
  for(i in 1:B2) {
    b2[i] ~ dnorm(0, sb2^-2) # Prior on effect 2
  }
  for(i in 1:B3) {
    b3[i] ~ dnorm(0, sb3^-2) # Prior on effect 2
  }

  # Interaction beta 1
  for(q in 1:Q) {
    # generate values of theta from a Normal distribution
    for(i in 1:B1) {
      thetaBeta1[i,q] ~ dnorm(0, 1^-2)
    }
    # calculate the mean of each column
    mBeta1[q] = sum(thetaBeta1[1:B1,q])/B1

    # centering
    for(i in 1:B1) {
      thetaBeta1New[i,q] = thetaBeta1[i,q] - mBeta1[q]
    }

    # get sd of each column
    sqrtThetaBeta1[q] = sqrt(sum(thetaBeta1New[1:B1,q]^2)) + 0.00000001

    # get the final beta
    for(i in 1:B1) {
      beta1[i,q] = (thetaBeta1New[i,q] + 10000*M1[q]) / (sqrtThetaBeta1[q] + 10000*M1[q])
    }
    M1[q] ~ dbern(p1[q]) # M1 = 0 means normal, M1 = 1 means all beta1 = 1
    p1[q] ~ dbeta(1, 10)
  }

  # Interaction beta 2
  for(q in 1:Q) {
    for(i in 1:B2) {
      thetaBeta2[i,q] ~ dnorm(0, 1^-2)
    }
    mBeta2[q] = sum(thetaBeta2[1:B2,q])/B2

    for(i in 1:B2) {
      thetaBeta2New[i,q] = thetaBeta2[i,q] - mBeta2[q]
    }
    sqrtThetaBeta2[q] = sqrt(sum(thetaBeta2New[1:B2,q]^2)) + 0.00000001
    for(i in 1:B2) {
      beta2[i,q] = (thetaBeta2New[i,q] + 10000*M2[q]) / (sqrtThetaBeta2[q] + 10000*M2[q])
    }
    M2[q] ~ dbern(p2[q])
    p2[q] ~ dbeta(1, 10)
  }

    # Interaction beta 3
  for(q in 1:Q) {
    for(i in 1:B3) {
      thetaBeta3[i,q] ~ dnorm(0, 1^-2)
    }
    mBeta3[q] = sum(thetaBeta3[1:B3,q])/B3

    for(i in 1:B3) {
      thetaBeta3New[i,q] = thetaBeta3[i,q] - mBeta3[q]
    }
    sqrtThetaBeta3[q] = sqrt(sum(thetaBeta3New[1:B3,q]^2)) + 0.00000001
    for(i in 1:B3) {
      beta3[i,q] = (thetaBeta3New[i,q] + 10000*M3[q]) / (sqrtThetaBeta3[q] + 10000*M3[q])
    }
    M3[q] ~ dbern(p3[q])
    p3[q] ~ dbeta(1, 10)
  }

  # Prior on lambda
  for(q in 1:Q) {
    lambda_raw[q] ~ dnorm(0, slambda^-2)T(0,)
  }
  lambda = sort(lambda_raw)

  # Priors
  sb1 ~ dt(0, 1^-2, 1)T(0,)
  sb2 ~ dt(0, 1^-2, 1)T(0,)
  sb3 ~ dt(0, 1^-2, 1)T(0,)
  slambda ~ dt(0, 1^-2, 1)T(0,)
  sy ~ dgamma(a, b) # Prior on residual sd
}"

  B1 <- Bv[1]
  B2 <- Bv[2]
  B3 <- Bv[3]

  x <- expand.grid(lapply(Bv, function(x) 1:x))
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
    a = a,
    b = b
  )

  # Choose the parameters to watch
  modelParameters <- c(
    "b1", "b2", "b3","lambda", "beta1", "beta2", "beta3", "sy", "muall",
    "int", "mu", "p1", "p2", "p3"
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

# Model 8 - Andrew's code (4 variables)

model_jags_8 <- function(data, Q, mmu, smu, a,b, nthin, nburnin, alpha_beta){

  Y <- data$y
  V <- length(data$Bv)
  Bv <- data$Bv
  N <- length(data$y)



  modelCode <- "
model{
  # Likelihood
  for (n in 1:N) {
    Y[n] ~ dnorm(mu[n], sy^-2)
    mu[n] = muall + b1[var1[n]] + b2[var2[n]] + b3[var3[n]] + b4[var4[n]] + int[n]
    int[n] = sum(lambda[1:Q] * beta1[var1[n],1:Q] * beta2[var2[n],1:Q] * beta3[var3[n],1:Q]* beta4[var4[n],1:Q])
  }

  # Priors
  muall ~ dnorm(mmu, smu^-2) # grand mean

  # Main effects
  for(i in 1:B1) {
    b1[i] ~ dnorm(0, sb1^-2) # Prior on effect 1
  }
  for(i in 1:B2) {
    b2[i] ~ dnorm(0, sb2^-2) # Prior on effect 2
  }
  for(i in 1:B3) {
    b3[i] ~ dnorm(0, sb3^-2) # Prior on effect 2
  }

  for(i in 1:B4) {
    b4[i] ~ dnorm(0, sb4^-2) # Prior on effect 2
  }

  # Interaction beta 1
  for(q in 1:Q) {
    # generate values of theta from a Normal distribution
    for(i in 1:B1) {
      thetaBeta1[i,q] ~ dnorm(0, 1^-2)
    }
    # calculate the mean of each column
    mBeta1[q] = sum(thetaBeta1[1:B1,q])/B1

    # centering
    for(i in 1:B1) {
      thetaBeta1New[i,q] = thetaBeta1[i,q] - mBeta1[q]
    }

    # get sd of each column
    sqrtThetaBeta1[q] = sqrt(sum(thetaBeta1New[1:B1,q]^2)) + 0.00000001

    # get the final beta
    for(i in 1:B1) {
      beta1[i,q] = (thetaBeta1New[i,q] + 10000*M1[q]) / (sqrtThetaBeta1[q] + 10000*M1[q])
    }
    M1[q] ~ dbern(p1[q]) # M1 = 0 means normal, M1 = 1 means all beta1 = 1
    p1[q] ~ dbeta(1, 10)
  }

  # Interaction beta 2
  for(q in 1:Q) {
    for(i in 1:B2) {
      thetaBeta2[i,q] ~ dnorm(0, 1^-2)
    }
    mBeta2[q] = sum(thetaBeta2[1:B2,q])/B2

    for(i in 1:B2) {
      thetaBeta2New[i,q] = thetaBeta2[i,q] - mBeta2[q]
    }
    sqrtThetaBeta2[q] = sqrt(sum(thetaBeta2New[1:B2,q]^2)) + 0.00000001
    for(i in 1:B2) {
      beta2[i,q] = (thetaBeta2New[i,q] + 10000*M2[q]) / (sqrtThetaBeta2[q] + 10000*M2[q])
    }
    M2[q] ~ dbern(p2[q])
    p2[q] ~ dbeta(1, 10)
  }

    # Interaction beta 3
  for(q in 1:Q) {
    for(i in 1:B3) {
      thetaBeta3[i,q] ~ dnorm(0, 1^-2)
    }
    mBeta3[q] = sum(thetaBeta3[1:B3,q])/B3

    for(i in 1:B3) {
      thetaBeta3New[i,q] = thetaBeta3[i,q] - mBeta3[q]
    }
    sqrtThetaBeta3[q] = sqrt(sum(thetaBeta3New[1:B3,q]^2)) + 0.00000001
    for(i in 1:B3) {
      beta3[i,q] = (thetaBeta3New[i,q] + 10000*M3[q]) / (sqrtThetaBeta3[q] + 10000*M3[q])
    }
    M3[q] ~ dbern(p3[q])
    p3[q] ~ dbeta(1, 10)
  }

   # Interaction beta 4
  for(q in 1:Q) {
    for(i in 1:B4) {
      thetaBeta4[i,q] ~ dnorm(0, 1^-2)
    }
    mBeta4[q] = sum(thetaBeta4[1:B4,q])/B4

    for(i in 1:B4) {
      thetaBeta4New[i,q] = thetaBeta4[i,q] - mBeta4[q]
    }
    sqrtThetaBeta4[q] = sqrt(sum(thetaBeta4New[1:B4,q]^2)) + 0.00000001
    for(i in 1:B4) {
      beta4[i,q] = (thetaBeta4New[i,q] + 10000*M4[q]) / (sqrtThetaBeta4[q] + 10000*M4[q])
    }
    M4[q] ~ dbern(p4[q])
    p4[q] ~ dbeta(1, 10)
  }

  # Prior on lambda
  for(q in 1:Q) {
    lambda_raw[q] ~ dnorm(0, slambda^-2)T(0,)
  }
  lambda = sort(lambda_raw)

  # Priors
  sb1 ~ dt(0, 1^-2, 1)T(0,)
  sb2 ~ dt(0, 1^-2, 1)T(0,)
  sb3 ~ dt(0, 1^-2, 1)T(0,)
  sb4 ~ dt(0, 1^-2, 1)T(0,)
  slambda ~ dt(0, 1^-2, 1)T(0,)
  sy ~ dgamma(a, b) # Prior on residual sd
}"

  B1 <- Bv[1]
  B2 <- Bv[2]
  B3 <- Bv[3]
  B4 <- Bv[4]

  x <- expand.grid(lapply(Bv, function(x) 1:x))
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
    a = a,
    b = b
  )

  # Choose the parameters to watch
  modelParameters <- c(
    "b1", "b2", "b3","lambda", "beta1", "beta2", "beta3", "sy", "muall",
    "int", "mu", "p1", "p2", "p3", 'p4'
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

# Model 9 - (4 variables)
model_jags_9 <- function(data, Q, mmu, smu, a,b, stheta = 1, nthin, nburnin){

  Y <- data$y
  V <- length(data$Bv)
  Bv <- data$Bv
  N <- length(data$y)


  modelCode <- "
    model{
      # Likelihood
      for (n in 1:N) {
          Y[n] ~ dnorm(mu[n], sy)
          mu[n] = muall + b1[var1[n]] + b2[var2[n]] + b3[var3[n]] + b4[var4[n]] + int[n]
          int[n] = sum(lambda[1:Q] * beta1[var1[n],1:Q] * beta2[var2[n],1:Q] * beta3[var3[n],1:Q])
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
          b4[i] ~ dnorm(0, sb4^-2) # Prior on effect 3
      }


      for(q in 1:Q){
          for(i in 1:B1){
            thetaBeta1[i,q] ~ dnorm(0, 1/stheta^2)
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
            thetaBeta2[i,q] ~ dnorm(0, 1/stheta^2)
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
            thetaBeta3[i,q] ~ dnorm(0, 1/stheta^2)
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
      sb4 ~ dt(0, 1, 1)T(0,)
      slambda ~ dt(0, 1, 1)T(0,)

      sy ~ dgamma(a, b) # Prior on residual standard deviation - inverse of sy

    }
    "

  B1 <- Bv[1]
  B2 <- Bv[2]
  B3 <- Bv[3]
  B4 <- Bv[4]

  x <- expand.grid(lapply(Bv, function(x) 1:x))
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
    "b1", "b2", "b3", 'b4', "lambda", "beta1", "beta2", "sy", "muall",
    "int", "mu"
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


# Function to simulate data -----------------------------------------------------------

# Generate betas
genInt <- function(index = 10, Q = 1, stheta = 1) {
  # Generate theta matrix
  theta <- matrix(rnorm(index * Q, 0, stheta), nrow = index, ncol = Q)
  # Calculate the column means of theta
  m <- colMeans(theta)
  # Center the columns of theta
  thetaN <- sweep(theta, 2, m)
  # Calculate sqrtTheta
  sqrtTheta <- sapply(1:Q, function(q) {
    sqrt(1 / sum(thetaN[, q]^2))
  })
  # Calculate variable
  variable <- sweep(thetaN, 2, sqrtTheta, "*")
  return(variable)
}

# Function to generate data from Model 1

genData_1 <- function(B1, B2, B3, sigma, sb, lambda){ #sb1, sb2, sb3, sb4){

  Bv <- c(B1, B2, B3)
  #Bv <- c(B1, B2)
  N <- Reduce("*", Bv)
  Q <- length(lambda)
  V <- length(Bv)
  mu <- 10

  x <- expand.grid(1:Bv[1], 1:Bv[2], 1:Bv[3]) |>
    dplyr::mutate_if(is.numeric,as.factor)

  colnames(x) <- c('var1', 'var2', 'var3')

  # x <- expand.grid(1:Bv[1], 1:Bv[2]) |>
  #   dplyr::mutate_if(is.numeric,as.factor)
  #
  # colnames(x) <- c('var1', 'var2')


  b1 <- rnorm(B1, 0, sb)
  b1 <- b1 - mean(b1)

  b2 <- rnorm(B2, 0, sb)
  b2 <- b2 - mean(b2)

  b3 <- rnorm(B3, 0, sb)
  b3 <- b3 - mean(b3)


  beta1 <- genInt(B1, Q)
  beta2 <- genInt(B2, Q)
  beta3 <- genInt(B3, Q)

  blin = rep(0, N)
  for (k in 1:Q) {
    blin <- blin + lambda[k]*beta1[x[,'var1'],k]*beta2[x[,'var2'],k]
    #print(blin)
  }

  m <-  mu + b1[x[,'var1']] + b2[x[,'var2']] + b3[x[,'var3']] + blin
  y <-  rnorm(N, m, sigma)

  return(list(mu = mu,
              b1 = b1,
              b2 = b2,
              b3 = b3,
              #b4 = b4,
              blin = blin,
              y = y,
              betas = list(beta1 = beta1, beta2 = beta2, beta3 = beta3),#, beta4 = beta4),
              Bv = Bv))

}

# Function to generate data from Model 2

genData_2 <- function(B1, B2, B3, sigma, sb, lambda){
  Bv <- c(B1, B2, B3)
  #Bv <- c(B1, B2)
  N <- Reduce("*", Bv)
  Q <- length(lambda)
  V <- length(Bv)
  mu <- 10

  x <- expand.grid(1:Bv[1], 1:Bv[2], 1:Bv[3]) |>
    dplyr::mutate_if(is.numeric,as.factor)

  colnames(x) <- c('var1', 'var2', 'var3')

  # x <- expand.grid(1:Bv[1], 1:Bv[2]) |>
  #   dplyr::mutate_if(is.numeric,as.factor)
  #
  # colnames(x) <- c('var1', 'var2')


  b1 <- rnorm(B1, 0, sb)
  b1 <- b1 - mean(b1)

  b2 <- rnorm(B2, 0, sb)
  b2 <- b2 - mean(b2)

  b3 <- rnorm(B3, 0, sb)
  b3 <- b3 - mean(b3)


  beta1 <- genInt(B1, Q)
  beta2 <- genInt(B2, Q)
  beta3 <- genInt(B3, Q)

  blin = rep(0, N)
  for (k in 1:Q) {
    blin <- blin + lambda[k]*beta1[x[,'var1'],k]*beta2[x[,'var2'],k]
    #print(blin)
  }


  blin2 = rep(0, N)
  for (k in 1:Q) {
    blin2 <- blin2 + lambda[k]*beta1[x[,'var1'],k]*beta2[x[,'var2'],k]*beta3[x[,'var3'],k]
    #print(blin)
  }

  m <-  mu + b1[x[,'var1']] + b2[x[,'var2']] + b3[x[,'var3']] + blin + blin2
  y <-  rnorm(N, m, sigma)

  return(list(mu = mu,
              b1 = b1,
              b2 = b2,
              b3 = b3,
              #b4 = b4,
              blin = blin,
              blin2 = blin2,
              y = y,
              betas = list(beta1 = beta1, beta2 = beta2, beta3 = beta3),#, beta4 = beta4),
              Bv = Bv))
}

# Function to generate data from Model 3

genData_3 <- function(B1, B2, B3, sigma, sb, lambda){
  Bv <- c(B1, B2, B3)
  #Bv <- c(B1, B2)
  N <- Reduce("*", Bv)
  Q <- length(lambda)
  V <- length(Bv)
  mu <- 10

  x <- expand.grid(1:Bv[1], 1:Bv[2], 1:Bv[3]) |>
    dplyr::mutate_if(is.numeric,as.factor)

  colnames(x) <- c('var1', 'var2', 'var3')

  # x <- expand.grid(1:Bv[1], 1:Bv[2]) |>
  #   dplyr::mutate_if(is.numeric,as.factor)
  #
  # colnames(x) <- c('var1', 'var2')


  b1 <- rnorm(B1, 0, sb)
  b1 <- b1 - mean(b1)

  b2 <- rnorm(B2, 0, sb)
  b2 <- b2 - mean(b2)

  b3 <- rnorm(B3, 0, sb)
  b3 <- b3 - mean(b3)


  beta1 <- genInt(B1, Q)
  beta2 <- genInt(B2, Q)
  beta3 <- genInt(B3, Q)

  blin = rep(0, N)
  for (k in 1:Q) {
    blin <- blin + lambda[k]*beta1[x[,'var1'],k]*beta2[x[,'var2'],k]*beta3[x[,'var3'],k]
    #print(blin)
  }

  m <-  mu + b1[x[,'var1']] + b2[x[,'var2']] + b3[x[,'var3']] + blin
  y <-  rnorm(N, m, sigma)

  return(list(mu = mu,
              b1 = b1,
              b2 = b2,
              b3 = b3,
              #b4 = b4,
              blin = blin,
              y = y,
              betas = list(beta1 = beta1, beta2 = beta2, beta3 = beta3),#, beta4 = beta4),
              Bv = Bv))

}

# Four variables

genData_4 <- function(B1, B2, B3, B4, sigma, sb, lambda){
  Bv <- c(B1, B2, B3, B4)
  #Bv <- c(B1, B2)
  N <- Reduce("*", Bv)
  Q <- length(lambda)
  V <- length(Bv)
  mu <- 10

  x <- expand.grid(1:Bv[1], 1:Bv[2], 1:Bv[3], 1:Bv[3]) |>
    dplyr::mutate_if(is.numeric,as.factor)

  colnames(x) <- c('var1', 'var2', 'var3', 'var4')

  # x <- expand.grid(1:Bv[1], 1:Bv[2]) |>
  #   dplyr::mutate_if(is.numeric,as.factor)
  #
  # colnames(x) <- c('var1', 'var2')


  b1 <- rnorm(B1, 0, sb)
  b1 <- b1 - mean(b1)

  b2 <- rnorm(B2, 0, sb)
  b2 <- b2 - mean(b2)

  b3 <- rnorm(B3, 0, sb)
  b3 <- b3 - mean(b3)

  b4 <- rnorm(B4, 0, sb)
  b4 <- b4 - mean(b4)


  beta1 <- genInt(B1, Q)
  beta2 <- genInt(B2, Q)
  beta3 <- genInt(B3, Q)
  beta4 <- genInt(B4, Q)

  blin = rep(0, N)
  for (k in 1:Q) {
    blin <- blin + lambda[k]*beta1[x[,'var1'],k]*beta2[x[,'var2'],k]*beta3[x[,'var3'],k]
    #print(blin)
  }

  m <-  mu + b1[x[,'var1']] + b2[x[,'var2']] + b3[x[,'var3']] + b4[x[,'var4']] + blin
  y <-  rnorm(N, m, sigma)

  return(list(mu = mu,
              b1 = b1,
              b2 = b2,
              b3 = b3,
              b4 = b4,
              #b4 = b4,
              blin = blin,
              y = y,
              betas = list(beta1 = beta1, beta2 = beta2, beta3 = beta3, beta4 = beta4),#, beta4 = beta4),
              Bv = Bv))

}



# Run models --------------------------------------------------------------

# gen from genData_1

B1 <- 6
B2 <- 3
B3 <- 2
#B4 <- 2
lambda <- 10#c(8,10, 12)
sb <- 1
sigma <- 1

set.seed(022)
#x <- .Random.seed
dat <- genData_1(B1 = B1, B2 = B2, B3 = B3,  sb = sb, sigma = sigma, lambda = lambda)


# run the two models

model_1 <- model_jags_1(data = dat,
                        Q = 1,
                        mmu = 10,
                        smu = 1,
                        a = 0.1,
                        b = 0.1,
                        nthin = 2,
                        nburnin = 2000)


model_7 <- model_jags_7(data = dat,
                        Q = 1,
                        mmu = 10,
                        smu = 1,
                        a = 0.1,
                        b = 0.1,
                        nthin = 2,
                        nburnin = 2000)


caret::RMSE(dat$blin, model_1$BUGSoutput$mean$int)
caret::RMSE(dat$blin, model_7$BUGSoutput$mean$int)

qplot(dat$blin, model_1$BUGSoutput$mean$int) + geom_abline() + theme_bw() +
  labs(x = 'True interaction', y = 'Predicted',
       title = quote(paste('Model 1: ', '\u03A3'[q]^Q, ' ', lambda[q], beta[iq]^(1),beta[jq]^(2), ' - RMSE = 0.3121'))) +
qplot(dat$blin, model_7$BUGSoutput$mean$int) + geom_abline() + theme_bw() +
  labs(x = 'True interaction', y = 'Predicted',
       title =  quote(paste('Model 3: ', '\u03A3'[q]^Q, ' ', lambda[q],
                            beta[iq]^(1),beta[jq]^(2),beta[kq]^(3), ' - RMSE = 0.3166')))




# gen from genData_2

B1 <- 6
B2 <- 3
B3 <- 2
#B4 <- 2
lambda <- 10#c(8,10, 12)
sb <- 1
sigma <- 1

set.seed(022)
#x <- .Random.seed
dat <- genData_2(B1 = B1, B2 = B2, B3 = B3,  sb = sb, sigma = sigma, lambda = lambda)


# run the two models

model_1 <- model_jags_1(data = dat,
                        Q = 1,
                        mmu = 10,
                        smu = 1,
                        a = 0.1,
                        b = 0.1,
                        nthin = 2,
                        nburnin = 2000)

model_2 <- model_jags_2(data = dat,
                        Q = 1,
                        mmu = 10,
                        smu = 1,
                        a = 0.1,
                        b = 0.1,
                        nthin = 2,
                        nburnin = 2000)


model_7 <- model_jags_7(data = dat,
                        Q = 2,
                        mmu = 10,
                        smu = 1,
                        a = 0.1,
                        b = 0.1,
                        nthin = 2,
                        nburnin = 2000)

true_int <- dat$blin + dat$blin2
int2 <-  model_2$BUGSoutput$mean$int_p1 + model_2$BUGSoutput$mean$int_p2

caret::RMSE(true_int, model_1$BUGSoutput$mean$int)
caret::RMSE(true_int, int2)
caret::RMSE(true_int, model_7$BUGSoutput$mean$int)



qplot(true_int, model_1$BUGSoutput$mean$int) + geom_abline() + theme_bw() +
  labs(x = 'True interaction', y = 'Predicted',
       title = quote(paste('Model 1: ', '\u03A3'[q]^Q, ' ', lambda[q], beta[iq]^(1),beta[jq]^(2), ' - RMSE = 1.9845'))) +
qplot(true_int, int2) + geom_abline() + theme_bw() +
  labs(x = 'True interaction', y = 'Predicted',
       title = quote(atop(paste('Model 2: ', '\u03A3'[q]^Q[1], ' ', lambda[q]^(1),
                           beta[iq]^(1),beta[jq]^(2), '+', '\u03A3'[q]^Q[2], ' ', lambda[q]^(2),
                           beta[iq]^(1),beta[jq]^(2),beta[kq]^(3)), ' RMSE = 0.4985'))) +
qplot(true_int, model_7$BUGSoutput$mean$int) + geom_abline() + theme_bw() +
  labs(x = 'True interaction', y = 'Predicted',
       title =  quote(paste('Model 3: ', '\u03A3'[q]^Q, ' ', lambda[q],
                            beta[iq]^(1),beta[jq]^(2),beta[kq]^(3), ' - RMSE = 0.5128')))






# Simulate and run with 4 variables


B1 <- 6
B2 <- 3
B3 <- 2
B4 <- 2
lambda <- 10#c(8,10, 12)
sb <- 1
sigma <- 1

set.seed(022)
#x <- .Random.seed
dat <- genData_4(B1 = B1, B2 = B2, B3 = B3, B4 = B4,  sb = sb, sigma = sigma, lambda = lambda)



model_1 <- model_jags_9(data = dat,
                        Q = 1,
                        mmu = 10,
                        smu = 1,
                        a = 0.1,
                        b = 0.1,
                        nthin = 2,
                        nburnin = 2000)


model_7 <- model_jags_8(data = dat,
                        Q = 1,
                        mmu = 10,
                        smu = 1,
                        a = 0.1,
                        b = 0.1,
                        nthin = 2,
                        nburnin = 2000)


caret::RMSE(dat$blin, model_1$BUGSoutput$mean$int)
caret::RMSE(dat$blin, model_7$BUGSoutput$mean$int)

qplot(dat$blin, model_1$BUGSoutput$mean$int) + geom_abline() + theme_bw() +
  labs(x = 'True interaction', y = 'Predicted',
       title = quote(paste('Model 1: ', '\u03A3'[q]^Q, ' ', lambda[q], beta[iq]^(1),beta[jq]^(2), beta[kq]^(3), ' - RMSE = 0.3795'))) +
  qplot(dat$blin, model_7$BUGSoutput$mean$int) + geom_abline() + theme_bw() +
  labs(x = 'True interaction', y = 'Predicted',
       title =  quote(paste('Model 3: ', '\u03A3'[q]^Q, ' ', lambda[q],
                            beta[iq]^(1),beta[jq]^(2),beta[kq]^(3),beta[lq]^(4), ' - RMSE = 0.3790')))








#attr(dat, "seed") <- x
dat2 <- genData_1(B1 = B1, B2 = B2, B3 = B3,  sb = sb, sigma = sigma, lambda = lambda)


# check the constraints
dat$betas$beta2[,1] |> sum()
dat$betas$beta2[,1]^2 |> sum()

# run jags models
model_1 <- model_jags_1(data = dat,
                    Q = 1,
                    mmu = 10,
                    smu = 1,
                    a = 0.1,
                    b = 0.1,
                    nthin = 2,
                    nburnin = 2000)

model_2 <- model_jags_2(data = dat,
                        Q = 1,
                        mmu = 10,
                        smu = 1,
                        a = 0.1,
                        b = 0.1,
                        nthin = 2,
                        nburnin = 2000)

model_3 <- model_jags_3(data = dat,
                        Q = 1,
                        mmu = 10,
                        smu = 1,
                        a = 0.1,
                        b = 0.1,
                        nthin = 2,
                        nburnin = 2000)

int2 <- model_2$BUGSoutput$mean$int_p1 + model_2$BUGSoutput$mean$int_p2

qplot(dat$blin, model_1$BUGSoutput$mean$int) + geom_abline() + theme_bw() +
  labs(x = 'True interaction', y = 'Predicted',
       title = quote(paste('Model 1: ', '\u03A3'[q]^Q, ' ', lambda[q],
                           beta[iq]^(1),beta[jq]^(2)))) +
qplot(dat$blin, model_2$BUGSoutput$mean$int_p1) + geom_abline() + theme_bw() +
  labs(x = 'True interaction', y = 'Predicted',
       title = quote(paste('Model 2: ', '\u03A3'[q]^Q[1], ' ', lambda[q]^(1),
                           beta[iq]^(1),beta[jq]^(2)))) +
qplot(dat$blin, model_2$BUGSoutput$mean$int_p2) + geom_abline() + theme_bw() +
  labs(x = 'True interaction', y = 'Predicted',
       title = quote(paste('Model 2: ', '\u03A3'[q]^Q[2], ' ', lambda[q]^(2),
                           beta[iq]^(1),beta[jq]^(2),beta[kq]^(3)))) +
  qplot(dat$blin, int2) + geom_abline() + theme_bw() +
  labs(x = 'True interaction', y = 'Predicted',
       title = quote(paste('Model 2: ', '\u03A3'[q]^Q[1], ' ', lambda[q]^(1),
                           beta[iq]^(1),beta[jq]^(2), '+', '\u03A3'[q]^Q[2], ' ', lambda[q]^(2),
                           beta[iq]^(1),beta[jq]^(2),beta[kq]^(3)))) +
qplot(dat$blin, model_3$BUGSoutput$mean$int) + geom_abline() + theme_bw() +
  labs(x = 'True interaction', y = 'Predicted',
       title = quote(paste('Model 3: ', '\u03A3'[q]^Q, ' ', lambda[q],
                           beta[iq]^(1),beta[jq]^(2),beta[kq]^(3))))


rmse <- function(true, pred){
  caret::RMSE(true, pred)
}


# RMSE interaction
rmse(dat$blin, model_1$BUGSoutput$mean$int)
rmse(dat$blin, model_2$BUGSoutput$mean$int_p1)
rmse(dat$blin, model_2$BUGSoutput$mean$int_p2)
rmse(dat$blin, int2)
rmse(dat$blin, model_3$BUGSoutput$mean$int)


# RMSE y
rmse(dat$y, model_1$BUGSoutput$mean$mu)
rmse(dat$y, model_2$BUGSoutput$mean$mu)
rmse(dat$y, model_3$BUGSoutput$mean$mu)

qplot(dat$y, model_3$BUGSoutput$mean$mu) + geom_abline() + theme_bw()


model_1$BUGSoutput$mean$beta1 |> sum()
model_1$BUGSoutput$mean$beta1^2 |> sum()

model_all$BUGSoutput$mean$beta1 |> sum()



model_4 <- model_jags_4(data = dat,
                        Q = 1,
                        mmu = 10,
                        smu = 1,
                        a = 0.1,
                        b = 0.1,
                        nthin = 2,
                        nburnin = 2000)

qplot(dat$blin, model_4$BUGSoutput$mean$int) + geom_abline() + theme_bw()

model_4$BUGSoutput$mean$beta1 |> sum()
sum(model_4$BUGSoutput$mean$beta1^2)


# Model 5
model_5 <- model_jags_5(data = dat,
                        Q = 1,
                        mmu = 10,
                        smu = 1,
                        a = 0.1,
                        b = 0.1,
                        nthin = 2,
                        nburnin = 2000)
plot(model_5)
qplot(dat$blin, model_5$BUGSoutput$mean$int) + geom_abline() + theme_bw()
caret::RMSE(dat$blin, model_5$BUGSoutput$mean$int)



# Model 6

model_6 <- model_jags_6(data = dat,
                        Q = 1,
                        mmu = 10,
                        smu = 1,
                        a = 0.1,
                        b = 0.1,
                        nthin = 2,
                        nburnin = 2000,
                        alpha_beta = list(alpha_beta1 = rep(1, B1),
                                          alpha_beta2 = rep(1, B2),
                                          alpha_beta3 = rep(1, B3)))
plot(model_6)
qplot(dat$blin, model_6$BUGSoutput$mean$int) + geom_abline() + theme_bw()


# Model 7

model_7 <- model_jags_7(data = dat,
                        Q = 1,
                        mmu = 10,
                        smu = 1,
                        a = 0.1,
                        b = 0.1,
                        nthin = 2,
                        nburnin = 2000)
plot(model_7)
qplot(dat$blin, model_7$BUGSoutput$mean$int) + geom_abline() + theme_bw()
caret::RMSE(dat$blin, model_7$BUGSoutput$mean$int)


# Two methods

qplot(dat$blin, model_5$BUGSoutput$mean$int) + geom_abline() +
  theme_bw() + labs(x = 'True interaction', y = 'Predicted', title = 'Alessandra - RMSE = 0.5346') +
qplot(dat$blin, model_7$BUGSoutput$mean$int) + geom_abline() +
  theme_bw() + labs(x = 'True interaction', y = 'Predicted', title = 'Andrew - RMSE = 0.3166')


# Simulate from model 2 ---------------------------------------------------


# generate data

genData2 <- function(B1, B2, B3, sigma, sb, lambda){ #sb1, sb2, sb3, sb4){

  Bv <- c(B1, B2, B3)
  #Bv <- c(B1, B2)
  N <- Reduce("*", Bv)
  Q <- length(lambda)
  V <- length(Bv)
  mu <- 10

  x <- expand.grid(1:Bv[1], 1:Bv[2], 1:Bv[3]) |>
    dplyr::mutate_if(is.numeric,as.factor)

  colnames(x) <- c('var1', 'var2', 'var3')

  # x <- expand.grid(1:Bv[1], 1:Bv[2]) |>
  #   dplyr::mutate_if(is.numeric,as.factor)
  #
  # colnames(x) <- c('var1', 'var2')


  b1 <- rnorm(B1, 0, sb)
  b1 <- b1 - mean(b1)

  b2 <- rnorm(B2, 0, sb)
  b2 <- b2 - mean(b2)

  b3 <- rnorm(B3, 0, sb)
  b3 <- b3 - mean(b3)


  beta1 <- genInt(B1, Q)
  beta2 <- genInt(B2, Q)
  beta3 <- genInt(B3, Q)

  blin = rep(0, N)
  for (k in 1:Q) {
    blin <- blin + lambda[k]*beta1[x[,'var1'],k]*beta2[x[,'var2'],k]
    #print(blin)
  }


  blin2 = rep(0, N)
  for (k in 1:Q) {
    blin2 <- blin2 + lambda[k]*beta1[x[,'var1'],k]*beta2[x[,'var2'],k]*beta3[x[,'var3'],k]
    #print(blin)
  }

  m <-  mu + b1[x[,'var1']] + b2[x[,'var2']] + b3[x[,'var3']] + blin + blin2
  y <-  rnorm(N, m, sigma)

  return(list(mu = mu,
              b1 = b1,
              b2 = b2,
              b3 = b3,
              #b4 = b4,
              blin = blin,
              blin2 = blin2,
              y = y,
              betas = list(beta1 = beta1, beta2 = beta2, beta3 = beta3),#, beta4 = beta4),
              Bv = Bv))

}


B1 <- 6
B2 <- 3
B3 <- 2
#B4 <- 2
lambda <- 10#c(8,10, 12)
sb <- 1
sigma <- 1

set.seed(022)
dat <- genData2(B1 = B1, B2 = B2, B3 = B3,  sb = sb, sigma = sigma, lambda = lambda)


model_2_1 <- model_jags_2(data = dat,
                        Q = 1,
                        mmu = 10,
                        smu = 1,
                        a = 0.1,
                        b = 0.1,
                        nthin = 2,
                        nburnin = 2000)

model_2_2 <- model_jags_3(data = dat,
                        Q = 1,
                        mmu = 10,
                        smu = 1,
                        a = 0.1,
                        b = 0.1,
                        nthin = 2,
                        nburnin = 2000)



int2_1 <- model_2_1$BUGSoutput$mean$int_p1 + model_2_1$BUGSoutput$mean$int_p2

qplot(dat$blin, model_2_1$BUGSoutput$mean$int_p1) + geom_abline() + theme_bw() +
  labs(x = 'True interaction', y = 'Predicted',
       title = quote(paste('Model 2: ', '\u03A3'[q]^Q[1], ' ', lambda[q]^(1),
                           beta[iq]^(1),beta[jq]^(2)))) +
  qplot(dat$blin2, model_2_1$BUGSoutput$mean$int_p2) + geom_abline() + theme_bw() +
  labs(x = 'True interaction', y = 'Predicted',
       title = quote(paste('Model 2: ', '\u03A3'[q]^Q[2], ' ', lambda[q]^(2),
                           beta[iq]^(1),beta[jq]^(2),beta[kq]^(3)))) +
  qplot(dat$blin + dat$blin2, int2_1) + geom_abline() + theme_bw() +
  labs(x = 'True interaction', y = 'Predicted',
       title = quote(paste('Model 2: ', '\u03A3'[q]^Q[1], ' ', lambda[q]^(1),
                           beta[iq]^(1),beta[jq]^(2), '+', '\u03A3'[q]^Q[2], ' ', lambda[q]^(2),
                           beta[iq]^(1),beta[jq]^(2),beta[kq]^(3)))) +
  qplot(dat$blin + dat$blin2, model_2_2$BUGSoutput$mean$int) + geom_abline() + theme_bw() +
  labs(x = 'True interaction', y = 'Predicted',
       title = quote(paste('Model 3: ', '\u03A3'[q]^Q, ' ', lambda[q],
                           beta[iq]^(1),beta[jq]^(2),beta[kq]^(3))))



model_5 <- model_jags_5(data = dat,
                        Q = 1,
                        mmu = 10,
                        smu = 1,
                        a = 0.1,
                        b = 0.1,
                        nthin = 2,
                        nburnin = 2000)
plot(model_5)
qplot(dat$blin, model_5$BUGSoutput$mean$int) + geom_abline() + theme_bw()




# Model 6

model_6 <- model_jags_6(data = dat,
                        Q = 1,
                        mmu = 10,
                        smu = 1,
                        a = 0.1,
                        b = 0.1,
                        nthin = 2,
                        nburnin = 2000,
                        alpha_beta = list(alpha_beta1 = rep(1, B1),
                                          alpha_beta2 = rep(1, B2),
                                          alpha_beta3 = rep(1, B3)))
plot(model_6)
qplot(dat$blin, model_6$BUGSoutput$mean$int) + geom_abline() + theme_bw()



# Simulate from the full model --------------------------------------------

# generate data
genData <- function(B1, B2, B3, sigma, sb, lambda){ #sb1, sb2, sb3, sb4){

  Bv <- c(B1, B2, B3)
  #Bv <- c(B1, B2)
  N <- Reduce("*", Bv)
  Q <- length(lambda)
  V <- length(Bv)
  mu <- 10

  x <- expand.grid(1:Bv[1], 1:Bv[2], 1:Bv[3]) |>
    dplyr::mutate_if(is.numeric,as.factor)

  colnames(x) <- c('var1', 'var2', 'var3')

  # x <- expand.grid(1:Bv[1], 1:Bv[2]) |>
  #   dplyr::mutate_if(is.numeric,as.factor)
  #
  # colnames(x) <- c('var1', 'var2')


  b1 <- rnorm(B1, 0, sb)
  b1 <- b1 - mean(b1)

  b2 <- rnorm(B2, 0, sb)
  b2 <- b2 - mean(b2)

  b3 <- rnorm(B3, 0, sb)
  b3 <- b3 - mean(b3)


  beta1 <- genInt(B1, Q)
  beta2 <- genInt(B2, Q)
  beta3 <- genInt(B3, Q)

  blin = rep(0, N)
  for (k in 1:Q) {
    blin <- blin + lambda[k]*beta1[x[,'var1'],k]*beta2[x[,'var2'],k]*beta3[x[,'var3'],k]
    #print(blin)
  }

  m <-  mu + b1[x[,'var1']] + b2[x[,'var2']] + b3[x[,'var3']] + blin
  y <-  rnorm(N, m, sigma)

  return(list(mu = mu,
              b1 = b1,
              b2 = b2,
              b3 = b3,
              #b4 = b4,
              blin = blin,
              y = y,
              betas = list(beta1 = beta1, beta2 = beta2, beta3 = beta3),#, beta4 = beta4),
              Bv = Bv))

}


B1 <- 6
B2 <- 3
B3 <- 2
#B4 <- 2
lambda <- 10#c(8,10, 12)
sb <- 1
sigma <- 1

set.seed(022)
#x <- .Random.seed
dat <- genData(B1 = B1, B2 = B2, B3 = B3,  sb = sb, sigma = sigma, lambda = lambda)


model_full_1 <- model_jags_1(data = dat,
                        Q = 1,
                        mmu = 10,
                        smu = 1,
                        a = 0.1,
                        b = 0.1,
                        nthin = 2,
                        nburnin = 2000)

model_full_2 <- model_jags_2(data = dat,
                        Q = 1,
                        mmu = 10,
                        smu = 1,
                        a = 0.1,
                        b = 0.1,
                        nthin = 2,
                        nburnin = 2000)

model_full_3 <- model_jags_3(data = dat,
                        Q = 1,
                        mmu = 10,
                        smu = 1,
                        a = 0.1,
                        b = 0.1,
                        nthin = 2,
                        nburnin = 2000)


for (i in 1:6000) {
  print(model_full_3$BUGSoutput$sims.list$beta1[i, ,]^2 |> sum())
}


int_full_2 <- model_full_2$BUGSoutput$mean$int_p1 + model_full_2$BUGSoutput$mean$int_p2

qplot(dat$blin, model_full_1$BUGSoutput$mean$int) + geom_abline() + theme_bw() +
  labs(x = 'True interaction', y = 'Predicted',
       title = quote(paste('Model 1: ', '\u03A3'[q]^Q, ' ', lambda[q],
                           beta[iq]^(1),beta[jq]^(2)))) +
  qplot(dat$blin, int_full_2) + geom_abline() + theme_bw() +
  labs(x = 'True interaction', y = 'Predicted',
       title = quote(paste('Model 2: ', '\u03A3'[q]^Q[1], ' ', lambda[q]^(1),
                           beta[iq]^(1),beta[jq]^(2), '+', '\u03A3'[q]^Q[2], ' ', lambda[q]^(2),
                           beta[iq]^(1),beta[jq]^(2),beta[kq]^(3)))) +
  qplot(dat$blin, model_full_3$BUGSoutput$mean$int) + geom_abline() + theme_bw() +
  labs(x = 'True interaction', y = 'Predicted',
       title = quote(paste('Model 3: ', '\u03A3'[q]^Q, ' ', lambda[q],
                           beta[iq]^(1),beta[jq]^(2),beta[kq]^(3))))

model_5 <- model_jags_5(data = dat,
                        Q = 1,
                        mmu = 10,
                        smu = 1,
                        a = 0.1,
                        b = 0.1,
                        nthin = 2,
                        nburnin = 2000)
plot(model_5)
qplot(dat$blin, model_5$BUGSoutput$mean$int) + geom_abline() + theme_bw()




# Model 6

model_6 <- model_jags_6(data = dat,
                        Q = 1,
                        mmu = 10,
                        smu = 1,
                        a = 0.1,
                        b = 0.1,
                        nthin = 2,
                        nburnin = 2000,
                        alpha_beta = list(alpha_beta1 = rep(1, B1),
                                          alpha_beta2 = rep(1, B2),
                                          alpha_beta3 = rep(1, B3)))
plot(model_6)
qplot(dat$blin, model_6$BUGSoutput$mean$int) + geom_abline() + theme_bw()



# -------------------------------------------------------------------------

model <- bammitJags(data = dat,
                    Q = 1,
                    mmu = 10,
                    smu = 1,
                    a = 0.1,
                    b = 0.1,
                    nthin = 2,
                    nburnin = 2000)

qplot(dat$y, model$BUGSoutput$mean$mu) + geom_abline() + theme_bw() +
  labs(x = 'y', y = expression(hat(y)), title = 'Model 1') +
  qplot(dat$y, model$BUGSoutput$mean$mu2) + geom_abline() + theme_bw() +
    labs(x = 'y', y = expression(hat(y)), title = 'Model 2') +
  qplot(dat$y, model$BUGSoutput$mean$mu3) + geom_abline() + theme_bw() +
  labs(x = 'y', y = expression(hat(y)), title = 'Model 3')


qplot(dat$blin, model$BUGSoutput$mean$int) + geom_abline() + theme_bw() +
  labs(x = 'True interaction', y = 'Predicted', title = 'Model 1') +
  qplot(dat$blin, model$BUGSoutput$mean$int2) + geom_abline() + theme_bw() +
  labs(x = 'True interaction', y = '', title = 'Model 2') +
  qplot(dat$blin, model$BUGSoutput$mean$int_p2) + geom_abline() + theme_bw() +
  labs(x = 'True interaction', y = '', title = 'Model 3')


caret::RMSE(dat$blin, model$BUGSoutput$mean$int)
caret::RMSE(dat$blin, model$BUGSoutput$mean$int2)

model$BUGSoutput$mean$beta1


# Second way --------------------------------------------------------------


simLObammit <- function(V = 3,
                        Q = 1,
                        Bv = c(3,2,2),
                        mu = 100,
                        VInt = c(1,2),
                        lambda = c(10),
                        sb = 1,
                        sB = 10,
                        sy = 1){

  N <- Reduce("*", Bv)
  BvComp <- Bv[which(!(seq(1:V) %in% VInt))]

  # generate main effects
  bv <- lapply(Bv, function(x) {
    bvaux <- rnorm(x, 0, sb)
    bvF <- bvaux - mean(bvaux)
    return(bvF)
  })
  #names(bv) <- paste0("b", 1:length(Bv))

  meff <- rowSums(rev(expand.grid(rev(bv))))

  # generate bilinear term
  Beta <- vector(mode = "list", length = length(VInt))
  for (i in 1:length(VInt)) {
    Beta[[i]] <- genInt(Bv[i], Q)
  }

  Ones <- vector(mode = "list", length = V-length(VInt))
  for (i in 1:(V-length(VInt))) {
    Ones[[i]] <- matrix(1, BvComp[i], Q)
  }

  k <- as.list(rep(1, Q))
  for (j in 1:length(k)) {
    for (i in 1:length(Beta)) {
      k[[j]] <- kronecker(k[[j]], Beta[[i]][,j])
    }
    for (i in 1:length(Ones)) {
      k[[j]] <- kronecker(k[[j]], Ones[[i]][,j])
    }
    k[[j]] <- k[[j]]*lambda[j]
  }

  int <- Reduce("+", k)

  # generate y
  m <- mu + meff + int
  y <- rnorm(N, m, sy)

  return(list(mu = mu,
              bv = bv,
              int = int,
              y = y,
              Bv = Bv))

}



V = 4
Q = 1
Bv = c(4,3,2,2)
mu = 10
lambda = c(10)
sb = 10
sB = 10
sy = 1.5
VInt = c(1,2)

dat <- simLObammit(V = V, Q = Q, Bv = Bv, mu = mu, lambda = lambda,
                   sb = sb, sB = sB, sy = sy, VInt = VInt)


model <- bammitJags(data = dat,
                    Q = 1,
                    mmu = 10,
                    smu = 1,
                    a = 0.1,
                    b = 0.1,
                    nthin = 2,
                    nburnin = 2000)

qplot(dat$y, model$BUGSoutput$mean$mu) + geom_abline() + theme_bw() +
  labs(x = 'y', y = expression(hat(y)), title = 'Model 1') +
  qplot(dat$y, model$BUGSoutput$mean$mu2) + geom_abline() + theme_bw() +
  labs(x = 'y', y = expression(hat(y)), title = 'Model 2')


qplot(dat$int, model$BUGSoutput$mean$int) + geom_abline() + theme_bw() +
  labs(x = 'True interaction', y = 'Predicted', title = 'Model 1') +
  qplot(dat$int, model$BUGSoutput$mean$int2) + geom_abline() + theme_bw() +
  labs(x = 'True interaction', y = 'Predicted', title = 'Model 2')



Beta[[1]] <- matrix(c(c(1,2,3), c(1,2,3)), ncol = 2)
Beta[[2]] <- matrix(c(c(1,2), c(1,2)), ncol = 2)

beta1 <- c(1,2,3)
beta2 <- c(1,2)
beta3 <- c(3,5)
k1 <- kronecker(beta1, beta2)
k1
k2 <- kronecker(k1, beta3)
k2

kronecker(k1, c(1,1))




# -------------------------------------------------------------------------


# Stan model --------------------------------------------------------------------

library(rstan)
ammiStan <- function(data, Q, mmu, smu, a,b, mtheta = 0, stheta = 1, nthin, nburnin){

  data = dat
  Q = 1
  mtheta = 0
  stheta = 1
  mmu = 100
  smu = 10
  a = 0.1
  b = 0.1


  Y <- data$y
  V <- length(data$Bv)
  Bv <- data$Bv
  N <- length(data$y)

  B1 <- Bv[1]
  B2 <- Bv[2]

  x <- expand.grid(lapply(Bv, function(x) 1:x))
  colnames(x) <- paste0('var', 1:2)
  var1 <- x$var1
  var2 <- x$var2

  # Set up the data
  stan_data <- list(
    N = N,
    B1 = B1,
    B2 = B2,
    Q = Q,
    Y = Y,
    var1 = var1,
    var2 = var2,
    mmu = mmu,
    smu = smu,
    a = a,
    b = b,
    mtheta = mtheta,
    stheta = stheta
  )

  # Run the model
  modelRun <- stan(
    file = '~/Documents/GitHub/bammit/Results/codes/ammiStan.stan',  # the model code
    data = stan_data,  # the data for the model
    iter = 500,  # the number of MCMC iterations
    chains = 2  # the number of MCMC chains
  )

  return(modelRun)

}









# Model with andrew correction and V = 4 ----------------------------------

library(R2jags)
library(ggplot2)
library(patchwork)

rm(list = ls())

# Generate betas
genInt <- function(index = 10, Q = 1, stheta = 1) {
  # Generate theta matrix
  theta <- matrix(rnorm(index * Q, 0, stheta), nrow = index, ncol = Q)
  # Calculate the column means of theta
  m <- colMeans(theta)
  # Center the columns of theta
  thetaN <- sweep(theta, 2, m)
  # Calculate sqrtTheta
  sqrtTheta <- sapply(1:Q, function(q) {
    sqrt(1 / sum(thetaN[, q]^2))
  })
  # Calculate variable
  variable <- sweep(thetaN, 2, sqrtTheta, "*")
  return(variable)
}


# Model 8 - Andrew's code
model_jags_8 <- function(data, Q, mmu, smu, a,b, nthin, nburnin, alpha_beta){

  Y <- data$y
  V <- length(data$Bv)
  Bv <- data$Bv
  N <- length(data$y)



  modelCode <- "
model{
  # Likelihood
  for (n in 1:N) {
    Y[n] ~ dnorm(mu[n], sy^-2)
    mu[n] = muall + b1[var1[n]] + b2[var2[n]] + b3[var3[n]] + b4[var4[n]] + int[n]
    int[n] = sum(lambda[1:Q] * beta1[var1[n],1:Q] * beta2[var2[n],1:Q] * beta3[var3[n],1:Q] * beta4[var4[n],1:Q])
  }

  # Priors
  muall ~ dnorm(mmu, smu^-2) # grand mean

  # Main effects
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

  # Interaction beta 1
  for(q in 1:Q) {
    # generate values of theta from a Normal distribution
    for(i in 1:B1) {
      thetaBeta1[i,q] ~ dnorm(0, 1^-2)
    }
    # calculate the mean of each column
    mBeta1[q] = sum(thetaBeta1[1:B1,q])/B1

    # centering
    for(i in 1:B1) {
      thetaBeta1New[i,q] = thetaBeta1[i,q] - mBeta1[q]
    }

    # get sd of each column
    sqrtThetaBeta1[q] = sqrt(sum(thetaBeta1New[1:B1,q]^2)) + 0.00000001

    # get the final beta
    for(i in 1:B1) {
      beta1[i,q] = (thetaBeta1New[i,q] + 10000*M1[q]) / (sqrtThetaBeta1[q] + 10000*M1[q])
    }
    M1[q] ~ dbern(p1[q]) # M1 = 0 means normal, M1 = 1 means all beta1 = 1
    p1[q] ~ dbeta(1, 10)
  }

  # Interaction beta 2
  for(q in 1:Q) {
    for(i in 1:B2) {
      thetaBeta2[i,q] ~ dnorm(0, 1^-2)
    }
    mBeta2[q] = sum(thetaBeta2[1:B2,q])/B2

    for(i in 1:B2) {
      thetaBeta2New[i,q] = thetaBeta2[i,q] - mBeta2[q]
    }
    sqrtThetaBeta2[q] = sqrt(sum(thetaBeta2New[1:B2,q]^2)) + 0.00000001
    for(i in 1:B2) {
      beta2[i,q] = (thetaBeta2New[i,q] + 10000*M2[q]) / (sqrtThetaBeta2[q] + 10000*M2[q])
    }
    M2[q] ~ dbern(p2[q])
    p2[q] ~ dbeta(1, 10)
  }

  # Interaction beta 3
  for(q in 1:Q) {
    for(i in 1:B3) {
      thetaBeta3[i,q] ~ dnorm(0, 1^-2)
    }
    mBeta3[q] = sum(thetaBeta3[1:B3,q])/B3

    for(i in 1:B3) {
      thetaBeta3New[i,q] = thetaBeta3[i,q] - mBeta3[q]
    }
    sqrtThetaBeta3[q] = sqrt(sum(thetaBeta3New[1:B3,q]^2)) + 0.00000001
    for(i in 1:B3) {
      beta3[i,q] = (thetaBeta3New[i,q] + 10000*M3[q]) / (sqrtThetaBeta3[q] + 10000*M3[q])
    }
    M3[q] ~ dbern(p3[q])
    p3[q] ~ dbeta(1, 10)
  }


    # Interaction beta 4
  for(q in 1:Q) {
    for(i in 1:B4) {
      thetaBeta4[i,q] ~ dnorm(0, 1^-2)
    }
    mBeta4[q] = sum(thetaBeta4[1:B4,q])/B4

    for(i in 1:B4) {
      thetaBeta4New[i,q] = thetaBeta4[i,q] - mBeta4[q]
    }
    sqrtThetaBeta4[q] = sqrt(sum(thetaBeta4New[1:B4,q]^2)) + 0.00000001
    for(i in 1:B4) {
      beta4[i,q] = (thetaBeta4New[i,q] + 10000*M4[q]) / (sqrtThetaBeta4[q] + 10000*M4[q])
    }
    M4[q] ~ dbern(p4[q])
    p4[q] ~ dbeta(1, 10)
  }


  # Prior on lambda
  for(q in 1:Q) {
    lambda_raw[q] ~ dnorm(0, slambda^-2)T(0,)
  }
  lambda = sort(lambda_raw)

  # Priors
  sb1 ~ dt(0, 1^-2, 1)T(0,)
  sb2 ~ dt(0, 1^-2, 1)T(0,)
  sb3 ~ dt(0, 1^-2, 1)T(0,)
  sb4 ~ dt(0, 1^-2, 1)T(0,)
  slambda ~ dt(0, 1^-2, 1)T(0,)
  sy ~ dgamma(a, b) # Prior on residual sd
}"

  B1 <- Bv[1]
  B2 <- Bv[2]
  B3 <- Bv[3]
  B4 <- Bv[4]

  x <- expand.grid(lapply(Bv, function(x) 1:x))
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
    a = a,
    b = b
  )

  # Choose the parameters to watch
  modelParameters <- c(
    "b1", "b2", "b3", 'b4', "lambda", "beta1", "beta2", "beta3", "beta4", "sy", "muall",
    "int", "mu", "p1", "p2", "p3", 'p4'
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

# Model 9

model_jags_9 <- function(data, Q, mmu, smu, a,b, stheta = 1, nthin, nburnin){

  Y <- data$y
  V <- length(data$Bv)
  Bv <- data$Bv
  N <- length(data$y)


  modelCode <- "
    model{
      # Likelihood
      for (n in 1:N) {
          Y[n] ~ dnorm(mu[n], sy)
          mu[n] = muall + b1[var1[n]] + b2[var2[n]] + b3[var3[n]] + int_p1[n] + int_p2[n] + int_p3[n]
          int_p1[n] = sum(lambda_1[1:Q] * beta1[var1[n],1:Q] * beta2[var2[n],1:Q])
          int_p2[n] = sum(lambda_2[1:Q] * beta1[var1[n],1:Q] * beta2[var2[n],1:Q]* beta3[var3[n],1:Q])
          int_p3[n] = sum(lambda_2[1:Q] * beta1[var1[n],1:Q] * beta2[var2[n],1:Q]* beta3[var3[n],1:Q]* beta4[var4[n],1:Q])
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
            thetaBeta1[i,q] ~ dnorm(0, 1/stheta^2)
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
            thetaBeta2[i,q] ~ dnorm(0, 1/stheta^2)
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
        lambda_raw_1[q] ~ dnorm(0, slambda^-2)T(0,)
      }
      lambda_1 = sort(lambda_raw_1)


      for(q in 1:Q) {
        lambda_raw_2[q] ~ dnorm(0, slambda^-2)T(0,)
      }
      lambda_2 = sort(lambda_raw_2)

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

  x <- expand.grid(lapply(Bv, function(x) 1:x))
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
    "b1", "b2", "b3", "b4", "lambda_1", 'lambda_2', "beta1", "beta2", 'beta3', 'beta4', "sy", "muall",
    "mu", 'int_p1', 'int_p2', 'int_p3'
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


# Function to generate data from Model

genData_3 <- function(B1, B2, B3, B4, sigma, sb, lambda){
  Bv <- c(B1, B2, B3, B4)
  #Bv <- c(B1, B2)
  N <- Reduce("*", Bv)
  Q <- length(lambda)
  V <- length(Bv)
  mu <- 10

  x <- expand.grid(1:Bv[1], 1:Bv[2], 1:Bv[3], 1:Bv[4]) |>
    dplyr::mutate_if(is.numeric,as.factor)

  colnames(x) <- c('var1', 'var2', 'var3', 'var4')

  # x <- expand.grid(1:Bv[1], 1:Bv[2]) |>
  #   dplyr::mutate_if(is.numeric,as.factor)
  #
  # colnames(x) <- c('var1', 'var2')


  b1 <- rnorm(B1, 0, sb)
  b1 <- b1 - mean(b1)

  b2 <- rnorm(B2, 0, sb)
  b2 <- b2 - mean(b2)

  b3 <- rnorm(B3, 0, sb)
  b3 <- b3 - mean(b3)

  b4 <- rnorm(B4, 0, sb)
  b4 <- b4 - mean(b4)


  beta1 <- genInt(B1, Q)
  beta2 <- genInt(B2, Q)
  beta3 <- genInt(B3, Q)
  beta4 <- genInt(B4, Q)

  blin = rep(0, N)
  for (k in 1:Q) {
    blin <- blin + lambda[k]*beta1[x[,'var1'],k]*beta2[x[,'var2'],k]
    #print(blin)
  }

  blin2 = rep(0, N)
  for (k in 1:Q) {
    blin2 <- blin2 + lambda[k]*beta1[x[,'var1'],k]*beta2[x[,'var2'],k]*beta3[x[,'var3'],k]
    #print(blin)
  }

  blin3 = rep(0, N)
  for (k in 1:Q) {
    blin3 <- blin3 + lambda[k]*beta1[x[,'var1'],k]*beta2[x[,'var2'],k]*beta3[x[,'var3'],k]*beta4[x[,'var4'],k]
    #print(blin)
  }

  m <-  mu + b1[x[,'var1']] + b2[x[,'var2']] + b3[x[,'var3']] + b4[x[,'var4']] + blin + blin2 + blin3
  y <-  rnorm(N, m, sigma)

  return(list(mu = mu,
              b1 = b1,
              b2 = b2,
              b3 = b3,
              b4 = b4,
              #b4 = b4,
              blin = blin,
              blin2 = blin2,
              blin3 = blin3,
              y = y,
              betas = list(beta1 = beta1, beta2 = beta2, beta3 = beta3, beta4 = beta4),
              Bv = Bv))
}



# Simulating with Q = 1 -------------------------------------------------------------------

B1 <- 6
B2 <- 4
B3 <- 3
B4 <- 2
lambda <- 10#c(8,10, 12)
sb <- 1
sigma <- 1

set.seed(022)
dat <- genData_3(B1 = B1, B2 = B2, B3 = B3, B4 = B4, sb = sb, sigma = sigma, lambda = lambda)

# run the two models
model_Q1 <- model_jags_8(data = dat,
                         Q = 1,
                         mmu = 10,
                         smu = 1,
                         a = 0.1,
                         b = 0.1,
                         nthin = 2,
                         nburnin = 2000)

model_9_Q1 <- model_jags_9(data = dat,
                           Q = 1,
                           mmu = 10,
                           smu = 1,
                           a = 0.1,
                           b = 0.1,
                           nthin = 2,
                           nburnin = 2000)

plot(model_Q1)
plot(model_9_Q1)

true_blin <- dat$blin + dat$blin2 + dat$blin3
int_Q1 <-  model_9_Q1$BUGSoutput$mean$int_p1 + model_9_Q1$BUGSoutput$mean$int_p2 + model_9_Q1$BUGSoutput$mean$int_p3

caret::RMSE(true_blin, model_Q1$BUGSoutput$mean$int)
caret::RMSE(true_blin, int_Q1)

qplot(true_blin, int_Q1) + geom_abline() + theme_bw() +
  labs(x = 'True interaction', y = 'Predicted', title = 'True interaction fitted') +
  qplot(true_blin, model_Q1$BUGSoutput$mean$int) + geom_abline() + theme_bw() +
  labs(x = 'True interaction', y = '', title = 'Our method')



