# Code to get the results of the paper

library(R2jags)
library(ggplot2)
library(patchwork)
library(reshape2)

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

# Models

# Model 1 - Interaction just between beta1 and beta 2
model_jags_one_way <- function(data, Q, mmu, smu, a, b, stheta = 1, nthin, nburnin){

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
model_jags_one_two_way <- function(data, Q, mmu, smu, a,b, stheta = 1, nthin, nburnin){

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

# Model BAMMIT (V = 2)
model_bammit_V2 <- function(data, Q, mmu, smu, a,b, nthin, nburnin, alpha_beta){

  Y <- data$y
  V <- length(data$Bv)
  Bv <- data$Bv
  N <- length(data$y)



  modelCode <- "
model{
  # Likelihood
  for (n in 1:N) {
    Y[n] ~ dnorm(mu[n], sy^-2)
    mu[n] = muall + b1[var1[n]] + b2[var2[n]] +  int[n]
    int[n] = sum(lambda[1:Q] * beta1[var1[n],1:Q] * beta2[var2[n],1:Q])
  }

  # Priors
  muall ~ dnorm(mmu, smu^-2) # grand mean

  # Main effects
  for(i in 1:B1) {
    b1_aux[i] ~ dnorm(0, sb1^-2) # Prior on effect 1
  }

  m_b1_aux <- sum(b1_aux)/B1
  for(i in 1:B1) {
    b1[i] = b1_aux[i] - m_b1_aux
  }

  for(i in 1:B2) {
    b2_aux[i] ~ dnorm(0, sb2^-2) # Prior on effect 2
  }

  m_b2_aux <- sum(b2_aux)/B2
  for(i in 1:B2) {
    b2[i] = b2_aux[i] - m_b2_aux
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
      beta1[i,q] = (thetaBeta1New[i,q]) / (sqrtThetaBeta1[q])
    }
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
      beta2[i,q] = (thetaBeta2New[i,q]) / (sqrtThetaBeta2[q])
    }
  }

  # Prior on lambda
  for(q in 1:Q) {
    lambda_raw[q] ~ dnorm(0, slambda^-2)T(0,)
  }
  lambda = sort(lambda_raw)

  # Priors
  sb1 ~ dt(0, 1^-2, 1)T(0,)
  sb2 ~ dt(0, 1^-2, 1)T(0,)
  slambda ~ dt(0, 1^-2, 1)T(0,)
  sy ~ dgamma(a, b) # Prior on residual sd
}"

  B1 <- Bv[1]
  B2 <- Bv[2]

  x <- expand.grid(lapply(Bv, function(x) 1:x))
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
    a = a,
    b = b
  )

  # Choose the parameters to watch
  modelParameters <- c(
    "b1", "b2", "lambda", "beta1", "beta2",  "sy", "muall",
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

# Model BAMMIT (V = 3)
model_bammit_V3 <- function(data, Q, mmu, smu, a,b, nthin, nburnin, alpha_beta){

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
    b1_aux[i] ~ dnorm(0, sb1^-2) # Prior on effect 1
  }

  m_b1_aux <- sum(b1_aux)/B1
  for(i in 1:B1) {
    b1[i] = b1_aux[i] - m_b1_aux
  }


  for(i in 1:B2) {
    b2_aux[i] ~ dnorm(0, sb2^-2) # Prior on effect 2
  }

  m_b2_aux <- sum(b2_aux)/B2
  for(i in 1:B2) {
    b2[i] = b2_aux[i] - m_b2_aux
  }


  for(i in 1:B3) {
    b3_aux[i] ~ dnorm(0, sb3^-2) # Prior on effect 2
  }

  m_b3_aux <- sum(b3_aux)/B3
  for(i in 1:B3) {
    b3[i] = b3_aux[i] - m_b3_aux
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

# Model BAMMIT (V = 4)
model_bammit_V4 <- function(data, Q, mmu, smu, a,b, nthin, nburnin, alpha_beta){

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
    b1_aux[i] ~ dnorm(0, sb1^-2) # Prior on effect 1
  }

  m_b1_aux <- sum(b1_aux)/B1
  for(i in 1:B1) {
    b1[i] = b1_aux[i] - m_b1_aux
  }

  for(i in 1:B2) {
    b2_aux[i] ~ dnorm(0, sb2^-2) # Prior on effect 2
  }

  m_b2_aux <- sum(b2_aux)/B2
  for(i in 1:B2) {
    b2[i] = b2_aux[i] - m_b2_aux
  }

  for(i in 1:B3) {
    b3_aux[i] ~ dnorm(0, sb3^-2) # Prior on effect 3
  }

  m_b3_aux <- sum(b3_aux)/B3
  for(i in 1:B3) {
    b3[i] = b3_aux[i] - m_b3_aux
  }

  for(i in 1:B4) {
    b4_aux[i] ~ dnorm(0, sb4^-2) # Prior on effect 4
  }

  m_b4_aux <- sum(b4_aux)/B4
  for(i in 1:B4) {
    b4[i] = b4_aux[i] - m_b4_aux
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
    "b1", "b2", "b3", "b4", "lambda", "beta1", "beta2", "beta3", "beta4", "sy", "muall",
    "int", "mu", "p1", "p2", "p3", "p4"
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

# Model BFM
model_bfm <- function(data, mmu, smu, a,b, nthin, nburnin){

  Y <- data$y
  V <- length(data$Bv)
  Bv <- data$Bv
  N <- length(data$y)


  modelCode <- "
   model
   {
   # Likelihood
    for (n in 1:N) {
     Y[n] ~ dnorm(mu[n], sy)
     mu[n] = muall + b1[var1[n]] + b2[var2[n]] + b3[var3[n]] + int[n]
     int[n] = b1[var1[n]]*b2[var2[n]] + b1[var1[n]]*b3[var3[n]] + b2[var2[n]]*b3[var3[n]] + b1[var1[n]]*b2[var2[n]]*b3[var3[n]]
    }

   # Priors
   # Prior on grand mean
    muall ~ dnorm(mmu, smu^-2)

   # Prior on genotype effect
   for(i in 1:B1) {
   b1[i] ~ dnorm(0, sb1^-2) # Prior on genotype effect
   }

   # Prior on environment effect
   for(j in 1:B2) {
   b2[j] ~ dnorm(0, sb2^-2) # Prior on environment effect
   }

   # Prior on time effect
   for(k in 1:B3) {
   b3[k] ~ dnorm(0, sb3^-2) # Prior on time effect
   }

   #Priors
   sb1 ~ dt(0, 1, 1)T(0,)
   sb2 ~ dt(0, 1, 1)T(0,)
   sb3 ~ dt(0, 1, 1)T(0,)
   slambda ~ dt(0, 1, 1)T(0,)

   # Prior on residual standard deviation
    sy ~ dgamma(a, b) # inverse of sy
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
    "b1", "b2", "b3", "sy", "muall", "int", "mu"
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


# Funtions to generate the data

# Function to generate data from Model 1
genData_V3_one_way <- function(B1, B2, B3, sigma, sb, lambda){ #sb1, sb2, sb3, sb4){

  Bv <- c(B1, B2, B3)
  #Bv <- c(B1, B2)
  N <- Reduce("*", Bv)
  Q <- length(lambda)
  V <- length(Bv)
  mu <- 100

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
              betas = list(beta1 = beta1, beta2 = beta2),#, beta4 = beta4),
              Bv = Bv))

}
# Function to generate data from Model 2
genData_V3_one_two_way <- function(B1, B2, B3, sigma, sb, lambda){
  Bv <- c(B1, B2, B3)
  #Bv <- c(B1, B2)
  N <- Reduce("*", Bv)
  Q <- length(lambda)
  V <- length(Bv)
  mu <- 100

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
              int = blin + blin2,
              y = y,
              betas = list(beta1 = beta1, beta2 = beta2, beta3 = beta3),#, beta4 = beta4),
              Bv = Bv))
}
# Function to generate data from BAMMIT model with V = 2
genData_V2 <- function(B1, B2, sigma, sb, lambda){
  Bv <- c(B1, B2)
  #Bv <- c(B1, B2)
  N <- Reduce("*", Bv)
  Q <- length(lambda)
  V <- length(Bv)
  mu <- 100

  x <- expand.grid(1:Bv[1], 1:Bv[2]) |>
    dplyr::mutate_if(is.numeric,as.factor)

  colnames(x) <- c('var1', 'var2')

  # x <- expand.grid(1:Bv[1], 1:Bv[2]) |>
  #   dplyr::mutate_if(is.numeric,as.factor)
  #
  # colnames(x) <- c('var1', 'var2')


  b1 <- rnorm(B1, 0, sb)
  b1 <- b1 - mean(b1)

  b2 <- rnorm(B2, 0, sb)
  b2 <- b2 - mean(b2)


  beta1 <- genInt(B1, Q)
  beta2 <- genInt(B2, Q)

  blin = rep(0, N)
  for (k in 1:Q) {
    blin <- blin + lambda[k]*beta1[x[,'var1'],k]*beta2[x[,'var2'],k]
  }


  m <-  mu + b1[x[,'var1']] + b2[x[,'var2']] +  blin
  y <-  rnorm(N, m, sigma)

  return(list(mu = mu,
              b1 = b1,
              b2 = b2,
              blin = blin,
              y = y,
              betas = list(beta1 = beta1, beta2 = beta2),#, beta4 = beta4),
              Bv = Bv))
}
# Function to generate data from BAMMIT model with V = 3
genData_V3 <- function(B1, B2, B3, sigma, sb, lambda){
  Bv <- c(B1, B2, B3)
  #Bv <- c(B1, B2)
  N <- Reduce("*", Bv)
  Q <- length(lambda)
  V <- length(Bv)
  mu <- 100

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
# Function to generate data from BAMMIT model with V = 4
genData_V4 <- function(B1, B2, B3, B4, sigma, sb, lambda){
  Bv <- c(B1, B2, B3, B4)
  #Bv <- c(B1, B2)
  N <- Reduce("*", Bv)
  Q <- length(lambda)
  V <- length(Bv)
  mu <- 100

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
  b4 <- b3 - mean(b4)


  beta1 <- genInt(B1, Q)
  beta2 <- genInt(B2, Q)
  beta3 <- genInt(B3, Q)
  beta4 <- genInt(B4, Q)

  blin = rep(0, N)
  for (k in 1:Q) {
    blin <- blin + lambda[k]*beta1[x[,'var1'],k]*beta2[x[,'var2'],k]*beta3[x[,'var3'],k]*beta4[x[,'var4'],k]
    #print(blin)
  }

  m <-  mu + b1[x[,'var1']] + b2[x[,'var2']] + b3[x[,'var3']] + b4[x[,'var4']] + blin
  y <-  rnorm(N, m, sigma)

  return(list(mu = mu,
              b1 = b1,
              b2 = b2,
              b3 = b3,
              b4 = b4,
              blin = blin,
              y = y,
              betas = list(beta1 = beta1, beta2 = beta2, beta3 = beta3, beta4 = beta4),
              Bv = Bv))
}


# Generate the results


# Comp times --------------------------------------------------------------

B1 <- 12
B2 <- 10
B3 <- 4
B4 <- 2
lambda <- 10 #c(8,10, 12)
sb <- 1
sigma <- 1

set.seed(022)
dat_V2 <- genData_V2(B1 = B1, B2 = B2, sb = sb, sigma = sigma, lambda = lambda)
dat_V3 <- genData_V3(B1 = B1, B2 = B2, B3 = B3, sb = sb, sigma = sigma, lambda = lambda)
dat_V4 <- genData_V4(B1 = B1, B2 = B2, B3 = B3, B4 = B4, sb = sb, sigma = sigma, lambda = lambda)

model_V4_Q1 <- model_bammit_V4(data = dat_V4,
                               Q = 1,
                               mmu = 100,
                               smu = 1,
                               a = 0.1,
                               b = 0.1,
                               nthin = 2,
                               nburnin = 2000)

plot(model_V4_Q1)
qplot(dat_V4$y, model_V4_Q1$BUGSoutput$mean$mu) + geom_abline()
qplot(dat_V4$blin, model_V4_Q1$BUGSoutput$mean$int) + geom_abline()
qplot(dat_V4$b1, model_V4_Q1$BUGSoutput$mean$b1) + geom_abline()
qplot(dat_V4$b2, model_V4_Q1$BUGSoutput$mean$b2) + geom_abline()
qplot(dat_V4$b3, model_V4_Q1$BUGSoutput$mean$b3) + geom_abline()
qplot(dat_V4$b4, model_V4_Q1$BUGSoutput$mean$b4) + geom_abline()

complexity <- microbenchmark::microbenchmark(m1 = model_bammit_V2(data = dat_V2,
                                                                  Q = 1,
                                                                  mmu = 100,
                                                                  smu = 1,
                                                                  a = 0.1,
                                                                  b = 0.1,
                                                                  nthin = 2,
                                                                  nburnin = 2000),
                                             m2 = model_bammit_V2(data = dat_V2,
                                                                  Q = 2,
                                                                  mmu = 100,
                                                                  smu = 1,
                                                                  a = 0.1,
                                                                  b = 0.1,
                                                                  nthin = 2,
                                                                  nburnin = 2000),
                                             m3 = model_bammit_V2(data = dat_V2,
                                                                  Q = 3,
                                                                  mmu = 100,
                                                                  smu = 1,
                                                                  a = 0.1,
                                                                  b = 0.1,
                                                                  nthin = 2,
                                                                  nburnin = 2000),
                                             m4 = model_bammit_V3(data = dat_V3,
                                                                  Q = 1,
                                                                  mmu = 100,
                                                                  smu = 1,
                                                                  a = 0.1,
                                                                  b = 0.1,
                                                                  nthin = 2,
                                                                  nburnin = 2000),
                                             m5 = model_bammit_V3(data = dat_V3,
                                                                  Q = 2,
                                                                  mmu = 100,
                                                                  smu = 1,
                                                                  a = 0.1,
                                                                  b = 0.1,
                                                                  nthin = 2,
                                                                  nburnin = 2000),
                                             m6 = model_bammit_V3(data = dat_V3,
                                                                  Q = 3,
                                                                  mmu = 100,
                                                                  smu = 1,
                                                                  a = 0.1,
                                                                  b = 0.1,
                                                                  nthin = 2,
                                                                  nburnin = 2000),
                                             m7 = model_bammit_V4(data = dat_V4,
                                                                  Q = 1,
                                                                  mmu = 100,
                                                                  smu = 1,
                                                                  a = 0.1,
                                                                  b = 0.1,
                                                                  nthin = 2,
                                                                  nburnin = 2000),
                                             m8 = model_bammit_V4(data = dat_V4,
                                                                  Q = 2,
                                                                  mmu = 100,
                                                                  smu = 1,
                                                                  a = 0.1,
                                                                  b = 0.1,
                                                                  nthin = 2,
                                                                  nburnin = 2000),
                                             m9 = model_bammit_V4(data = dat_V4,
                                                                  Q = 3,
                                                                  mmu = 100,
                                                                  smu = 1,
                                                                  a = 0.1,
                                                                  b = 0.1,
                                                                  nthin = 2,
                                                                  nburnin = 2000),
                                             times = 1)

df_comp <- data.frame(V = as.factor(c(2,2,2,3,3,3,4,4,4)),
                      Q = as.numeric(c(1,2,3,1,2,3,1,2,3)),
                      time = c(37.65172, 119.49784, 119.49784,
                               186.98443, 602.86332, 1046.84481,
                               404.10346, 1607.02512, 2026.73302)/60)

df_comp |> ggplot(aes(x = Q, y = time, colour = V)) +
  geom_point(size = 3) +
  geom_line(linewidth = 1) +
  ylab('Time (minutes)') +
  theme_bw()

plot(complexity, type = 'p')


# save(complexity, file = "Running models/complexity.RData")



# Figure 1 ----------------------------------------------------------------

B1 <- 12
B2 <- 10
B3 <- 4
lambda <- 10 #c(8,10, 12)
sb <- 1
sigma <- 1.5

set.seed(022)
dat_V3_Q1_one_way <- genData_V3_one_way(B1 = B1, B2 = B2, B3 = B3,  sb = sb, sigma = sigma, lambda = lambda)


model_1_F1_Q1 <- model_jags_one_way(data = dat_V3_Q1_one_way,
                                    Q = 1,
                                    mmu = 10,
                                    smu = 1,
                                    a = 0.1,
                                    b = 0.1,
                                    nthin = 2,
                                    nburnin = 2000)

model_7_F1_Q1 <- model_bammit_V3(data = dat_V3_Q1_one_way,
                                 Q = 1,
                                 mmu = 10,
                                 smu = 1,
                                 a = 0.1,
                                 b = 0.1,
                                 nthin = 2,
                                 nburnin = 2000)


plot(model_7_F1_Q1)
qplot(dat_V3_Q1_one_way$y, model_7_F1_Q1$BUGSoutput$mean$mu) + geom_abline()
qplot(dat_V3_Q1_one_way$blin, model_7_F1_Q1$BUGSoutput$mean$int) + geom_abline()
qplot(dat_V3_Q1_one_way$b1, model_7_F1_Q1$BUGSoutput$mean$b1) + geom_abline()
qplot(dat_V3_Q1_one_way$b2, model_7_F1_Q1$BUGSoutput$mean$b2) + geom_abline()
qplot(dat_V3_Q1_one_way$b3, model_7_F1_Q1$BUGSoutput$mean$b3) + geom_abline()

# fig1_2_way <- model_1_F1_Q1
# fig2_bammit <- model_7_F1_Q1
#
# save(fig1_2_way, file = "Running models/fig1_2_way.RData")
# save(fig2_bammit, file = "Running models/fig2_bammit.RData")

# Create a data frame
fig1_df <- data.frame(
  model = c(rep("2-way interaction", length(model_1_F1_Q1$BUGSoutput$mean$int)),
            rep("BAMMIT", length(model_7_F1_Q1$BUGSoutput$mean$int))),
  est = c(model_1_F1_Q1$BUGSoutput$mean$int, model_7_F1_Q1$BUGSoutput$mean$int),
  q1 = c(apply(model_1_F1_Q1$BUGSoutput$sims.list$int, 2, function(x) quantile(x, 0.025)),
         apply(model_7_F1_Q1$BUGSoutput$sims.list$int, 2, function(x) quantile(x, 0.025))),
  q2 = c(apply(model_1_F1_Q1$BUGSoutput$sims.list$int, 2, function(x) quantile(x, 0.975)),
         apply(model_7_F1_Q1$BUGSoutput$sims.list$int, 2, function(x) quantile(x, 0.975))),
  true = c(dat_V3_Q1_one_way$blin, dat_V3_Q1_one_way$blin)
)


fig1_df |>
  ggplot(aes(x = true, y = est)) +
  geom_abline(size = 0.3, col = "gray") +
  geom_linerange(aes(ymin =  q1, ymax = q2), alpha = 0.5, size = 0.4) +
  geom_point(colour =  "steelblue", size = 1) +
  facet_wrap(~model) +
  labs(x = "True", y = "Estimated") +
  theme_bw(base_size = 16) +
  theme(strip.background = element_blank(),
        panel.spacing.x = unit(5,"mm"),
        legend.position = c(0.75, 0.04))



# Figure 2 - behavior of the BAMMIT interaction when increase the --------

# Simulating with Q = 1 ---

B1 <- 10
B2 <- 5
B3 <- 2
lambda <- 10 #c(8,10, 12)
sb <- 1
sigma <- 1

set.seed(022)
dat <- genData_2(B1 = B1, B2 = B2, B3 = B3,  sb = sb, sigma = sigma, lambda = lambda)


# run the two models

model_7_Q1 <- model_bammit_V3(data = dat,
                           Q = 1,
                           mmu = 10,
                           smu = 1,
                           a = 0.1,
                           b = 0.1,
                           nthin = 2,
                           nburnin = 2000)


qplot(dat$y, model_7_Q1$BUGSoutput$median$mu) +
  geom_abline() +
  geom_linerange(aes(ymin =  apply(model_7_Q1$BUGSoutput$sims.list$mu, 2, function(x) quantile(x, 0.05)),
                     ymax = apply(model_7_Q1$BUGSoutput$sims.list$mu, 2, function(x) quantile(x, 0.95))), alpha = 0.5, size = 0.4) +
  theme_bw() +
  labs(x = 'True', y = 'Estimated')


model_7_Q2 <- model_bammit(data = dat,
                           Q = 2,
                           mmu = 10,
                           smu = 1,
                           a = 0.1,
                           b = 0.1,
                           nthin = 2,
                           nburnin = 2000)

model_7_Q3 <- model_bammit(data = dat,
                           Q = 3,
                           mmu = 10,
                           smu = 1,
                           a = 0.1,
                           b = 0.1,
                           nthin = 2,
                           nburnin = 2000)

model_7_Q4 <- model_bammit(data = dat,
                           Q = 4,
                           mmu = 10,
                           smu = 1,
                           a = 0.1,
                           b = 0.1,
                           nthin = 2,
                           nburnin = 2000)


fig2_bammit_Q1 <- model_7_Q1
fig2_bammit_Q2 <- model_7_Q2
fig2_bammit_Q3 <- model_7_Q3
fig2_bammit_Q4 <- model_7_Q4

# save(fig2_bammit_Q1, file = "Running models/fig2_bammit_Q1.RData")
# save(fig2_bammit_Q2, file = "Running models/fig2_bammit_Q4.RData")
# save(fig2_bammit_Q3, file = "Running models/fig2_bammit_Q3.RData")
# save(fig2_bammit_Q4, file = "Running models/fig2_bammit_Q4.RData")
#
#
# caret::RMSE(dat$blin + dat$blin2, model_7_Q1$BUGSoutput$mean$int)
# caret::RMSE(dat$blin + dat$blin2, model_7_Q2$BUGSoutput$mean$int)
# caret::RMSE(dat$blin + dat$blin2, model_7_Q3$BUGSoutput$mean$int)
# caret::RMSE(dat$blin + dat$blin2, model_7_Q4$BUGSoutput$mean$int)
#
# qplot(dat$blin + dat$blin2, model_7_Q1$BUGSoutput$mean$int) + geom_abline() +
# qplot(dat$blin + dat$blin2, model_7_Q2$BUGSoutput$mean$int) + geom_abline() +
# qplot(dat$blin + dat$blin2, model_7_Q3$BUGSoutput$mean$int) + geom_abline() +
# qplot(dat$blin + dat$blin2, model_7_Q4$BUGSoutput$mean$int) + geom_abline()

# Create a data frame
fig2_df <- data.frame(
  model = c(rep("Q = 1", length(model_7_Q1$BUGSoutput$mean$int)),
            rep("Q = 2", length(model_7_Q2$BUGSoutput$mean$int)),
            rep("Q = 3", length(model_7_Q3$BUGSoutput$mean$int)),
            rep("Q = 4", length(model_7_Q4$BUGSoutput$mean$int))),
  est = c(model_7_Q1$BUGSoutput$mean$int,
          model_7_Q2$BUGSoutput$mean$int,
          model_7_Q3$BUGSoutput$mean$int,
          model_7_Q4$BUGSoutput$mean$int),
  q1 = c(apply(model_7_Q1$BUGSoutput$sims.list$int, 2, function(x) quantile(x, 0.025)),
         apply(model_7_Q2$BUGSoutput$sims.list$int, 2, function(x) quantile(x, 0.025)),
         apply(model_7_Q3$BUGSoutput$sims.list$int, 2, function(x) quantile(x, 0.025)),
         apply(model_7_Q4$BUGSoutput$sims.list$int, 2, function(x) quantile(x, 0.025))),
  q2 = c(apply(model_7_Q1$BUGSoutput$sims.list$int, 2, function(x) quantile(x, 0.975)),
         apply(model_7_Q2$BUGSoutput$sims.list$int, 2, function(x) quantile(x, 0.975)),
         apply(model_7_Q3$BUGSoutput$sims.list$int, 2, function(x) quantile(x, 0.975)),
         apply(model_7_Q4$BUGSoutput$sims.list$int, 2, function(x) quantile(x, 0.975))),
  true = c(dat$blin + dat$blin2, dat$blin + dat$blin2, dat$blin + dat$blin2, dat$blin + dat$blin2)
)



fig2_df |>
  ggplot(aes(x = true, y = est)) +
  geom_abline(size = 0.3, col = "gray") +
  geom_linerange(aes(ymin =  q1, ymax = q2), alpha = 0.5, size = 0.4) +
  geom_point(colour =  "steelblue", size = 1) +
  facet_wrap(~model, nrow = 1) +
  labs(x = "True", y = "Estimated") +
  theme_bw(base_size = 16) +
  theme(strip.background = element_blank(),
        panel.spacing.x = unit(5,"mm"),
        legend.position = c(0.75, 0.04))




# Table 2 -----------------------------------------------------------------

# RMSE of the BAMMIT interaction and the 2-way+3-way interaction when increase
# the number of components (Q) and simulate with different values of Q.

B1 <- 12
B2 <- 10
B3 <- 4
sb <- 1
sigma <- 1.5


# Generate data with Q = 1

lambda <- 10 #c(8,10, 12)

set.seed(022)
dat_V3_Q1_one_two_way <- genData_V3_one_two_way(B1 = B1, B2 = B2, B3 = B3,  sb = sb, sigma = sigma, lambda = lambda)

model_V3_Q1_Qsim1 <- model_bammit_V3(data = dat_V3_Q1_one_two_way,
                                     Q = 1,
                                     mmu = 100,
                                     smu = 1,
                                     a = 0.1,
                                     b = 0.1,
                                     nthin = 2,
                                     nburnin = 2000)

model_V3_Q2_Qsim1 <- model_bammit_V3(data = dat_V3_Q1_one_two_way,
                                     Q = 2,
                                     mmu = 100,
                                     smu = 1,
                                     a = 0.1,
                                     b = 0.1,
                                     nthin = 2,
                                     nburnin = 2000)

model_V3_Q4_Qsim1 <- model_bammit_V3(data = dat_V3_Q1_one_two_way,
                                     Q = 4,
                                     mmu = 100,
                                     smu = 1,
                                     a = 0.1,
                                     b = 0.1,
                                     nthin = 2,
                                     nburnin = 2000)

model_V3_Q6_Qsim1 <- model_bammit_V3(data = dat_V3_Q1_one_two_way,
                                     Q = 6,
                                     mmu = 100,
                                     smu = 1,
                                     a = 0.1,
                                     b = 0.1,
                                     nthin = 2,
                                     nburnin = 2000)


set.seed(022)
dat_V3_Q1_one_two_way_test <- genData_V3_one_two_way(B1 = B1, B2 = B2, B3 = B3,  sb = sb, sigma = sigma, lambda = lambda)
caret::RMSE(model_V3_Q1_Qsim1$BUGSoutput$mean$mu, dat_V3_Q1_one_two_way_test$y)
caret::RMSE(model_V3_Q2_Qsim1$BUGSoutput$mean$mu, dat_V3_Q1_one_two_way_test$y)
caret::RMSE(model_V3_Q4_Qsim1$BUGSoutput$mean$mu, dat_V3_Q1_one_two_way_test$y)
caret::RMSE(model_V3_Q6_Qsim1$BUGSoutput$mean$int, dat_V3_Q1_one_two_way$int)



# Generate data with Q = 2

lambda <- c(8,10)

set.seed(022)
dat_V3_Q2_one_two_way <- genData_V3_one_two_way(B1 = B1, B2 = B2, B3 = B3,  sb = sb, sigma = sigma, lambda = lambda)

model_V3_Q1_Qsim2 <- model_bammit_V3(data = dat_V3_Q2_one_two_way,
                                     Q = 1,
                                     mmu = 100,
                                     smu = 1,
                                     a = 0.1,
                                     b = 0.1,
                                     nthin = 2,
                                     nburnin = 2000)

model_V3_Q2_Qsim2 <- model_bammit_V3(data = dat_V3_Q2_one_two_way,
                                     Q = 2,
                                     mmu = 100,
                                     smu = 1,
                                     a = 0.1,
                                     b = 0.1,
                                     nthin = 2,
                                     nburnin = 2000)

model_V3_Q4_Qsim2 <- model_bammit_V3(data = dat_V3_Q2_one_two_way,
                                     Q = 4,
                                     mmu = 100,
                                     smu = 1,
                                     a = 0.1,
                                     b = 0.1,
                                     nthin = 2,
                                     nburnin = 2000)

model_V3_Q6_Qsim2 <- model_bammit_V3(data = dat_V3_Q2_one_two_way,
                                     Q = 6,
                                     mmu = 100,
                                     smu = 1,
                                     a = 0.1,
                                     b = 0.1,
                                     nthin = 2,
                                     nburnin = 2000)

set.seed(023)
dat_V3_Q2_one_two_way_test <- genData_V3_one_two_way(B1 = B1, B2 = B2, B3 = B3,  sb = sb, sigma = sigma, lambda = lambda)
caret::RMSE(model_V3_Q1_Qsim2$BUGSoutput$mean$int, dat_V3_Q2_one_two_way_test$int)
caret::RMSE(model_V3_Q2_Qsim2$BUGSoutput$mean$mu, dat_V3_Q2_one_two_way_test$y)
caret::RMSE(model_V3_Q4_Qsim2$BUGSoutput$mean$mu, dat_V3_Q2_one_two_way_test$y)
caret::RMSE(model_V3_Q6_Qsim2$BUGSoutput$mean$mu, dat_V3_Q2_one_two_way_test$y)




# Generate data with Q = 3

lambda <- c(8,10,12)

set.seed(022)
dat_V3_Q3_one_two_way <- genData_V3_one_two_way(B1 = B1, B2 = B2, B3 = B3,  sb = sb, sigma = sigma, lambda = lambda)

model_V3_Q1_Qsim3 <- model_bammit_V3(data = dat_V3_Q3_one_two_way,
                                     Q = 1,
                                     mmu = 100,
                                     smu = 1,
                                     a = 0.1,
                                     b = 0.1,
                                     nthin = 2,
                                     nburnin = 2000)

model_V3_Q2_Qsim3 <- model_bammit_V3(data = dat_V3_Q3_one_two_way,
                                     Q = 2,
                                     mmu = 100,
                                     smu = 1,
                                     a = 0.1,
                                     b = 0.1,
                                     nthin = 2,
                                     nburnin = 2000)

model_V3_Q4_Qsim3 <- model_bammit_V3(data = dat_V3_Q3_one_two_way,
                                     Q = 4,
                                     mmu = 100,
                                     smu = 1,
                                     a = 0.1,
                                     b = 0.1,
                                     nthin = 2,
                                     nburnin = 2000)

model_V3_Q6_Qsim3 <- model_bammit_V3(data = dat_V3_Q3_one_two_way,
                                     Q = 6,
                                     mmu = 100,
                                     smu = 1,
                                     a = 0.1,
                                     b = 0.1,
                                     nthin = 2,
                                     nburnin = 2000)




set.seed(023)
dat_V3_Q3_one_two_way_test <- genData_V3_one_two_way(B1 = B1, B2 = B2, B3 = B3,  sb = sb, sigma = sigma, lambda = lambda)
caret::RMSE(model_V3_Q1_Qsim3$BUGSoutput$mean$int, dat_V3_Q3_one_two_way$int)
caret::RMSE(model_V3_Q2_Qsim3$BUGSoutput$mean$mu, dat_V3_Q3_one_two_way_test$y)
caret::RMSE(model_V3_Q4_Qsim3$BUGSoutput$mean$mu, dat_V3_Q3_one_two_way_test$y)
caret::RMSE(model_V3_Q6_Qsim3$BUGSoutput$mean$int, dat_V3_Q3_one_two_way$int)


  # Figure 3 ----------------------------------------------------------------



# Figure 4 ----------------------------------------------------------------



# Table 3 -----------------------------------------------------------------



# Table 4 -----------------------------------------------------------------




# Figure 5 ----------------------------------------------------------------




model <- model_Q4


plot(model$BUGSoutput$mean$p1)
plot(model$BUGSoutput$mean$p2)
plot(model$BUGSoutput$mean$p3)
plot(model$BUGSoutput$mean$p4)


df_gammas <- data.frame(Q = c(1,2,3,4),
                        gamma1 = model$BUGSoutput$mean$p1,
                        gamma2 = model$BUGSoutput$mean$p2,
                        gamma3 = model$BUGSoutput$mean$p3)

df_gammas <- reshape2::melt(df_gammas, id = 'Q')


gammas_ic <- list(
  gamma1 = data.frame(Q = c(1,2,3,4),
                      q1 = apply(model$BUGSoutput$sims.list$p1, 2, function(x) quantile(x, 0.05)),
                      q2 = apply(model$BUGSoutput$sims.list$p1, 2, function(x) quantile(x, 0.95))),
  gamma2 = data.frame(Q = c(1,2,3,4),
                      q1 = apply(model$BUGSoutput$sims.list$p2, 2, function(x) quantile(x, 0.05)),
                      q2 = apply(model$BUGSoutput$sims.list$p2, 2, function(x) quantile(x, 0.95))),
  gamma3 = data.frame(Q = c(1,2,3,4),
                      q1 = apply(model$BUGSoutput$sims.list$p3, 2, function(x) quantile(x, 0.05)),
                      q2 = apply(model$BUGSoutput$sims.list$p3, 2, function(x) quantile(x, 0.95)))
) |> plyr::ldply()



gammas_df <- data.frame(Q = df_gammas$Q, gammas = df_gammas$variable, true = df_gammas$value,
                        q1 = gammas_ic$q1, q2 = gammas_ic$q2)
#
# gammas_df$variable <- factor(df_gammas$variable,
#                              labels = c(bquote(p^(1)),
#                                         bquote(p^(2)),
#                                         bquote(p^(3)),
#                                         bquote(p^(4))))


gammas_df$variable <- factor(df_gammas$variable,
                             labels = c(bquote(Genotype),
                                        bquote(Environment),
                                        bquote(Year)))


gammas_df |>
  ggplot(aes(x = Q, y = true)) +
  geom_linerange(aes(ymin =  q1, ymax = q2), alpha = 1, size = 0.5) +
  geom_point(colour =  "steelblue", size = 3) +
  facet_wrap(~variable, nrow = 1,
             labeller = label_parsed) +
  labs(x = "q", y = bquote(hat(p)[q])) +
  theme_bw(base_size = 16) +
  theme(strip.background = element_blank(),
        panel.spacing.x = unit(5,"mm"),
        legend.position = c(0.75, 0.04))





# Old stuff  (I need to delete it) ---------------------------------------------------------------




plot(model_V3_Q1_Qsim1)


set.seed(023)
dat_V3_Q1_one_two_way_test <- genData_V3_one_two_way(B1 = B1, B2 = B2, B3 = B3,  sb = sb, sigma = sigma, lambda = lambda)

caret::RMSE(dat_V3_Q1_one_two_way_test$y, model_V3_Q1_Qsim1$BUGSoutput$median$mu)
caret::RMSE(dat_V3_Q1_one_two_way$int, model_V3_Q1_Qsim1$BUGSoutput$median$int)
caret::RMSE(dat_V3_Q1_one_two_way$blin + dat_V3_Q1_one_two_way$blin2, model_V3_Q6_Qsim1$BUGSoutput$median$int)

qplot(dat_V3_Q1_one_two_way_test$y, model_V3_Q1_Qsim1$BUGSoutput$median$mu) + geom_abline()
caret::RMSE(dat_V3_Q1_one_two_way$y, model_V3_Q1_Qsim1$BUGSoutput$median$mu)



# plot with intervals

V3_Q1_Dataint <-  list(Q1 = dat$blin + dat$blin2,
                       Q2 =  dat$blin + dat$blin2,
                       Q3 = dat$blin + dat$blin2,
                       Q4 = dat$blin + dat$blin2) |>
  melt() |>
  select(value, L1)

V3_Q1_int <- list(Q1 = model_7_Q1$BUGSoutput$mean$int,
                  Q2 = model_7_Q2$BUGSoutput$mean$int,
                  Q3 = model_7_Q3$BUGSoutput$mean$int,
                  Q4 = model_7_Q4$BUGSoutput$mean$int) |>
  melt() |>
  select(value, L1)

quant_V3_Q1_int <- list(
  Q1 = data.frame(q1 = apply(model_7_Q1$BUGSoutput$sims.list$int, 2, function(x) quantile(x, 0.05)),
                    q2 = apply(model_7_Q1$BUGSoutput$sims.list$int, 2, function(x) quantile(x, 0.95))),
  Q2 = data.frame(q1 = apply(model_7_Q2$BUGSoutput$sims.list$int, 2, function(x) quantile(x, 0.05)),
                      q2 = apply(model_7_Q2$BUGSoutput$sims.list$int, 2, function(x) quantile(x, 0.95))),
  Q3 = data.frame(q1 = apply(model_7_Q3$BUGSoutput$sims.list$int, 2, function(x) quantile(x, 0.05)),
                 q2 = apply(model_7_Q3$BUGSoutput$sims.list$int, 2, function(x) quantile(x, 0.95))),
  Q4 = data.frame(q1 = apply(model_7_Q4$BUGSoutput$sims.list$int, 2, function(x) quantile(x, 0.05)),
                  q2 = apply(model_7_Q4$BUGSoutput$sims.list$int, 2, function(x) quantile(x, 0.95)))
) |> plyr::ldply()


V3_Q1_int <- cbind(V3_Q1_int, quant_V3_Q1_int$q1, quant_V3_Q1_int$q2)
colnames(V3_Q1_int) <- c("value", "L1", "q1", "q2")

df_v3_Q1_int <- data.frame(v = V3_Q1_int$L1, est = V3_Q1_int$value, true = V3_Q1_int$value,
                           q1 = V3_Q1_int$q1, q2 = V3_Q1_int$q2)

lab <- function(x) c("Q = 1", "Q = 2", "Q = 3", "Q = 4")
df_v3_Q1_int |>
  ggplot(aes(x = true, y = est)) +
  geom_abline(size = 0.3, col = "gray") +
  geom_linerange(aes(ymin =  q1, ymax = q2), alpha = 0.5, size = 0.4) +
  geom_point(colour =  "steelblue", size = 1) +
  facet_wrap(~v, scales = 'free',
             labeller = as_labeller(lab)) +
  labs(x = "True", y = "Estimated") +
  theme_bw(base_size = 16) +
  theme(strip.background = element_blank(),
        panel.spacing.x = unit(5,"mm"),
        legend.position = c(0.75, 0.04))



# Figure 3 - Main effects with V = 4


# Table 2

# V = 3

B1 <- 12
B2 <- 10
B3 <- 4
lambda <- c(8,10)
sb <- 1.5
sigma <- 1.5

set.seed(022)
dat_train_V3_Q2 <- genData_1(B1 = B1, B2 = B2, B3 = B3,  sb = sb, sigma = sigma, lambda = lambda)
set.seed(023)
dat_test_V3_Q2 <- genData_1(B1 = B1, B2 = B2, B3 = B3,  sb = sb, sigma = sigma, lambda = lambda)

# Run models

bammit_V3_Q2 <- model_bammit(data = dat_train_V3_Q2,
                             Q = 2,
                             mmu = 10,
                             smu = 1,
                             a = 0.1,
                             b = 0.1,
                             nthin = 2,
                             nburnin = 2000)


bfm_V3_Q2 <- model_bfm(data = dat_train_V3_Q2,
                       mmu = 10,
                       smu = 1,
                       a = 0.1,
                       b = 0.1,
                       nthin = 2,
                       nburnin = 2000)




pred_bammit <- function(model, data){
  Bv <- data$Bv
  x <- expand.grid(1:Bv[1], 1:Bv[2], 1:Bv[3]) |>
    dplyr::mutate_if(is.numeric,as.factor)
  colnames(x) <- paste0('var', 1:3)
  var1 <- x$var1
  var2 <- x$var2
  var3 <- x$var3

  muhat <- model$mean$muall
  b1 <- model$mean$b1
  b2 <- model$mean$b2
  b3 <- model$mean$b3
  inthat <- model$mean$int

  N <- length(data$y)
  yhat <- rep(muhat, N) + b1[x$var1] + b2[x$var2] + b3[x$var3] +  inthat
  return(yhat)
}



caret::RMSE(bammit_V3_Q2$BUGSoutput$mean$mu, dat_test_V3_Q2$y)
caret::RMSE(yhat, dat_test_V3_Q2$y)
caret::RMSE(bfm_V3_Q2$BUGSoutput$mean$mu, dat_test_V3_Q2$y)


# V = 4



# Real data

library(caret)

predictionBAMMITReal <- function(model, data) {
  genotype <- data$Genotype
  environment <- data$Environment
  time <- data$Year
  bloc <- data$Bloc

  muhat <- model$mean$muall
  g <- model$mean$b1
  e <- model$mean$b2
  t <- model$mean$b3
  bl <- model$mean$b4
  inthat <- model$mean$int

  N <- length(data$Yield)
  yhat <- rep(muhat, N) + g[genotype] + e[environment] + t[time] + bl[bloc]  +  inthat

  return(yhat)
}


load("~/Documents/GitHub/bammit/Real data/train_ireland.RData")
load("~/Documents/GitHub/bammit/Real data/test_ireland.RData")

train$Bloc <- gsub('.*_(\\d)$', 'b\\1', train$Bloc) |>  as.factor()
test$Bloc <- gsub('.*_(\\d)$', 'b\\1', test$Bloc) |>  as.factor()
train$Bloc|>  levels()
test$Bloc |> levels()

load("~/Documents/GitHub/bammit/Running models/model_Q6.RData")
plot(model_Q6)
plot(model_Q6$BUGSoutput$mean$M1)
plot(model_Q6$BUGSoutput$mean$M2)
plot(model_Q6$BUGSoutput$mean$M3)
plot(model_Q6$BUGSoutput$mean$M4)

df_gammas <- data.frame(Q = c(1,2,3,4,5,6),
                        gamma1 = 1-model_Q6$BUGSoutput$mean$M1,
                        gamma2 = 1-model_Q6$BUGSoutput$mean$M2,
                        gamma3 = 1-model_Q6$BUGSoutput$mean$M3,
                        gamma4 = 1-model_Q6$BUGSoutput$mean$M4)


df_gammas <- reshape2::melt(df_gammas, id = 'Q')

df_gammas |> ggplot(aes(x = Q, y = value, colour = variable)) +
  geom_point()



labs <- function(x) c(bquote(gamma[1]),
                     bquote(gamma[2]),
                     bquote(gamma[3]),
                     bquote(gamma[4]))

df_gammas$variable <- factor(df_gammas$variable,
                             labels = c(bquote(gamma^(1)),
                                        bquote(gamma^(2)),
                                        bquote(gamma^(3)),
                                        bquote(gamma^(4))))

df_gammas |> ggplot(aes(x = Q, y = value)) +
  geom_point(size = 3, colour = 'black')+
 facet_wrap(~variable, nrow = 1,
            labeller = label_parsed) +
  theme_bw(base_size = 16) +
  ylab('Estimated value')+
  theme(strip.background = element_blank(),
        panel.spacing.x = unit(5,"mm"),
        legend.position = c(0.75, 0.04))


# Generate stratified sample for validation set
set.seed(023)
validation_indices <- createDataPartition(test$Yield, p = 0.5, list = FALSE)
val_set <- test[validation_indices, ]
val_set <- subset(test, Bloc %in% c('b3'))

# The rest of the data goes into the test set
test_set <- test[-validation_indices, ]
test_set <- subset(test, Bloc %in% c('b4'))

predicted_yield = predictionBAMMITReal(model_Q6$BUGSoutput, val_set)


model_Q6$BUGSoutput$mean$b4
beta1 <- model_Q6$BUGSoutput$mean$beta1
beta2 <- model_Q6$BUGSoutput$mean$beta2
beta3 <- model_Q6$BUGSoutput$mean$beta3
beta4 <- model_Q6$BUGSoutput$mean$beta4
lambda <- rnorm(6,0,1)


Q <- Q
Y <- data$Yield
N <- length(data$Yield)
B1 <- length(levels(unique(data$Genotype)))
B2 <- length(levels(unique(data$Environment)))
B3 <- length(levels(unique(data$Year)))
B4 <- length(levels(unique(data$Bloc)))

var1 <- data$Genotype
var2 <- data$Environment
var3 <- data$Year
var4 <- data$Bloc

# int <- c()
# for(n in 1:N){
#   int[n] <- sum(lambda[1:Q]*beta1[var1[n],1:Q]*beta2[var2[n],1:Q]*beta3[var3[n], 1:Q]*beta4[var4[n], 1:Q])
# }

int = rep(0, N)
for (k in 1:Q) {
  int <- int + lambda[k]*beta1[var1,k]*beta2[var2,k]*beta3[var3,k]*beta4[var4,k]
}


length(int)

muhat <- model$mean$muall
g <- model$mean$b1
e <- model$mean$b2
t <- model$mean$b3
bl <- model$mean$b4

yhat <- rep(muhat, N) + g[var1] + e[var2] + t[var3] + bl[var4]  +  int
yhat

caret::RMSE(yhat, test$Yield)



model_aov <- aov(Yield~ Genotype*Year*Environment + Error(Environment:Bloc), data = train)

library(lme4)
model_flmm <- lmer(Yield ~ Genotype * Environment * Year + (1|Bloc), data = train)
predictions <- predict(model_flmm, newdata = test)#, re.form = NA) # The re.form = NA argument is used to exclude the random effects (in this case, the block effect) from the predictions.




## with the new data
library(dplyr)
library(reshape2)

model_Q6_test <- model_Q6


df_gammas <- data.frame(Q = c(1,2,3,4,5,6),
                        gamma1 = model_Q6_test$BUGSoutput$mean$M1,
                        gamma2 = model_Q6_test$BUGSoutput$mean$M2,
                        gamma3 = model_Q6_test$BUGSoutput$mean$M3,
                        gamma4 = model_Q6_test$BUGSoutput$mean$M4)

df_gammas <- reshape2::melt(df_gammas, id = 'Q')


gammas_ic <- list(
  gamma1 = data.frame(Q = c(1,2,3,4,5,6),
                    q1 = apply(model_Q6_test$BUGSoutput$sims.list$M1, 2, function(x) quantile(x, 0.05)),
                    q2 = apply(model_Q6_test$BUGSoutput$sims.list$M1, 2, function(x) quantile(x, 0.95))),
  gamma2 = data.frame(Q = c(1,2,3,4,5,6),
                      q1 = apply(model_Q6_test$BUGSoutput$sims.list$M2, 2, function(x) quantile(x, 0.05)),
                      q2 = apply(model_Q6_test$BUGSoutput$sims.list$M2, 2, function(x) quantile(x, 0.95))),
  gamma3 = data.frame(Q = c(1,2,3,4,5,6),
                      q1 = apply(model_Q6_test$BUGSoutput$sims.list$M3, 2, function(x) quantile(x, 0.05)),
                      q2 = apply(model_Q6_test$BUGSoutput$sims.list$M3, 2, function(x) quantile(x, 0.95))),
  gamma4 = data.frame(Q = c(1,2,3,4,5,6),
                      q1 = apply(model_Q6_test$BUGSoutput$sims.list$M4, 2, function(x) quantile(x, 0.05)),
                      q2 = apply(model_Q6_test$BUGSoutput$sims.list$M4, 2, function(x) quantile(x, 0.95)))
) |> plyr::ldply()



gammas_df <- data.frame(Q = df_gammas$Q, gammas = df_gammas$variable, true = df_gammas$value,
                               q1 = gammas_ic$q1, q2 = gammas_ic$q2)

gammas_df$variable <- factor(df_gammas$variable,
                             labels = c(bquote(gamma^(1)),
                                        bquote(gamma^(2)),
                                        bquote(gamma^(3)),
                                        bquote(gamma^(4))))
gammas_df |>
  ggplot(aes(x = Q, y = true)) +
  geom_linerange(aes(ymin =  q1, ymax = q2), alpha = 0.5, size = 0.4) +
  geom_point(colour =  "steelblue", size = 1) +
  facet_wrap(~variable, nrow = 1,
             labeller = label_parsed) +
  labs(x = "True", y = "Estimated") +
  theme_bw(base_size = 16) +
  theme(strip.background = element_blank(),
        panel.spacing.x = unit(5,"mm"),
        legend.position = c(0.75, 0.04))


ggplot(aes(x = Q, y = value)) +
  geom_point(size = 3, colour = 'black')+
  facet_wrap(~variable, nrow = 1,
             labeller = label_parsed) +
  theme_bw(base_size = 16) +
  ylab('Estimated value')+
  theme(strip.background = element_blank(),
        panel.spacing.x = unit(5,"mm"),
        legend.position = c(0.75, 0.04))







