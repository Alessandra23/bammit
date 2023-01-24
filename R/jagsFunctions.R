#' Jags code BAMMIT model
#' @description A Jags code to estimated the parameters of BAMMIT model
#' @param data A list containing the variables y and Bv.
#' @param mmu Mean of the grand mean
#' @param smu sd of grand mean
#' @param Q number of components of interaction term
#' @param a shape parameter of a Gamma distribution
#' @param b scale parameter of a Gamma distribution
#' @param nthin thinning rate
#' @param nburnin length of burn in
#' @return The output of Jags function
#' @export
#' @importFrom R2jags 'jags'
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
      "b1", "b2",  "lambda", "beta1", "beta2", "sy", "muall", "int", "mu"
    )


  }

  if(V == 3){

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
      "b1", "b2", "b3" ,"lambda", "beta1", "beta2", "beta3", "sy", "muall", "int", "mu"
    )

  }

  if(V == 4){

    modelCode <- "
    model{
      # Likelihood
      for (n in 1:N) {
          Y[n] ~ dnorm(mu[n], sy)
          mu[n] = muall + b1[var1[n]] + b2[var2[n]] + b3[var3[n]] + b4[var4[n]]  + int[n]
          int[n] = sum(lambda[1:Q] * beta1[var1[n],1:Q] * beta2[var2[n],1:Q] * beta3[var3[n],1:Q] * beta4[var4[n],1:Q])
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
      "b1", "b2", "b3", "b4", "lambda", "beta1", "beta2", "beta3", "beta4", "sy", "muall", "int", "mu"
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

#' Jags code AR-BAMMIT model
#' @description A Jags code to estimated the parameters of AR-BAMMIT model
#' @param data A list containing the variables y and Bv.
#' @param mmu Mean of the grand mean
#' @param smu sd of grand mean
#' @param Q number of components of interaction term
#' @param a shape parameter of a Gamma distribution
#' @param b scale parameter of a Gamma distribution
#' @param nthin thinning rate
#' @param nburnin length of burn in
#' @return The output of Jags function
#' @export
#' @importFrom R2jags 'jags'
ARbammitJags <- function(data, Q, mmu, smu, a,b, stheta = 1, somega = 1, nthin, nburnin){

  Y <- data$y
  V <- length(data$Bv)
  Bv <- data$Bv
  N <- length(data$y)


  if(V == 3){

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

      # Prior on effect 3 - auto-regressive
      bb3[1] ~ dnorm(0,1)
      phib ~ dunif(-1,1)
      for (k in 2:B3) {
        muk[k] = alphab + phib * bb3[k-1]
        bb3[k] ~ dnorm(muk[k], sb3^-2)
      }
      muB3 = sum(bb3)/B3
      for(k in 1:B3){
        b3[k] = bb3[k] - muB3
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
        thetaB3[1,q] ~ dnorm(0,stheta)
       }

       phiB3 ~ dunif(-1,1)
       for(q in 1:Q){
        for(k in 2:B3){
          thetaB3[k,q] ~ dnorm(alphaBeta + phiB3 * thetaB3[k-1,q], somega)
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
      alphab ~ dnorm(0, 10^-2)
      alphaBeta ~ dnorm(0, 10^-2)

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
      somega = somega,
      a = a,
      b = b
    )

    # Choose the parameters to watch
    modelParameters <- c(
      "b1", "b2", "b3" ,"lambda", "beta1", "beta2", "beta3", "sy", "muall", "int", "mu"
    )

  }

  if(V == 4){

    modelCode <- "
    model{
      # Likelihood
      for (n in 1:N) {
          Y[n] ~ dnorm(mu[n], sy)
          mu[n] = muall + b1[var1[n]] + b2[var2[n]] + b3[var3[n]] + b4[var4[n]]  + int[n]
          int[n] = sum(lambda[1:Q] * beta1[var1[n],1:Q] * beta2[var2[n],1:Q] * beta3[var3[n],1:Q] * beta4[var4[n],1:Q])
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

      # Prior on effect 4 - auto-regressive
      bb4[1] ~ dnorm(0,1)
      phib ~ dunif(-1,1)
      for (k in 2:B4) {
        muk[k] = alphab + phib * bb4[k-1]
        bb4[k] ~ dnorm(muk[k], sb4^-2)
      }
      muB4 = sum(bb4)/B4
      for(k in 1:B4){
        b4[k] = bb4[k] - muB4
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
        thetaB4[1,q] ~ dnorm(0,stheta)
       }
       phiB4 ~ dunif(-1,1)
       for(q in 1:Q){
        for(k in 2:B4){
          thetaB4[k,q] ~ dnorm(alphaBeta + phiB4 * thetaB4[k-1,q], somega)
        }
        mB4[q] = sum(thetaB4[1:B4,q])/B4
        for(k in 1:B4){
        thetaB4New[k,q] = thetaB4[k,q] - mB4[q]
        }
        sqrtThetaB4[q] = sqrt(1/(sum(thetaB4New[1:B4,q]^2+ 0.000001)))
        for(k in 1:B4){
          beta4[k,q] = thetaB4New[k,q]*sqrtThetaB4[q]
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
      alphab ~ dnorm(0, 10^-2)
      alphaBeta ~ dnorm(0, 10^-2)

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
      somega = somega,
      a = a,
      b = b
    )

    # Choose the parameters to watch
    modelParameters <- c(
      "b1", "b2", "b3", "b4", "lambda", "beta1", "beta2", "beta3", "beta4", "sy", "muall", "int", "mu"
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

#' Jags code BAMMIT model Real Data
#' @description A Jags code to estimated the parameters of BAMMIT model
#' @param data A data frame containing the columns Yield, Genotype, Environment, Year and Block
#' @param mmu Mean of the grand mean
#' @param smu sd of grand mean
#' @param Q number of components of interaction term
#' @param a shape parameter of a Gamma distribution
#' @param b scale parameter of a Gamma distribution
#' @param nthin thinning rate
#' @param nburnin length of burn in
#' @return The output of Jags function
#' @export
#' @importFrom R2jags 'jags'
#'
bammitJagsRealData <- function(data, Q = 1, mmu = 10, smu = 2, stheta = 1, a = 0.1, b = 0.1,
                               nthin = 1, nburnin = 2000){

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
    "g", "e", "t", "bl", "lambda", "gamma", "delta", "rho", "kappa","sy", "muall", "blin", "mu"
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

#' Jags code AR BAMMIT model
#' @description A Jags code to estimated the parameters of AR BAMMIT model
#' @param data A data frame containing the columns Yield, Genotype, Environment, Year and Block
#' @param mmu Mean of the grand mean
#' @param smu sd of grand mean
#' @param Q number of components of interaction term
#' @param a shape parameter of a Gamma distribution
#' @param b scale parameter of a Gamma distribution
#' @param nthin thinning rate
#' @param nburnin length of burn in
#' @return The output of Jags function
#' @export
#' @importFrom R2jags 'jags'
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
    "g", "e", "t", "bl", "lambda", "beta1", "beta2", "beta3", "beta4","sy", "muall", "blin", "mu"
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

#' Jags code AMMI model
#' @description A Jags code to estimated the parameters of AMMITmodel
#' @param data A data frame containing the columns Yield, Genotype and Environment
#' @param mmu Mean of the grand mean
#' @param smu sd of grand mean
#' @param Q number of components of interaction term
#' @param a shape parameter of a Gamma distribution
#' @param b scale parameter of a Gamma distribution
#' @param nthin thinning rate
#' @param nburnin length of burn in
#' @return The output of Jags function
#' @export
#' @importFrom R2jags 'jags'
AMMIJagsRealData <- function(data, Q = 1, mmu = 10, smu = 2, stheta = 1, a = 0.1, b = 0.1,
                             nthin = 1, nburnin = 2000){

  Q <- Q
  Y <- data$Yield
  I <- length(levels(unique(data$Genotype)))
  J <- length(levels(unique(data$Environment)))
  N <- length(data$Yield)

  genotype <- data$Genotype
  environment <- data$Environment

  modelCode <- "
  model
  {
  # Likelihood
   for (n in 1:N) {
    Y[n] ~ dnorm(mu[n], sy)
    mu[n] = muall + g[genotype[n]] + e[environment[n]] + blin[n]
    blin[n] = sum(lambda[1:Q] * gamma[genotype[n],1:Q] * delta[environment[n],1:Q])
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


  # Prior on eigenvalues
  for(q in 1:Q) {
    lambda_raw[q] ~ dnorm(0, slambda^-2)T(0,)
  }
  lambda = sort(lambda_raw)

  #Priors
  sg ~ dt(0, 1, 1)T(0,)
  se ~ dt(0, 1, 1)T(0,)
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
    Q = Q,
    genotype = genotype,
    environment = environment,
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
    "g", "e",  "lambda", "gamma", "delta","sy", "muall", "blin"
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
