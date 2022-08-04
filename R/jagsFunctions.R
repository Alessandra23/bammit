#' Jags code BAMMIT model
#' @description A Jags code to estimated the parameters of AMMI model followuing Josse's approach
#' @param data A list containing the simulated values of y, the values of I, J, and Q and a mtrix od g and e.
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
      "b1", "b2",  "lambda", "beta1", "beta2", "sy", "muall", "int"
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
      "b1", "b2", "b3" ,"lambda", "beta1", "beta2", "beta3", "sy", "muall", "int"
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
      "b1", "b2", "b3", "b4", "lambda", "beta1", "beta2", "beta3", "beta4", "sy", "muall", "int"
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
