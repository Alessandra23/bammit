# Example Andrew

library(R2jags)
library(ggplot2)

rm(list = ls())

# Simulate the data -------------------------------------------------------

# Function to impose the restrictions
genInt <- function(index = 10, Q = 2, stheta = 1) {
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

  # Get variable
  variable <- sweep(thetaN, 2, sqrtTheta, "*")

  return(variable)
}

# Function to simulate the data
genData <- function(B1, B2, mu, sigma, sb, lambda) {
  Bv <- c(B1, B2)
  N <- Reduce("*", Bv)
  Q <- length(lambda)
  V <- length(Bv)

  x <- expand.grid(1:Bv[1], 1:Bv[2]) |>
    dplyr::mutate_if(is.numeric, as.factor)

  colnames(x) <- c("var1", "var2")


  b1 <- rnorm(B1, 0, sb)
  b1 <- b1 - mean(b1)

  b2 <- rnorm(B2, 0, sb)
  b2 <- b2 - mean(b2)


  beta1 <- genInt(B1, Q)
  beta2 <- genInt(B2, Q)

  int <- rep(0, N)
  for (k in 1:Q) {
    int <- int + lambda[k] * beta1[x[, "var1"], k] * beta2[x[, "var2"], k]
  }

  # m <-  mu + b1[x[,'var1']] + b2[x[,'var2']] +  int
  m <- int
  y <- rnorm(N, m, sigma)

  return(list(
    mu = mu,
    b1 = b1,
    b2 = b2,
    int = int,
    y = y,
    betas = list(beta1 = beta1, beta2 = beta2),
    Bv = Bv
  ))
}

# Generate the data
B1 <- 12
B2 <- 10
lambda <- 10 # c(8,10, 12)
sb <- 1
sigma <- 1
mu <- 100

set.seed(020)
dat <- genData(B1 = B1, B2 = B2, mu = mu, sb = sb, sigma = sigma, lambda = lambda)

# Check constraints
dat$betas$beta1[, 1] |> sum()
dat$betas$beta1[, 1]^2 |> sum()

dat$betas$beta2[, 1] |> sum()
dat$betas$beta2[, 1]^2 |> sum()


# Jags model --------------------------------------------------------------------

ammiJags <- function(data, Q, mmu, smu, a, b, mtheta = 0, stheta = 1, nthin, nburnin) {
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

      # Main effects
      for(i in 1:B1) {
          b1[i] ~ dnorm(0, sb1^-2) # Prior on effect 1
      }

      for(i in 1:B2) {
          b2[i] ~ dnorm(0, sb2^-2) # Prior on effect 2
      }

      # Interaction beta 1
      for(q in 1:Q){
          # generate values of theta from a Normal distribution
          for(i in 1:B1){
            thetaBeta1[i,q] ~ dnorm(mtheta,1/stheta^2)
          }
          # calculate the mean of each column
          mBeta1[q] = sum(thetaBeta1[1:B1,q])/B1

          # centering
          for(i in 1:B1){
            thetaBeta1New[i,q] = thetaBeta1[i,q] - mBeta1[q]
          }

          # get sd of each column
          sqrtThetaBeta1[q] = sqrt(sum(thetaBeta1New[1:B1,q]^2)) + 0.00000001

          # get the final beta
          for(i in 1:B1){
            beta1[i,q] = thetaBeta1New[i,q]/sqrtThetaBeta1[q]
          }
      }

       # Interaction beta 2
       for(q in 1:Q){
          for(i in 1:B2){
            thetaBeta2[i,q] ~ dnorm(mtheta,1/stheta^2)
          }
          mBeta2[q] = sum(thetaBeta2[1:B2,q])/B2

          for(i in 1:B2){
            thetaBeta2New[i,q] = thetaBeta2[i,q] - mBeta2[q]
          }
          sqrtThetaBeta2[q] = sqrt(sum(thetaBeta2New[1:B2,q]^2)) + 0.00000001
          for(i in 1:B2){
            beta2[i,q] = thetaBeta2New[i,q]/sqrtThetaBeta2[q]
          }
      }

      # Prior on lambda
      for(q in 1:Q) {
        lambda_raw[q] ~ dnorm(0, slambda^-2)T(0,)
      }
      lambda = sort(lambda_raw)

      # Priors
      sb1 ~ dt(0, 1, 1)T(0,)
      sb2 ~ dt(0, 1, 1)T(0,)
      slambda ~ dt(0, 1, 1)T(0,)

      sy ~ dgamma(a, b) # Prior precision

    }
    "

  B1 <- Bv[1]
  B2 <- Bv[2]

  x <- expand.grid(lapply(Bv, function(x) 1:x))
  colnames(x) <- paste0("var", 1:2)
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
    mtheta = mtheta,
    a = a,
    b = b
  )

  # Choose the parameters to watch
  modelParameters <- c(
    "b1", "b2", "lambda", "beta1", "beta2", "sy", "muall", "int", "mu",
    "thetaBeta1", "thetaBeta1New", "sqrtThetaBeta1", "mBeta1"
  )

  # Run the model
  modelRun <- jags(
    data = modelData,
    parameters.to.save = modelParameters,
    model.file = textConnection(modelCode),
    # progress.bar = "none",
    n.thin = nthin,
    n.burnin = nburnin,
    n.iter = 2000 * nthin + nburnin
  )

  return(modelRun)
}


# Run Jags model ----------------------------------------------------------

model <- ammiJags(
  data = dat,
  Q = 1,
  mtheta = 0,
  stheta = 1,
  mmu = 100,
  smu = 10,
  a = 0.1,
  b = 0.1,
  nthin = 2,
  nburnin = 2000
)

plot(model)

# Check betas
model$BUGSoutput$mean$beta1 |> sum()
model$BUGSoutput$mean$beta1^2 |> sum()

model$BUGSoutput$mean$beta2 |> sum()
model$BUGSoutput$mean$beta2^2 |> sum()


# Check thetas
model$BUGSoutput$mean$thetaBeta1
model$BUGSoutput$mean$thetaBeta1New
model$BUGSoutput$mean$sqrtThetaBeta1
model$BUGSoutput$mean$mBeta1


# Plots

# interaction
qplot(dat$int, model$BUGSoutput$mean$int) + geom_abline() + theme_bw() + labs(x = "True int", y = "Pred int")
# y
qplot(dat$y, model$BUGSoutput$mean$mu) + geom_abline() + theme_bw() + labs(x = "y", y = expression(hat(y)))
