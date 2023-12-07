testJags <- function(data, nthin = 1, nburnin = 1000){

  # from data
  Y <- data$y
  N <- length(Y)
  J <- data$J

  modelCode <- "
  model
  {
  # Likelihood
   for (k in 1:N) {
    Y[k] ~ dnorm(mu[k], sy^-2)
    mu[k] = sum(a[k,1:J])
   }

  for(j in 1:J){
    for(i in 1:(N-1)){
      aux[i,j] ~ dnorm(0, 1)
    }
    aux[N,j] = -sum(aux[1:(N-1),j])
    auxSum[j] = sqrt(sum(aux[1:N,j]^2)) + 0.000001
    for(i in 1:N){
      a[i,j] = aux[i,j]/auxSum[j]
    }
  }

   sy ~ dgamma(0.1, 0.1)
  }
  "

  # Set up the data
  modelData <- list(
    N = N,
    Y = Y,
    J = J
  )

  # Choose the parameters to watch
  modelParameters <- c(
    "a", "mu", "sy"
  )

  # Run the model
  modelRun <- jags(
    data = modelData,
    parameters.to.save = modelParameters,
    model.file = textConnection(modelCode),
    progress.bar = "none",
    n.thin = nthin,
    n.burnin = nburnin,
    n.iter = 2000 * nthin + nburnin
  )

  return(modelRun)
}



# simulate data

simTest <- function(N, J, sy) {

  a <- matrix(rnorm(N*J,0,1), ncol = J)
  mu <- apply(a, 1, sum)
  y <- rnorm(N, mu, sy)

  return(list(
    y = y,
    a = a,
    N = N,
    J = J
  ))
}

# test
data <-  simTest(6,3,2)
modelJags <- testJags(data = data)

modelJags$BUGSoutput$sims.list$a[,,1]
qqplot(data$a[,1], modelJags$BUGSoutput$sims.list$a[,,1]);
abline(0,1)


modelJags$BUGSoutput$mean$a[,1]^2 |> sum()

dat <- data.frame(
  x = c(1,2,3,4,5),
  y = c(10,9,8,7,6),
  z = c(2,4,6,8,10)
)

testVec <- c(2,4,6,8, 10)

which(dat$z %in% testVec)



theta = matrix(rnorm(27, 50, 1), ncol = 3, nrow = 9)
theta = rbind(theta, -apply(theta,2,sum))
apply(theta, 2,sum)

gamma <- matrix(NA, ncol = 3, nrow = 10)
thetaSum <- sqrt(sum(theta[,1]^2))
gamma[1,1] <- theta[1,1]/thetaSum
aa = theta[,1]/sqrt(sum(theta[,1]^2))
sum(aa^2)


