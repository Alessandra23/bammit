library(R2jags)
library(reshape2)

square_root_matrix <- function(x) {

  # When Q = 1, x will be a scalar
  if (nrow(x) == 1) {
    return(sqrt(x))
  }

  # When Q > 1, then x will be a matrix
  if (nrow(x) > 1) {
    # Jordan normal form
    X <- eigen(x)
    P <- X$vectors
    A <- diag(X$values)

    A_sqrt <- diag(sqrt(X$values))
    P_inv <- solve(P)
    x_sqrt <- P %*% A_sqrt %*% P_inv
    return(x_sqrt)
  }
}

generate_gamma_delta <- function(INDEX, Q, stheta = 1) {
  first_row <- TRUE

  while (first_row) {
    raw_par <- matrix(rnorm(INDEX * Q, sd = stheta), ncol = Q)
    par_mean <- matrix(rep(apply(raw_par, 2, mean), each = nrow(raw_par)), ncol = Q)
    par_aux <- raw_par - par_mean

    # Constraints ----
    # apply(par_aux,2,sum)
    parTpar <- solve(t(par_aux) %*% (par_aux))
    A <- square_root_matrix(parTpar)
    samples <- par_aux %*% A

    # Force the first to be positive
    for (i in 1:nrow(samples)) {
      row1 <- samples[1, ]
      if (all(samples[i, ] > 0)) {
        aux <- samples[i, ]
        samples[1, ] <- aux
        samples[i, ] <- row1
        return(samples)
      }
    }
    # t(samples)%*%samples == 0
    # apply(samples,2,sum) == diag(Q)
  }
}


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

  g <- g - mean(g)
  e <- e - mean(e)

  # generate lambda, gamma, delta and kappa
  gamma <- generate_gamma_delta(I, Q)
  delta <- generate_gamma_delta(J, Q)

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


data <- simulateDataAmmi(I = 12, J = 10, mu = 100, sg = 10, se = 10, sy = 1.5, lambda = 10, stheta = 1)

josseModel <- josseJags(
  data = data, mmu = 90, smu = 10,
  sg = 10, se = 10, slambda = 1,
  a = 0.1, b = 0.1, nthin = 1, nburnin = 3000
)

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
