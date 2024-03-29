
#' Auxiliar function to generate the interaction term
#' @export
genInt <- function(index, Q){
  theta <- matrix(NA, nrow = index, ncol = Q)
  variable <- matrix(NA, nrow = index, ncol = Q)
  sqrtTheta <-  vector()
  for(q in 1:Q){
    for (i in 1:index) {
      theta[i,q] <- rnorm(1)
    }
    m <- apply(theta, 2, mean)
    thetaN <- apply(theta, 1, function(x){x-m})  |>  as.matrix()

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

#' Simulate data BAMMIT
#' @export
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
    bvaux <- rnorm(x, 0, sb)
    bvF <- bvaux - mean(bvaux)
    return(bvF)
  })
  #names(bv) <- paste0("b", 1:length(Bv))

  meff <- rowSums(rev(expand.grid(rev(bv))))

  # generate bilinear term
  Beta <- vector(mode = "list", length = V)
  for (i in 1:V) {
    Beta[[i]] <- genInt(Bv[i], Q)
  }

  k <- as.list(rep(1, Q))
  for (j in 1:length(k)) {
    for (i in 1:length(Beta)) {
      k[[j]] <- kronecker(k[[j]], Beta[[i]][,j])
    }
    k[[j]] <- k[[j]]*lambda[j]
  }

  int <- Reduce("+", k)

  # generate y
  m <- mu + meff + int
  y <- rnorm(N, m, sy[])

  return(list(mu = mu,
              bv = bv,
              int = int,
              y = y,
              Bv = Bv))

}



#' Simulate data BAMMIT
#' @export
simBAMMIT_li <- function(V = 3,
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
    bvaux <- rnorm(x, 0, sb)
    bvF <- bvaux - mean(bvaux)
    return(bvF)
  })
  #names(bv) <- paste0("b", 1:length(Bv))

  meff <- rowSums(rev(expand.grid(rev(bv))))

  # generate bilinear term
  Beta <- vector(mode = "list", length = V)
  for (i in 1:V) {
    Beta[[i]] <- genInt(Bv[i], Q)
  }

  k <- as.list(rep(1, Q))
  for (j in 1:length(k)) {
    for (i in 1:length(Beta)) {
      k[[j]] <- kronecker(k[[j]], Beta[[i]][,j])
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


#' Auxiliar function to generate the interaction term for the AR-BAMMIT
#' @export
genARInt <- function(index, Q, somega = 1){

  theta <- matrix(NA, nrow = index, ncol = Q)
  variable <- matrix(NA, nrow = index, ncol = Q)
  sqrtTheta <-  vector()

  alphaTheta <- rnorm(1,0,100)
  phiTheta <- runif(1,-1,1)
  for(q in 1:Q){
    theta[1,q] <- rnorm(1,0,1)
    for (i in 2:index) {
      theta[i,q] <- rnorm(1, alphaTheta + phiTheta * theta[i-1], somega)
    }

    m <- apply(theta, 2, mean)
    thetaN <- as.matrix(apply(theta, 1, function(x){x-m}))

    if(Q>1){
      thetaN <- t(thetaN)
    }

    sqrtTheta[q] <- sqrt(1/sum(thetaN[,q]^2))

    for (i in 1:index) {
      variable[i,q] <- (thetaN[i,q])*sqrtTheta[q]
    }
  }

  return(list(variable = variable,
              phiTheta = phiTheta,
              alphaTheta = alphaTheta))
}

#' Simulate data AR-BAMMIT
#' @example data <- simArBAMMIT(V = 3, Q = 1, Bv = c(6,4,6), mu = 100, lambda = c(10), sb = 1, sB = 10, seta = 1, somega = 1, sy = 1)
#' @export
simArBAMMIT <- function(V = 3,
                        Q = 1,
                        Bv = c(3,2,6),
                        mu = 100,
                        lambda = c(10),
                        sb = 1,
                        sB = 10,
                        seta = 1,
                        somega = 1,
                        sy = 1){

  N <- Reduce("*", Bv)

  # generate main effects
  bv <- lapply(Bv[-length(Bv)], function(x) {
    bvaux <- rnorm(x, 0, sb)
    bvF <- bvaux - mean(bvaux)
    return(bvF)
  })
  #names(bv) <- paste0("b", 1:length(Bv))
  bv[[length(Bv)]] <- vector(mode = "numeric", length = Bv[length(Bv)])
  bv[[length(Bv)]][1] <- rnorm(1,0,1)

  alphab <- rnorm(1,0,sqrt(10))
  phib <- runif(1,-1,1)
  for (k in 2:Bv[length(Bv)]) {
    muk <-  alphab + phib * bv[[length(Bv)]][k-1]
    bv[[length(Bv)]][k] <-  rnorm(1, muk, seta)
  }

  bv[[length(Bv)]] <- bv[[length(Bv)]] - mean(bv[[length(Bv)]])

  meff <- rowSums(rev(expand.grid(rev(bv))))

  # generate bilinear term
  Beta <- vector(mode = "list", length = V)
  for (i in 1:(V-1)) {
    Beta[[i]] <- genInt(Bv[i], Q)
  }

  Beta[[V]] <- genARInt(Bv[V], Q)

  alphaPhiVectors <- list(Beta[[V]]$alphaTheta, Beta[[V]]$phiTheta)
  Beta[[V]] <- Beta[[V]]$variable


  k <- as.list(rep(1, Q))
  for (j in 1:length(k)) {
    for (i in 1:length(Beta)) {
      k[[j]] <- kronecker(k[[j]], Beta[[i]][,j])
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
              Bv = Bv,
              alphaPhiVectors,
              sy = sy))

}


#' Generate data sets
#' @export
genDataSets <- function(Bv,Q, base_dir = '~/Documents/GitHub/bammit/Simulated data/test/'){

  if(Q == 1){
    lambda = 10
  }

  if(Q == 2){
    lambda = c(8,10)
  }

  if(Q == 3){
    lambda = c(8,10, 12)
  }

  V <- length(Bv)
  N <- Reduce('*', Bv)


  data <- simBAMMIT(V = V,
                    Q = Q,
                    Bv = Bv,
                    mu = 100,
                    lambda = lambda,
                    sb = 1,
                    sB = 1,
                    sy = 1)
  fileName <- paste0(base_dir, 'V', V, 'N', N, 'Q' , Q, '.rds')
  saveRDS(data, file = fileName)

}

#' Run the simulated data
#' @export
runModel <- function(data, Q, base_dir = '~/Documents/GitHub/bammit/Running models/'){

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
  fileName <- paste0(base_dir, 'model', 'V', V, 'N', N, 'Q' , Q, '.rds')
  saveRDS(model, file = fileName)

}

