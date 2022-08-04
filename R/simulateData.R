
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
    thetaN <- apply(theta, 1, function(x){x-m}) %>% as.matrix()

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
    bv <- rnorm(x, 0, sb)
    bv <- bv - mean(bv)
    return(bv)
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
