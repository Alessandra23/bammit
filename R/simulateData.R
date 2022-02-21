#' Simulate data from AMMI model
#'
#' @description A funtion to generate values of a AMMI model
#' @param I Number of genotypes
#' @param J Number of environments
#' @param mu Grand mean
#' @param sg standard deviation of genotypes
#' @param se standard deviation of environments
#' @param sy standard deviation of the model
#' @param lambda vector of values of $\lambda$
#' @return A list containing the response y, genotypes, environments, bilinear term, Q, and a matrix of g and e.
#' @example
#' data <- simulateDataAmmi(I = 6, J = 4, mu = 100, sg = 10, se = 10, sy = 2, lambda = c(10,12))
#' @export
#'
simulateDataAmmi <- function(I, J, mu, sg, se, sy, lambda) {
  # number of obs
  N <- I * J
  Q <- length(lambda)

  # generate main effects

  g <- rnorm(I, 0, sg)
  e <- rnorm(J, 0, se)

  # generate lambda, gamma, delta and kappa
  gamma <- generate_gamma_delta(I, Q)
  delta <- generate_gamma_delta(J, Q)

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
    blin = blin,
    I = I,
    J = J,
    Q = Q,
    x = x
  ))
}


#' Simulate data from Bayesian AMMIT model
#'
#' @description A funtion to generate values of the BAMMIT model
#' @param I Number of genotypes.
#' @param J Number of environments.
#' @param K Number of terms of time.
#' @param mu Grand mean.
#' @param sg standard deviation of genotypes.
#' @param se standard deviation of environments.
#' @param st standard deviation of time.
#' @param sy standard deviation of the model.
#' @param lambda vector of values of $\lambda$.
#' @return A list containing the response y, genotypes, environments, time,
#' bilinear term, Q, and a matrix of g, e and t.
#' @example data <- simulateDataBammit(I = 6, J = 4, K = 2, mu = 100, sg = 10, se = 10, st = 10, sy = 2, lambda = c(10,12))
#' @export
#'
simulateDataBammit <- function(I, J, K, mu, sg, se, st, sy, lambda){
  N <- I * J * K
  Q <- length(lambda)

  # generate main effects

  g <- rnorm(I, 0, sg)
  e <- rnorm(J, 0, se)
  t <- rnorm(K, 0, st)

  # generate lambda, gamma, delta and kappa
  gamma <- generate_gamma_delta(I, Q)
  delta <- generate_gamma_delta(J, Q)
  kappa <- generate_gamma_delta(K, Q)

  x <- expand.grid(1:I, 1:J, 1:K)
  names(x) <- c("g", "e", "t")
  x$g <- as.factor(x$g)
  x$e <- as.factor(x$e)
  x$t <- as.factor(x$t)

  # generate bilinear term
  blin <- rep(0, N)
  for (k in 1:length(lambda)) {
    blin <- blin + lambda[k] * gamma[x[, "g"], k] * delta[x[, "e"], k] *kappa[x[, "t"], k]
  }

  mu_ij <- mu + g[x[, "g"]] + e[x[, "e"]] + blin

  # generate y
  y <- rnorm(N, mu_ij, sy)

  return(list(
    y = y,
    g = g,
    e = e,
    t = t,
    blin = blin,
    I = I,
    J = J,
    K = K,
    Q = Q,
    x = x
  ))
}

# Function from ambarti package (https://github.com/ebprado/AMBARTI)

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

generate_gamma_delta <- function(INDEX, Q) {
  first_row <- TRUE

  while (first_row) {
    raw_par <- matrix(rnorm(INDEX * Q), ncol = Q)
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
