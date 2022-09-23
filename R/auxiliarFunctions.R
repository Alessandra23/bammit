#' Root Mean Square Error
#' @description A function to calculate the root mean square error.
#' @param true A vector of the true values of the response variable.
#' @param predicted A vector of the predicted values of the response variable.
#' @return A number.
#' @export
#'
RMSE <- function(true, predicted){
  sqrt(mean(((true - predicted))^2))
}

#' R2
#' @description A function to calculate the R2.
#' @param true A vector of the true values of the response variable.
#' @param predicted A vector of the predicted values of the response variable.
#' @return A number.
#' @export
#'
R2 <- function(pred, obs){
  yobshat <- mean(obs)
  R2 <- 1 - (sum((obs - pred)^2)/sum((obs - yobshat)^2))
  return(R2)
}

#' Predicted values of BAMMIT model
#' @description A function to calculate the predicted values of a BAMMIT model from a JAGS object
#' @param model A object from Jags model
#' @param data A object from simulated data.
#' @return A vector of predicted values of the response variable.
#' @export
#'
predictionBAMMIT <- function(model, data) {

  V <- length(data$Bv)

  if(V == 2){
    Bv <- data$Bv
    x <- rev(expand.grid(lapply(rev(Bv), function(x) 1:x)))
    colnames(x) <- paste0('var', 1:2)
    var1 <- x$var1
    var2 <- x$var2

    muhat <- model$mean$muall
    b1 <- model$mean$b1
    b2 <- model$mean$b2
    inthat <- model$mean$blin

    N <- length(data$y)
    yhat <- rep(muhat, N) + b1[x$var1] + b2[x$var2] +   inthat
  }

  if(V == 3){
    Bv <- data$Bv
    x <- rev(expand.grid(lapply(rev(Bv), function(x) 1:x)))
    colnames(x) <- paste0('var', 1:3)
    var1 <- x$var1
    var2 <- x$var2
    var3 <- x$var3

    muhat <- model$mean$muall
    b1 <- model$mean$b1
    b2 <- model$mean$b2
    b3 <- model$mean$b3
    inthat <- model$mean$blin

    N <- length(data$y)
    yhat <- rep(muhat, N) + b1[x$var1] + b2[x$var2] + b3[x$var3] +  inthat
  }

  if(V == 4){
    Bv <- data$Bv
    x <- rev(expand.grid(lapply(rev(Bv), function(x) 1:x)))
    colnames(x) <- paste0('var', 1:4)
    var1 <- x$var1
    var2 <- x$var2
    var3 <- x$var3
    var4 <- x$var4

    muhat <- model$mean$muall
    b1 <- model$mean$b1
    b2 <- model$mean$b2
    b3 <- model$mean$b3
    b4 <- model$mean$b4
    inthat <- model$mean$blin

    N <- length(data$y)
    yhat <- rep(muhat, N) + b1[x$var1] + b2[x$var2] + b3[x$var3] + b4[x$var4]  +  inthat
  }

  return(yhat)
}

#' Predicted values of BAMMIT model for real data
#' @description A function to calculate the predicted values of a BAMMIT model from a JAGS object
#' @param model A object from Jags model
#' @param data A object from simulated data.
#' @return A vector of predicted values of the response variable.
#' @export
#'
predictionBAMMITReal <- function(model, data) {

  genotype <- data$Genotype
  environment <- data$Environment
  time <- data$Year
  bloc <- data$Bloc

  muhat <- model$mean$muall
  g <- model$mean$g
  e <- model$mean$e
  t <- model$mean$t
  bl <- model$mean$bl
  inthat <- model$mean$blin

  N <- length(data$Yield)
  yhat <- rep(muhat, N) + g[genotype] + e[environment] + t[time] + bl[bloc]  +  inthat


  return(yhat)
}

#' Predicted values of AMMI model for real data
#' @description A function to calculate the predicted values of a AMMI model from a JAGS object
#' @param model A object from Jags model
#' @param data A object from simulated data.
#' @return A vector of predicted values of the response variable.
#' @export
#'
predictionAMMIReal <- function(model, data) {


  genotype <- data$Genotype
  environment <- data$Environment

  muhat <- model$mean$muall
  g <- model$mean$g
  e <- model$mean$e
  inthat <- model$mean$blin

  N <- length(data$Yield)
  yhat <- rep(muhat, N) + g[genotype] + e[environment]  +  inthat


  return(yhat)
}

#' @export
#'
loadRData <- function(fileName){
  load(fileName)
  get(ls()[ls() != "fileName"])
}


#' Tranform list of data into data set
#' @export
#' @param data A list
transdat <- function(data){

  Y <- data$y
  V <- length(data$Bv)
  Bv <- data$Bv
  N <- length(data$y)

  if(V == 2){

    B1 <- Bv[1]
    B2 <- Bv[2]

    x <- rev(expand.grid(lapply(rev(Bv), function(x) 1:x)))
    colnames(x) <- paste0('var', 1:2)
    x$var1 <- as.factor(x$var1)
    x$var2 <- as.factor(x$var2)
    dat <- cbind(x, Y)
  }

  if(V == 3){

    B1 <- Bv[1]
    B2 <- Bv[2]
    B3 <- Bv[3]

    x <- rev(expand.grid(lapply(rev(Bv), function(x) 1:x)))
    colnames(x) <- paste0('var', 1:3)
    x$var1 <- as.factor(x$var1)
    x$var2 <- as.factor(x$var2)
    x$var3 <- as.factor(x$var3)
    dat <- cbind(x, Y)

  }


  if(V == 4){

    B1 <- Bv[1]
    B2 <- Bv[2]
    B3 <- Bv[3]
    B4 <- Bv[4]

    x <- rev(expand.grid(lapply(rev(Bv), function(x) 1:x)))
    colnames(x) <- paste0('var', 1:4)
    x$var1 <- as.factor(x$var1)
    x$var2 <- as.factor(x$var2)
    x$var3 <- as.factor(x$var3)
    x$var4 <- as.factor(x$var4)
    dat <- cbind(x, Y)

  }

  return(dat)

}
