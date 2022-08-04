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

#' Predicted values of BAMMIT model
#' @description A function to calculate the predicted values of a BAMMIT model from a JAGS
#' @param model A object from Jags model
#' @param data A object from simulated data.
#' @return A vector of predicted values of the response variable.
#' @export
#'
predictionBAMMIT <- function(model, data) {
  muhat <- model$BUGSoutput$mean$muall
  ghat <- model$BUGSoutput$mean$g
  ehat <- model$BUGSoutput$mean$e
  that <- model$BUGSoutput$mean$t
  blinhat <- model$BUGSoutput$mean$blin

  N <- length(data$y)

  yhat <- rep(muhat, N) + ghat[data$x$g] + ehat[data$x$e] + that[data$x$t] +  blinhat

  return(yhat)
}

#' Predicted values of BAMMIT model
#' @description A function to calculate the predicted values of a BAMMIT model from a JAGS
#' @param model A object from Jags model
#' @param data A object from simulated data.
#' @return A vector of predicted values of the response variable.
#' @export
#'
predictionBAMMITRealData <- function(model, data) {
  muhat <- model$BUGSoutput$mean$muall
  ghat <- model$BUGSoutput$mean$g
  ehat <- model$BUGSoutput$mean$e
  that <- model$BUGSoutput$mean$t
  blinhat <- model$BUGSoutput$mean$blin

  N <- length(data$Mean)

  yhat <- rep(muhat, N) + ghat[data[, "Genotype"]$Genotype] +
    ehat[data[, "Environment"]$Environment] + that[data[, "Year"]$Year] +  blinhat[1:168]

  return(yhat)
}


#' @export
#'
loadRData <- function(fileName){
  load(fileName)
  get(ls()[ls() != "fileName"])
}
