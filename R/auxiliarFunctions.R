#' @export
RMSE <- function(true, predicted){
  sqrt(mean(((true - predicted))^2))
}

#' @export
#'
predictionAMMI <- function(model, data) {
  muhat <- model$BUGSoutput$mean$muall
  ghat <- model$BUGSoutput$mean$g
  ehat <- model$BUGSoutput$mean$e
  blinhat <- model$BUGSoutput$mean$blin

  N <- length(data$y)

  yhat <- rep(muhat, N) + ghat[data$x$g] + ehat[data$x$e] +  blinhat

  return(yhat)
}

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
