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
  g <- model$mean$b1
  e <- model$mean$b2
  t <- model$mean$b3
  bl <- model$mean$b4
  inthat <- model$mean$int

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
predictionAMMIReal <- function(model, data){#train_data, test_data) {

  # x <- test_data[,c('Genotype', 'Environment')]
  # colnames(x) <- c('genotype', 'environment')
  #
  # muhat <- model$BUGSoutput$mean$muall
  # g_aux <- model$BUGSoutput$mean$g
  # e_aux <- model$BUGSoutput$mean$e
  #
  # names(g_aux) <- levels(train_data$Genotype)
  # names(e_aux) <- levels(train_data$Environment)
  #
  #
  # g <- rep(0,length(levels(x$genotype)))
  # names(g) <- levels(x$genotype)
  # g[names(g) %in% names(g_aux)] <- g_aux
  #
  # e <- rep(0,length(levels(x$environment)))
  # names(e) <- levels(x$environment)
  # e[names(e) %in% names(e_aux)] <- e_aux
  #
  # inthat_aux <- model$BUGSoutput$mean$blin
  # inthat_aux <- matrix(inthat_aux, ncol = length(g_aux), nrow = length(e_aux))
  # colnames(inthat_aux) <- levels(train_data$Genotype)
  # rownames(inthat_aux) <- levels(train_data$Environment)
  #
  # inthat <- matrix(0, ncol = length(levels(x$genotype)), nrow = length(levels(x$environment)))
  # colnames(inthat) <- levels(x$genotype)
  # rownames(inthat) <- levels(x$environment)
  # inthat[rownames(inthat_aux), colnames(inthat_aux)] <- inthat_aux
  #
  # df <- reshape2::melt(inthat)
  # names(df) <- c('environment', 'genotype', 'int')
  # df <- plyr::join(x, df)
  #
  # N <- length(test_data$Yield)
  # yhat <- rep(muhat, N) + g[df[,'genotype']]  + e[df[,'environment']]  +  df[,'int']
  #
  # return(as.vector(yhat))


  genotype <- data$Genotype
  environment <- data$Environment

  muhat <- model$mean$muall
  g <- model$mean$g
  e <- model$mean$e
  inthat <- model$mean$blin

  N <- length(data$Yield)
  yhat <- rep(muhat, N) + g[genotype] + e[environment] +  inthat

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


#' @export
classical_AMMI <- function(data, Q){

  y <- data$y

  g = data$g
  e = data$e


  # Fit the linear model
  linear_mod = aov(y ~ g + e + g:e)
  # linear_mod = lm(y_train ~ g + e)

  # Get the residuals for the interaction g:e
  #interaction_tab = matrix(residuals(linear_mod), ncol = length(unique(g)), length(unique(env)))
  interaction_tab = model.tables(linear_mod, type='effects', cterms = 'g:e')
  interaction_tab = interaction_tab$tables$`g:e`

  # Get the number of PCs
  if (is.null(data$Q) == FALSE) {Q = data$Q}

  # Run the Singular Value Decomposition (SVD) to compute lambda, gamma, and delta
  sv_dec <- svd(interaction_tab, nu = Q, nv = Q)

  # Get parameter estimates
  # mu_hat     = linear_mod$coefficients[1] # slightly biased compared to mean(y_train)
  mu_hat     = mean(y)
  g_hat      = aggregate(x = y - mu_hat, by = list(g), FUN = "mean")[,2]
  e_hat      = aggregate(x = y - mu_hat, by = list(e), FUN = "mean")[,2]
  lambda_hat = sv_dec$d[1:Q]
  gamma_hat  = -1*sv_dec$u
  delta_hat  = -1*sv_dec$v

  return(list(mu_hat     = mu_hat,
              g_hat      = g_hat,
              e_hat      = e_hat,
              lambda_hat = lambda_hat,
              gamma_hat  = gamma_hat,
              delta_hat  = delta_hat))
}
