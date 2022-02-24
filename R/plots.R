# sigmaComp <- function(sigma, lambda = 10, ntimes = 10, sizes = c(12, 6)){
#   I <- sizes[1]
#   J <- sizes[2]
#
#   blinSigma <- list(blinReal = list(),
#                     blinPredJosse = list(),
#                     blinPredCrossa = list())
#   blinSigma <- vector(mode = "list", length = length(sigma))
#
#   for (j in 1:length(sigma)) {
#     for(i in 1:ntimes){
#       data <- simulateDataAmmi(I = I, J = J, mu = 100, sg = 10, se = 10, sy = sigma[j], lambda = lambda)
#       josseModel <- josseJags(data = data,  mmu = 90, smu = 10,
#                               sg = 10, se = 10, slambda = 1,
#                               a = 0.1, b = 0.1, nthin = 1, nburnin = 500)
#
#       crossaModel <- crossaJags(data = data, mmu = 90, smu = 10,
#                                 mug = 0, sg = 10, mue = 0, se = 10,
#                                 mulambda = 10, slambda = 1, a = 0.1, b = 0.1, stheta = 1,
#                                 nthin = 1, nburnin = 500)
#
#       blinSigma[[j]]$blinReal[i] <- mean(data$blin)
#       blinSigma[[j]]$blinPredJosse[i]<- mean(josseModel$BUGSoutput$mean$blin)
#       blinSigma[[j]]$blinPredCrossa[i] <- mean(crossaModel$BUGSoutput$mean$blin)
#     }
#   }
#
#   df <- melt(blinSigma)
#   names(df) <- c("blin", "model", "sigma")
#
#   lbs <- paste0('\u03c3', " = ", sigma)
#   custLab <- function (x){
#     lbs
#   }
#
#   p1 <- df %>% ggplot(aes(x = model, y = blin, fill = model)) +
#     geom_boxplot() + scale_fill_brewer(palette="RdBu", labels = c("Crossa", "Josse", "Real")) +
#     facet_wrap(~sigma, labeller =  as_labeller(custLab)) +
#     labs(x = " ", y = "Bilinear term", fill = "Values") +
#     theme_minimal(base_size = 18) +
#     theme(axis.ticks.x = element_blank(),
#           axis.text.x = element_blank())
#   return(list(df = df,
#               plot = p1))
#
# }


#' @export
sigmaComp <- function(sigma, lambda = 10, sizes = c(12, 6), slambda = 1, ncol = ncol){
  I <- sizes[1]
  J <- sizes[2]

  blinSigma <- vector(mode = "list", length = length(sigma))

  for (j in 1:length(sigma)) {
      data <- simulateDataAmmi(I = I, J = J, mu = 100, sg = 10, se = 10, sy = sigma[j], lambda = lambda)
      josseModel <- josseJags(data = data,  mmu = 90, smu = 10,
                              sg = 10, se = 10, slambda = slambda,
                              a = 0.1, b = 0.1, nthin = 1, nburnin = 500)

      crossaModel <- crossaJags(data = data, mmu = 90, smu = 10,
                                mug = 0, mue = 0,
                                a = 0.1, b = 0.1, stheta = 1,
                                nthin = 1, nburnin = 500)

      blinSigma[[j]]$blinReal <- data$blin
      blinSigma[[j]]$blinPredJosse <- josseModel$BUGSoutput$mean$blin
      blinSigma[[j]]$blinPredCrossa <- crossaModel$BUGSoutput$mean$blin
  }

  blinSigma <- lapply(blinSigma, as.vector)
  df <- melt(blinSigma)
  names(df) <- c("blin", "NA", "model", "sigma")

  lbs <- paste0('\u03c3', " = ", sigma)
  custLab <- function (x){
    lbs
  }

  p1 <- df %>% ggplot(aes(x = model, y = blin, fill = model)) +
    geom_boxplot() + scale_fill_brewer(palette="RdBu", labels = c("Crossa", "Josse", "Real")) +
    facet_wrap(~sigma, labeller =  as_labeller(custLab)) +
    labs(x = " ", y = "Bilinear term", fill = "Values") +
    theme_minimal(base_size = 18) +
    theme(axis.ticks.x = element_blank(),
          axis.text.x = element_blank())

  p2 <- df %>%  filter(model != "blinReal") %>% ggplot(aes(x = model, y = blin)) +
    geom_boxplot(aes(fill = model)) +
    facet_wrap(~sigma, labeller =  as_labeller(custLab)) +
    labs(x = " ", y = "Bilinear term", fill = "Values") +
    scale_fill_brewer(palette="RdBu", labels = c("Crossa", "Josse", "Real")) +
    theme_minimal(base_size = 18) +
    geom_hline(yintercept = mean(df$blin[df$model=="blinReal"]), linetype="dashed", color = "#B2182B")+
    theme(legend.position="bottom",
          axis.ticks.x = element_blank(),
          axis.text.x = element_blank())

  pJosse <- lapply(blinSigma, function(x){
    qplot(x$blinReal, x$blinPredJosse, ylab = "Estimated bilinear - Josse", xlab = "True bilinear") + geom_abline() +
      geom_abline() + theme_minimal(base_size = 18) + geom_point(colour = "#F4A582", size = 1.2)
  } )

  pCrossa <- lapply(blinSigma, function(x){
    qplot(x$blinReal, x$blinPredCrossa, ylab = "Estimated bilinear - Crossa", xlab = "True bilinear") + geom_abline() +
      geom_abline() + theme_minimal(base_size = 18) + geom_point(colour = "#F4A582", size = 1.2)
  } )

  pJosseCrossa <- lapply(blinSigma, function(x){
    qplot(x$blinPredJosse, x$blinPredCrossa, ylab = "Estimated bilinear - Crossa", xlab = "Estimated bilinear - Josse") + geom_abline() +
      geom_abline() + theme_minimal(base_size = 18) + geom_point(colour = "#F4A582", size = 1.2)
  } )


  plotsJosse <- cowplot::plot_grid(arrangeGrob(grobs = pJosse, ncol = ncol))
  plotsCrossa <- cowplot::plot_grid(arrangeGrob(grobs = pCrossa, ncol = ncol))
  plotsJosseCrossa <- cowplot::plot_grid(arrangeGrob(grobs = pJosseCrossa, ncol = ncol))


  return(list(df = df,
              plot = list(p1, p2,
                          plotsJosse,
                          plotsCrossa,
                          plotsJosseCrossa)))

}

#' @export
#'

vizSimAmmi <- function(data){
  indexQ <- 1:data$Q
  indexQ <- paste0("Q", indexQ)
  gamma <- data.frame(data$gamma)
  delta <- data.frame(data$delta)
  colnames(gamma) <- colnames(delta) <- indexQ
  gamma <- melt(gamma)
  delta <- melt(delta)
  gammaDelta <- merge(gamma, delta, by = "variable")
  colnames(gammaDelta) <- c("Q", "gamma", "delta")
  gammaDeltaMelt <- melt(gammaDelta)

  mainEffects <- merge(data$g, data$e)
  colnames(mainEffects) <- c("g", "e")
  mainEffects <- melt(mainEffects)

  geBoxplot <- mainEffects %>% ggplot() +
    geom_boxplot(aes(x = variable, y = value, fill = variable)) +
    theme_minimal(base_size = 16) +
    scale_fill_brewer(palette="RdBu")+
    labs( y = "Value", x = "Parameter") +
    theme(legend.position = "none")

  boxplotGamma <- gammaDelta %>% ggplot() +
    geom_boxplot(aes(x = Q, y = gamma, fill = Q)) +
    scale_fill_brewer(palette="RdBu") +
    theme_minimal(base_size = 16) +
    labs(x = " ", y = expression(gamma))+
    theme(legend.position = "none")

  boxplotDelta <- gammaDelta %>% ggplot() +
    geom_boxplot(aes(x = Q, y = delta, fill = Q)) +
    scale_fill_brewer(palette="RdBu") +
    theme_minimal(base_size = 16) +
    labs(x = " ", y = expression(delta)) +
    theme(legend.position = "none")


  boxplotGammaDelta <- gammaDeltaMelt %>% ggplot() +
    geom_boxplot(aes(x = Q, y = value, fill = variable)) +
    theme_minimal(base_size = 16) +
    scale_fill_brewer(palette="RdBu", labels = c(expression(gamma), expression(delta)))+
    labs(x = " ", y = "Value", fill = "Parameter")

  return(list(boxplotGamma = boxplotGamma,
              boxplotDelta = boxplotDelta,
              boxplotGammaDelta = boxplotGammaDelta,
              geBoxplot = geBoxplot))
}



#' @export
#'

vizSimBammit <- function(data){
  indexQ <- 1:data$Q
  indexQ <- paste0("Q", indexQ)
  gamma <- data.frame(data$gamma)
  delta <- data.frame(data$delta)
  kappa <- data.frame(data$kappa)
  colnames(gamma) <- colnames(delta) <- colnames(kappa) <- indexQ
  gamma <- melt(gamma)
  delta <- melt(delta)
  kappa <- melt(kappa)
  gammaDeltaKappa <- list(gamma, delta, kappa) %>% purrr::reduce(inner_join, by = "variable")
  colnames(gammaDeltaKappa) <- c("Q", "gamma", "delta", "kappa")
  gammaDeltaKappaMelt <- melt(gammaDeltaKappa)

  mainEffects <- list(data.frame(g = data$g),
                      data.frame(e = data$e),
                      data.frame(t = data$t)) %>% melt()
  colnames(mainEffects) <- c("variable", "value", "Q")

  getBoxplot <- mainEffects %>% ggplot() +
    geom_boxplot(aes(x = variable, y = value, fill = variable)) +
    theme_minimal(base_size = 16) +
    scale_fill_brewer(palette="RdBu")+
    labs( y = "Value", x = "Parameter") +
    theme(legend.position = "none")

  boxplotGamma <- gammaDeltaKappa %>% ggplot() +
    geom_boxplot(aes(x = Q, y = gamma, fill = Q)) +
    scale_fill_brewer(palette="RdBu") +
    theme_minimal(base_size = 16) +
    labs(x = " ", y = expression(gamma))+
    theme(legend.position = "none")

  boxplotDelta <- gammaDeltaKappa %>% ggplot() +
    geom_boxplot(aes(x = Q, y = delta, fill = Q)) +
    scale_fill_brewer(palette="RdBu") +
    theme_minimal(base_size = 16) +
    labs(x = " ", y = expression(gamma)) +
    theme(legend.position = "none")

  boxplotKappa <- gammaDeltaKappa %>% ggplot() +
    geom_boxplot(aes(x = Q, y = kappa, fill = Q)) +
    scale_fill_brewer(palette="RdBu") +
    theme_minimal(base_size = 16) +
    labs(x = " ", y = expression(kappa)) +
    theme(legend.position = "none")


  boxplotGammaDeltaKappa <- gammaDeltaKappaMelt %>% ggplot() +
    geom_boxplot(aes(x = Q, y = value, fill = variable)) +
    theme_minimal(base_size = 16) +
    scale_fill_brewer(palette="RdBu", labels = c(expression(gamma), expression(delta), expression(kappa)))+
    labs(x = " ", y = "Value", fill = "Parameter")

  return(list(boxplotGamma = boxplotGamma,
              boxplotDelta = boxplotDelta,
              boxplotKappa = boxplotKappa,
              boxplotGammaDeltaKappa = boxplotGammaDeltaKappa,
              getBoxplot = getBoxplot))
}


