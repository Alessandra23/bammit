## Get a small data set of data2
load("Real data/data2.RData")
years <- c("2005")
data2small <- data2 |> dplyr::filter(YEAR %in% years) |>
  rename(Genotype = GENOTYPE,
         Environment = LOCATION,
         Year = YEAR,
         Mean = MEAN)


modelData2 <- bammitJagsRealData(
  data = data2small, mmu = 90, smu = 10, stheta = 1/100, a = 0.1, b = 0.1,
  nthin = 1, nburnin = 10, Q = 1
)


fitcmcm <- as.mcmc(modelData2)
fitmat <- as.matrix(fitcmcm)
fitdf <- as.data.frame(fitmat)

fitg <- fitdf[, grep(x = colnames(fitdf), pattern = "g[", fixed = TRUE)]
fite <- fitdf[, grep(x = colnames(fitdf), pattern = "e[", fixed = TRUE)]
fitt <- fitdf[, grep(x = colnames(fitdf), pattern = "t[", fixed = TRUE)]
fitblin <- fitdf[, grep(x = colnames(fitdf), pattern = "blin[", fixed = TRUE)]
fitgamma <- fitdf[, grep(x = colnames(fitdf), pattern = "gamma[", fixed = TRUE)]
#get all
#gPostCoeff <- select(fitdf,paste0("g[", seq(1, ncol(fitg), by = 1), "]"))
gPostCoeff <- select(fitdf,paste0("g[", seq(10, 20, by = 1), "]"))
names(gPostCoeff) <- paste0("g", 1:ncol(gPostCoeff))
gPostCoeffLong <- gather(gPostCoeff)

gPostCoeffLong <- gather(fitg)

gennames <- c("g[1]", 'g[20]', 'g[24]')
gPostCoeffLong <- gPostCoeffLong |> filter(key %in% gennames)


g1 <- ggplot(data = gPostCoeffLong, aes(x = value, y = key)) +
  stat_density_ridges(quantile_lines = TRUE,
                      quantiles = c(0.025, 0.5, 0.975),
                      alpha = 0.7,
                      colour = "steelblue",
                      fill = "steelblue") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  theme_bw() + xlab("Posterior estimate") + ylab("")
g1



# frequestist -------------------------------------------------------------

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

data <- data2small |> group_by(Year) |>
  mutate(y = mean(Mean)) |>
  select(Genotype, Environment, y)
colnames(data) <- c("year", "g", "e", "y")



