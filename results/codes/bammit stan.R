# Test stan code

library(rstan)


# Stan Data ---------------------------------------------------------------

data <- simBAMMIT(V = 2,
                  Q = 1,
                  Bv = c(12,10),
                  mu = 100,
                  lambda = c(12),
                  sb = 1,
                  sB = 1,
                  sy = 1)

x <- rev(expand.grid(lapply(rev(data$Bv), function(x) 1:x)))
colnames(x) <- paste0('var', 1:2)
var1 <- x$var1
var2 <- x$var2

stanData <- list(N = length(data$y),
             B1 = data$Bv[1],
             B2 = data$Bv[2],
             Q = 1,
             y = data$y,
             var1 = var1,
             var2 = var2,
             mmu = 100,
             smu = 10 ,
             stheta = 1,
             a = 0.01,
             b = 0.01)


# Stan model --------------------------------------------------------------

stanModelBammit <- stan(file = 'Results/codes/bammit.stan',
                        data = stanData,
                        chains = 2,
                        warmup = 100,
                        iter = 500,
                        thin = 2)
traceplot(stanModelBammit)
