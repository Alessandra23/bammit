# --------------------- Running AMMI (each year) ----------------- #

# packages

library(bammit)
library(ggplot2)
library(stringr)
library(R2jags)
library(dplyr)

load("Real data/train_ireland.RData")
load("Real data/test_ireland.RData")


# Run AMMI models

# base_dir = '/Volumes/ALESSA HD/PhD/bammit/Real data/'
# year <- '2010'
# fileName <- paste0(base_dir, 'Ireland_VCU_', year, '_data.RData')
# load(fileName)

# all years

model <- AMMIJagsRealData(train, Q = 1, mmu = 10, smu = 2, stheta = 1, a = 0.1, b = 0.1,
                          nthin = 2, nburnin = 2000)

saveRDS(model, file = "Running models/modelAmmiAllYears.rds")

pred <- predictionAMMIReal(model = model$BUGSoutput, test)

caret::RMSE(pred, test$Yield)
caret::R2(pred, test$Yield)



# by year


train_data <- split(train, train$Year)
test_data <- split(test, test$Year)
model <- lapply(train_data, function(x){
  AMMIJagsRealData(x, Q = 1, mmu = 10, smu = 2, stheta = 1, a = 0.1, b = 0.1,
                   nthin = 2, nburnin = 2000)
})

modelAmmi <- list()
for(i in 1:10){
  modelAmmi[[i]] <- model[[i]]$BUGSoutput
}
saveRDS(modelAmmi, file = "Running models/modelAmmi.rds")
pred <- list()
for(i in 1:10){
  pred[[i]] <- predictionAMMIReal(model = model[[i]]$BUGSoutput, test_data[[i]])
}

caret::RMSE(unlist(pred), test$Yield)
caret::R2(unlist(pred), test$Yield)


caret::RMSE(pred[[3]], test_data[[3]]$Yield)
caret::R2(unlist(pred), test$Yield)

train_data <- train |>
  filter(Year == '2012') |>
  # group_by(Genotype, Environment) |>
  # summarise(Yield = mean(Yield)) |>
  # ungroup() |>
  mutate(Genotype = factor(Genotype, levels = unique(Genotype))) |>
  mutate(Environment = factor(Environment, levels = unique(Environment)))

model <- AMMIJagsRealData(train_data, Q = Q, mmu = 10, smu = 2, stheta = 1, a = 0.1, b = 0.1,
                          nthin = 2, nburnin = 2000)

# predictions
pred <- predictionAMMIReal(model = model, test_data = test, train_data = train_data)
caret::RMSE(pred, test$Yield)
caret::R2(pred, test$Yield)





