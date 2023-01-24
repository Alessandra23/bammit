library(R2jags)

load("~/Documents/GitHub/bammit/Real data/train_ireland.RData")
load("~/Documents/GitHub/bammit/Real data/test_ireland.RData")

# BAMMIT Real data
trainQ1 <- bammitJagsRealData(data = train, Q = 1)
saveRDS(trainQ1, file = "~/pineR/bammit/Running models/trainQ1.rds")
trainQ2 <- bammitJagsRealData(data = train, Q = 2)
saveRDS(trainQ2, file = "~/pineR/bammit/Running models/trainQ2.rds")
trainQ3 <- bammitJagsRealData(data = train, Q = 3)
saveRDS(trainQ3, file = "~/pineR/bammit/Running models/trainQ3.rds")

predR <- predictionBAMMITReal(trainQ3$BUGSoutput, test)
RMSE(predR, test$Yield)
R2(predR, test$Yield)



# AR BAMMIT real data -----------------------------------------------------

trainQ1ar <- arbammitJagsRealData(data = train, Q = 1)
saveRDS(trainQ1ar$BUGSoutput, file = "~/pineR/bammit/Running models/trainQ1ar.rds")

trainQ1ar$BUGSoutput$median$blin
predRar <- predictionBAMMITReal(trainQ1ar$BUGSoutput, test)
RMSE(predRar, test$Yield)



# Read the data -----------------------------------------------------------

rd <- readRDS("/Users/alessalemos/Documents/GitHub/bammit/Running models/trainQ1.rds")
pred_rd <- predictionBAMMITReal(rd$BUGSoutput, test)
caret::RMSE(pred_rd, test$Yield)
caret::R2(pred_rd, test$Yield)

