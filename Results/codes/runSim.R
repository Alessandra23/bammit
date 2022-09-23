# --------------------- Running simulations ----------------- #

# packages

library(bammit)
library(ggplot2)
library(dplyr)


# examples simulating data and running model ------------------------------

data <- simBAMMIT(V = 3,
                  Q = 1,
                  Bv = c(12,10, 2),
                  mu = 100,
                  lambda = c(12),
                  sb = 1,
                  sB = 1,
                  sy = 1)

model <- bammitJags(data = data,
                    Q = 1,
                    mmu = 100,
                    smu = 10,
                    a = 0.1,
                    b = 0.1,
                    nthin = 2,
                    nburnin = 2000)

# simulate data -----------------------------------------------------------

# Generate training data

# Bv <- list(c(12,10),
#            c(12,10,4),
#            c(12,10,4,2),
#            c(100,10,5))
# Q = c(1,2,3)

# create data over the Bv and Q
# for (q in Q) {
#   lapply(Bv, genDataSets, Q = q)
# }


# Reading the train data --------------------------------------------------------

base_dir_train <- '~/Documents/GitHub/bammit/Simulated data/train'
all_files_train <- list.files(base_dir_train, pattern="*.rds", full.names = TRUE)
train <- lapply(all_files_train, readRDS)
names(train) <- sub(pattern = "(.*)\\..*$", replacement = "\\1", basename(all_files_train))


# Reading the test data --------------------------------------------------------

base_dir_test <- '~/Documents/GitHub/bammit/Simulated data/test'
all_files_test <- list.files(base_dir_test, pattern="*.rds", full.names = TRUE)
test <- lapply(all_files_test, readRDS)
names(test) <- sub(pattern = "(.*)\\..*$", replacement = "\\1", basename(all_files_test))

# Running models for simulated data ---------------------------------------

runModel <- function(data, Q, name, base_dir = '~/pineR/bammit/sim run/'){

  V <- length(data$Bv)
  N <- Reduce('*', data$Bv)

  modelAux <- bammitJags(data = data,
                         Q = Q,
                         mmu = 100,
                         smu = 10,
                         a = 0.1,
                         b = 0.1,
                         nthin = 2,
                         nburnin = 2000)

  model <- modelAux$BUGSoutput
  fileName <- paste0(base_dir, 'model', name,'QM' , Q, '.rds')
  saveRDS(model, file = fileName)

}

# Q = c(1,2,3)
# for (q in Q) {
#   for (i in 1:length(train)) {
#     runModel(train[[i]], Q = q, name = names(train)[i])
#   }
# }


# Read jags output --------------------------------------------------------

# all outputs
base_dir_out <- '~/Documents/GitHub/bammit/Running models'
all_files_out <- list.files(base_dir_out, pattern="*.rds", full.names = TRUE)
out <- lapply(all_files_out, readRDS)
names(out) <- sub(pattern = "(.*)\\..*$", replacement = "\\1", basename(all_files_out))


# just the models used on the paper
modelV4N960Q1QM1 <- readRDS("/Volumes/Alessa HD/PhD/bammit/Running models/modelV4N960Q1QM1.rds")

# Fitting ML models -------------------------------------------------------

library(neuralnet)
library(randomForest)

datV2train <- transdat(data = train[[2]])
datV3train <- transdat(data = train[[5]]) # N = 480
datV4train <- transdat(data = train[[11]])

datV2test <- transdat(data = test[[2]])
datV3test <- transdat(data = test[[5]]) # N = 480
datV4test <- transdat(data = test[[11]])


# ----- Random Forest

# models

mrfV2 <- randomForest(Y ~ ., data=datV2train, mtry = 2, importance=TRUE, na.action=na.omit)
mrfV3 <- randomForest(Y ~ ., data=datV3train, mtry = 2, importance=TRUE, na.action=na.omit)
mrfV4 <- randomForest(Y ~ ., data=datV4train, mtry = 2, importance=TRUE, na.action=na.omit)

# predictions

predrfV2 <- predict(mrfV2, data = datV2test[,-3])
predrfV3 <- predict(mrfV3, data = datV3test[,-4])
predrfV4 <- predict(mrfV4, data = datV4test[,-5])

# rmse

rmserfV2 <- RMSE(predrfV2, datV2test[,3])
rmserfV3 <- RMSE(predrfV3, datV3test[,4])
rmserfV4 <- RMSE(predrfV4, datV4test[,5])

c(rmserfV2, rmserfV3, rmserfV4)

# r2

r2rfV2 <- R2(predrfV2, datV2test[,3])
r2rfV3 <- R2(predrfV3, datV3test[,4])
r2rfV4 <- R2(predrfV4, datV4test[,5])

c(r2rfV2, r2rfV3, r2rfV4)


# ----- NN

# models

mrfV2 <- randomForest(Y ~ ., data=datV2train, mtry = 2, importance=TRUE, na.action=na.omit)
mrfV3 <- randomForest(Y ~ ., data=datV3train, mtry = 2, importance=TRUE, na.action=na.omit)
mrfV4 <- randomForest(Y ~ ., data=datV4train, mtry = 2, importance=TRUE, na.action=na.omit)

# predictions

predrfV2 <- predict(mrfV2, data = datV2test[,-3])
predrfV3 <- predict(mrfV3, data = datV3test[,-3])
predrfV4 <- predict(mrfV4, data = datV4test[,-3])

# rmse

rmserfV2 <- RMSE(predrfV2, datV2test[,3])
rmserfV3 <- RMSE(predrfV3, datV3test[,3])
rmserfV4 <- RMSE(predrfV4, datV4test[,3])






n <- names(datV2)
f <- as.formula(paste("Y ~", paste(n[!n %in% "Y"], collapse = " + ")))
nn <- neuralnet(f, data=datV2)

nn <- neuralnet(Y ~ var1 + var2,
                data=datV2, hidden=c(2,1), linear.output=FALSE, threshold=0.01)
nn$result.matrix
plot(nn)


m <- model.matrix(
  ~ Y + var1 + var2,
  data = datV2
)

nam <- colnames(m)[-c(1:2)]
f <- as.formula(paste("Y ~", paste(nam, collapse = " + ")))
r <- neuralnet(f,data=m)



m <- randomForest(Y ~ ., data=datV2, mtry = 2,
             importance=TRUE, na.action=na.omit)

m$ntree


library(caret)

set.seed(022)
modelrfV2 <- train(Y~.,
                   data = datV2,
                   method = "rf",
                   ntree = 500,
                   prox=TRUE,
                   metric='RMSE',
                   allowParallel=TRUE)
modelrfV2$results




set.seed(022)
modelbrnnV2 <- train(Y~.,
                     data = datV2,
                     method = "nnet",
                     #ntree = 100,
                     prox=TRUE,
                     metric='RMSE',
                     allowParallel=TRUE)
modelbrnnV2




# Make training and validation data
set.seed(1)
train <- sample(nrow(iris), nrow(iris)*0.5)
valid <- seq(nrow(iris))[-train]
iristrain <- iris[train,]
irisvalid <- iris[valid,]

# Binarize the categorical output
iristrain <- cbind(iristrain, iristrain$Species == 'setosa')
iristrain <- cbind(iristrain, iristrain$Species == 'versicolor')
iristrain <- cbind(iristrain, iristrain$Species == 'virginica')
names(iristrain)[6:8] <- c('setosa', 'versicolor', 'virginica')

# Fit model
nn <- neuralnet(
  setosa+versicolor+virginica ~ Sepal.Length + Sepal.Width + Petal.Length + Petal.Width,
  data=iristrain,
  hidden=c(3)
)
plot(nn)


nn <- neuralnet(
  setosa+versicolor+virginica ~ Sepal.Length + Sepal.Width + Petal.Length + Petal.Width,
  data=Boston,
  hidden=c(3)
)
plot(nn)


library(neuralnet)
n <- names(Boston)
f <- as.formula(paste("medv ~", paste(n[!n %in% "medv"], collapse = " + ")))
nn <- neuralnet(f,data=Boston,hidden=c(5,3), linear.output=T)

# Predictions -------------------------------------------------------------

datV2test <- transdat(data = V2N120Q2)

m <- randomForest(Y ~ ., data=datV2, mtry = 2,
                  importance=TRUE, na.action=na.omit)
predrfV2 <- predict(m, data = datV2test[,-3])
RMSE(predrfV2, datV2test[,3])

