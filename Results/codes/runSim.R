# --------------------- Running simulations ----------------- #

# packages

library(bammit)
library(ggplot2)
library(dplyr)
library(reshape2)
library(randomForest) # to run RF
library(xgboost) # to run xboost
library(arm) # to run bayesian glm
library(AMBARTI)

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


# Reading just some of the datasets ---------------------------------------------

V2N120Q2_train <- train[[2]]
V3N480Q2_train <- train[[5]]
V4N960Q2_train <- train[[11]]
V4N960Q1_train <- train[[10]]
#V3N5000Q2_train <- train[[8]]

V2N120Q2_test <- test[[2]]
V3N480Q2_test <- test[[5]]
V4N960Q2_test <- test[[11]]
V4N960Q1_test <- test[[10]]
#V3N5000Q2_test <- test[[8]]


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
# base_dir_out <- '~/Documents/GitHub/bammit/Running models'
# all_files_out <- list.files(base_dir_out, pattern="*.rds", full.names = TRUE)
# out <- lapply(all_files_out, readRDS)
# names(out) <- sub(pattern = "(.*)\\..*$", replacement = "\\1", basename(all_files_out))

# just the models used on the paper
modelV4N960Q1QM1 <- readRDS("/Volumes/Alessa HD/PhD/bammit/Running models/modelV4N960Q1QM1.rds")
modelV2N120Q2QM2 <- readRDS("/Volumes/Alessa HD/PhD/bammit/Running models/modelV2N120Q2QM2.rds")
modelV3N480Q2QM2 <- readRDS("/Volumes/Alessa HD/PhD/bammit/Running models/modelV3N480Q2QM2.rds")
modelV4N960Q2QM2 <- readRDS("/Volumes/Alessa HD/PhD/bammit/Running models/modelV4N960Q2QM2.rds")
#modelV3N5000Q2QM2 <- readRDS("/Volumes/Alessa HD/PhD/bammit/Running models/modelV3N5000Q2QM2.rds")


# Plots -------------------------------------------------------------------

# Fig1 : true vs estimated additive term scenario (iii), V = 4, Q = 2

# df estimated values
dfEstV4 <- list(b1 = modelV4N960Q2QM2$mean$b1,
                b2 = modelV4N960Q2QM2$mean$b2,
                b3 = modelV4N960Q2QM2$mean$b3,
                b4 = modelV4N960Q2QM2$mean$b4) |>
  melt() |> select(L1, value)

# df simulated data
dfDataV4 <- list(b1 = V4N960Q2_train$bv[[1]],
                 b2 = V4N960Q2_train$bv[[2]],
                 b3 = V4N960Q2_train$bv[[3]],
                 b4 = V4N960Q2_train$bv[[4]]) |>
  melt() |> select(L1, value)

dfMainEffV4 <- cbind(dfEstV4, dfDataV4)[,c(1,2,4)]
colnames(dfMainEffV4) <- c("var", "est", "true")

quantV4mainEff <- list(
  b1 = data.frame(q1 = apply(modelV4N960Q2QM2$sims.list$b1, 2, function(x) quantile(x, 0.05)),
                  q2 = apply(modelV4N960Q2QM2$sims.list$b1, 2, function(x) quantile(x, 0.95))),
  b2 = data.frame(q1 = apply(modelV4N960Q2QM2$sims.list$b2, 2, function(x) quantile(x, 0.05)),
                  q2 = apply(modelV4N960Q2QM2$sims.list$b2, 2, function(x) quantile(x, 0.95))),
  b3 = data.frame(q1 = apply(modelV4N960Q2QM2$sims.list$b3, 2, function(x) quantile(x, 0.05)),
                  q2 = apply(modelV4N960Q2QM2$sims.list$b3, 2, function(x) quantile(x, 0.95))),
  b4 = data.frame(q1 = apply(modelV4N960Q2QM2$sims.list$b4, 2, function(x) quantile(x, 0.05)),
                  q2 = apply(modelV4N960Q2QM2$sims.list$b4, 2, function(x) quantile(x, 0.95)))
) |> plyr::ldply()

dfMainEffV4 <- cbind(dfMainEffV4, quantV4mainEff$q1, quantV4mainEff$q2 )
colnames(dfMainEffV4) <- c("var", "est", "true", "q1", "q2")


dfMainEffV4$facet <- factor(dfMainEffV4$var,
                            labels = c(expression(bold(b)^{(1)}),
                                       expression(bold(b)^{(2)}),
                                       expression(bold(b)^{(3)}),
                                       expression(bold(b)^{(4)})))

dfMainEffV4 |>
  ggplot(aes(x = true, y = est)) +
  geom_abline(size = 0.3, col = "gray") +
  geom_linerange(aes(ymin =  q1, ymax = q2), alpha = 0.5, size = 0.4) +
  geom_point(colour =  "steelblue", size = 1) +
  facet_wrap(~facet, ncol = 4, labeller = label_parsed) +
  labs(x = "True", y = "Estimated") +
  theme_bw(base_size = 16) +
  xlim(-2,2)+
  theme(strip.background = element_blank(),
        panel.spacing.x = unit(5,"mm"),
        legend.position = c(0.75, 0.04))




# Fig2 : true vs estimated interaction term scenarios (i), (ii) and (iii), V = 2,3,4, Q = 2

# simulated data
VQ2Dataint <- list(V2 = V2N120Q2_train$int,
                   V3 = V3N480Q2_train$int,
                   V4 = V4N960Q2_train$int) |>
  melt() |>
  select(value, L1)

VQ2int <- list(V2 = modelV2N120Q2QM2$mean$blin,
               V3 = modelV3N480Q2QM2$mean$blin,
               V4 = modelV4N960Q2QM2$mean$blin) |>
  melt() |>
  select(value, L1)

quantVQ2int <- list(
  V2 = data.frame(q1 = apply(modelV2N120Q2QM2$sims.list$blin, 2, function(x) quantile(x, 0.05)),
                  q2 = apply(modelV2N120Q2QM2$sims.list$blin, 2, function(x) quantile(x, 0.95))),
  V3 = data.frame(q1 = apply(modelV3N480Q2QM2$sims.list$blin, 2, function(x) quantile(x, 0.05)),
                  q2 = apply(modelV3N480Q2QM2$sims.list$blin, 2, function(x) quantile(x, 0.95))),
  V4 = data.frame(q1 = apply(modelV4N960Q2QM2$sims.list$blin, 2, function(x) quantile(x, 0.05)),
                  q2 = apply(modelV4N960Q2QM2$sims.list$blin, 2, function(x) quantile(x, 0.95)))
) |> plyr::ldply()

VQ2int <- cbind(VQ2int, quantVQ2int$q1, quantVQ2int$q2)
colnames(VQ2int) <- c("value", "L1", "q1", "q2")

dfVQ2int <- data.frame(v = VQ2int$L1, est = VQ2int$value, true = VQ2Dataint$value,
                       q1 = VQ2int$q1, q2 = VQ2int$q2)

lab <- function(x) c("V = 2", "V = 3", "V = 4")
dfVQ2int |>
  ggplot(aes(x = true, y = est)) +
  geom_abline(size = 0.3, col = "gray") +
  geom_linerange(aes(ymin =  q1, ymax = q2), alpha = 0.2, size = 0.1) +
  geom_point(colour =  "steelblue", size = 0.2) +
  facet_wrap(~v,
             labeller = as_labeller(lab)) +
  labs(x = "True", y = "Estimated") +
  theme_bw(base_size = 14) +
  theme(strip.background = element_blank(),
        panel.spacing.x = unit(5,"mm"),
        legend.position = c(0.75, 0.04))


# Fitting ML models -------------------------------------------------------

datV2train <- transdat(data = V2N120Q2_train)
datV3train <- transdat(data = V3N480Q2_train) # N = 480
datV4train <- transdat(data = V4N960Q2_train)

datV2test <- transdat(data = V2N120Q2_test)
datV3test <- transdat(data = V3N480Q2_test) # N = 480
datV4test <- transdat(data = V4N960Q2_test)


# ----- Random Forest

# models

set.seed(02)
mrfV2 <- randomForest(Y ~ ., data=datV2train, mtry = 2, importance=TRUE, na.action=na.omit)
mrfV3 <- randomForest(Y ~ ., data=datV3train, mtry = 2, importance=TRUE, na.action=na.omit)
mrfV4 <- randomForest(Y ~ ., data=datV4train, mtry = 2, importance=TRUE, na.action=na.omit)

# predictions

predrfV2 <- predict(mrfV2, data = datV2test[,-3])
predrfV3 <- predict(mrfV3, data = datV3test[,-4])
predrfV4 <- predict(mrfV4, data = datV4test[,-5])

# rmse

rmserfV2 <- caret::RMSE(predrfV2, datV2test[,3])
rmserfV3 <- caret::RMSE(predrfV3, datV3test[,4])
rmserfV4 <- caret::RMSE(predrfV4, datV4test[,5])

caret::RMSE(predrfV2, datV2test[,3])
caret::R2(predrfV2, datV2test[,3])

c(rmserfV2, rmserfV3, rmserfV4)

# r2

r2rfV2 <- caret::R2(predrfV2, datV2test[,3])
r2rfV3 <- caret::R2(predrfV3, datV3test[,4])
r2rfV4 <- caret::R2(predrfV4, datV4test[,5])

c(r2rfV2, r2rfV3, r2rfV4)



# xboost ------------------------------------------------------------------

# x and y for train and test  data

# V2
train_x_V2 = data.matrix(datV2train[, -3])
train_y_V2 = datV2train[,3]

test_x_V2 = data.matrix(datV2test[, -3])
test_y_V2 = datV2test[, 3]

xgb_train_V2 = xgb.DMatrix(data = train_x_V2, label = train_y_V2)
xgb_test_V2 = xgb.DMatrix(data = test_x_V2, label = test_y_V2)

# V3
train_x_V3 = data.matrix(datV3train[, -4])
train_y_V3 = datV3train[,4]

test_x_V3 = data.matrix(datV3test[, -4])
test_y_V3 = datV3test[, 4]

xgb_train_V3 = xgb.DMatrix(data = train_x_V3, label = train_y_V3)
xgb_test_V3 = xgb.DMatrix(data = test_x_V3, label = test_y_V3)


# V4
train_x_V4 = data.matrix(datV4train[, -5])
train_y_V4 = datV4train[,5]

test_x_V4 = data.matrix(datV4test[, -5])
test_y_V4 = datV4test[, 5]

xgb_train_V4 = xgb.DMatrix(data = train_x_V4, label = train_y_V4)
xgb_test_V4 = xgb.DMatrix(data = test_x_V4, label = test_y_V4)

# models

set.seed(02)
mxbV2 <- xgboost(data = xgb_train_V2, nrounds = 50)
mxbV3 <- xgboost(data = xgb_train_V3, nrounds = 50)
mxbV4 <- xgboost(data = xgb_train_V4, nrounds = 50)

# predictions

predxbV2 <- predict(mxbV2, xgb_test_V2)
predxbV3 <- predict(mxbV3, xgb_test_V3)
predxbV4 <- predict(mxbV4, xgb_test_V4)

# rmse

rmsexbV2 <- caret::RMSE(predxbV2, test_y_V2)
rmsexbV3 <- caret::RMSE(predxbV3, test_y_V3)
rmsexbV4 <- caret::RMSE(predxbV4, test_y_V4)

cat("V2: ", rmsexbV2,
    "V3: ", rmsexbV3,
    "V4: ", rmsexbV4)

# r2

r2rfV2 <- caret::R2(predxbV2, test_y_V2)
r2rfV3 <- caret::R2(predxbV3, test_y_V3)
r2rfV4 <- caret::R2(predxbV4, test_y_V4)

cat("V2: ", r2rfV2,
    "V3: ", r2rfV3,
    "V4: ", r2rfV4)



# bayesian glm ------------------------------------------------------------


# models

set.seed(02)
mbglmV2 <- bayesglm(Y ~ var1 + var2 , data = datV2train)
mbglmV3 <- bayesglm(Y ~ var1 + var2 + var3 , data = datV3train)
mbglmV4 <- bayesglm(Y ~ var1 + var2 + var4, data = datV4train)

# predictions

predbglmV2 <- predict(mbglmV2, data = datV2test[,-3])
predbglmV3 <- predict(mbglmV3, data = datV3test[,-4])
predbglmV4 <- predict(mbglmV4, data = datV4test[,-5])

# rmse

rmsebglmV2 <- caret::RMSE(predbglmV2, datV2test[,3])
rmsebglmV3 <- caret::RMSE(predbglmV3, datV3test[,4])
rmsebglmV4 <- caret::RMSE(predbglmV4, datV4test[,5])

cat("V2: ", rmsebglmV2,
    "V3: ", rmsebglmV3,
    "V4: ", rmsebglmV4)

# r2

r2bglmV2 <- caret::R2(predbglmV2, datV2test[,3])
r2bglmV3 <- caret::R2(predbglmV3, datV3test[,4])
r2bglmV4 <- caret::R2(predbglmV4, datV4test[,5])

cat("V2: ", r2bglmV2,
    "V3: ", r2bglmV3,
    "V4: ", r2bglmV4)




# AMMI --------------------------------------------------------------------


colnames(datV3train)[c(1,2,4)] <- c('Genotype', 'Environment','Yield')
colnames(datV4train)[c(1,2,5)] <- c('Genotype', 'Environment','Yield')

colnames(datV3test)[c(1,2,4)] <- c('Genotype', 'Environment','Yield')
colnames(datV4test)[c(1,2,5)] <- c('Genotype', 'Environment','Yield')


ammi_datV3train <- AMMIJagsRealData(datV3train, Q = 1, mmu = 100, smu = 10, stheta = 1, a = 0.1, b = 0.1,
                                    nthin = 2, nburnin = 2000)

ammi_datV4train <- AMMIJagsRealData(datV4train, Q = 1, mmu = 100, smu = 10, stheta = 1, a = 0.1, b = 0.1,
                                    nthin = 2, nburnin = 2000)

saveRDS(ammi_datV3train, "/Volumes/Alessa HD/PhD/bammit/Running models/ammi_datV3train.rds")
saveRDS(ammi_datV4train, "/Volumes/Alessa HD/PhD/bammit/Running models/ammi_datV4train.rds")

pred_datV3train <- predictionAMMIReal(model = ammi_datV3train$BUGSoutput, datV3test)
pred_datV4train <- predictionAMMIReal(model = ammi_datV4train$BUGSoutput, datV4test)

caret::RMSE(pred_datV3train, datV3test$Yield)
caret::R2(pred_datV3train, datV3test$Yield)


caret::RMSE(pred_datV4train, datV4test$Yield)
caret::R2(pred_datV4train, datV4test$Yield)

# AMBARTI -----------------------------------------------------------------


# V = 3

df = datV3train
df$g = gsub('g','', df$Genotype)
df$e = gsub('e','', df$Environment)
df$y = df$Yield
df = df[,-which(colnames(df) %in% c('Genotype','Environment','var3','Yield'))]
df = as.data.frame(df)
y = df$y
x = df[,-which(colnames(df) == 'y')]
fit.ambarti.ammi = ambarti(x, y, ntrees = 50, nburn = 500, npost = 1000, nsteps = 1)
#saveRDS(fit.ambarti.ammi, "~/Documents/GitHub/bammit/Running models/ambarti_datV3train.RData")

qq = var_used_trees(fit.ambarti.ammi)
df2 = datV3test
df2$g = gsub('g','', df2$Genotype)
df2$e = gsub('e','', df2$Environment)
df2$y = df2$Yield
df2 = df2[,-which(colnames(df2) %in% c('Genotype','Environment','var3','Yield'))]
df2 = as.data.frame(df2)
y_test = df2$y
x_test = df2[,-which(colnames(df2) == 'y')]
yhat_ambarti2_ammi = predict_ambarti_alessa(object=fit.ambarti.ammi, newdata = x_test, type = 'mean')
saveRDS(yhat_ambarti2_ammi, "~/Documents/GitHub/bammit/Running models/yhat_ambarti_datV3test.RData")

caret::RMSE(y_test, yhat_ambarti2_ammi[,1])
caret::R2(y_test, yhat_ambarti2_ammi[,1])


# V = 4

df = datV4train
df$g = gsub('g','', df$Genotype)
df$e = gsub('e','', df$Environment)
df$y = df$Yield
df = df[,-which(colnames(df) %in% c('Genotype','Environment','var3', 'var4' ,'Yield'))]
df = as.data.frame(df)
y = df$y
x = df[,-which(colnames(df) == 'y')]
fit.ambarti.ammi = ambarti(x, y, ntrees = 50, nburn = 500, npost = 1000, nsteps = 1)
saveRDS(fit.ambarti.ammi, "~/Documents/GitHub/bammit/Running models/ambarti_datV4train.RData")

qq = var_used_trees(fit.ambarti.ammi)
df2 = datV4test
df2$g = gsub('g','', df2$Genotype)
df2$e = gsub('e','', df2$Environment)
df2$y = df2$Yield
df2 = df2[,-which(colnames(df2) %in% c('Genotype','Environment','var3', 'var4','Yield'))]
df2 = as.data.frame(df2)
y_test = df2$y
x_test = df2[,-which(colnames(df2) == 'y')]
yhat_ambarti2_ammi = predict_ambarti_alessa(object=fit.ambarti.ammi, newdata = x_test, type = 'mean')
saveRDS(yhat_ambarti2_ammi, "~/Documents/GitHub/bammit/Running models/yhat_ambarti_datV4test.RData")

caret::RMSE(y_test, yhat_ambarti2_ammi[,1])
caret::R2(y_test, yhat_ambarti2_ammi[,1])





# ----- NN

# tidy formula (factor to dummies)

# V2
# datV2train_matrix <- model.matrix(~var1 + var2 + Y, data=datV2train)
# col_list_datV2train <- paste(c(colnames(datV2train_matrix[, -c(1,ncol(datV2train_matrix))])),collapse="+")
# col_list_datV2train <- paste(c("Y~",col_list_datV2train),collapse="")
# f_datV2train <- formula(col_list_datV2train)
# datV2test_matrix <- model.matrix(~var1 + var2 + Y, data=datV2test)
# col_list_datV2test <- paste(c(colnames(datV2test_matrix[, -c(1,ncol(datV2test_matrix))])),collapse="+")
# col_list_datV2test <- paste(c("Y~",col_list_datV2test),collapse="")
# f_datV2test <- formula(col_list_datV2test)
# # models
# set.seed(02)
# mnnV2 <- neuralnet(f_datV2train, data = datV2train_matrix, hidden = 3)
# # predictions
# prednnV2 <- predict(mnnV2, newdata = datV2train_matrix)
# # trying nnet
# library(nnet)
# library(NeuralNetTools)
# model <- nnet(datV2test[,c('var1', 'var2')], datV2test[,'Y'], size = 3)
# model <- nnet(Y ~ var1 + var2, size = 3, data = datV2test)
# model <- neuralnet(Y ~ var1 + var2, data = datV2test)





# library(caret)
# set.seed(022)
# modelrfV2 <- train(Y~.,
#                    data = datV2,
#                    method = "rf",
#                    ntree = 500,
#                    prox=TRUE,
#                    metric='RMSE',
#                    allowParallel=TRUE)
# set.seed(02)
# modelbrnnV2 <- train(Y~.,
#                      data = datV2train,
#                      method = "nnet",
#                      metric='RMSE',
#                      allowParallel=TRUE)


