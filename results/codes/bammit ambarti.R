#devtools::install_github("ebprado/AMBARTI/R package", ref='main')
library(AMBARTI)
library(bammit)

load("~/Documents/GitHub/bammit/Real data/train_ireland.RData")
load("~/Documents/GitHub/bammit/Real data/test_ireland.RData")
df = train
df$g = gsub('g','', df$Genotype)
df$e = gsub('e','', df$Environment)
df$t = as.factor(df$Year)
df$y = df$Yield
df = df[,-which(colnames(df) %in% c('Genotype','Environment','Year','Yield', 'Bloc'))]
df = as.data.frame(df)
y = df$y
x = df[,-which(colnames(df) == 'y')]
fit.ambarti = alessa_ambarti(x, y, ntrees = 50, nburn = 1000, npost = 2000, nsteps = 1)
saveRDS(fit.ambarti, "~/Documents/GitHub/bammit/Running models/fit.ambarti.RData")

qq = var_used_trees(fit.ambarti)
df2 = test
df2$g = gsub('g','', df2$Genotype)
df2$e = gsub('e','', df2$Environment)
df2$t = as.factor(df2$Year)
df2$y = df2$Yield
df2 = df2[,-which(colnames(df2) %in% c('Genotype','Environment','Year','Yield', 'Bloc'))]
df2 = as.data.frame(df2)
y_test = df2$y
x_test = df2[,-which(colnames(df2) == 'y')]
yhat_ambarti2 = predict_ambarti_alessa(object=fit.ambarti, newdata = x_test, type = 'mean')
saveRDS(yhat_ambarti2, "~/Documents/GitHub/bammit/Running models/yhat_ambarti2.RData")

RMSE(y_test, yhat_ambarti2)
R2(y_test, yhat_ambarti2)


# AMBARTI 2 variables -----------------------------------------------------

# all years
df = train
df$g = gsub('g','', df$Genotype)
df$e = gsub('e','', df$Environment)
df$y = df$Yield
df = df[,-which(colnames(df) %in% c('Genotype','Environment','Year','Yield', 'Bloc'))]
df = as.data.frame(df)
y = df$y
x = df[,-which(colnames(df) == 'y')]
fit.ambarti.ammi = ambarti(x, y, ntrees = 50, nburn = 1000, npost = 2000, nsteps = 1)
saveRDS(fit.ambarti.ammi, "~/Documents/GitHub/bammit/Running models/fit.ambarti.ammi.RData")

qq = var_used_trees(fit.ambarti.ammi)
df2 = test
df2$g = gsub('g','', df2$Genotype)
df2$e = gsub('e','', df2$Environment)
df2$y = df2$Yield
df2 = df2[,-which(colnames(df2) %in% c('Genotype','Environment','Year','Yield', 'Bloc'))]
df2 = as.data.frame(df2)
y_test = df2$y
x_test = df2[,-which(colnames(df2) == 'y')]
yhat_ambarti2_ammi = predict_ambarti_alessa(object=fit.ambarti.ammi, newdata = x_test, type = 'mean')
saveRDS(yhat_ambarti2_ammi, "~/Documents/GitHub/bammit/Running models/yhat_ambarti2.RData")

caret::RMSE(y_test, yhat_ambarti2_ammi[,1])
caret::R2(y_test, yhat_ambarti2_ammi[,1])

## by year

train_data <- split(train, train$Year)
test_data <- split(test, test$Year)

# run the model

model <- lapply(train_data, function(i){
  df = i
  df$g = gsub('g','', df$Genotype)
  df$e = gsub('e','', df$Environment)
  df$y = df$Yield
  df = df[,-which(colnames(df) %in% c('Genotype','Environment','Year','Yield', 'Bloc'))]
  df = as.data.frame(df)
  y = df$y
  x = df[,-which(colnames(df) == 'y')]
  fit.ambarti.ammi = ambarti(x, y, ntrees = 50, nburn = 1000, npost = 1000, nsteps = 1)
  return(fit.ambarti.ammi)
})

#saveRDS(model, file = "Running models/modelAmbartiByYear.rds")
modelAmbartiByYear <- readRDS("~/Documents/GitHub/bammit/Running models/modelAmbartiByYear.rds")

# get predictions

tidy_pred <- function(mod, data){
  qq = var_used_trees(mod)
  df2 = data
  df2$g = gsub('g','', df2$Genotype)
  df2$e = gsub('e','', df2$Environment)
  df2$y = df2$Yield
  df2 = df2[,-which(colnames(df2) %in% c('Genotype','Environment','Year','Yield', 'Bloc'))]
  df2 = as.data.frame(df2)
  y_test = df2$y
  x_test = df2[,-which(colnames(df2) == 'y')]
  return(x_test)
}

pred <- list()
x_test <- list()
for(i in 1:10){
  x_test[[i]] <- tidy_pred(mod = modelAmbartiByYear[[i]], data = test_data[[i]])
  pred[[i]] = predict_ambarti_alessa(object = modelAmbartiByYear[[i]], newdata = x_test[[i]], type = 'mean')
  print(paste0('pred', i))
}


qq = var_used_trees(model$`2011`)
df2 = test_data$`2011`
df2$g = gsub('g','', df2$Genotype)
df2$e = gsub('e','', df2$Environment)
df2$y = df2$Yield
df2 = df2[,-which(colnames(df2) %in% c('Genotype','Environment','Year','Yield', 'Bloc'))]
df2 = as.data.frame(df2)
y_test = df2$y
x_test = df2[,-which(colnames(df2) == 'y')]
yhat_ambarti2_ammi = predict_ambarti_alessa(object=model$`2011`, newdata = x_test, type = 'mean')


newPred <- readRDS("~/Documents/GitHub/bammit/Running models/newPred.rds")
caret::RMSE(unlist(newPred), test$Yield)
caret::R2(unlist(newPred), test$Yield)






