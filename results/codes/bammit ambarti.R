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
