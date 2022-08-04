devtools::install_github("ebprado/AMBARTI/R package", ref='main')
library(AMBARTI)

load("Real data/ireland.RData")
df = ireland
df$g = gsub('g','', df$Genotype)
df$e = gsub('e','', df$Environment)
df$t = as.factor(df$Year)
df$y = df$Mean
df = df[,-which(colnames(df) %in% c('Genotype','Environment','Year','Mean'))]
df = as.data.frame(df)
y = df$y
x = df[,-which(colnames(df) == 'y')]
fit.ambarti = alessa_ambarti(x, y, ntrees = 50, nburn = 200, npost = 100, nsteps = 1)
qq = var_used_trees(fit.ambarti)

fit.ambarti$y_hat

