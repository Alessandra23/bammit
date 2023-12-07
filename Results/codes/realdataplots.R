library(stringr)
library(ggplot2)
library(bammit)
# load("Real data/train_ireland.RData")
# load("Real data/test_ireland.RData")

# modeltrainQ1 <- readRDS("~/Documents/GitHub/bammit/Running models/trainQ1.rds")
# model <- modeltrainQ1
# data <- test


load("~/Documents/GitHub/bammit/Real data/ireland.RData")
load("~/Documents/GitHub/bammit/Running models/model_Q4.RData")

ireland$Bloc <- gsub('.*_(\\d)$', 'b\\1', ireland$Bloc)

# separate the data set in train and test
train <- subset(ireland,  grepl('1$|2$|3$', Bloc))
test <- subset(ireland,  grepl('4$', Bloc))

train$Bloc <-train$Bloc |> as.factor()
test$Bloc <- test$Bloc |>  as.factor()

levels(train$Bloc)
levels(test$Bloc)

test_original <- test

train$Yield <- scale(train$Yield) |> as.numeric()
test$Yield <- scale(test$Yield) |> as.numeric()



model <- model_Q4
data <- test_original

# pred <- predictionBAMMITReal(model = model$BUGSoutput, data = data)
# caret::RMSE(pred,test$Yield)
# caret::R2(pred,test$Yield)
# model$BUGSoutput

genotype <- data$Genotype
environment <- data$Environment
time <- data$Year
bloc <- data$Bloc
muhat <- model$BUGSoutput$mean$muall
g <- model$BUGSoutput$sims.list$g
e <- model$BUGSoutput$sims.list$e
t <- model$BUGSoutput$sims.list$t
bl <- model$BUGSoutput$sims.list$bl
inthat <- model$BUGSoutput$sims.list$blin
N <- length(data$Yield)

# yhat <- matrix(0, nrow = nrow(g), ncol = N)
# for(i in 1:nrow(g)){
#   yhat[i,] <- rep(muhat, N) + g[i, genotype] + e[i, environment] + t[i, time] + bl[i, bloc]  +  inthat[i,]
# }
#
# inthat_eff <- matrix(0, nrow = nrow(g), ncol = N)
# for(i in 1:nrow(g)){
#   inthat_eff[i,] <- g[i, genotype] + inthat[i,]
# }

#
# yhat <- model$BUGSoutput$median$mu_test
# inthat_eff <- model$BUGSoutput$median$int_test


# set up names

g <- data$Genotype |> as.character()
e <- data$Environment |> as.character()
t <- data$Year |> as.character()
b <- data$Bloc |> as.character()

# get yhats
# pred <- model$BUGSoutput$sims.list$mu_test
# yhat <- apply(pred, 2, median) * sd(data$Yield) + mean(data$Yield)# predict values for yhat in the original scale
# # get sd
# yhatSD <-  apply(pred, 2, sd) * sd(data$Yield)


years_pred <- list()
for(i in 1:ncol(model$BUGSoutput$sims.list$mu_test)){
  years_pred[[i]] <- model$BUGSoutput$sims.list$mu_test[,i]*sd(data$Yield) + mean(data$Yield)
}
length(years_pred)

yhat <- lapply(years_pred, median) |> unlist()
yhatSD <- lapply(years_pred, sd) |> unlist()

caret::RMSE(yhat, data$Yield)
caret::R2(yhat, data$Yield)

# create a df
df.ambarti <- data.frame(yhat = yhat, e = e, g = g, t = t, b = b, sd = yhatSD)
#saveRDS(df.ambarti, "Running models/bammitVsupOut.rds")
#bammitVsupOut <- readRDS("~/Documents/GitHub/bammit/Running models/bammitVsupOut.rds")
#df.ambarti <- bammitVsupOut
df.ambarti <- subset(df.ambarti,  t == c("2019"))
#df.ambarti <- subset(df.ambarti,  grepl('3$', b))
df.ambarti <- df.ambarti |> dplyr::select("yhat", "e", "g", "t", "sd")
df.ambarti$e <- factor(df.ambarti$e, levels = unique(df.ambarti$e))
df.ambarti$g <- factor(df.ambarti$g, levels = unique(df.ambarti$g))
df.ambarti$t <- factor(df.ambarti$t, levels = unique(df.ambarti$t))
#df.ambarti$b <- factor(df.ambarti$b, levels = unique(df.ambarti$b))



fill <- df.ambarti$yhat
breakPattern_g <- as.numeric(str_extract_all(levels(reorder(df.ambarti$g, fill)), "[0-9]+"))
breakPattern_e <- as.numeric(str_extract_all(levels(reorder(df.ambarti$e, fill)), "[0-9]+"))
breakPattern_t <- as.numeric(str_extract_all(levels(reorder(df.ambarti$t, fill)), "[0-9]+"))
# breakPattern_b <- as.numeric(str_extract_all(levels(b), "[0-9]+"))

breakPattern_e <- rev(breakPattern_e)



# vsup plots --------------------------------------------------------------

pal = rev(colorspace::diverging_hcl(palette = "Blue-Red 3", n = 8, alpha = 4))

# colors <- scales::colour_ramp(
#   colors = c(blue = 'steelblue', red = 'firebrick')
# )((0:7)/7)
# wheel <- function(col, radius = 1)
#   pie(rep(1, length(col)), col = col, radius = radius)
# wheel(colors)
# pal <- colors

# r_yhat <- range(df.ambarti$yhat)
# r_yhat[1] <- floor(r_yhat[1])
# r_yhat[2] <- ceiling(r_yhat[2])
#
# r_sd <- range(df.ambarti$sd)
# # r_sd[1] <- floor(r_sd[1])
# # r_sd[2] <- ceiling(r_sd[2])
# r_sd[1] <- floor(r_sd[1] * 10) / 10
# r_sd[2] <- ceiling(r_sd[2] * 100) / 100
# r_sd


range(df.ambarti$yhat)
r_yhat <- c(9,14)
range(df.ambarti$sd)
r_sd <- c(0,0.6)

pp <- ggplot(df.ambarti) +
  #geom_raster(aes(e, g, fill = zip(yhat, sd))) +
  geom_raster(aes(reorder(e, -fill), reorder(g, fill), fill = zip(yhat, sd))) +
  bivariate_scale(
    name = c("\U0177", "sd"),
    aesthetics = "fill",
    limits = list(r_yhat, r_sd),
    palette = pal_vsup(
      values = pal,
      unc_levels = 4,
      max_desat = 0.6,
      pow_desat = 0.2,
      max_light = 0.6,
      pow_light = 1
    ),
    guide = "colorfan"
  ) +
  theme_bw(base_size = 11) +
  labs(x = "Environment", y = "Genotype") #+
#facet_wrap(~t,scales = "free")

pp


ggsave("~/Documents/GitHub/bammit/Results/figures paper/real data/vsup2019.pdf", width = 7.55, height = 5.04, device = cairo_pdf)
dev.off()
# dim: 8.55 and 6.04



# intercation -------------------------------------------------------------



# set up names

g <- data$Genotype |> as.character()
e <- data$Environment |> as.character()
t <- data$Year |> as.character()
b <- data$Bloc |> as.character()

# get yhats
pred <- inthat_eff
yhat <- apply(pred, 2, median) # predict values for yhat
# get sd
yhatSD <-  apply(pred, 2, sd)

# create a df
df.ambarti <- data.frame(yhat = yhat, e = e, g = g, t = t, b = b, sd = yhatSD)
# saveRDS(df.ambarti, "Running models/bammitVsupOut.rds")
# bammitVsupOut <- readRDS("~/Documents/GitHub/bammit/Running models/bammitVsupOut.rds")
# df.ambarti <- bammitVsupOut
df.ambarti <- subset(df.ambarti,  t == c("2015"))
df.ambarti <- subset(df.ambarti,  grepl('3$', b))
df.ambarti <- df.ambarti |> dplyr::select("yhat", "e", "g", "t", "sd")
df.ambarti$e <- factor(df.ambarti$e, levels = unique(df.ambarti$e))
df.ambarti$g <- factor(df.ambarti$g, levels = unique(df.ambarti$g))
df.ambarti$t <- factor(df.ambarti$t, levels = unique(df.ambarti$t))
#df.ambarti$b <- factor(df.ambarti$b, levels = unique(df.ambarti$b))



fill <- df.ambarti$yhat
breakPattern_g <- as.numeric(str_extract_all(levels(reorder(df.ambarti$g, fill)), "[0-9]+"))
breakPattern_e <- as.numeric(str_extract_all(levels(reorder(df.ambarti$e, fill)), "[0-9]+"))
breakPattern_t <- as.numeric(str_extract_all(levels(reorder(df.ambarti$t, fill)), "[0-9]+"))
# breakPattern_b <- as.numeric(str_extract_all(levels(b), "[0-9]+"))

breakPattern_e <- rev(breakPattern_e)



# vsup plots --------------------------------------------------------------

pal = rev(colorspace::diverging_hcl(palette = "Blue-Red 3", n = 8, alpha = 4))

# colors <- scales::colour_ramp(
#   colors = c(blue = 'steelblue', red = 'firebrick')
# )((0:7)/7)
# wheel <- function(col, radius = 1)
#   pie(rep(1, length(col)), col = col, radius = radius)
# wheel(colors)
# pal <- colors

# r_yhat <- range(df.ambarti$yhat)
# r_yhat[1] <- floor(r_yhat[1])
# r_yhat[2] <- ceiling(r_yhat[2])
#
# r_sd <- range(df.ambarti$sd)
# # r_sd[1] <- floor(r_sd[1])
# # r_sd[2] <- ceiling(r_sd[2])
# r_sd[1] <- floor(r_sd[1] * 10) / 10
# r_sd[2] <- ceiling(r_sd[2] * 100) / 100
# r_sd


r_yhat <- range(df.ambarti$yhat)
r_yhat <- c(8,16)
r_sd <- range(df.ambarti$sd)
r_sd <- c(0.48,0.6)

pp <- ggplot(df.ambarti) +
  #geom_raster(aes(e, g, fill = zip(yhat, sd))) +
  geom_raster(aes(reorder(e, -fill), reorder(g, fill), fill = zip(yhat, sd))) +
  bivariate_scale(
    name = c(expression(hat(g) + hat(g x e)), "sd"),
    aesthetics = "fill",
    limits = list(r_yhat, r_sd),
    palette = pal_vsup(
      values = pal,
      unc_levels = 4,
      max_desat = 0.6,
      pow_desat = 0.2,
      max_light = 0.6,
      pow_light = 1
    ),
    guide = "colorfan"
  ) +
  theme_bw(base_size = 14) +
  labs(x = "Environment", y = "Genotype") #+
#facet_wrap(~t,scales = "free")

pp

# dim: 8.55 and 6.04

# -------------------------------------------------------------------------


name <- "y_hat"
intLims <- range(yhat)
limitsInt <- range(labeling::rpretty(intLims[1], intLims[2]))

library(ggplot2)

p_main <- ggplot(dfNew, aes(e, g)) +
  geom_tile(aes(fill = yhat))+
  facet_wrap(~t)
p_main

#scale_y_discrete(breaks = breakPattern_g, labels = levels(reorder(g, fill))) +
#scale_x_discrete(breaks =  breakPattern_e, labels = levels(reorder(e, -fill)))

intPal = rev(colorspace::diverging_hcl(palette = "Blue-Red 3", n = 100))
lims <- intLims





p_main <- p_main +
  scale_fill_gradientn(
    limits = lims,
    colors = intPal,
    guide = guide_colorbar(
      frame.colour = "black",
      ticks.colour = "black"
    ),
    name = name
  ) +
  theme_classic()

p_main <- p_main + labs(x = "Environment", y = "Genotype")
p_main # original plot

# -------------------------------------------------------------------------
# plot vsup

# library(scales)
# library(tibble)
# library(dplyr)
# library(purrr)
# library(gtable)
# library(scales)
# library(purrr)
# library(grid)

# you can change the max_desat, pow_desat, max_light, pow_light settings
# to get something you like. but i think the default settings are good

pal = rev(colorspace::diverging_hcl(palette = "Blue-Red 3", n = 8))

pp <- ggplot(df.ambarti) +
  geom_raster(aes(e, g, fill = zip(yhat, sd))) +
  bivariate_scale(
    name = c("Value", "sd"),
    aesthetics = "fill",
    limits = list(c(4, 12), c(0.5, 0.55)),
    palette = pal_vsup(
      values = pal,
      unc_levels = 4,
      max_desat = 0.6,
      pow_desat = 0.2,
      max_light = 0.6,
      pow_light = 1
    ),
    guide = "colorfan"
  ) +
  theme_bw(base_size = 14) +
  labs(x = "Environment", y = "Genotype") #+
  #facet_wrap(~t,scales = "free")

pp

df.ambarti$sd |> range()
df.ambarti$yhat |> range()

pp


load("~/Documents/GitHub/bammit/Real data/test_ireland.RData")
test

RMSE(test$Yield, pred)




# plot year ---------------------------------------------------------------

# bammitVsupOut <- readRDS("~/Documents/GitHub/bammit/Running models/bammitVsupOut.rds")
# df.ambarti <- bammitVsupOut

df.ambarti <- data.frame(yhat = yhat, e = e, g = g, t = t, b = b, sd = yhatSD)

mean_real <- data |> dplyr::group_by(Year) |> dplyr::summarise(mean = mean(Yield),
                                                               sd = sd(Yield))

df.ambarti |> ggplot(aes(x = factor(t), y = yhat)) +
  geom_boxplot(size = 0.2) +
  geom_point(data = mean_real, aes(x = factor(Year), y = mean),
             colour = "steelblue", size = 2.5) +
  theme_bw(base_size = 16) +
  labs(x = 'Year', y = 'Yield')




# df_years <- list()
# df_years <- lapply(unique(df.ambarti$t), function(x){
#   df <- subset(df.ambarti,  t == x)
#   #df <- subset(df,  grepl('3$', b))
# })
#
# dat <- dplyr::bind_rows(df_years)
#
# means_year <- lapply(df_years, function(x){
#   data.frame(m = median(x[, 'yhat']),
#              sd = sd(x[, 'yhat']))
# }) |>
#   dplyr::bind_rows() |>
#   dplyr::mutate(t = unique(dat$t),
#                 lower = m - 2*sd,
#                 upper = m + 2*sd,
#                 true = mean_real$mean)
# #pred = true - m)
#
#
# df_long <- means_year |> tidyr::pivot_longer(cols = c('true', 'm'),
#                                              names_to = 'key',
#                                              values_to = 'value')
#
# df_long |> ggplot(aes(x = t, y = value, group = key, colour = key)) +
#   geom_point() +  # Draw lines
#   geom_pointrange(aes(ymin = lower, ymax = upper)) +
#   labs(
#     x = "Year",
#     y = "Yield"
#   )+
#   theme_bw(base_size = 16)

# p1 <- ggplot(data = means_year, aes(t, pred))+
#     #geom_point() +
#     labs(x = 'Year', y = '\U0177') +
#     geom_pointrange(aes(ymin = lower, ymax = upper), colour = "steelblue") +
#     theme_bw(base_size = 16)
# p1







