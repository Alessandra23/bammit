library(stringr)
load("~/Documents/GitHub/bammit/Real data/train_ireland.RData")
modeltrainQ1 <- readRDS("~/Documents/GitHub/bammit/Real data/modeltrainQ1.rds")
model <- modeltrainQ1
data <- train

pred <- predictionBAMMITReal(model = model, data = data)

# set up names

g <- data$Genotype |> as.character()
e <- data$Environment |> as.character()
t <- data$Year |> as.character()
b <- data$Bloc |> as.character()

# get yhats
yhat <- apply(pred, 2, median) # predict values for yhat
# get sd
yhatSD <-  apply(pred, 2, sd)

# create a df
df.ambarti <- data.frame(yhat = yhat, e = e, g = g, t = t, b = b, sd = yhatSD)
df.ambarti <- subset(df.ambarti,  t == c("2019"))
df.ambarti <- subset(df.ambarti,  grepl('1$', b))
df.ambarti <- df.ambarti |> dplyr::select("yhat", "e", "g", "t", "sd")
df.ambarti$e <- factor(df.ambarti$e, levels = unique(df.ambarti$e))
df.ambarti$g <- factor(df.ambarti$g, levels = unique(df.ambarti$g))
df.ambarti$t <- factor(df.ambarti$t, levels = unique(df.ambarti$t))
#df.ambarti$b <- factor(df.ambarti$b, levels = unique(df.ambarti$b))



fill <- yhat
breakPattern_g <- as.numeric(str_extract_all(levels(reorder(g, fill)), "[0-9]+"))
breakPattern_e <- as.numeric(str_extract_all(levels(reorder(e, fill)), "[0-9]+"))
breakPattern_t <- as.numeric(str_extract_all(levels(reorder(t, fill)), "[0-9]+"))
# breakPattern_b <- as.numeric(str_extract_all(levels(b), "[0-9]+"))

breakPattern_e <- rev(breakPattern_e)



# vsup plots --------------------------------------------------------------

pal = rev(colorspace::diverging_hcl(palette = "Blue-Red 3", n = 8))

range(df.ambarti$yhat)
range(df.ambarti$sd)

pp <- ggplot(df.ambarti) +
  geom_raster(aes(e, g, fill = zip(yhat, sd))) +
  bivariate_scale(
    name = c("Value", "sd"),
    aesthetics = "fill",
    limits = list(c(9, 14), c(0.5, 0.55)),
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






