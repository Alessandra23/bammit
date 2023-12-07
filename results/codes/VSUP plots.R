library(stringr)
load("~/Documents/GitHub/bammit/Real data/Ireland_VCU_2015_AMBARTI.RData")

data <- ambarti


# set up names
features <- attributes(ambarti$x)
g_names <- colnames(features$contrasts$g)
e_names <- colnames(features$contrasts$e)
e <- rep(e_names, times = length(g_names))
g <- rep(g_names, each = length(e_names))

# get yhats
yhat <- apply(ambarti$y_hat, 2, median) # predict values for yhat
yhatBart <- apply(ambarti$y_hat_bart, 2, median) # predict values for bart

# get sd
yhatSD <-  apply(ambarti$y_hat, 2, sd)

# create a df
df.ambarti <- data.frame(yhatBart = yhatBart, yhat = yhat, e = e, g = g, sd = yhatSD)
df.ambarti$e <- factor(df.ambarti$e, levels = unique(df.ambarti$e))
df.ambarti$g <- factor(df.ambarti$g, levels = unique(df.ambarti$g))

fill <- yhat
breakPattern_g <- as.numeric(str_extract_all(levels(reorder(g, fill)), "[0-9]+"))
breakPattern_e <- as.numeric(str_extract_all(levels(reorder(e, fill)), "[0-9]+"))

breakPattern_e <- rev(breakPattern_e)


# -------------------------------------------------------------------------

name <- "y_hat"
intLims <- range(yhat)
limitsInt <- range(labeling::rpretty(intLims[1], intLims[2]))

library(ggplot2)

p_main <- ggplot(df.ambarti, aes(reorder(e, -fill), reorder(g, fill))) +
  geom_tile(aes(fill = fill)) +
scale_y_discrete(breaks = breakPattern_g, labels = levels(reorder(g, fill))) +
scale_x_discrete(breaks =  breakPattern_e, labels = levels(reorder(e, -fill)))

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
  geom_raster(aes(reorder(e, -fill), reorder(g, fill), fill = zip(yhat, sd))) +
  bivariate_scale(
    name = c("Value", "sd"),
    aesthetics = "fill",
    limits = list(c(10, 15), c(0.2, 0.31)),
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
  labs(x = "Environment", y = "Genotype")


pp






