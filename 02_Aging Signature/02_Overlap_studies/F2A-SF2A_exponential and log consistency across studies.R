## Distribution of number of genes and their consistency
#Data extracted from gene_in_AS_pyramid ("generate freqs for pyramid")
#A. Svendsen Feb 2020

#Adapted from: 
#https://stackoverflow.com/questions/37329074/geom-smooth-and-exponential-fits
#https://stackoverflow.com/questions/47822911/how-to-add-linear-model-results-adj-r-squared-slope-and-p-value-onto-regressi

setwd("/Users/Art/Drive/PhD/Experiments/Aging Signature/000_FINAL/")
library(tidyverse)
library(patchwork)
source("/Users/Art/Drive/PhD/Experiments/Aging Signature/plotting_theme.R")

mydata <- read.csv("00_sourceFiles/07_Summary_overlaps.csv", header = T, stringsAsFactors = F, na.strings = "NA")

#subset meta and reanalysis
df.meta <-data.frame(x = mydata$overlaps, y = mydata$meta)
df.re <-data.frame(x = mydata$overlaps, y = mydata$reanalysis)

#lm fit model for individual sets
fit.meta <- lm(y ~ exp(-x), data = df.meta)
fit.re <- lm(y ~ exp(-x), data = df.re)
summary(fit.meta)
summary(fit.re)

fit.meta.lm <- lm(log(y) ~ x, data = df.meta[1:9,])
fit.re.lm <- lm(log(y) ~ x, data = df.re[1:11,])
summary(fit.meta.lm)
summary(fit.re.lm)

#The regression equations
# y = a*b + b
# a = Estimate (Intercept)
# b = Estimate (x)
# meta <- -1.10405*x + 9.07949
# reanl <- y = -0.7281*x + 7.9921

#Predicted
df.predicted <- data.frame()
for (i in 1:12){
  y.meta = -1.10388 * i + 9.07689
  df.predicted <- rbind.data.frame(df.predicted, cbind.data.frame(consistency = i, num = exp(y.meta), analysis = "meta", model = "predicted"))
}

for (i in 1:12){
  y.reanl = -0.72134*i + 7.72582
  df.predicted <- rbind.data.frame(df.predicted, cbind.data.frame(consistency = i, num = exp(y.reanl), analysis = "reanl", model = "predicted"))
}
for (i in 1:12){
  if (i <= 9){
    y.meta = filter(df.meta, x == i)[,2]
  } else {
    y.meta = c(NA)
  }
  df.predicted <- rbind.data.frame(df.predicted, cbind.data.frame(consistency = i, num =  y.meta, analysis = "meta", model = "observed"))
}
for (i in 1:12){
  if (i <= 11){
    y.reanl = filter(df.re, x == i)[,2]
  } else {
    y.reanl = c(NA)
  }
  df.predicted <- rbind.data.frame(df.predicted, cbind.data.frame(consistency = i, num = y.reanl, analysis = "reanl", model = "observed"))
}

#Exponential merge curves
p_exp <- ggplot() + 
  geom_point(data = df.meta, aes(x=x, y=y), color = '#984ea3', shape = 1, size = 7) +
  geom_smooth(data = df.meta, aes(x=x, y=y), method="glm", formula= (y ~ exp(-x)), se=FALSE, linetype = 1, color = '#984ea3', alpha = 0.5) +
  geom_point(data = df.re, aes(x=x, y=y), color = '#4daf4a', shape = 1, size = 7) +
  geom_smooth(data = df.re, aes(x=x, y=y), method="glm", formula= (y ~ exp(-x)), se=FALSE, linetype = 1, color = '#4daf4a', alpha = 0.5) +
  scale_x_continuous(breaks = c(1:12)) +
  labs(y = "# of DE genes",
       x = "number of overlaps across studies") +
  theme_pb() +
  theme(aspect.ratio = 1)

#Individual log transformed plots with stats
#Meta
p_meta <- ggplot() +
  geom_point(data = filter(df.predicted, analysis == "meta"), aes(x = consistency, y = log(num), color = model), shape = 1, size = 7) +
  geom_smooth(data = filter(df.predicted, analysis == "meta", model == "predicted"), aes(x= consistency, y= log(num)), method="lm", formula= (y ~ x), se=FALSE, linetype = 1, color = 'grey', alpha = 0.5) +
  geom_smooth(data = filter(df.predicted, analysis == "meta", model == "observed"), aes(x= consistency, y= log(num)), method="lm", formula= (y ~ x), se=FALSE, linetype = 1, color = '#984ea3', alpha = 0.5) +
  scale_color_manual(values = c('grey', '#984ea3')) +
  scale_y_continuous(limits = c(-5,10)) +
  scale_x_continuous(limits = c(0,15)) +
  labs(title = "meta", y = "log10 - Reported DE genes") +
  geom_label(
    aes(x = 1.5, y = 0), hjust = 0,
    color = "#984ea3",
    label = paste("metanalysis",
                  "\nAdj R2 = ",signif(summary(fit.meta)$adj.r.squared, 5),
    size = 5)) +
  theme_pb() +
  theme(aspect.ratio = 1)
p_meta

#ReAn
p_reanl <- ggplot() +
  geom_point(data = filter(df.predicted, analysis == "reanl"), aes(x = consistency, y = log(num), color = model), shape = 1, size = 7) +
  geom_smooth(data = filter(df.predicted, analysis == "reanl", model == "predicted"), aes(x= consistency, y= log(num)), method="lm", formula= (y ~ x), se=FALSE, linetype = 1, color = 'grey', alpha = 0.5) +
  geom_smooth(data = filter(df.predicted, analysis == "reanl", model == "observed"), aes(x= consistency, y= log(num)), method="lm", formula= (y ~ x), se=FALSE, linetype = 1, color = '#4daf4a', alpha = 0.5) +
  scale_y_continuous(limits = c(-5,10)) +
  scale_x_continuous(limits = c(0,15)) +
  scale_color_manual(values = c('grey', '#4daf4a')) +
  labs(title = "reanl", y = "log10 - Reported DE genes") +
  geom_label(
    aes(x = 1.5, y = 0), hjust = 0,
    color = "#4daf4a",
    label = paste("reanl",
                  "\nAdj R2 = ",signif(summary(fit.re)$adj.r.squared, 5),
                  size = 5)) +
  theme_pb() +
  theme(aspect.ratio = 1)
p_reanl

p_exp
