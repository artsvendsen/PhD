##Global Life Expectancy 
##Data from https://ourworldindata.org/life-expectancy (orginaly from World Bank Data)
setwd("/Users/Art/Drive/PhD/Manuscripts/2020_Hoffman 8th edition/data/")
library(tidyverse)

mydata <-read.csv("Life expentancy/life-expectancy.csv", header = T, stringsAsFactors = T)[,-2]
colnames(mydata) <- c("country", "year", "expectancy")



ggplot(data = filter(mydata, year > 1950)) +
  #geom_smooth(geom='line', alpha=0.05, se=FALSE, aes(x = year, y = expectancy, group = country), color = 'grey')
  geom_line(aes(x = year, y = expectancy, group = country),
            stat = "smooth", method = 'loess', formula = y ~ x,
            color = 'grey',
            alpha = 0.25) +
  geom_smooth(aes(x = year, y = expectancy)) +
  labs(y = "Life expectancy (in years)") +
  scale_y_continuous(limits = c(0,100)) +
  theme_bw() +
  theme(aspect.ratio = 1)