#Mortality and case rates by leukemmia types and age
#Data was extracted from Cancer Research UK
#(https://www.cancerresearchuk.org/health-professional/cancer-statistics/statistics-by-cancer-type/leukaemia/incidence#heading-One)
#A. Svenden March 2020
setwd("/Users/Art/Drive/PhD/Manuscripts/2020_Hoffman 8th edition/data/")
library(ggplot2)
library(tidyverse)
library(patchwork)


theme_pb <- function(){
  theme_bw() +
    theme(
      panel.grid = element_blank(),
      panel.border = element_rect(size = 0.8, color = "black"),
      #axis.line = element_line(size = 0.5, colour = "black"),
      axis.ticks.length = unit(0.25, "cm"),
      axis.ticks = element_line(size = 1, color = "black"),
      # axis.title=element_text(size = 18, color = 'black', vjust = 1),
      # axis.text.x = element_text(size = 15, color = 'black', vjust = 0.5),
      # axis.text.y = element_text(size = 15, color = 'black', vjust = 0.5),
      plot.title = element_text(size = 18, color = 'black', hjust = 0.5)
    )
}

#Mortality
mydata <- read.csv("leukemia incidence/Leukemias_mortality_by_Age_UK_combined.csv", header = T, stringsAsFactors = T)

mortality <- ggplot(mydata, aes(x = Age, y = avarege, group = kind)) + 
  #geom_line(aes(color = kind)) +
  #geom_point(alpha = 0.5) +
  geom_smooth(se = F, aes(color = kind)) +
  scale_color_manual(values = c('#173F5F', '#20639B', '#3CAEA3', '#F6D55C', '#ED553B')) +
  labs(x = "", y = "Mortality rate per 100,000") +
  scale_y_continuous(limits = c(0,100))+
  theme_bw() +
  theme(aspect.ratio = 1,
        panel.grid = element_blank(),
        panel.border = element_rect(size = 0.8, color = "black"),
        axis.ticks = element_line(size = 1, color = "black"),
        axis.title=element_text(color = 'black', vjust = 1),
        axis.text.x = element_text(color = 'black', hjust = 1, angle = 90),
        axis.text.y = element_text(color = 'black', vjust = 0.5)
        #axis.text.x = element_text(angle = 90, hjust = 1)
        ) +
  scale_fill_discrete(name="Condition",
                      breaks = c("ALL", "AML", "CLL", "CML", "leukemia"),
                      labels = c("ALL", "AML", "CLL", "CML", "total leukemia"))

mortality

mydata2 <- read.csv("leukemia cases/Leukemias_cases_by_Age_UK_combined.csv", header = T, stringsAsFactors = T)

cases <- ggplot(mydata2, aes(x = Age, y = avarege, group = kind)) + 
  #geom_line(aes(color = kind)) +
  #geom_point(alpha = 0.5) +
  geom_smooth(se = F, aes(color = kind)) +
  scale_color_manual(values = c('#173F5F', '#20639B', '#3CAEA3', '#F6D55C', '#ED553B')) +
  labs(x = "", y = "Incidence rate per 100,000") +
  scale_y_continuous(limits = c(0,100))+
  theme_bw() +
  theme(aspect.ratio = 1,
        panel.grid = element_blank(),
        panel.border = element_rect(size = 0.8, color = "black"),
        axis.ticks = element_line(size = 1, color = "black"),
        axis.title=element_text(color = 'black', vjust = 1),
        axis.text.x = element_text(color = 'black', hjust = 1, angle = 90),
        axis.text.y = element_text(color = 'black', vjust = 0.5)
        #axis.text.x = element_text(angle = 90, hjust = 1)
  ) +
  scale_fill_discrete(name="Condition",
                      breaks = c("ALL", "AML", "CLL", "CML", "leukemia"),
                      labels = c("ALL", "AML", "CLL", "CML", "total leukemia"))

cases

cases | mortality

