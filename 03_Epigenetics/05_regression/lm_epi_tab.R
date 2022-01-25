#try to apply lm() to epigenetic summary table
#table compiled by A.Svendsen
#by LVB, Oct 2021

setwd("~/Documents/students_colleagues/2021/Arthur")
mytab<-read.csv("table_epigenetics.csv", row.names = 1)
mytab[mytab=="nd"]<-"absent"
mylist<-colnames(mytab)
#lm() function consumes table easily, no adjustments required
fit=lm(Log2FC~0+ Freq_group + accessibility+ H3K4m3+ H3K27m3+
         H3K36m3+ DNA,mytab)
anova(fit)
summary(fit)
#quite amazingly data explain at least 50% of the expressions.

#this part performs poorly, not enough variation, too much redundancy
#check predicting quality as ROC
library(pROC)
g<-roc(response=mytab$Log2FC, predictor=predict(fit), auc=T)
plot(g, 
     xlim=c(1,0), #does not work well
     col="red", 
     legacy.axes=TRUE, #choose type of x coordinate
     print.thres=TRUE) #show optimal threshold
#suggested bent at 0.336
#clumsy plot, too many identical values

# by hand
plot(mytab$Log2FC,predict(fit))
abline(h=0.336, col="blue")
#does not make sense, ignore this part