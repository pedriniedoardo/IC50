
# load the package --------------------------------------------------------
library(drc)
library(tidyverse)
library(sandwich)
library(lmtest)
library(multcomp)

# load the data -----------------------------------------------------------
ic_test <- read_csv("IC50.csv")

str(ic_test)

# data manipulation -------------------------------------------------------
# fix the data for data analysis
df_long <- ic_test%>%
  gather(key = conc_label,value = colony)%>%
  mutate(conc_label = factor(conc_label,levels = c("Cont","10uM","50uM","100uM","1mM")))

str(df_long)

# translate the labels in quantities unsing a lut
label <- unique(df_long$conc_label)
d_conc <- c(0,10,50,100,1000)
lut <- data.frame(conc_label = label,conc_n = d_conc)

# copy the dataset
df_long <- left_join(df_long,lut,by = "conc_label")


# build a model -----------------------------------------------------------
model <- drm(colony ~ conc_n,data = df_long,fct = LL.3())
summary(model)

# use the R packages lmtest and sandwich to obtain robust standard errors
# to address the fact that some variance heterogeneity is present.
coeftest(model, vcov = sandwich)

# Simultaneous inference is also possible through the use of the function glht() in the R
# package multcomp:
summary(glht(model))

# Estimating effective doses ED5, ED10, and ED50 is accomplished using ED():
ed<-ED(model, c(5, 10, 50), interval = "delta")

# use the model -----------------------------------------------------------
# buidl a table of values to collect the prediction ofthe model
# change the canvas based on the starting and final point
newdata <- expand.grid(conc_n=exp(seq(log(0.1), log(1000), length=100)))
# produce the predictions and estimte the confidence intervals
pm <- predict(model, newdata=newdata, interval="confidence")
# collect the result of the predictions
newdata$p <- pm[,1]
newdata$pmin <- pm[,2]
newdata$pmax <- pm[,3]

# shift conc == 0 a bit up toavoid messing with the log
# other paper suggest a more strict shift of the position of the 0, using 1/100 of the lowest drug concentration
# df_long_lm$conc0[df_long_lm$conc0 == 0] <- 0.5

# lowest non zero conc
low_conc <- sort(unique(df_long$conc_n))[2]
zero_conc <- low_conc/100

df_long <- df_long%>%
  mutate(conc0 = ifelse(conc_n == 0 ,yes = zero_conc,no = conc_n))


# plot data ---------------------------------------------------------------
df_long%>%
  ggplot(aes(y=colony,x=conc_label))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(width = 0.1,alpha = 0.5)+
  xlab("")+
  ggsave("test_IC50_boxplot.pdf",height = 3,width = 4)

# plotting the data
df_long%>%
  ggplot(aes(x = conc0, y = colony)) +
  geom_jitter(width = 0.1) +
  # stat_summary(aes(y = colonies,group=1), fun.y=mean, colour="red", geom="line",group=1)+
  geom_ribbon(data=newdata, aes(x=conc_n, y=p, ymin=pmin, ymax=pmax), alpha=0.2)+
  geom_line(data=newdata, aes(x=conc_n, y=p))+
  scale_x_log10()+
  xlab("uM drug") + ylab("number of colonies")+
  # add a red line in correspondence of the IC50
  geom_vline(xintercept = ed[3,1],col="red")+
  ggsave("test_IC50_plot.pdf",height = 3,width = 3)
