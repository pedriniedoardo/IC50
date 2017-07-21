library(readr)
df <- read_csv("Copy of Toshiki4IC50.csv")
# structure of the input dataset
str(df)
# fix the data for data analysis
library(dplyr)
library(tidyr)
df_long <- df %>% mutate(sample=letters[1:nrow(df)])%>%
  # gather the column by treatment
  gather(key = conc,value = colonies,-sample)
#structure of the output dataset
str(df_long)
# make the boxplot (treatment as categorical variable)
library(ggplot2)
ggplot(df_long,aes(y=colonies,x=factor(conc,levels = c("Cont","10uM","50uM","100uM","1mM"))))+
  geom_boxplot()+
  xlab("")
#########################
#########################
## IMPLEMENT IT BETTER ##
#########################
#########################

# make the conc as number using a lookup table
label <- unique(df_long$conc)
d_conc <- c(0,10,50,100,1000)
# change the character with the numeric variable

# copy the dataset
df_long_lm <- df_long
# loop to update the velues according to the lookup table
for(i in 1:nrow(df_long_lm)){
  id <- match(df_long_lm$conc[i],label)
  df_long_lm$conc[i] <- d_conc[id]
}
# make the value as numeric
df_long_lm$conc <- as.numeric(df_long_lm$conc)

# build the model base on the data

model <- drm(colonies ~ conc, 
             data = df_long_lm, 
             fct = LL.3())
summary(model)

# use the R packages lmtest and sandwich to obtain robust standard errors
# to address the fact that some variance heterogeneity is present.
library(sandwich)
library(lmtest)
coeftest(model, vcov = sandwich)

# Simultaneous inference is also possible through the use of the function glht() in the R
# package multcomp:

library(multcomp)
summary(glht(model))

# Estimating effective doses ED5, ED10, and ED50 is accomplished using ED():
ed<-ED(model, c(5, 10, 50), interval = "delta")

# proudce the plot
# buidl a table of values to collect the prediction ofthe model
newdata<-{}
# change the canvas based on the starting and final point
newdata <- expand.grid(conc=exp(seq(log(0.1), log(1000), length=100)))
# produce the predictions and estimte the confidence intervals
pm <- predict(model, newdata=newdata, interval="confidence")
# collect the result of the predictions
newdata$p <- pm[,1]
newdata$pmin <- pm[,2]
newdata$pmax <- pm[,3]
# plot curve
library(ggplot2)
# shift conc == 0 a bit up toavoid messing with the log
df_long_lm$conc0 <- df_long_lm$conc
# other paper suggest a more strict shift of the position of the 0, using 1/100 of the lowest drug concentration
# df_long_lm$conc0[df_long_lm$conc0 == 0] <- 0.5
df_long_lm$conc0[df_long_lm$conc0 == 0] <- (1/100*sort(unique(df_long_lm$conc))[2])
# plotting the data
ggplot(df_long_lm, aes(x = conc0, y = colonies)) +
  geom_point() +
  # stat_summary(aes(y = colonies,group=1), fun.y=mean, colour="red", geom="line",group=1)+
  geom_ribbon(data=newdata, aes(x=conc, y=p, ymin=pmin, ymax=pmax), alpha=0.2)+
  geom_line(data=newdata, aes(x=conc, y=p))+
  coord_trans(x="log")+
  xlab("uM drug") + ylab("number of colonies")+
  # add a red line in correspondence of the IC50
  geom_vline(xintercept = ed[3,1],col="red") +
  # fix the scale
  scale_x_continuous(breaks = c(1,3,10,20,50,100,500,1000))
