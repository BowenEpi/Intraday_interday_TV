###############################################################################
#
#
# Title: A new method to separate the impacts of interday and intraday 
#        temperature variability on mortality
# 
# Author: Bo Wen, Yao Wu, Yuming Guo, Shanshan Li
#
# Affiliation: Monash University
#
# Example code for definition, and model
#
#
###############################################################################
library(dlnm);library(tidyverse);library(splines);library(mixmeta);library(car);
library(tsModel)

###############################################################################  
# Define interday and intraday temperature variability (TV) 0-7
# Note: definition should be conducted within the same city or province

duration <- 7

## INTERDAY AND INTRADAY TV
data$mean <- apply(data[,c("tmin","tmax")],1,mean,na.rm=T)
data$tv_inter <- sqrt(apply(with(data,rbind(Lag(mean,0:duration))),1,var,na.rm=T)*duration*2/((duration+1)*2-1))
data$tv_intra <- sqrt(apply(Lag(apply(data[,c("tmin","tmax")],1,var,na.rm=T),0:duration),1,sum,na.rm=T)/((duration+1)*2-1))

## TOTAL TV
data$tv <- apply(with(data,cbind(Lag(tmin,0:duration),Lag(tmax,0:duration))),1,sd,na.rm=T)

###############################################################################  
## GENERATE STRATA VARS
data$year <- as.numeric(format(data$date,"%Y"))
data$month <- as.numeric(format(data$date,"%m"))
data$ymdow <- format(data$date,"%Y%m%w")
data$pymdow <- as.factor(paste(data$citycode, data$ymdow,sep="."))

iqr_inter <- IQR(dt.fr$tv_inter,na.rm=T)
iqr_intra <- IQR(dt.fr$tv_intra,na.rm=T)
iqr_tv <- IQR(dt.fr$tv,na.rm=T)

## PARAMETER SETTING
temp.lag = 21
df.lag = 4

## CONSTRUCT CROSSBASIS FUNCTIONS
basistemp<-crossbasis(data$tmean,lag=temp.lag, argvar=list(fun="ns",df=df.lag),
                      arglag=list(fun="ns",knots=logknots(temp.lag, fun="ns", df=df.lag)),
                      group=data$cityname)

basisrh<-crossbasis(data$RH,lag=temp.lag, argvar=list(fun="ns",df=df.lag),
                    arglag=list(fun="ns",knots=logknots(temp.lag, fun="ns", df=df.lag)),
                    group=data$cityname)

## FIT MODEL WITH TV
fit1 <- gnm(y~I(tv/iqr_tv)+basistemp+basisrh,data=data,family=quasipoisson,eliminate = pymdow)

## FIT MODEL WITH INTERDAY AND INTRADAY TV
fit2 <- gnm(y~I(tv_inter/iqr_inter)+I(tv_intra/iqr_intra)+basistemp+basisrh,data=data,family=quasipoisson,eliminate = pymdow)

## EXTRACT COEEFICIENTS
dt.result <- data.frame(index=c(paste0("TV0-",duration),paste0("Inter-TV0-",duration),
                                paste0("Intra-TV0-",duration)),
                        coef=NA,se=NA,coef_low=NA,coef_high=NA)

## MODEL WITH TV
dt.result[1,"coef"] <- summary(fit1)$coefficients["I(tv/iqr_tv)","Estimate"]
dt.result[1,"se"] <- summary(fit1)$coefficients["I(tv/iqr_tv)","Std. Error"]

## MODEL WITH INTERDAY AND INTRADAY TV
dt.result[2,"coef"] <- summary(fit2)$coefficients["I(tv_inter/iqr_inter)","Estimate"]
dt.result[2,"se"] <- summary(fit2)$coefficients["I(tv_inter/iqr_inter)","Std. Error"]

dt.result[3,"coef"] <- summary(fit2)$coefficients["I(tv_intra/iqr_intra)","Estimate"]
dt.result[3,"se"] <- summary(fit2)$coefficients["I(tv_intra/iqr_intra)","Std. Error"]

## CALCULATE 95% CI
dt.result$coef_low <- dt.result$coef-1.96*dt.result$se
dt.result$coef_high <- dt.result$coef+1.96*dt.result$se








