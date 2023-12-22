# THE OCG #
##################################################################
# Chlorophyll-a prediction based on level 3 MODIS Aqua data      #
##################################################################

#a data frame with the bands: Rrs_412,Rrs_443,Rrs_488,Rrs_547,Rrs_555,Rrs_667
myData<-data.frame(Rrs_412=0.009470001,
                   Rrs_443=0.007502001,
                   Rrs_469=0.006674001,
                   Rrs_488=0.005760001,
                   Rrs_531=0.004334001,
                   Rrs_547=0.003742001,
                   Rrs_555=0.003296001,
                   Rrs_645=0.000664001,
                   Rrs_667=0.000388001,
                   Rrs_678=0.000500001)

# or use csv: myData<-read.csv("satellite-matchups.csv")
########################################
# how to calculate a selected quantile #
########################################

library(gamlss)
library(gamlss.dist)

# mu has log link
mu_prediction<-exp(-19.157+28.326*sqrt(myData$Rrs_412)-240.12*myData$Rrs_443-2.360*log(myData$Rrs_488)-333.163*myData$Rrs_547+114.507*sqrt(myData$Rrs_555)+6.768*sqrt(myData$Rrs_667))

# sigma has log link
sigma_prediction<-exp(0.7915-32.5579*myData$Rrs_443+0.2316*log(myData$Rrs_555))

# nu has identity link
nu_prediction<- 0.1957-62.8881*myData$Rrs_412

# tau has log link
tau_prediction<-rep(exp(1.626),nrow(myData))

# quantile I want to predict
my_quantile<-0.5

quantilePrediction<-qBCTo(my_quantile,mu=mu_prediction,sigma=sigma_prediction,nu=nu_prediction,tau=tau_prediction)

############################################
# how to calculate the qCV  #
############################################

IQRpred<-qBCTo(0.75,mu=mu_prediction,sigma=sigma_prediction,nu=nu_prediction,tau=tau_prediction)-qBCTo(0.25,mu=mu_prediction,sigma=sigma_prediction,nu=nu_prediction,tau=tau_prediction)
medpred<-qBCTo(0.5,mu=mu_prediction,sigma=sigma_prediction,nu=nu_prediction,tau=tau_prediction)
qCV<-IQRpred/medpred

################################################
# how to calculate p > threshold e.g., 5 ug/L  #
################################################

threshold<-5 

pPrediction<-1-pBCTo(threshold,mu=mu_prediction,sigma=sigma_prediction,nu=nu_prediction,tau=tau_prediction)
