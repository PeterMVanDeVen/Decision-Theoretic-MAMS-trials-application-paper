# SCRIPTS TO EVALUATE TYPE I ERROR IN THE "PICK ALL TREATMENTS SUPERIOR TO CONTROL" SETTING
# FOR DELTA = 0.01 AND DELTA = 0.001

rm(list=ls())
set.seed(26062012)

# load packages

library(snow)
library(doSNOW)
library(doParallel)
library(foreign)

# load functions
setwd("C:/wd") #set to own working directory where functions are stored
source("Functions PickAllTreatmentsSuperiorToControl.R")


# type I error when success rates in all arms equal are 0.65 (success rates UC)

resp = 0.65
  
prior = c(1,1)
maxpt = 3*240 
burn=12
batch=12
K=3
Ntrials = 10000
Ydata <- array(dim=c(maxpt, K, Ntrials))
for(i in 1:Ntrials){
   Ydata[,,i]=matrix(rbinom(3*maxpt, 1, resp), ncol= 3)
}


delta = 0.1
dropctrl = FALSE


# gamma of 0.001

gamma = 0.001


trial.drop          <- prdrop(Ydata, delta, gamma, burn=12, Ntrials, batch=12, prior, maxpt=maxpt, dropctrl)
res.drop            <- pcalc(trial.drop, prior, delta, nr_dec=4)
res.drop

trial.nodrop        <- prnodrop(Ydata, delta, gamma, burn=12, Ntrials, batch=12, prior, maxpt=maxpt)
res.nodrop          <- pcalc(trial.nodrop, prior, delta, nr_dec=4)
res.nodrop

# gamma of 0.01

newgamma = 0.01

trial.drop.newgamma   <- reevaluate(trial.drop, Ydata, newgamma)
res.drop.newgamma     <- pcalc(trial.drop.newgamma, prior, delta, nr_dec=4)
res.drop.newgamma 

trial.nodrop.newgamma <- reevaluate(trial.nodrop, Ydata, newgamma)
res.nodrop.newgamma   <- pcalc(trial.nodrop.newgamma, prior, delta, nr_dec=4)
res.nodrop.newgamma 

# type I error and mean N for delta = 0.001
# with dropping

res.drop$dec
mean(res.drop$N)

# type I error: 16.3%
# mean N: 174

# without dropping

res.nodrop$dec
mean(res.nodrop$N)

# type I error: 5.3%
# mean N: 195


# with delta = 0.01
# with dropping

res.drop.newgamma$dec
mean(res.drop.newgamma$N)

# type I error 18.4 %
# mean N: 125.8

# without dropping

res.nodrop.newgamma$dec
mean(res.nodrop.newgamma$N)

# type I error 9.7%
# mean N: 138.3




save.image("Results PACES Type I error.Rdata")
