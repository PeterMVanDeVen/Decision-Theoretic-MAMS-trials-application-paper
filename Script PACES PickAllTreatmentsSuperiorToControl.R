rm(list=ls())

# load packages

library(foreign)

# load functions

setwd("C:/wd") #set to own working directory where functions are stored
source("Functions PickAllTreatmentsSuperiorToControl.R")

# load PACES dataset

data = read.spss("PACES minimal.sav", use.value.labels = FALSE)

# load PACES dataset

data = read.spss("PACES minimal.sav", use.value.labels = FALSE)

# Outcome 
# success = 1: no dose adjustment
# success = 0: dose adjustment
#
# Treatment 
# arm 0: Usual care
# arm 1: OncoMove
# arm 2: OnTrack

arm = data$Group
success = 1-data$doseadj

# Cross Table of Results in all 230 subjects

all230 = table(arm, success)
all230

# Obtain number of subjects per arm (n) and number of successes per arm (x) 
n = all230[,2]+all230[,1]
x = all230[,2]
maxpatients = min(n)

# Make matrix with outcomes for each arm in different columns (keeping ordering of inclusion in arms)

Y=matrix(c(
  success[which(arm==0)][1:maxpatients],
  success[which(arm==1)][1:maxpatients],
  success[which(arm==2)][1:maxpatients]),ncol=3, byrow=FALSE)


#########################################
# TRIALS WITH FIXED SAMPLE SIZE OF 230  #
# BAYESIAN ANALYSIS AT THE END          #
#########################################

# SET G: HIGHER G MORE PRECISE ESTIMATES

G = 1000000
G_Tdata = 10000000

# MARGIN FOR DECLARING SUPERIORITY OF EXPERIMENTAL ARM OVER CONTROL ARM: DELTA = 0.1

set.seed(18082021)
prior = c(1,1)
delta = 0.1

Tdata(x, n, prior, delta, G = G_Tdata, print=TRUE)
phat = x/n
rbind(x,n, phat)

# Best decision: Only On Track superior to Usual Care
# Posterior probability that decision minimizes loss: 86.8%
#
# Output Tdata: 0.040834 0.000070 0.867911 0.091185
#
# 0.868: Only On Track preferred over Usual Care
# 0.0912: Both On Track and OncoMove preferred over Usual Care
# Marginally: 0.959 On Track preferred over Usual Care
# Marginally: 0.0912 OncoMove preferred over Usual Care

# KEEP CONTROL ARM (ANY TIME)

dropctrl = FALSE


# Trials that continue until maximum number of 76 in an arm is reached 

set.seed(19082021)

K=3
Ydata <- array(dim=c(maxpatients, K, 1))
Ydata[,,1] <- Y

#################################################
# EVALUATE TRIAL USING BAYESIAN-ADAPTIVE METHOD #
# BOTH WITH DROPPING (trial.out)                #
# AND WITHOUT DROPPING (trial_nodrop.out)       #
# GAMMA PARAMETER SET TO 0 -> No EARLY STOPPING #
#################################################

set.seed(05012022)

trial.out        <- prdrop(Ydata, delta, gamma=0, burn=12, Ntrials=1, batch=12, prior, maxpt=1*12+2*76, dropctrl, G)
trial_nodrop.out <- prnodrop(Ydata, delta, gamma=0, burn=12, Ntrials=1, batch=12, prior, maxpt=3*76, G)

trial.out
trial_nodrop.out

##################################################################################
# SUMMARY                                                                        #
# RESULTS FOR TRIAL ALLOWING EARLY DROPPING OF EXPERIMENAL ARMS (MAX 76 PER ARM) #
##################################################################################

n_drop = c(trial.out[[1]],trial.out[[2]],trial.out[[3]])
x_drop = c(trial.out[[4]],trial.out[[5]],trial.out[[6]])
phat_drop = x_drop/n_drop

rbind(n_drop, x_drop, phat_drop)
sum(n_drop)

# ONCOMOVE DROPPED AT FIRST INTERIM ANALYSIS
# 45/66 (UC), 8/12 (OM), 61/66 (OT) SUCCESSES

set.seed(20082021)

Tdata(x_drop, n_drop, prior, delta, G = G_Tdata, print=TRUE)

# Best decision: Only On Track superior to Usual Care
# Posterior probability that decision minimizes loss: 81.3%
#
# Output Tdata: 0.019248 0.000487 0.813285 0.166980
#
# 0.813: Only On Track preferred over Usual Care
# 0.167: Both On Track and OncoMove preferred over Usual Care
# Marginally: 0.980 On Track preferred over Usual Care
# Marginally: 0.167 OncoMove preferred over Usual Care


##############################################################################
# RESULTS TRIAL WITHOUT EARLY DROPPING OF EXPERIMENAL ARMS (MAX 76 PER ARM)  #
##############################################################################

n_nodrop = c(trial_nodrop.out[[1]],trial_nodrop.out[[2]],trial_nodrop.out[[3]])
x_nodrop = c(trial_nodrop.out[[4]],trial_nodrop.out[[5]],trial_nodrop.out[[6]])
phat_nodrop = x_nodrop/n_nodrop

rbind(n_nodrop, x_nodrop, phat_nodrop)

# 47/72 (UC), 49/72 (OM), 65/72 (OT) SUCCESSES

set.seed(21082021)

Tdata(x_nodrop, n_nodrop, prior, delta, G = G_Tdata, print=TRUE)

# Best decision: Only On Track superior to Usual Care
# Posterior probability that decision minimizes loss: 81.2%
#
# Output Tdata: 0.014818 0.000049 0.811759 0.173374
#
# 0.812: Only On Track preferred over Usual Care
# 0.173: Both On Track and OncoMove preferred over Usual Care
# Marginally: 0.985 On Track preferred over Usual Care
# Marginally: 0.173 OncoMove preferred over Usual Care




#################################################
# EVALUATE TRIAL USING BAYESIAN-ADAPTIVE METHOD #
# BOTH WITH AND WITHOUT DROPPING                #
# GAMMA = 0.01 FOR EARLY STOPPING OF TRIAL      #
#################################################

nstages001.out = min(which(trial.out$delta_benefit < 0.01))
nstages001_nodrop.out = min(which(trial_nodrop.out$delta_benefit < 0.01))

# WITH DROPPING

n001_drop <- apply(trial.out$n[1:nstages001.out,],2,sum)
x001_drop <- c(sum(Y[(1:n001_drop[1]),1]),sum(Y[(1:n001_drop[2]),2]),sum(Y[(1:n001_drop[3]),3]))
phat001_drop <- x001_drop/n001_drop

rbind(n001_drop, x001_drop, phat001_drop)

# ONCOMOVE DROPPED AT FIRST INTERIM ANALYSIS
# 20/30 (UC), 8/12 (OM), 28/30 (OT) SUCCESSES
# 72 patients included

set.seed(20220501)

Tdata(x001_drop, n001_drop, prior, delta, G = G_Tdata, print=TRUE)


# Best decision: Only On Track superior to Usual Care
# Posterior probability that decision minimizes loss: 71.3%
#
# Output Tdata: 0.057857 0.002143 0.712979 0.227021
#
# 0.713: Only On Track preferred over Usual Care
# 0.227: Both On Track and OncoMove preferred over Usual Care
# Marginally: 0.940 On Track preferred over Usual Care
# Marginally: 0.229 OncoMove preferred over Usual Care

# WITHOUT DROPPING

n001_nodrop <- apply(trial_nodrop.out$n[1:nstages001_nodrop.out,],2,sum)
x001_nodrop <- c(sum(Y[(1:n001_nodrop[1]),1]),sum(Y[(1:n001_nodrop[2]),2]),sum(Y[(1:n001_nodrop[3]),3]))
phat001_nodrop <- x001_nodrop/n001_nodrop

rbind(n001_nodrop, x001_nodrop, phat001_nodrop)

# 25/36 (UC), 24/36 (OM), 33/36 (OT) SUCCESSES
# 108 patients

set.seed(20220601)

Tdata(x001_nodrop, n001_nodrop, prior, delta, G = G_Tdata, print=TRUE)

# Best decision: Only On Track superior to Usual Care
# Posterior probability that decision minimizes loss: 77.6%
#
# Output Tdata: 0.1054479 0.0006099 0.7763859 0.1175563
#
# 0.776: Only On Track preferred over Usual Care
# 0.118: Both On Track and OncoMove preferred over Usual Care
# Marginally: 0.894 On Track preferred over Usual Care
# Marginally: 0.118 OncoMove preferred over Usual Care




#################################################
# EVALUATE TRIAL USING BAYESIAN-ADAPTIVE METHOD #
# BOTH WITH AND WITHOUT DROPPING                #
# GAMMA = 0.001 FOR EARLY STOPPING OF TRIAL     #
#################################################

nstages0001.out = min(which(trial.out$delta_benefit < 0.001))
nstages0001_nodrop.out = min(which(trial_nodrop.out$delta_benefit < 0.001))

# WITH DROPPING

n0001_drop <- apply(trial.out$n[1:nstages0001.out,],2,sum)
x0001_drop <- c(sum(Y[(1:n0001_drop[1]),1]),sum(Y[(1:n0001_drop[2]),2]),sum(Y[(1:n0001_drop[3]),3]))
phat0001_drop <- x0001_drop/n0001_drop

rbind(n0001_drop, x0001_drop, phat0001_drop)

# ONCOMOVE DROPPED AT FIRST INTERIM ANALYSIS
# 45/66 (UC), 8/12 (OM), 61/66 (OT) SUCCESSES
# 144 patienten
  
set.seed(20220321)

Tdata(x0001_drop, n0001_drop, prior, delta, G = G_Tdata, print=TRUE)


# Best decision: Only On Track superior to Usual Care
# Posterior probability that decision minimizes loss: 81.4%
#
# Output Tdata: 0.0191110 0.0005049 0.8135657 0.1668184
#
# 0.814: Only On Track preferred over Usual Care
# 0.167: Both On Track and OncoMove preferred over Usual Care
# Marginally: 0.990 On Track preferred over Usual Care
# Marginally: 0.167 OncoMove preferred over Usual Care


# WITHOUT DROPPING 

n0001_nodrop <- apply(trial_nodrop.out$n[1:nstages0001_nodrop.out,],2,sum)
x0001_nodrop <- c(sum(Y[(1:n0001_nodrop[1]),1]),sum(Y[(1:n0001_nodrop[2]),2]),sum(Y[(1:n0001_nodrop[3]),3]))
phat0001_nodrop <- x0001_nodrop/n0001_nodrop

rbind(n0001_nodrop, x0001_nodrop, phat0001_nodrop)

# 34/48 (UC), 33/48 (OM), 44/48 (OT) SUCCESSES
# 144 patienten

set.seed(20220601)

Tdata(x0001_nodrop, n0001_nodrop, prior, delta, G = G_Tdata, print=TRUE)


# Best decision: Only On Track superior to Usual Care
# Posterior probability that decision minimizes loss: 77.6%
#
# Output Tdata: 0.0942760 0.0003086 0.8104322 0.0949832
#
# 0.776: Only On Track preferred over Usual Care
# 0.118: Both On Track and OncoMove preferred over Usual Care
# Marginally: 0.894 On Track preferred over Usual Care
# Marginally: 0.118 OncoMove preferred over Usual Care

save.image("Results PickAllTreatmentsSuperiorToControl.Rdata")



