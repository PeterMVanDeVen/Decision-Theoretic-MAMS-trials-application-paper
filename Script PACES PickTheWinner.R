rm(list=ls())

# Load packages

library(foreign)

# Load functions

setwd("C:/wd") #set to own working directory where functions are stored
source("Functions PickTheWinner.R")

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

# MARGIN FOR DECLARING SUPERIORITY OF EXPERIMENTAL ARM OVER CONTROL ARM: DELTA = 0

set.seed(18082021)
prior = c(1,1)
delta = 0

phat = x/n
rbind(x,n, phat)

Tdata(x, n, prior, delta, G = 1000000, print=TRUE)

# Best decision: Only On Track superior to Usual Care
# Posterior probability that decision minimizes loss: 99.9%
#
# Output of Tdata: 0.000581 0.000562 0.998857

# Allow dropping of all arms (including usual care)

dropctrl = TRUE

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

# ONCOMOVE DROPPED AT FIRST INTERIM ANALYSIS
# 45/66 (UC), 8/12 (OM), 61/66 (OT) SUCCESSES

set.seed(20082021)

Tdata(x_drop, n_drop, prior, delta, G = G_Tdata, print=TRUE)

# Best decision: Only On Track superior to Usual Care
# Posterior probability that decision minimizes loss: 97.3%
#
# Output of Tdata: 0.0190316 0.0080356 0.9729328


##############################################################################
# SUMMARY                                                                    #
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
# Posterior probability that decision minimizes loss: 98.5%
#
# [1] 0.0148843 0.0003671 0.9847486


#################################################
# EVALUATE TRIAL USING BAYESIAN-ADAPTIVE METHOD #
# BOTH WITH AND WITHOUT DROPPING                #
# GAMMA = 0.01: EARLY STOPPING OF TRIAL         #
#################################################

nstages001.out = min(which(trial.out$delta_benefit < 0.01))
nstages001_nodrop.out = min(which(trial_nodrop.out$delta_benefit < 0.01))


# setting with dropping and delta = 0.01

n001_drop <- apply(trial.out$n[1:nstages001.out,],2,sum)
x001_drop <- c(sum(Y[(1:n001_drop[1]),1]),sum(Y[(1:n001_drop[2]),2]),sum(Y[(1:n001_drop[3]),3]))
phat001_drop <- x001_drop/n001_drop

rbind(n001_drop, x001_drop, phat001_drop)

# ONCOMOVE DROPPED AT FIRST INTERIM ANALYSIS
# 20/30 (UC), 8/12 (OM), 28/30 (OT) SUCCESSES
# 72 patients 

set.seed(20220501)

Tdata(x001_drop, n001_drop, prior, delta, G = G_Tdata, print=TRUE)


# Best decision: Only On Track superior to Usual Care
# Posterior probability that decision minimizes loss: 97.9%
# Output Tdata: 0.046352 0.0166673 0.9786975


# setting without dropping and delta = 0.01

n001_nodrop <- apply(trial_nodrop.out$n[1:nstages001_nodrop.out,],2,sum)
x001_nodrop <- c(sum(Y[(1:n001_nodrop[1]),1]),sum(Y[(1:n001_nodrop[2]),2]),sum(Y[(1:n001_nodrop[3]),3]))
phat001_nodrop <- x001_nodrop/n001_nodrop

rbind(n001_nodrop, x001_nodrop, phat001_nodrop)

# ONCOMOVE DROPPED AT FIRST INTERIM ANALYSIS
# 16/24 (UC), 16/24 (OM), 22/24 (OT) SUCCESSES
# 72 patients 

set.seed(20220601)

Tdata(x001_nodrop, n001_nodrop, prior, delta, G = G_Tdata, print=TRUE)


# Best decision: Only On Track superior to Usual Care
# Posterior probability that decision minimizes loss: 96.6%
#
# Output Tdata: 0.0168189 0.0168255 0.9663556



#################################################
# EVALUATE TRIAL USING BAYESIAN-ADAPTIVE METHOD #
# BOTH WITH AND WITHOUT DROPPING                #
# GAMMA = 0.001: EARLY STOPPING OF TRIAL        #
#################################################

nstages0001.out = min(which(trial.out$delta_benefit < 0.001))
nstages0001_nodrop.out = min(which(trial_nodrop.out$delta_benefit < 0.001))

# Setting with dropping and delta = 0.001

n0001_drop <- apply(trial.out$n[1:nstages0001.out,],2,sum)
x0001_drop <- c(sum(Y[(1:n0001_drop[1]),1]),sum(Y[(1:n0001_drop[2]),2]),sum(Y[(1:n0001_drop[3]),3]))
phat0001_drop <- x0001_drop/n0001_drop

rbind(n0001_drop, x0001_drop, phat0001_drop)

# ONCOMOVE DROPPED AT FIRST INTERIM ANALYSIS
# 20/30 (UC), 8/12 (OM), 28/30 (OT) SUCCESSES

set.seed(20220321)

Tdata(x0001_drop, n0001_drop, prior, delta, G = G_Tdata, print=TRUE)

# Best decision: On Track best
# Posterior probability that decision minimizes loss: 97.9%
#
# Output Tdata: 0.0046303 0.0167122 0.9786575


# Setting without dropping and delta = 0.001

n0001_nodrop <- apply(trial_nodrop.out$n[1:nstages0001_nodrop.out,],2,sum)
x0001_nodrop <- c(sum(Y[(1:n0001_nodrop[1]),1]),sum(Y[(1:n0001_nodrop[2]),2]),sum(Y[(1:n0001_nodrop[3]),3]))
phat0001_nodrop <- x0001_nodrop/n0001_nodrop

rbind(n0001_nodrop, x0001_nodrop, phat0001_nodrop)

# ONCOMOVE DROPPED AT FIRST INTERIM ANALYSIS
# 16/24 (UC), 16/24 (OM), 22/24 (OT) SUCCESSES

set.seed(20220601)

Tdata(x0001_nodrop, n0001_nodrop, prior, delta, G = G_Tdata, print=TRUE)

# Best decision: On Track best
# Posterior probability that decision minimizes loss: 96.6%
#
# Output Tdata: 0.0168189 0.0168255 0.9663556

save.image("Results PACES PickTheWinner.Rdata")


