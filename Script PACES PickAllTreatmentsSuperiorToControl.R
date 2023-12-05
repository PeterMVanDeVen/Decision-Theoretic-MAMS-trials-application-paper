# clear workspace

rm(list=ls())

# load packages

library(foreign)
library(snow)
library(doSNOW)
library(doParallel)

# Load required functions (R script "Functions PickAllTreatmentsSuperiorToControl.R")

# setwd("C:/wd") #set to own working directory where functions are stored

setwd("C:/Users/pven3/OneDrive - UMC Utrecht/Documenten/Onderzoek/Decision-theoretic design/Paper Laurien/Github/Improved") #set to own working directory where functions are stored

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

# OBTAIN ESTIMATES OF SUCCES RATES

phat = x/n
rbind(x,n, phat)

# Function Tdata calculates the posterior probability of each final decision being the correct decision
# entry 1: Posterior probability of neither arm 2 (OncoMove) nor arm 3 (OnTarget) being superior to arm 1 (Usual Care) 
# entry 2: Posterior probability of only arm 2 (OncoMove) being superior to arm 1 (Usual Care)
# entry 3: Posterior probability of only arm 3 (OnTarget) being superior to arm 1 (Usual Care)
# entry 4: Posterior probability of both arm 2 being (OncoMove) and arm 3 (OnTarget) being superior to arm 1 (Usual Care)

Tdata(x, n, prior, delta, G = G_Tdata, print=TRUE)

######################################################################
# Results                                                            #
#                                                                    #
# Bayesian decision-theoretic analysis of original trial             #
#                                                                    #
# Success rate arm 1: 66.2%                                          #
# Success rate arm 2: 66.2%                                          #
# Success rate arm 3: 88.2%                                          #
#                                                                    #
# Best decision: Only OnTrack is superior to Usual care              #
# Posterior probability that decision minimizes loss: 86.8%          #
#                                                                    #
# Output of Tdata: 0.041 0.00007 0.868 0.091                         #
#                                                                    #
######################################################################


######################################################################
# ANALYSIS 2                                                         #
#                                                                    #
# Reanalyses as Bayesian decision-theoretic MAMS trial               #
######################################################################

# We do not allow dropping of the usual care arm (set dropctrl = FALSE)

dropctrl = FALSE

set.seed(19082021)

K=3
Ydata <- array(dim=c(maxpatients, K, 1))
Ydata[,,1] <- Y


# To efficiently evaluate the trial under a Bayesian-adaptive decision-theoretic approach 
# with interim analyses after 36 patients for different thresholds for continuation, we 
# first calculate expected increases in the proportion of decisions for all interim analyses 
# (up to maximum of 76 is reached for two arms)
#
# We do this both for a setting where arms can be dropped (stored in trial.out) and a
# a setting where all arms remain in the trial (stored in trial_nodrop.out)
#
# To force continuation up to the maximum trial size and avoid early stopping,
# we set gamma parameter at 0


set.seed(05012022)

trial.out        <- prdrop(Ydata, delta, gamma=0, burn=12, Ntrials=1, batch=12, prior, maxpt=1*12+2*76, dropctrl, G)
trial_nodrop.out <- prnodrop(Ydata, delta, gamma=0, burn=12, Ntrials=1, batch=12, prior, maxpt=3*76,G)

trial.out
trial_nodrop.out


#################################################
# EVALUATE TRIAL USING BAYESIAN-ADAPTIVE METHOD #
# BOTH WITH AND WITHOUT DROPPING                #
# GAMMA = 0.01 FOR EARLY STOPPING OF TRIAL      #
#################################################

gamma <- 0.01
nstages001.out = min(which(trial.out$delta_benefit < gamma))
nstages001_nodrop.out = min(which(trial_nodrop.out$delta_benefit < gamma))

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


save.image("Results PickAllTreatmentsSuperiorToControl.Rdata")



