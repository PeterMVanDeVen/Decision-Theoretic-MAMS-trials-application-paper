# clear workspace

rm(list=ls())

# Load packages

library(foreign)
library(snow)
library(doSNOW)
library(doParallel)

# Load required functions (R script "Functions PickTheWinner.R")

# setwd("C:/wd") #set to own working directory where functions are stored

setwd("C:/Users/pven3/OneDrive - UMC Utrecht/Documenten/Onderzoek/Decision-theoretic design/Paper Laurien/Github/Improved") #set to own working directory where functions are stored
source("Functions PickTheWinner.R")

# Load PACES dataset

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

# Cross Table of results for all 230 subjects

all230 = table(arm, success)
all230

# Obtain the number of subjects per arm (n) and number of successes per arm (x) 

n = all230[,2]+all230[,1]
x = all230[,2]
maxpatients = min(n)

# Make matrix with outcomes for each arm in different columns (keeping ordering of inclusion in arms)

Y=matrix(c(
success[which(arm==0)][1:maxpatients],
success[which(arm==1)][1:maxpatients],
success[which(arm==2)][1:maxpatients]),ncol=3, byrow=FALSE)

# SET G: HIGHER G MORE PRECISE ESTIMATES

G = 10000000

######################################################################
# ANALYSIS 1                                                         #
#                                                                    #
# Bayesian decision-theoretic analysis of original trial             #
######################################################################

# PICK THE WINNER SETTING WITH DELTA = 0

set.seed(18082021)
prior = c(1,1)
delta = 0

# OBTAIN ESTIMATES OF SUCCES RATES

phat = x/n
rbind(x,n, phat)

# Function Tdata calculates the posterior probability of each final decision being the correct decision
# entry 1: Posterior probability of arm 1 (Usual Care) having highest success rate
# entry 2: Posterior probability of arm 2 (OncoMove) having highest success rate
# entry 3: Posterior probability of arm 3 (OnTarget) having highest success rate

Tdata(x, n, prior, delta, G = G, print=TRUE)

######################################################################
# Results                                                            #
#                                                                    #
# Bayesian decision-theoretic analysis of original trial             #
#                                                                    #
# Success rate arm 1: 66.2%                                          #
# Success rate arm 2: 66.2%                                          #
# Success rate arm 3: 88.2%                                          #
#                                                                    #
# Best decision: On Track is the treatment with highest success rate #
# Posterior probability that decision minimizes loss: 99.9%          #
#                                                                    #
# Output of Tdata: 0.0006 0.0006 0.9988                              #
#                                                                    #
######################################################################


######################################################################
# ANALYSIS 2                                                         #
#                                                                    #
# Reanalyses as Bayesian decision-theoretic MAMS trial               #
######################################################################

# In settings where early dropping of arms is allowed, we allow also 
# the usual care arm to be dropped (set dropctrl = TRUE)

dropctrl = TRUE

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

trial.out        <- prdrop(Ydata, delta, gamma=0, burn=12, Ntrials=1, batch=12, prior, maxpt=1*12+2*76, dropctrl)
trial_nodrop.out <- prnodrop(Ydata, delta, gamma=0, burn=12, Ntrials=1, batch=12, prior, maxpt=3*76)

trial.out
trial_nodrop.out

#################################################
# EVALUATE TRIAL USING BAYESIAN-ADAPTIVE METHOD #
# BOTH WITH AND WITHOUT DROPPING                #
# GAMMA = 0.01: EARLY STOPPING OF TRIAL         #
#################################################

gamma <- 0.01
nstages001.out = min(which(trial.out$delta_benefit < gamma))
nstages001_nodrop.out = min(which(trial_nodrop.out$delta_benefit < gamma))

#########################
# setting with dropping #
#########################

n001_drop <- apply(trial.out$n[1:nstages001.out,],2,sum)
x001_drop <- c(sum(Y[(1:n001_drop[1]),1]),sum(Y[(1:n001_drop[2]),2]),sum(Y[(1:n001_drop[3]),3]))
phat001_drop <- x001_drop/n001_drop

rbind(n001_drop, x001_drop, phat001_drop)

# ONCOMOVE DROPPED AT FIRST INTERIM ANALYSIS
# 20/30 (UC), 8/12 (OM), 28/30 (OT) SUCCESSES
# 72 patients 

set.seed(20220501)

Tdata(x001_drop, n001_drop, prior, delta, G = G, print=TRUE)

# Best decision: On Track best
# Posterior probability that decision minimizes loss: 97.9%
# Output Tdata: 0.046352 0.0166673 0.9786975

############################
# setting without dropping #
############################

n001_nodrop <- apply(trial_nodrop.out$n[1:nstages001_nodrop.out,],2,sum)
x001_nodrop <- c(sum(Y[(1:n001_nodrop[1]),1]),sum(Y[(1:n001_nodrop[2]),2]),sum(Y[(1:n001_nodrop[3]),3]))
phat001_nodrop <- x001_nodrop/n001_nodrop

rbind(n001_nodrop, x001_nodrop, phat001_nodrop)

# ONCOMOVE DROPPED AT FIRST INTERIM ANALYSIS
# 16/24 (UC), 16/24 (OM), 22/24 (OT) SUCCESSES
# 72 patients 

set.seed(20220601)

Tdata(x001_nodrop, n001_nodrop, prior, delta, G = G, print=TRUE)


# Best decision: On Track best treatment
# Posterior probability that decision minimizes loss: 96.6%
#
# Output Tdata: 0.0168189 0.0168255 0.9663556



#################################################
# EVALUATE TRIAL USING BAYESIAN-ADAPTIVE METHOD #
# BOTH WITH AND WITHOUT DROPPING                #
# GAMMA = 0.001: EARLY STOPPING OF TRIAL        #
#################################################

gamma <- 0.001
nstages0001.out = min(which(trial.out$delta_benefit < gamma))
nstages0001_nodrop.out = min(which(trial_nodrop.out$delta_benefit < gamma))

n0001_drop <- apply(trial.out$n[1:nstages0001.out,],2,sum)
x0001_drop <- c(sum(Y[(1:n0001_drop[1]),1]),sum(Y[(1:n0001_drop[2]),2]),sum(Y[(1:n0001_drop[3]),3]))
phat0001_drop <- x0001_drop/n0001_drop

rbind(n0001_drop, x0001_drop, phat0001_drop)

# ONCOMOVE DROPPED AT FIRST INTERIM ANALYSIS
# 20/30 (UC), 8/12 (OM), 28/30 (OT) SUCCESSES

set.seed(20220321)

Tdata(x0001_drop, n0001_drop, prior, delta, G = G, print=TRUE)

# Best decision: On Track best
# Posterior probability that decision minimizes loss: 97.9%
#
# Output Tdata: 0.0046303 0.0167122 0.9786575

############################
# Setting without dropping #
############################

n0001_nodrop <- apply(trial_nodrop.out$n[1:nstages0001_nodrop.out,],2,sum)
x0001_nodrop <- c(sum(Y[(1:n0001_nodrop[1]),1]),sum(Y[(1:n0001_nodrop[2]),2]),sum(Y[(1:n0001_nodrop[3]),3]))
phat0001_nodrop <- x0001_nodrop/n0001_nodrop

rbind(n0001_nodrop, x0001_nodrop, phat0001_nodrop)

# ONCOMOVE DROPPED AT FIRST INTERIM ANALYSIS
# 16/24 (UC), 16/24 (OM), 22/24 (OT) SUCCESSES

set.seed(20220601)

Tdata(x0001_nodrop, n0001_nodrop, prior, delta, G = G, print=TRUE)

# Best decision: On Track best
# Posterior probability that decision minimizes loss: 96.6%
#
# Output Tdata: 0.0168189 0.0168255 0.9663556

save.image("Results PACES PickTheWinner.Rdata")


