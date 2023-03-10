# Supplementary R code for manuscript "Bayesian adaptive decision-theoretic 
# designs for multi-arm multi stage clinical trials" published in 
# Statistical Methods in Medical Research
#
# Authors: A. Bassi, J. Berkhof, D. de Jong & P.M. van de Ven
#
# Correspondence/enquiries: P.M.vandeVen-3@umcutrecht.nl



# This script contains all functions required for the setting with
# multiple (up to 5) experimental arms.
#
# The following functions are included:
# create.d  : Simulates trial data
# prdrop    : Evaluates trials with adaptive stopping when arms can be dropped
# prnodrop  : Evaluates trials with adaptive stopping without dropping of arms
# reevaluate: Reevaluates trials using output from prdrop and prnodrop
#             for higher thresholds gamma (=C/Q) 
# pcalc     : Calculates operating characteristics (samples sizes and final decisions) 
#             using output from prdrop, prnodrop and reevaluate
# prsinglestage: Evaluates single-stage trials with fixed sample size and 
#                final decision based on minimization of expected loss
# evaldrop  : Calculates the proportion of correct decisions in case of 
#		  stopping and under available options for continuation for 
#             an additional stage when for setting with dropping of arms
#             (given the observed data)   
# evalnodrop: Calculates proportion of correct decisions in case of 
#		  stopping and continuation for an additional stage 
#             for setting without dropping of arms
#             (given the observed data)  
# Tdata     : Determines the decision that minimizes the expected loss given the 
#             observed data (specific for setting with only experimental arms)
# Ttheta    : Determines the correct decision given a true parameter vector theta
#             (specific for setting with only experimental arms)



# Function: create.d
#
# This function simulates binary outcome data for [Ntrials] trials 
# with [Nperarm] subjects per arm and response vector [resp]
#
# Required input:
# resp    : Vector containing the true response rates (success probabilities)
#         : in each of the arms
# Nperarm : Total number of patients per arm
# Ntrials : Total number of trials to be simulated (simulation size)
#
# Returns: array with (binary) outcomes for all trials


create.d <- function(resp, Nperarm, Ntrials){
  K <- length(resp) #no. arms
  Y <- array(dim=c(Nperarm, K, Ntrials))
  for(l in 1:Ntrials){
    for(k in 1:K){
      Y[,k,l] <- rbinom(n = Nperarm, size = 1, prob = resp[k])
    }
  }
  return(Y)
}



# Function: prdrop
#
# This function evaluates trials simulated using the create.d function
# using the decision-theoretic criteria for adaptive stopping and dropping of arms
#
# Required input:
# Y       : Trial data simulated using the create.d function
# delta   : Margin (delta) above which experimental arms (arms 2 and 3) are considered   
#           superior to the control arm (arm 1).   
#		Set to 0 for symmetric case with experimental arms only.
# gamma   : Minimum increase in the probability of a correct decision 
#           that is required for continuation of the trial (=C/Q in paper) 
# burn    : Batch size per arm for stage 1 
#         : If specified as 0 then batch size specified for 
#         : subsequent stages (batch) is also used for stage 1
# batch   : Number of subjects for stage 2,3,4,... divided by the number of
#           arms at the start
#           The total number of subjects in each stage is determined as number 
# 	      of arms at the start of the trial multiplied by the value assigned 
#           to batch (batch corresponds to n/K in paper)
#	      After arms have been dropped patients are divided as balanced 
#           as possible over the active arms
# Ntrials : Number of simulated trials contained in array Y
# prior   : Parameters c(a,b) for Beta(a,b) prior 
# maxpt   : Maximum for total number of patients that can be 
#           included in the trial (cap)
# dropctrl: Specify as TRUE if dropping of the control arm is allowed
#           or if only experimental arms are considered
#
# Uses functions: evaldrop, Ttheta, Tdata
#
# Returns: vector containing for each simulated trial:
#          - the total number of patients allocated to each arm (first K entries)
#          - the number of responders per arm (entries K+1:2K)
#          plus list with
#          - reduction in probability of an incorrect decision expected in next stage
#            (for each interim evaluation)
#          - arms dropped in each stage (0 if no arm is dropped) 
#          - number patients allocated to each the arm in each stage


prdrop <- function(Y, delta = 0, gamma, burn = 0, Ntrials, batch, prior, maxpt, dropctrl = TRUE){
  cl<-makeCluster(8) #change the 8 to your number of CPU cores
  registerDoParallel(cl)
  LL <- array(dim=c(Ntrials,(2*(dim(Y)[2]) + 3)))
  print(Sys.time())
  LL <- foreach(i=1:Ntrials, .combine = rbind, .export=c('evaldrop', 'Ttheta', 'Tdata')) %dopar% {
    n <- batch  # batch size at start (PER ARM)
    arms <- dim(Y)[2]
    y <- Y[,,i]   
    TR <- 1
    dat <- rep(0, arms)
    len <- rep(max(burn-n, 0),arms)
    active=rep(TRUE, arms)
    drop <- 0
    ben <- NULL
    isdr <- NULL
    ns<-c();
    first=TRUE;
    cdt <- (sum(len)+arms*n <= maxpt)
    while((gamma==0 | TR > gamma) & cdt){
      if(drop[1] == 0){
        len <- len+n
          dat <- NULL
        for(j in 1:length(len)){
          dat <- c(dat, sum(y[1:len[j], j]))
        }
        pr <- evaldrop(y = dat, n1 = len, n2 = rep(n, arms), prior, delta, gamma, dropctrl)
        TR <- max(pr$ben.cont - pr$ben.stop, na.rm = T)
        drop <- which.max(pr$ben.cont - pr$ben.stop)
        drop <- drop - 1
        ben <- c(ben, TR)
        isdr <- c(isdr, drop)
        if(first) {nstage <- len; first=FALSE} else {nstage <- rep(n, arms)};
        ns<-rbind(ns, nstage)
        active[drop]=FALSE
      } else {
        old_len=len;
        n_active=sum(active)
        active_arms = (1:arms)[active]
        len[active_arms] <- len[active_arms]+ floor((n*arms)/(arms-length(drop)))
        n_remain = n*arms - floor((n*arms)/(arms-length(drop)))*n_active
        if(n_remain > 0) {arms_with_min = active_arms[which(min(len[active_arms])==len[active_arms])];
        n_active_arms_with_minimum = length(arms_with_min)
          if(n_remain == n_active_arms_with_minimum) {len[arms_with_min]=len[arms_with_min]+1} else 
            if(n_remain < n_active_arms_with_minimum) {randsamp = sample(arms_with_min, n_remain); len[randsamp]=len[randsamp]+1} else
              {len[arms_with_min]=len[arms_with_min]+1; 
               n_remain = n_remain - n_active_arms_with_minimum; 
               randsamp = sample(active_arms, n_remain); len[randsamp]=len[randsamp]+1}}
        dat <- NULL
        for(j in 1:length(len)){
          dat <- c(dat, sum(y[1:len[j], j]))
        }
        nf <- len-old_len
        pr <- evaldrop(dat, len, nf, prior, delta, gamma, dropctrl)
        TR <- max(pr$ben.cont - pr$ben.stop, na.rm = T)
        dp <- which.max(pr$ben.cont - pr$ben.stop)
        nstage <- nf;
        ns <- rbind(ns, nstage)
        if(dp != 1){
          drop <- c(drop,dp-1)
          isdr <- c(isdr,dp-1);
          active[drop]=FALSE
        } else(isdr <- c(isdr, 0))
        ben <- c(ben, TR)
      }
      cdt <- (sum(len)+arms*n <= maxpt)
    }
    c(len,dat, list(benefit = ben, drop=isdr, n = ns))
  }
  print(Sys.time())
  stopCluster(cl)
  rm(cl)
  return(LL)
 }

# Function: prnodrop
#
# This function evaluates trials simulated using the create.d function 
# using the decision-theoretic criteria for adaptive stopping but without dropping of arms
#
# Required input:
# Y       : Trial data simulated using the create.d function
# delta   : Margin (delta) above which experimental arms (arms 2 and 3) are considered   
#           superior to the control arm (arm 1).   
#		Set to 0 for symmetric case with experimental arms only.
# gamma   : Minimum increase in the probability of a correct decision 
#           that is required for continuation of the trial (=C/Q in paper) 
# burn    : Batch size per arm for stage 1 
#         : If specified as 0 then batch size specified for 
#         : subsequent stages (batch) is also used for stage 1
# batch   : Number of subjects for stage 2,3,4,... divided by the number of
#           arms at the start (batch corresponds to n/K in paper)
#           The total number of subjects in each stage is determined as number 
# 	      of arms at the start of the trial multiplied by the value assigned 
#           to batch 
#	      After arms have been dropped patients are divided as balanced 
#           as possible over the active arms
# Ntrials : Number of simulated trials contained in array Y
# prior   : Parameters c(a,b) for Beta(a,b) prior 
# maxpt   : Maximum for total number of patients that can be 
#           included in the trial (cap)
#
# Uses functions: evalnodrop, Ttheta, Tdata
#
# Returns: vector containing for each simulated trial:
#          - the total number of patients allocated to each arm (first K entries)
#          - the number of responders per arm (entries K+1:2K)
#          plus list with
#          - reduction in probability of an incorrect decision expected in next stage
#            (for each interim evaluation)
#          - arms dropped in each stage (all entries 0 as dropping is not considered)
#          - number patients allocated to each the arm in each stage


prnodrop <- function(Y, delta = 0, gamma, burn = 0, Ntrials, batch, prior, maxpt){
  cl<-makeCluster(8) #change the 8 to your number of CPU cores
  registerDoParallel(cl)
  LL <- array(dim=c(Ntrials,(2*(dim(Y)[2]) + 3)))
  print(Sys.time())
  LL <- foreach(i=1:Ntrials, .combine = rbind, .export=c('evalnodrop', 'Ttheta', 'Tdata')) %dopar% {
    n <- batch  # batch size at start (PER ARM)
    arms <- dim(Y)[2]
    y <- Y[,,i]
    TR <- 1
    dat <- rep(0, arms)
    len <- rep(max(burn-n, 0),arms)
    ben <- NULL
    isdr <- NULL
    ns<-c();
    first=TRUE;
    cdt <- (sum(len)+arms*n <= maxpt)
    while((gamma==0 | TR > gamma) & (sum(len) + arms*n <= maxpt) != 0){
      len <- len+n
      dat <- NULL
      for(j in 1:length(len)){
        dat <- c(dat, sum(y[1:len[j], j]))
      }
      pr <- evalnodrop(y = dat, n1 = len, n2 = rep(n, arms), prior, delta, gamma)
      TR <- max(pr$ben.cont - pr$ben.stop)
      ben <- c(ben, TR)
      isdr <- c(isdr, 0)      
      if(first) {nstage <- len; first=FALSE} else {nstage=rep(n, arms)};ns <- rbind(ns,nstage);
    }  
    c(len,dat, list(benefit = ben, drop=isdr, n = ns))
  }
  print(Sys.time())
  stopCluster(cl)
  rm(cl)
  return(LL)}

# Function: reevaluate
#
# This function evaluates simulated trials that have already been evaluated using
# prdrop or prnodrop function using a higher threshold gamma (=C/Q) for continuation
# This procedure is much faster than running prdrop or prnodrop with the new threshold
#
# Required input:
# trialres: List generated as output of prdrop of prnodrop functions
# Y       : Simulated trial data used to obtain trialres
# newgamma: New minimally required increase (C/Q) in probability  
#           of correct decision required for continuation
#
# Returns : List containing for each simulated trial the number of 
#           subjects allocated to each arm and the number of
#           successes per arm


reevaluate <- function(trialres, Y, newgamma){
 Q <- dim(trialres)[2]
 L <- dim(trialres)[1]
   arms <- (Q-3)/2
   Yrec <- array(dim = c(L,Q-3))
   for(i in 1:L){
      ben <- unlist(trialres[i,Q-2])
      drp <- unlist(trialres[i,Q-1])
      nstage <- matrix(unlist(trialres[i,Q]), ncol=arms, byrow=FALSE)
      if(length(which(ben <= newgamma))==0) {st<-length(ben)} else {st<-min(which(ben <= newgamma))}
      # with new gamma stopped after st stages
      if(st>1) {npt<-apply(nstage[1:st,], 2, sum)} else {npt<-nstage[1,]}  
      dat<-c()
      for(j in 1:arms) {
       dat[j]=sum(Y[1:npt[j],j,i])}  
       Yrec[i,] <- c(npt, dat)
    } #endfor
return(Yrec)}


# Function pcalc
#
# This functions calculates the frequentist operating characteristics of trials 
# evaluated using prdrop, prnodrop or reevaluate
# 
# Required input:
# LL      : List generated as output of prdrop, prnodrop or reevaluate 
#           function
# prior   : Parameters c(a,b) for Beta(a,b) prior 
# delta   : The margin (delta) above which experimental arms (arms 2,3,..) are considered   
#           superior to the control arm (arm 1). This should be the same value used as input 
#           for prdrop and prnodrop functions
# nr_dec  : Number of possible final decisions
#
# Uses function: Tdata
#
# Returns: list with vector 'dec' containing proportion of trials in 
#          which each of the final decisions is made and vector 'N' 
#          with total trial sizes


pcalc <- function(LL, prior, delta = 0, nr_dec){
	if(!is.numeric(LL)){
		Yrec <- matrix(as.numeric(LL[,1:(dim(LL)[2]-3)]), ncol = (dim(LL)[2]-3)) 
	} else{
		Yrec <- LL
	}   
	L <- dim(Yrec)[1]
	K <- dim(Yrec)[2]
	y <- Yrec[,((K/2)+1):K]
	n <- Yrec[,1:(K/2)]
	Pr <- array()
	for(j in 1:L){
		Pr[j] <- Tdata(y = y[j,], n = n[j,], prior, delta)
	}
	tbl <- table(factor(Pr, levels = 1:nr_dec))/L
	return (list(dec = tbl, N = apply(n,1,sum)))} 


# Function: prsinglestage
#
# This function calculates the proportion of each final decision for single-stage trials
# with fixed trial size using the decision-theoretic criterion
#
# Y       : Trial data simulated using the create.d function
# delta   : Margin (delta) above which experimental arms (arms 2,3,..) are considered   
#           superior to the control arm (arm 1).   
#		Set to 0 for symmetric case with experimental arms only.
# burn    : Total sample size per arm 
# Ntrials : Number of simulated trials contained in array Y
# prior   : Parameters c(a,b) for Beta(a,b) prior 
# nr_dec  : Number of possible final decisions
#
# Uses function: Tdata
#
# Returns: list with vector 'dec' containing proportion of trials in 
#          which each of the final decisions is made and 'N' 
#          which is the total trial size


prsinglestage <- function(Y, delta = 0, burn, Ntrials, prior, nr_dec){
  arms <- dim(Y)[2]
  cl<-makeCluster(2) #change the 2 to your number of CPU cores
  registerDoParallel(cl)
  LL <- array(dim=c(Ntrials,1))
  print(Sys.time())
  LL <- foreach(i=1:Ntrials, .combine = rbind, .export=c('Tdata')) %dopar% {
    n <- rep(burn, arms) # batch size at start (PER ARM)
    y <- apply(Y[,,i][1:burn,],2,sum)
    Ty = Tdata(y, n, prior, delta)
    c(decisions=Ty)
    }
   print(Sys.time())
   stopCluster(cl)
   rm(cl)
   tbl <- table(factor(LL[,1], levels = 1:nr_dec))/Ntrials
   return(list(dec = tbl, N = arms*burn))
}


# Function: evaldrop
#
# This function computes posterior (predictive) probabilities of a correct decision in case the 
# trial is stopped and in case the  trial is continued for an additional stage in settings where  
# dropping of arms is allowed.
#
# Required input:
# y       : Vector with number of successes observed so far in the trial
#           (number of responders/successes in each arm)
# n1      : Vector with total number of subjects for which outcomes
#           have been observed in each arm of the trial
# n2      : Vector with number of subjects to be included in 
#           each arm in the additional stage of the trial
#           Note 1: 0-entries correspond to arms that have been dropped
#           Note 2: The function currently redistributes the total number of patients 
#			  over the active arms to achieve a balanced distribution   
# prior   : Parameters c(a,b) for Beta(a,b) prior 
# delta   : Margin (delta) above which experimental arms (arms 2 and 3) are considered   
#           superior to the control arm (arm 1).   
#		Set to 0 for symmetric case with experimental arms only.
# gamma   : Minimum increase in probability of a correct decision that is 
#           required for continuation of the trial (C/Q in paper) 
# dropctrl: Specify as TRUE if dropping of the control arm is allowed
#           or if only experimental arms are considered
#
# Uses functions: Ttheta, Tdata
#
# Returns: List with posterior probability of a correct decision in case of stopping (ben.stop) 
#          and vector with predicted probabilities of correct decisions after continuing for an 
#          additional stage under different options for continuation (ben.cont)


evaldrop <- function(y, n1, n2, prior, delta = 0, gamma, dropctrl = TRUE){
   if(delta < 0){ 
     cat('ERROR: delta must be a non-negative number')
     return(NULL)}
   ndr <- sum(n2 == 0)
   dr <- which(n2 == 0)
   J <- length(y)          # total number of arms
   active <- rep(TRUE, J)
   active[dr]=FALSE
   active_arms = (1:J)[active]
   Jnew <- length(y) - ndr # number of active arms
   if(gamma>0) {G <- ceiling(2.5*(1/gamma))} else {G<-10000} # number of draws per arm
   a <- prior[1]
   b <- prior[2]
   benefit <- NA
   theta <- matrix(rbeta(G*J, shape1 = a + y, shape2 = b + n1 - y), ncol = J, byrow = T)
   Tt <- apply(theta, 1, Ttheta, delta)
   Ty <- Tdata(y = y, n = n1, prior, delta = delta)
   n2new=rep(0,J)
   n2new[active_arms] = floor(sum(n2)/Jnew)
   n_remain = sum(n2)-sum(n2new)
   if(n_remain > 0) {len=n1+n2new;arms_with_min = active_arms[which(min(len[active_arms])==len[active_arms])];
        n_active_arms_with_minimum = length(arms_with_min)
          if(n_remain == n_active_arms_with_minimum) {len[arms_with_min]=len[arms_with_min]+1; n2new=len-n1} else 
            if(n_remain < n_active_arms_with_minimum) {randsamp = sample(arms_with_min, n_remain); len[randsamp]=len[randsamp]+1; n2new=len-n1} else
              {len[arms_with_min]=len[arms_with_min]+1; 
               n_remain = n_remain - n_active_arms_with_minimum; 
               randsamp = sample(active_arms, n_remain); len[randsamp]=len[randsamp]+1}; n2new=len-n1}
   Ynew <- matrix(rbinom(n = G*J, size = n2new, prob = t(theta)), ncol = J, byrow = T)
   Tynew <- apply( t(t(Ynew) + y), 1, Tdata, n1+n2new, prior, delta)
   b0 <- sum(Ty == Tt)/G; 
   benefit[1] <- sum(Tynew == Tt)/G;
   if(Jnew > 2){
       if(ndr == 0){
         ifelse(dropctrl, idx <- (1:J), idx <- (2:J))  
            } else{ifelse(dropctrl, idx <- (1:J)[-(dr)], idx <- (2:J)[-(dr-1)])}
      for(k in idx){
         n2drop=rep(0,J)
         n2drop[active_arms] = floor(sum(n2)/(Jnew-1))
         n2drop[k]=0
         n_remain = sum(n2)-sum(n2drop) 
          if(n_remain > 0) {len=n1+n2drop;active_after_drop=active; active_after_drop[k]=FALSE; active_arms_after_drop = (1:J)[active_after_drop];
                arms_not_dropped_with_min = active_arms_after_drop[which(min(len[active_arms_after_drop])==len[active_arms_after_drop])];
                n_active_arms_with_minimum = length(arms_not_dropped_with_min)
                if(n_remain == n_active_arms_with_minimum) {len[arms_not_dropped_with_min]=len[arms_not_dropped_with_min]+1; n2drop=len-n1} else 
                if(n_remain < n_active_arms_with_minimum) {randsamp = sample(arms_not_dropped_with_min, n_remain); len[randsamp]=len[randsamp]+1; n2drop=len-n1} else
              {len[arms_not_dropped_with_min]=len[arms_not_dropped_with_min]+1; 
               n_remain = n_remain - n_active_arms_with_minimum; 
               randsamp = sample(active_arms_after_drop, n_remain); len[randsamp]=len[randsamp]+1}; n2drop=len-n1}
        Ynew1d <- matrix(rbinom(n = G*J, size = n2drop, prob = t(theta)), ncol = J, byrow = T)
        Tynew <- apply( t(t(Ynew1d) + y), 1, Tdata, n1 + n2drop, prior,  delta)
        benefit[k+1] <- sum(Tynew == Tt)/G;      }
   }
return(list(ben.stop = b0, ben.cont = benefit))}



# Function: evalnodrop
#
# This function computes posterior (predictive) probabilities of a correct decision in case the 
# trial is stopped and in case the  trial is continued for an additional stage in settings where
# dropping of arms is not considered.
#
# Required input:
# y       : Vector with number of successes observed so far in the trial
#           (number of responders/successes in each arm)
# n1      : Vector with total number of subjects for which outcomes
#           have been observed in each arm of the trial
# n2      : Vector with number of subjects to be included in 
#           each arm in the additional stage of the trial
# prior   : Parameters c(a,b) for Beta(a,b) prior 
# delta   : Margin (delta) above which experimental arms (arms 2,3) are considered   
#           superior to the control arm (arm 1).   
#		Set to 0 for symmetric case with experimental arms only.
# gamma   : Minimum increase in probability of a correct decision that is 
#           required for continuation of the trial (C/Q in paper) 
#
# Uses functions: Ttheta, Tdata
#
# Returns: list with posterior probability of correct decisions in case of stopping (ben.stop) 
#          and predicted probability of a correct decision in case of continuing for an addtional 
#          stage (ben.cont)


evalnodrop <- function(y, n1, n2, prior, delta = 0, gamma){
  if(delta < 0){ 
    cat("ERROR: delta must be a non-negative number \n")
    return(NULL)
  }
  J <- length(y)
  if(gamma>0) {G <- ceiling(2.5*(1/gamma))} else {G<-10000} # number of draws per arm
  a <- prior[1]
  b <- prior[2]
  benefit <- NA
  theta <- matrix(rbeta(G*J, shape1 = a + y, shape2 = b + n1 - y), ncol = J, byrow = T)
  Tt <- apply(theta, 1, Ttheta, delta)
  Ty <- Tdata(y = y, n = n1, prior, delta = delta)
  Ynew <- matrix(rbinom(n = G*J, size = n2, prob = t(theta)), ncol = J, byrow = T)
  Tynew <- apply( t(t(Ynew) + y), 1, Tdata, n1+n2, prior, delta)
  b0 <- sum(Ty == Tt)/G
  benefit[1] <- sum(Tynew == Tt)/G
  return(list(ben.stop = b0, ben.cont = benefit))}


# Function: Tdata
#
# This function returns the decision that minimizes the (expected) loss given the data. 
# The function included here is for the loss function that appears as equation (2) in the manuscript 
# and can be used when selecting from a set of multiple experimental arms the one with highest the response rate.
# Currently implemented for 3, 4 and 5 experimental arms
#
# Required input:
# y      : Vector with number of successes observed so far in the trial
#          (number of responders/successes in each arm)
# n      : Vector with total number of subjects for which outcomes
#          have been observed in each arm of the trial
# prior  : Parameters c(a,b) for Beta(a,b) prior 
# delta  : Not used here
#
# Returns: 1 if experimental arm 1 has highest posterior probability
#            of having the highest response rate
#          2 if experimental arm 2 has highest posterior probability
#            of having the highest response rate
#          3 if experimental arm 3 has highest posterior probability
#            of having the highest response rate
#          etc

Tdata <- function(y, n, prior, delta = 0, G = 1000, print=FALSE){
  if(length(y) != length(n)){
    cat("ERROR: dimension of y and n must coincide")
    return(NULL)
  }
  a <- prior[1]
  b <- prior[2]
  K <- length(y)
  theta <- matrix(rbeta(G*K, shape1 = a + y, shape2 = b + n - y), ncol = K, byrow = T)
  if(K==3) {prob <- c(mean(theta[,1]-theta[,2] >= -delta & theta[,1]-theta[,3] >  -delta), 
                      mean(theta[,2]-theta[,1] >  delta & theta[,2]-theta[,3] >= 0),
                      mean(theta[,3]-theta[,1] >= delta & theta[,3]-theta[,2] >  0))  } 
  
  if(K==4) {prob <- c(mean(theta[,1]-theta[,2] >= -delta & theta[,1]-theta[,3] >  -delta  & theta[,1]-theta[,4] >= -delta), 
                      mean(theta[,2]-theta[,1] >  delta & theta[,2]-theta[,3] >= 0  & theta[,2]-theta[,4] >  0),
                      mean(theta[,3]-theta[,1] >= delta & theta[,3]-theta[,2] >  0  & theta[,3]-theta[,4] >= 0),
                      mean(theta[,4]-theta[,1] >  delta & theta[,4]-theta[,2] >= 0  & theta[,4]-theta[,3] >  0)) }
  
  if(K==5) {prob <- c(mean(theta[,1]-theta[,2] >= -delta & theta[,1]-theta[,3] > -delta  & theta[,1]-theta[,4] >= -delta & theta[,1]-theta[,5] >  -delta), 
                      mean(theta[,2]-theta[,1] >  delta & theta[,2]-theta[,3] >= 0  & theta[,2]-theta[,4] >  0 & theta[,2]-theta[,5] >= 0),
                      mean(theta[,3]-theta[,1] >= delta & theta[,3]-theta[,2] >  0  & theta[,3]-theta[,4] >= 0 & theta[,3]-theta[,5] >  0),
                      mean(theta[,4]-theta[,1] >  delta & theta[,4]-theta[,2] >= 0  & theta[,4]-theta[,3] >  0 & theta[,4]-theta[,5] >= 0),
                      mean(theta[,5]-theta[,1] >= delta & theta[,5]-theta[,2] >  0  & theta[,5]-theta[,3] >= 0 & theta[,5]-theta[,4] >  0)) }
  
  best <- max(prob); indices_max = which(prob==best)
  if (length(indices_max)>1) {tm=sample(indices_max,1)} else {tm = indices_max[1]};if(print){print(prob)} 
  return(tm)}


# Function: Ttheta
#
# This function returns the decision that minimizes the (expected) loss given the true 
# parameter vector theta.
# The function included here is for the loss function that appears as equation (2) in the manuscript 
# and can be used when selecting from a set of multiple experimental arms the one with highest the response rate.
#
# Required input:
# theta  : Vector containing response/success probabilities for the  
#        : three arms 
# delta  : Not used here
#
# Returns: 1 if experimental arm 1 has highest response rate
#          2 if experimental arm 2 has highest response rate
#          3 if experimental arm 3 has highest response rate
#          etc
         
Ttheta <- function(theta, delta = 0){
	pos <- which.max(theta -c(0, rep(delta, length(theta)-1)))
  return(pos)}


