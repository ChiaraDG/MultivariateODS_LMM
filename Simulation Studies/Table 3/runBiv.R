########################################################################################################################
#### Run Simulation ####################################################################################################
########################################################################################################################

# load all packages needed
library(nlme)
library(mvtnorm)
library(mitools)
library(MASS)
library(tidyverse)
library(Matrix)
library(numDeriv)

# load R files with function
source("ACML.R")
source("DatGenFunBiv.R")
source("setupBiv.R")
source("FitImputeMI.R")

fisherZ <- function(x) log((1+x)/(1-x))
n       <- 1

########################################################################################################################
#### Function to return summary we need (from Shawn Garbett's code) ###################################################
########################################################################################################################

## The estimates returned from list of fits
ests     <- function(fits)
{
  result <- t(sapply(fits, function(f) f$coefficients[1:12]))
  colnames(result) <- c("est.int.Y1", "est.grp.Y1", "est.conf.Y1", "est.time.Y1", "est.tg.Y1", "est.tc.Y1",
                        "est.int.Y2", "est.grp.Y2", "est.conf.Y2", "est.time.Y2", "est.tg.Y2", "est.tc.Y2")
  as.data.frame(result)
}

Ds    <- function(fits)
{
  result <- t(sapply(fits, function(f) f$D[1:16]))
  colnames(result) <- c("D1D1", "D1t1", "D1D2", "D1t2",
                        "t1D1", "t1t1", "t1D2", "t1t2",
                        "D2D1", "D2t1", "D2D2", "D2t2",
                        "t2D1", "t2t1", "t2D2", "t2t2")
  as.data.frame(result)
}

Ss    <- function(fits)
{
  result <- t(sapply(fits, function(f) f$Sigma[1:2]))
  colnames(result) <- c("S1", "S2")
  as.data.frame(result)
}

## Variance of estimates from list of fits (from Rubin's rule)
vars     <- function(fits)
{
  result <- t(sapply(fits, function(f) f$covariance[1:12]))
  colnames(result) <- c("var.int.Y1", "var.grp.Y1", "var.conf.Y1", "var.time.Y1", "var.tg.Y1", "var.tc.Y1",
                        "var.int.Y2", "var.grp.Y2", "var.conf.Y2", "var.time.Y2", "var.tg.Y2", "var.tc.Y2")
  as.data.frame(result)
}


## Coverage of fit
covered <- function(fit, truth)
{
  rng   <- 1:12
  ses   <- sqrt(fit$covariance)[rng]
  lci   <- fit$coefficients[rng] - qnorm(.975)*ses
  uci   <- fit$coefficients[rng] + qnorm(.975)*ses
  
  (truth[rng] >= lci) & (truth[rng] <= uci)
}

## Coverage of a list of fits
coverage <- function(fits, truth)
{
  result <- t(sapply(fits, function(f) covered(f, truth)))
  colnames(result) <- c("cover.int.Y1", "cover.grp.Y1", "cover.conf.Y1", "cover.time.Y1", "cover.tg.Y1", "cover.tc.Y1",
                        "cover.int.Y2", "cover.grp.Y2", "cover.conf.Y2", "cover.time.Y2", "cover.tg.Y2", "cover.tc.Y2")
  as.data.frame(result)
}

########################################################################################################################
#### Run Model #########################################################################################################
########################################################################################################################

simulation <- function(run, count)
{
  
  cat("Generating population cutoff")
  ## Generate cutoffs using a bigger population of 25000 for both outcomes
  dat.tmp <- datGenBiv(N = 25000, ni = ni.sim, prob = prev.grp.sim, c.parm = conf.param.sim,
                          beta_y1 = beta_y1.sim, beta_y2 = beta_y2.sim,
                          sigma2.a1 = sigma.a1.sim, sigma2.a2 = sigma.a2.sim, sigma2.b1 = sigma.b1.sim,
                          sigma.a1a2 = sigma.a1a2.sim, sigma.a1b1 = sigma.a1b1.sim,
                          sigma.a1b2 = sigma.a1b2.sim, sigma.a2b1 = sigma.a2b1.sim,
                          sigma.a2b2 = sigma.a2b2.sim, sigma.b1b2 = sigma.b1b2.sim,
                          sigma2.b2 = sigma.b2.sim, sigma2.e1 = sigma.e1.sim, sigma2.e2 = sigma.e2.sim,
                          rho.e = rho.e.sim, type = type)
  # estimate cutoffs for ODS and BLUP sampling
  cutoffs.Y1.ODS   <- est.cutoffsODS(Y = dat.tmp$Y1, time = dat.tmp$time, id = dat.tmp$id,
                                       PropInCentralRegion = p.central)
  cutoffs.Y2.ODS   <- est.cutoffsODS(Y = dat.tmp$Y2, time = dat.tmp$time,
                                       id = dat.tmp$id, PropInCentralRegion = p.central)
  cutoffs.Y1.BDS  <- est.cutoffsBDS(Y = dat.tmp$Y1, time = dat.tmp$time, conf = dat.tmp$conf,
                                        id = dat.tmp$id, PropInCentralRegion = p.central)
  cutoffs.Y2.BDS  <- est.cutoffsBDS(Y = dat.tmp$Y2, time = dat.tmp$time, conf = dat.tmp$conf,
                                        id = dat.tmp$id, PropInCentralRegion = p.central)
  
  ##############################################################################################################
  cat("Generating random data from known parameters")
  dat <- datGenBiv(N = N.sim, ni = ni.sim, prob = prev.grp.sim, c.parm = conf.param.sim,
                  beta_y1 = beta_y1.sim, beta_y2 = beta_y2.sim,
                  sigma2.a1 = sigma.a1.sim, sigma2.a2 = sigma.a2.sim, sigma2.b1 = sigma.b1.sim, 
                  sigma.a1a2 = sigma.a1a2.sim, sigma.a1b1 = sigma.a1b1.sim,
                  sigma.a1b2 = sigma.a1b2.sim, sigma.a2b1 = sigma.a2b1.sim, 
                  sigma.a2b2 = sigma.a2b2.sim, sigma.b1b2 = sigma.b1b2.sim, 
                  sigma2.b2 = sigma.b2.sim, sigma2.e1 = sigma.e1.sim, sigma2.e2 = sigma.e2.sim, 
                  rho.e = rho.e.sim, type = type)
  
  
  # generate sampling indicators, and sampling probabilities
  cat("Generating sampling indicator")
  
  #### For re-using data collected in previous phase 1 sample ####
  dat.samp <- SampledFunY1(d = dat, NsPerStratumUniv = NsPerStratumUniv.sim, 
                           Y1.cutoffsODS =  cutoffs.Y1.ODS,
                           Y1.cutoffsBDS =  cutoffs.Y1.BDS)
  fulldat <- dat.samp$dat
  
  ##########################################################
  #### Re-using Data: fit MLE with Y2 only to show bias ####
  lme.Y2 <- lme(Y2 ~ grp + conf + time + timegrp + timeconf, 
                 data = subset(fulldat, SampledSlpODS == 1), 
                 random = ~ 1 + time | id,  
                 control = lmeControl(maxIter = 1e9, opt = c("optim")))
  Y2only.MLE  <- list(coefficients = c(rep(NA, 6), fixef(lme.Y2)), 
                       covariance = c(rep(NA, 6), diag(vcov(lme.Y2))))
  ##########################################################
      
  # #### ODS Slope Sampling with ACML ####
  params <- c(beta_y1.sim, beta_y2.sim)
  truevals.biv.ods  <- c(params, log(sigma.a1.sim), 
                       log(sigma.b1.sim), log(sigma.a2.sim), log(sigma.b1.sim),
                       fisherZ(rep(0.1, 6)),log(sigma.e1.sim), log(sigma.e2.sim))
  
  cat("ODS Sampling")
  
  dat.long <- fulldat %>%
           gather(., key = "y_name", value = "Y", Y1, Y2) %>%
           mutate(D1 = as.integer(y_name == "Y1"), D2 = as.integer(y_name == "Y2")) %>%
           arrange(id)
    
  # prepare the data for acml fitting
  dat.long$NoWeighting <- 1
  dat.long             <- dat.long %>% filter(SampledSlpODS == 1)
  dat.long$Slp.w       <- "slope1"
  w.functionSlp        <- dat.long$Slp.w
    
  # get the data in the right format
  dat.long$t1  <- dat.long$D1*dat.long$time
  dat.long$g1  <- dat.long$D1*dat.long$grp
  dat.long$c1  <- dat.long$D1*dat.long$conf
  dat.long$gt1 <- dat.long$D1*dat.long$timegrp
  dat.long$ct1 <- dat.long$D1*dat.long$timeconf
   
  dat.long$t2  <- dat.long$D2*dat.long$time
  dat.long$g2  <- dat.long$D2*dat.long$grp
  dat.long$c2  <- dat.long$D2*dat.long$conf
  dat.long$gt2 <- dat.long$D2*dat.long$timegrp
  dat.long$ct2 <- dat.long$D2*dat.long$timeconf
   
  dat.long$SlpProbLow <- dat.samp$SampProbsODS[[1]][1]
  dat.long$SlpProbMed <- dat.samp$SampProbsODS[[1]][2]
  dat.long$SlpProbHig <- dat.samp$SampProbsODS[[1]][3]
  SampProbSlp <- cbind(dat.long$SlpProbLow, dat.long$SlpProbMed, dat.long$SlpProbHig)
    
  dat.long$SlpCutoff1 <- cutoffs.Y1.ODS$SlpCutUniv[1]
  dat.long$SlpCutoff2 <- cutoffs.Y1.ODS$SlpCutUniv[2]
  cutpointsSlp <- cbind(dat.long$SlpCutoff1, dat.long$SlpCutoff2)
  # fit model
  acml.ods <- acml.lmem2(formula.fixed = Y ~ 0 + D1 + g1 + c1 + t1 + gt1 + ct1 + D2 + g2 + c2 + t2 + gt2 + ct2,
                              formula.random = ~ 0 + D1 + t1 + D2 + t2,
                              data = dat.long, id = id, InitVals = truevals.biv.ods,
                              SampProb = SampProbSlp, cutpoints = cutpointsSlp, 
                          Weights = NoWeighting, w.function = w.functionSlp)
    
  acml.ods  <- list(coefficients = acml.ods$coefficients[1:12], 
                     covariance = diag(acml.ods$covariance)[1:12])
   
  ##########################################################
  
  #### Fit imputation iteratively ####
  cat("Fit Multiple Imputation Iteratively")
   
  dat.long <- fulldat %>%
     gather(., key = "y_name", value = "Y", Y1, Y2) %>%
     mutate(D1 = as.integer(y_name == "Y1"), D2 = as.integer(y_name == "Y2")) %>%
     arrange(id)
   
  dat.long$t1  <- dat.long$D1*dat.long$time
  dat.long$g1  <- dat.long$D1*dat.long$grp
  dat.long$c1  <- dat.long$D1*dat.long$conf
  dat.long$gt1 <- dat.long$D1*dat.long$timegrp
  dat.long$ct1 <- dat.long$D1*dat.long$timeconf
   
  dat.long$t2  <- dat.long$D2*dat.long$time
  dat.long$g2  <- dat.long$D2*dat.long$grp
  dat.long$c2  <- dat.long$D2*dat.long$conf
  dat.long$gt2 <- dat.long$D2*dat.long$timegrp
  dat.long$ct2 <- dat.long$D2*dat.long$timeconf

  dat.ods         <- dat.long
  dat.ods$grp     <- ifelse(dat.ods$SampledSlpODS == 1, dat.ods$grp, NA)
  dat.ods$timegrp <- ifelse(dat.ods$SampledSlpODS == 1, dat.ods$timegrp, NA)
  dat.ods$g1      <- ifelse(dat.ods$SampledSlpODS == 1, dat.ods$g1, NA)
  dat.ods$g2      <- ifelse(dat.ods$SampledSlpODS == 1, dat.ods$g2, NA)
  dat.ods$gt1     <- ifelse(dat.ods$SampledSlpODS == 1, dat.ods$gt1, NA)
  dat.ods$gt2     <- ifelse(dat.ods$SampledSlpODS == 1, dat.ods$gt2, NA)
  
  model.ods       <- MI(formula.fixed = Y ~ 0 + D1 + g1 + c1 + t1 + gt1 + ct1 + D2 + g2 + c2 + t2 + gt2 + ct2, 
                         formula.random = ~0 + D1 + t1 + D2 + t2 | id, 
                         formula.imp = grp ~ conf, data = dat.ods, n.imp = 85,
                         iter = 10, id = "id", grp = c("g1", "g2"), 
                         timegrp = c("gt1", "gt2"), time = c("t1", "t2"), yname = "y_name")
   
  # get the results
  cat("Saving Results")
  fits <- list(Y2only.MLE, acml.ods, model.ods)
  results <- cbind(Sampling = c("Y2 Only", "ACML ODS", "MI ODS"), 
                   coverage(fits, params), ests(fits), vars(fits))
  
  # save the results
  save(results, file=paste0("outputIIM/run-", run, "-", count, ".RData"))
}
