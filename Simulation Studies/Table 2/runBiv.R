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
source("FitImputeMI.R")
source("setupBiv.R")

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
  
  #### For Independent Sampling ####
  dat.samp <- SampledFun(d = dat, Nsample = Nsample.sim, NsPerStratumUniv = NsPerStratumUniv.sim,
                          Y1.cutoffsODS =  cutoffs.Y1.ODS, Y2.cutoffsODS = cutoffs.Y2.ODS,
                          Y1.cutoffsBDS =  cutoffs.Y1.BDS, Y2.cutoffsBDS = cutoffs.Y2.BDS)
  fulldat <- dat.samp$dat
  
  ##########################################################
  #### Independent sampling ################################
  ##########################################################
  
  # fit ML and ACML with nlme and Functions6.R
  cat("Fit Complete Data Only: ML and ACML")
  
  #### Random sampling with ML ####
  dat.longRan <- fulldat %>% filter(SampledRan == 1) %>% 
            gather(., key = "y_name", value = "value", Y1, Y2) %>%
            mutate(D1 = as.integer(y_name == "Y1"), D2 = as.integer(y_name == "Y2")) %>%
            arrange(id)
    
  lme.ML   <- lme(value ~ 0 + D1 + D1:grp + D1:time + D1:conf + D1:timegrp + D1:timeconf +
                             D2 + D2:grp + D2:time + D2:conf + D2:timegrp + D2:timeconf, data = dat.longRan,
                           random = ~ 0 + D1 + D1:time + D2 + D2:time | id,
                           na.action = na.omit, weights = varIdent(form = ~1 | y_name),
                           control = lmeControl(maxIter = 1e9, opt = c("optim")))
  CompleteCase.MLE  <- list(
                coefficients = fixef(lme.ML)[c("D1", "D1:grp", "D1:conf", "D1:time", "D1:timegrp","D1:timeconf",
                                               "D2", "grp:D2", "conf:D2", "time:D2", "timegrp:D2", "timeconf:D2")],
                covariance = diag(vcov(lme.ML)[c("D1", "D1:grp", "D1:conf", "D1:time", "D1:timegrp", "D1:timeconf",
                                                 "D2", "grp:D2", "conf:D2", "time:D2", "timegrp:D2", "timeconf:D2"),
                                               c("D1", "D1:grp", "D1:conf", "D1:time", "D1:timegrp", "D1:timeconf",
                                               "D2", "grp:D2", "conf:D2", "time:D2", "timegrp:D2", "timeconf:D2")]))
   
  #### ODS and BDS Slope Sampling with ACML ####
  params <- c(beta_y1.sim, beta_y2.sim)
  truevals.biv.ods  <- c(params, log(sigma.a1.sim), 
                     log(sigma.b1.sim), log(sigma.a2.sim), log(sigma.b1.sim),
                     fisherZ(rep(0.1, 6)),log(sigma.e1.sim), log(sigma.e2.sim))
  
  cat("ODS Sampling")
  dat.longODS <- fulldat %>% filter(SampledSlpODS == 1) %>%
    gather(., key = "y_name", value = "Y", Y1, Y2) %>%
    mutate(D1 = as.integer(y_name == "Y1"), D2 = as.integer(y_name == "Y2")) %>%
    arrange(id)
 
  dat.longODS$t1  <- dat.longODS$D1*dat.longODS$time
  dat.longODS$g1  <- dat.longODS$D1*dat.longODS$grp
  dat.longODS$c1  <- dat.longODS$D1*dat.longODS$conf
  dat.longODS$gt1 <- dat.longODS$D1*dat.longODS$timegrp
  dat.longODS$ct1 <- dat.longODS$D1*dat.longODS$timeconf
  dat.longODS$t2  <- dat.longODS$D2*dat.longODS$time
  dat.longODS$g2  <- dat.longODS$D2*dat.longODS$grp
  dat.longODS$c2  <- dat.longODS$D2*dat.longODS$conf
  dat.longODS$gt2 <- dat.longODS$D2*dat.longODS$timegrp
  dat.longODS$ct2 <- dat.longODS$D2*dat.longODS$timeconf
   
  dat.longODS$NoWeighting <- 1
  
  dat.longODS$MixProbLow  <- dat.samp$SampProbsODS[[1]][1]*(dat.longODS$var.sampled == "Y1") +
    dat.samp$SampProbsODS[[2]][1]*(dat.longODS$var.sampled == "Y2")
  dat.longODS$MixProbMid  <- dat.samp$SampProbsODS[[1]][2]*(dat.longODS$var.sampled == "Y1") +
    dat.samp$SampProbsODS[[2]][2]*(dat.longODS$var.sampled == "Y2")
  dat.longODS$MixProbHigh <- dat.samp$SampProbsODS[[1]][3]*(dat.longODS$var.sampled == "Y1") +
    dat.samp$SampProbsODS[[2]][3]*(dat.longODS$var.sampled == "Y2")
  SampProbMixODS <- cbind(dat.longODS$MixProbLow, dat.longODS$MixProbMid, dat.longODS$MixProbHigh)
  
  dat.longODS$MixCutoff1  <- cutoffs.Y1.ODS$SlpCutUniv[1]*(dat.longODS$var.sampled == "Y1") +
    cutoffs.Y2.ODS$SlpCutUniv[1]*(dat.longODS$var.sampled == "Y2")
  dat.longODS$MixCutoff2 <- cutoffs.Y1.ODS$SlpCutUniv[2]*(dat.longODS$var.sampled == "Y1") +
    cutoffs.Y2.ODS$SlpCutUniv[2]*(dat.longODS$var.sampled == "Y2")
  cutpointsMixODS <- cbind(dat.longODS$MixCutoff1, dat.longODS$MixCutoff2)
  
  dat.longODS$MixSlp  <- ifelse(dat.longODS$var.sampled == "Y1", "slope1", "slope2")
  w.functionMixSlpODS    <- dat.longODS$MixSlp
  
  acml.ods <- acml.lmem2(formula.fixed = Y ~ 0 + D1 + g1 + c1 + t1 + gt1 + ct1 + D2 + g2 + 
                           c2 + t2 + gt2 + ct2, 
                         formula.random = ~ 0 + D1 + t1 + D2 + t2,
                         data = dat.longODS, id = id, InitVals = truevals.biv.ods,
                         SampProb = SampProbMixODS, cutpoints = cutpointsMixODS,
                         Weights = NoWeighting, w.function = w.functionMixSlpODS)
  acml.ods  <- list(coefficients = acml.ods$coefficients[1:12], 
                    covariance = diag(acml.ods$covariance)[1:12])
  
  cat("BDS Sampling")
  
  # get the phase 1 coefficients for fitting BDS sampling with the ACML software
  dat.longFC <- fulldat %>%
    gather(., key = "y_name", value = "Y", Y1, Y2) %>%
    mutate(D1 = as.integer(y_name == "Y1"), D2 = as.integer(y_name == "Y2")) %>%
    arrange(id)
  
  dat.longFC$t1  <- dat.longFC$D1*dat.longFC$time
  dat.longFC$g1  <- dat.longFC$D1*dat.longFC$grp
  dat.longFC$c1  <- dat.longFC$D1*dat.longFC$conf
  dat.longFC$gt1 <- dat.longFC$D1*dat.longFC$timegrp
  dat.longFC$ct1 <- dat.longFC$D1*dat.longFC$timeconf
  dat.longFC$t2  <- dat.longFC$D2*dat.longFC$time
  dat.longFC$g2  <- dat.longFC$D2*dat.longFC$grp
  dat.longFC$c2  <- dat.longFC$D2*dat.longFC$conf
  dat.longFC$gt2 <- dat.longFC$D2*dat.longFC$timegrp
  dat.longFC$ct2 <- dat.longFC$D2*dat.longFC$timeconf
  
  dat.longFC$NoWeighting <- 1
  truevals.fc            <- c(params[c(1, 3, 4, 6, 7, 9, 10, 12)],
                              log(sigma.a1.sim), log(sigma.b1.sim), log(sigma.a2.sim), log(sigma.b1.sim),
                              fisherZ(rep(0.1, 6)),log(sigma.e1.sim), log(sigma.e2.sim))
  cutpointsFC  <- matrix(c(-70, 80), ncol = 2, byrow = TRUE, nrow = nrow(dat.longFC))
  SampProbFC   <- matrix(1, ncol=3, nrow = nrow(dat.longFC))
  w.functionFC <- rep("intercept1",  nrow(dat.longFC))   
  
  fit.fc2 <- acml.lmem2(formula.fixed = Y ~ 0 + D1 + c1 + t1 + ct1 + D2 + c2 + t2 + ct2,
                        formula.random = ~ 0 + D1 + t1 + D2 + t2,
                        data = dat.longFC, id=id,
                        InitVals=truevals.fc,
                        SampProb=SampProbFC, 
                        cutpoints=cutpointsFC, 
                        Weights = NoWeighting, 
                        w.function=w.functionFC)
  CoefPhase1 <- coef(fit.fc2)
  
  
  dat.longBDS <- fulldat %>% filter(SampledSlpBDS == 1) %>%
    gather(., key = "y_name", value = "Y", Y1, Y2) %>%
    mutate(D1 = as.integer(y_name == "Y1"), D2 = as.integer(y_name == "Y2")) %>%
    arrange(id)
  
  dat.longBDS$t1  <- dat.longBDS$D1*dat.longBDS$time
  dat.longBDS$g1  <- dat.longBDS$D1*dat.longBDS$grp
  dat.longBDS$c1  <- dat.longBDS$D1*dat.longBDS$conf
  dat.longBDS$gt1 <- dat.longBDS$D1*dat.longBDS$timegrp
  dat.longBDS$ct1 <- dat.longBDS$D1*dat.longBDS$timeconf
  dat.longBDS$t2  <- dat.longBDS$D2*dat.longBDS$time
  dat.longBDS$g2  <- dat.longBDS$D2*dat.longBDS$grp
  dat.longBDS$c2  <- dat.longBDS$D2*dat.longBDS$conf
  dat.longBDS$gt2 <- dat.longBDS$D2*dat.longBDS$timegrp
  dat.longBDS$ct2 <- dat.longBDS$D2*dat.longBDS$timeconf
  
  dat.longBDS$NoWeighting <- 1
  
  dat.longBDS$MixProbLow  <- dat.samp$SampProbsBDS[[1]][1]*(dat.longBDS$var.sampled == "Y1") +
    dat.samp$SampProbsBDS[[2]][1]*(dat.longBDS$var.sampled == "Y2")
  dat.longBDS$MixProbMid  <- dat.samp$SampProbsBDS[[1]][2]*(dat.longBDS$var.sampled == "Y1") +
    dat.samp$SampProbsBDS[[2]][2]*(dat.longBDS$var.sampled == "Y2")
  dat.longBDS$MixProbHigh <- dat.samp$SampProbsBDS[[1]][3]*(dat.longBDS$var.sampled == "Y1") +
    dat.samp$SampProbsBDS[[2]][3]*(dat.longBDS$var.sampled == "Y2")
  SampProbMixBDS <- cbind(dat.longBDS$MixProbLow, dat.longBDS$MixProbMid, dat.longBDS$MixProbHigh)
  
  
  dat.longBDS$MixCutoff1  <- cutoffs.Y1.BDS$SlpCutUniv[1]*(dat.longBDS$var.sampled == "Y1") +
    cutoffs.Y2.BDS$SlpCutUniv[1]*(dat.longBDS$var.sampled == "Y2")
  dat.longBDS$MixCutoff2 <- cutoffs.Y1.BDS$SlpCutUniv[2]*(dat.longBDS$var.sampled == "Y1") +
    cutoffs.Y2.BDS$SlpCutUniv[2]*(dat.longBDS$var.sampled == "Y2")
  cutpointsMixBDS <- cbind(dat.longBDS$MixCutoff1, dat.longBDS$MixCutoff2)
  
  dat.longBDS$MixSlp  <- ifelse(dat.longBDS$var.sampled == "Y1", "blup.slope1", "blup.slope2")
  w.functionMixSlpBDS <- dat.longBDS$MixSlp
  # need to change this
  acml.bds  <- acml.lmem2(formula.fixed = Y ~ 0 + D1 + g1 + c1 + t1 + gt1 + ct1 + D2 + g2 + c2 + t2 + gt2 + ct2,
                           formula.random = ~ 0 + D1 + t1 + D2 + t2,
                           data = dat.longBDS, id = id, InitVals = truevals.biv.ods,
                           SampProb = SampProbMixBDS, cutpoints = cutpointsMixBDS,
                           Weights = NoWeighting, 
                          w.function = w.functionMixSlpBDS,
                           xcol.phase1 = c(1, 3, 4, 6, 7, 9, 10, 12), ests.phase1 = CoefPhase1)
  acml.bds  <- list(coefficients = acml.bds$coefficients[1:12], 
                    covariance = diag(acml.bds$covariance)[1:12])

  
  #### Fit imputation iteratively for each study design ####
  
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
  
  dat.rs         <- dat.long
  dat.rs$grp     <- ifelse(dat.rs$SampledRan == 1, dat.rs$grp, NA)
  dat.rs$timegrp <- ifelse(dat.rs$SampledRan == 1, dat.rs$timegrp, NA)
  dat.rs$g1      <- ifelse(dat.rs$SampledRan == 1, dat.rs$g1, NA)
  dat.rs$g2      <- ifelse(dat.rs$SampledRan == 1, dat.rs$g2, NA)
  dat.rs$gt1     <- ifelse(dat.rs$SampledRan == 1, dat.rs$gt1, NA)
  dat.rs$gt2     <- ifelse(dat.rs$SampledRan == 1, dat.rs$gt2, NA)
   
  dat.ods         <- dat.long
  dat.ods$grp     <- ifelse(dat.ods$SampledSlpODS == 1, dat.ods$grp, NA)
  dat.ods$timegrp <- ifelse(dat.ods$SampledSlpODS == 1, dat.ods$timegrp, NA)
  dat.ods$g1      <- ifelse(dat.ods$SampledSlpODS == 1, dat.ods$g1, NA)
  dat.ods$g2      <- ifelse(dat.ods$SampledSlpODS == 1, dat.ods$g2, NA)
  dat.ods$gt1     <- ifelse(dat.ods$SampledSlpODS == 1, dat.ods$gt1, NA)
  dat.ods$gt2     <- ifelse(dat.ods$SampledSlpODS == 1, dat.ods$gt2, NA)
   
  dat.bds         <- dat.long
  dat.bds$grp     <- ifelse(dat.bds$SampledSlpBLUP == 1, dat.bds$grp, NA)
  dat.bds$timegrp <- ifelse(dat.bds$SampledSlpBLUP == 1, dat.bds$timegrp, NA)
  dat.bds$g1      <- ifelse(dat.bds$SampledSlpBLUP == 1, dat.bds$g1, NA)
  dat.bds$g2      <- ifelse(dat.bds$SampledSlpBLUP == 1, dat.bds$g2, NA)
  dat.bds$gt1     <- ifelse(dat.bds$SampledSlpBLUP == 1, dat.bds$gt1, NA)
  dat.bds$gt2     <- ifelse(dat.bds$SampledSlpBLUP == 1, dat.bds$gt2, NA)
   
  model.rs   <- MI(formula.fixed = Y ~ 0 + D1 + g1 + c1 + t1 + gt1 + ct1 + D2 + g2 + c2 + t2 + gt2 + ct2, 
                    formula.random = ~0 + D1 + t1 + D2 + t2 | id, 
                    formula.imp = grp ~ conf, data = dat.rs, 
                    n.imp = nrep.imp, n.burn = burn, id = "id", grp = c("g1", "g2"), 
                    timegrp = c("gt1", "gt2"), time = c("t1", "t2"), yname = "y_name")
  
  model.ods   <- MI(formula.fixed = Y ~ 0 + D1 + g1 + c1 + t1 + gt1 + ct1 + D2 + g2 + c2 + t2 + gt2 + ct2, 
                     formula.random = ~0 + D1 + t1 + D2 + t2 | id, 
                     formula.imp = grp ~ conf, data = dat.ods, 
                     n.imp = nrep.imp, n.burn = burn, id = "id", grp = c("g1", "g2"), 
                     timegrp = c("gt1", "gt2"), time = c("t1", "t2"), yname = "y_name")
   
  model.bds   <- MI(formula.fixed = Y ~ 0 + D1 + g1 + c1 + t1 + gt1 + ct1 + D2 + g2 + c2 + t2 + gt2 + ct2, 
                      formula.random = ~0 + D1 + t1 + D2 + t2 | id, 
                      formula.imp = grp ~ conf, data = dat.bds, 
                      n.imp = nrep.imp, n.burn = burn, id = "id", grp = c("g1", "g2"), 
                      timegrp = c("gt1", "gt2"), time = c("t1", "t2"), yname = "y_name")
   
  # get the results
  cat("Saving Results")
  fits    <- list(CompleteCase.MLE, acml.ods, acml.bds, model.rs, model.ods, model.bds)
  results <- cbind(Sampling = c("CompleteCaseMLE", "ACML ODS", "ACML BDS", 
                                "IIA RS", "IIA ODS", "IIA BDS"), 
                   coverage(fits, params), ests(fits), vars(fits))
  # save the results
  save(results, file=paste0("outputIIM/run-", run, "-", count, ".RData"))
}
