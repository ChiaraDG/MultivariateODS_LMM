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


# load R files with function
#source("./OffsetImputation/LHS/Bivariate/FitImputeFunsOffsetBiv.R")
source("./OffsetImputation/BivariateLMM/FitImputeMI.R")
#source("./OffsetImputation/BivariateLMM/Functions6.R")

########################################################################################################################
#### Function to return summary we need (from Shawn Garbett's code) ###################################################
########################################################################################################################

## The estimates returned from list of fits
ests     <- function(fits)
{
  result <- t(sapply(fits, function(f) f$coefficients[1:42]))
  colnames(result) <- c("est.int.Y1", "est.grp.Y1", "est.time.Y1", "est.cigs0.10.Y1", "est.pack0.20.Y1",
                        "est.female.Y1", "est.age.10.Y1", "est.bmi.Y1", "est.bmi.change.Y1", "est.tf.Y1",
                        "est.ta.Y1", "est.tg.Y1", "est.site2.Y1", "est.site3.Y1", "est.site4.Y1", "est.site5.Y1",
                        "est.site6.Y1", "est.site7.Y1", "est.site8.Y1", "est.site9.Y1", "est.site10.Y1",
                        "est.int.Y2", "est.grp.Y2", "est.time.Y2", "est.cigs0.10.Y2", "est.pack0.20.Y2",
                        "est.female.Y2", "est.age.10.Y2", "est.bmi.Y2", "est.bmi.change.Y2", "est.tf.Y2",
                        "est.ta.Y2", "est.tg.Y2", "est.site2.Y2", "est.site3.Y2", "est.site4.Y2", "est.site5.Y2",
                        "est.site6.Y2", "est.site7.Y2", "est.site8.Y2", "est.site9.Y2", "est.site10.Y2")
  as.data.frame(result)
}


## Variance of estimates from list of fits
vars     <- function(fits)
{
  result <- t(sapply(fits, function(f) f$covariance[1:42]))
  colnames(result) <- c("var.int.Y1", "var.grp.Y1", "var.time.Y1", "var.cigs0.10.Y1", "var.pack0.20.Y1",
                        "var.female.Y1", "var.age.10.Y1", "var.bmi.Y1", "var.bmi.change.Y1", "var.tf.Y1",
                        "var.ta.Y1", "var.tg.Y1", "var.site2.Y1", "var.site3.Y1", "var.site4.Y1", "var.site5.Y1",
                        "var.site6.Y1", "var.site7.Y1", "var.site8.Y1", "var.site9.Y1", "var.site10.Y1",
                        "var.int.Y2", "var.grp.Y2", "var.time.Y2", "var.cigs0.10.Y2", "var.pack0.20.Y2",
                        "var.female.Y2", "var.age.10.Y2", "var.bmi.Y2", "var.bmi.change.Y2", "var.tf.Y2",
                        "var.ta.Y2", "var.tg.Y2", "var.site2.Y2", "var.site3.Y2", "var.site4.Y2", "var.site5.Y2",
                        "var.site6.Y2", "var.site7.Y2", "var.site8.Y2", "var.site9.Y2", "var.site10.Y2")
  as.data.frame(result)
}


## covariance of estimates from list of fits
covars     <- function(fits)
{
  result <- t(sapply(fits, function(f) f$covariance[1:6]))
  colnames(result) <- c("timegrp.cov.Y1", "timefemale.cov.Y1", "timeage.cov.Y1",
                        "timegrp.cov.Y2", "timefemale.cov.Y2", "timeage.cov.Y2")
  as.data.frame(result)
}
# 

fisherZ <- function(x) log((1+x)/(1-x))
n       <- 1

########################################################################################################################
#### Prepare the dataset for the analysis (from Jonathan's Schildcrout code) ###########################################
########################################################################################################################

SampledFun <- function(d){
  
  NsPerStratumUniv      <- c(150, 100, 150)
  Nsample.sim           <- 2*sum(NsPerStratumUniv)
  
  cat("Sample at Random")
  SampledRan            <- random.sampling(id.long = d$id, n = Nsample.sim)
  d$SampledRan          <- SampledRan
  
  ##############################
  #### Compute the cut-offs ####
  ##############################
  
  #### ODS Sampling ####
  cat("ODS Sample")
  IntSlps            <- CalcSSIntSlp(Y = d$Y1, time = d$time, id = d$id)
  Slp.Y1             <- IntSlps[[2]]
  cutpoints.Y1       <- quantile(Slp.Y1, probs = c(0.15, 0.85))
  
  
  IntSlps            <- CalcSSIntSlp(Y = d$Y2, time = d$time, id = d$id)
  Slp.Y2             <- IntSlps[[2]]
  cutpoints.Y2       <- quantile(Slp.Y2, probs = c(0.15, 0.85))
  
  #### BLUP Sampling ####
  mod <- lme(Y1 ~ time*(age.10 + grp + female) + cigs0.10 + pack0.20 + 
               site2 + site3 + site4 + site5 + site6 + site7 + site8 + site9 + site10 +
               bmi0.5 + bmi.change.5, random = ~ 1 + time | id, data = d,
             control = lmeControl(maxIter = 1e9, opt = c("optim")))
  random.slope.Y1   <- ranef(mod)[,2]
  cutpoints.BDS.Y1 <- quantile(random.slope.Y1, probs = c(0.15, 0.85), na.rm = TRUE)
  
  mod <- lme(Y2 ~ time*(age.10 + grp + female) + cigs0.10 + pack0.20 + 
               site2 + site3 + site4 + site5 + site6 + site7 + site8 + site9 + site10 +
               bmi0.5 + bmi.change.5, random = ~ 1 + time | id, data = d,
             control = lmeControl(maxIter = 1e9, opt = c("optim")))
  random.slope.Y2   <- ranef(mod)[,2]
  cutpoints.BDS.Y2 <- quantile(random.slope.Y2, probs = c(0.15, 0.85), na.rm = TRUE)
  
  
  #####################################
  #### Identify stratum membership ####
  #####################################
  
  #### ODS Sampling ####
  d$StratSlp.Y1 <- identify.stratum(w.function="slope", cutpoints = cutpoints.Y1, Slp = Slp.Y1)
  d$StratSlp.Y2 <- identify.stratum(w.function="slope", cutpoints = cutpoints.Y2, Slp = Slp.Y2)
  
  #### BLUP Sampling ####
  StratSlpBLUP       <- identify.stratum(w.function="slope", cutpoints = cutpoints.BDS.Y1, 
                                         Slp = random.slope.Y1)
  nobs               <- d %>% group_by(id) %>% summarise(nobs = n()) %>% pull(nobs)
  d$StratSlpBLUP.Y1  <- rep(StratSlpBLUP, times = nobs)
  
  StratSlpBLUP       <- identify.stratum(w.function="slope", cutpoints = cutpoints.BDS.Y2, 
                                         Slp = random.slope.Y2)
  d$StratSlpBLUP.Y2  <- rep(StratSlpBLUP, times = nobs)
  
  ##################################################################################
  #### divide N/2 to be sampled based on Y_1 and N/2 to be sampled based on Y_2 ####
  ##################################################################################
  
  nobs             <- d %>% group_by(id) %>% summarise(nobs = n()) %>% pull(nobs)
  d$var.sampled    <- rep(sample(c("Y1", "Y2"), size = length(unique(d$id)), replace = TRUE,
                                 prob = c(0.5, 0.5)), times = nobs)
  SampProbsODS     <- SampProbsBDS <- CutoffsODS <- CutoffsBDS <- list()
  
  #### First outcome ####
  dY1              <- subset(d, d$var.sampled == "Y1")
  # identify those sampled along with individual sampling probs and stratum sampling probs
  # ODS sampling
  SampledSlp <- ods.sampling(id.long = dY1$id, stratum.long = dY1$StratSlp.Y1, 
                             SamplingStrategy = "IndepODS", NsPerStratum = NsPerStratumUniv)
  dY1$SampledSlpODS       <- SampledSlp[[1]]
  SampProbsODS[[1]]       <- SampledSlp$SampProbs
  CutoffsODS[[1]]         <- cutpoints.Y1
  
  
  # BLUP sampling
  SampledSlpBDS      <- ods.sampling(id.long = dY1$id, stratum.long = dY1$StratSlpBLUP.Y1, 
                                     SamplingStrategy = "IndepODS", NsPerStratum = NsPerStratumUniv)
  dY1$SampledSlpBLUP <- SampledSlpBDS$Sampled
  SampProbsBDS[[1]]  <- SampledSlpBDS$SampProbs
  CutoffsBDS[[1]]    <- cutpoints.BDS.Y1
  
  #### Second outcome ####
  dY2              <- subset(d, d$var.sampled == "Y2")
  # ODS sampling
  SampledSlp <- ods.sampling(id.long = dY2$id, stratum.long = dY2$StratSlp.Y2, 
                             SamplingStrategy = "IndepODS", NsPerStratum = NsPerStratumUniv)
  dY2$SampledSlpODS       <- SampledSlp[[1]]
  SampProbsODS[[2]]       <- SampledSlp$SampProbs
  CutoffsODS[[2]]         <- cutpoints.Y2
  
  # BLUP sampling
  SampledSlpBDS   <- ods.sampling(id.long = dY2$id, stratum.long = dY2$StratSlpBLUP.Y2, 
                                  SamplingStrategy = "IndepODS", NsPerStratum = NsPerStratumUniv)
  dY2$SampledSlpBLUP <- SampledSlp$Sampled
  SampProbsBDS[[2]]          <- SampledSlpBDS$SampProbs
  CutoffsBDS[[2]]         <- cutpoints.BDS.Y2
  
  ###################################################
  #### Put the data back together and sort by ID ####
  ###################################################
  
  d            <- bind_rows(dY1, dY2) %>% arrange(id)
  
  ################
  #### output ####
  ################
  out <- list(dat = d, SampProbsODS = SampProbsODS, CutoffsODS = CutoffsODS,
              SampProbsBDS = SampProbsBDS, CutoffsBDS = CutoffsBDS)
  return(out)
}

random.sampling <- function(id.long, n){
  s <- sample(unique(id.long), n)
  Sampled <- as.integer(id.long %in% s)
  return(Sampled)
}

LinRegFn <- function(data){  
  X  <- cbind(1, data$time)
  Xt <- t(X)
  solve(Xt %*% X) %*% Xt %*% data$Y
}

## calculate subject specific intercepts and slopes and output them with the same length as the longitudinal data
# INPUT: outcome, time and ID
# OUTPUT: subject specific intercept and slopes
CalcSSIntSlp <- function(Y, time, id){
  data.tmp  <- data.frame(id = id, Y = Y, time = time)
  data.list <- split(data.tmp, id)
  L.id      <- c(unlist(tapply(id,id,length)))
  mtx       <- matrix(unlist(lapply(data.list, LinRegFn)), byrow=TRUE, ncol=2)
  out       <- list(Int = rep(mtx[,1], L.id), Slp = rep(mtx[,2], L.id))
  return(out)
}


identify.stratum <- function(w.function, cutpoints, Int, Slp){
  ## under bivar sampling cutpoints should be of the form (int.lower, int.upper, slp.lower, slp.upper)
  if(w.function == "intercept"){
    stratum <- ifelse( Int <  cutpoints[1], 1, ifelse( Int >= cutpoints[2], 3, 2))}
  if(w.function == "slope"){
    stratum <- ifelse( Slp <  cutpoints[1], 1, ifelse( Slp >= cutpoints[2], 3, 2))}
  if(w.function=="bivar"){
    stratum <- ifelse(Int >= cutpoints[1] & Int < cutpoints[2] & Slp >=cutpoints[3] & 
                        Slp <cutpoints[4], 1, 2)}
  return(stratum)
}

ods.sampling <- function(id.long, stratum.long, SamplingStrategy, NsPerStratum){    
  # This should be of length 3 for univariate sampling and of length 2 for bivariate sampling
  
  strat.1       <- c(unlist(tapply(stratum.long, id.long, unique)))
  id.1          <- c(unlist(tapply(id.long, id.long, unique)))
  ni            <- c(unlist(tapply(id.long, id.long, length)))
  N             <- length(id.1)
  NPerStratum   <- c(unlist(tapply(id.1, strat.1, length)))
  SampleTooMany <- any(NsPerStratum>NPerStratum)
  if (SampleTooMany){ 
    print("Warning: You want to sample more people than you have in one of the strata.
          Sampling from that stratum with probability 1")
    WhichStratum <- which(NsPerStratum>NPerStratum)
    NsPerStratum[WhichStratum] <- NPerStratum[WhichStratum]}
  SampProbs <- NsPerStratum/NPerStratum
  SampProb.1 <- ifelse(strat.1==1, SampProbs[1],
                       ifelse(strat.1==2, SampProbs[2], SampProbs[3]))
  if (SamplingStrategy=="IndepODS"){Samp <- rbinom(N, 1, SampProb.1)}
  if (SamplingStrategy=="DepODS"){  Sampled.ids <- NULL
  for (mm in 1:length(NPerStratum)){
    Sampled.ids <- c( Sampled.ids, c(sample(id.1[strat.1==mm], NsPerStratum[mm], replace=FALSE)))}
  Samp <- ifelse(id.1 %in% Sampled.ids, 1, 0)}
  Sampled <- list(Sampled=rep(Samp, ni),SampProbi=rep(SampProb.1, ni), SampProbs=SampProbs)
  return(Sampled)
}

# # load dataset
load("./OffsetImputation/LHS/lhs4fev2.Rdata")
expit <- function(x) exp(x)/(1+exp(x))

dat              <- lhs4fev2
dat$fev10        <- dat$fev*10
dat$fvc10        <- dat$fvc*10
dat              <- dat[order(dat$id, dat$visit),]
dat$Y1           <- dat$fev10
dat$Y2           <- dat$fvc10
dat$time         <- dat$visit
dat$grp          <- dat$rs177852
dat$timegrp      <- dat$time*dat$grp
dat$bmi.change.5 <- dat$bmi5.change*10
dat$timeage      <- dat$time*dat$age.10
dat$timefemale   <- dat$time*dat$female


########################################################################################################################
#### Run Model #########################################################################################################
########################################################################################################################

simulation <- function(dataset = dat, run, count)
{
  cat("Prepare the data for ACL and IIM algorithm")
  dat.samp   <- SampledFun(dataset)
  fulldat    <- dat.samp[["dat"]]
  
  #####################################################################################################################
  ############################### random sampling #####################################################################
  #####################################################################################################################
  
  # cat("Fit Random Sampling with Maximum Likelihood")
  # dat.RS <- fulldat %>% filter(SampledRan == 1) %>%
  #   gather(., key = "y_name", value = "value", Y1, Y2) %>%
  #   mutate(D1 = as.integer(y_name == "Y1"), D2 = as.integer(y_name == "Y2")) %>%
  #   arrange(id)
  # 
  # lme.ML <- lme(value ~ -1 + D1 + D1:grp + D1:time + D1:cigs0.10 + D1:pack0.20 +
  #                 D1:female + D1:age.10 + D1:bmi0.5 +  D1:bmi.change.5 +
  #                 D1:timefemale + D1:timeage + D1:timegrp +
  #                 D1:site2 + D1:site3 + D1:site4 + D1:site5 + D1:site6 + D1:site7 +
  #                 D1:site8 + D1:site9 + D1:site10 +
  #                 D2 + D2:grp + D2:time + D2:cigs0.10 + D2:pack0.20 +
  #                 D2:female + D2:age.10 + D2:bmi0.5 + D2:bmi.change.5 +
  #                 D2:timefemale + D2:timeage + D2:timegrp +
  #                 D2:site2 + D2:site3 + D2:site4 + D2:site5 + D2:site6 + D2:site7 +
  #                 D2:site8 + D2:site9 + D2:site10,
  #               data = dat.RS, random = ~ 0 + D1 + D1:time + D2 + D2:time | id,
  #               weights = varIdent(form = ~1 | y_name), na.action = na.omit,
  #               control = lmeControl(maxIter = 1e9, opt = c("optim")))
  # # save coefficients and covariance
  # LME.b.y1 <- fixef(lme.ML)[c("D1", "D1:grp", "D1:time", "D1:cigs0.10", "D1:pack0.20",
  #                             "D1:female", "D1:age.10", "D1:bmi0.5", "D1:bmi.change.5",
  #                             "D1:timefemale", "D1:timeage", "D1:timegrp", "D1:site2", "D1:site3",
  #                             "D1:site4", "D1:site5", "D1:site6", "D1:site7", "D1:site8", "D1:site9", "D1:site10")]
  # 
  # LME.b.y2 <- fixef(lme.ML)[c("D2", "grp:D2", "time:D2", "cigs0.10:D2", "pack0.20:D2",
  #                             "female:D2", "age.10:D2", "bmi0.5:D2", "bmi.change.5:D2",
  #                             "timefemale:D2", "timeage:D2", "timegrp:D2", "site2:D2", "site3:D2",
  #                             "site4:D2", "site5:D2", "site6:D2", "site7:D2", "site8:D2", "site9:D2", "site10:D2")]
  # 
  # 
  # LME.var  <- diag(vcov(lme.ML)[c("D1", "D1:grp", "D1:time", "D1:cigs0.10", "D1:pack0.20",
  #                                 "D1:female", "D1:age.10", "D1:bmi0.5", "D1:bmi.change.5",
  #                                 "D1:timefemale", "D1:timeage", "D1:timegrp", "D1:site2", "D1:site3",
  #                                 "D1:site4", "D1:site5", "D1:site6", "D1:site7", "D1:site8", "D1:site9", "D1:site10",
  #                                 "D2", "grp:D2", "time:D2", "cigs0.10:D2", "pack0.20:D2",
  #                                 "female:D2", "age.10:D2", "bmi0.5:D2", "bmi.change.5:D2",
  #                                 "timefemale:D2", "timeage:D2", "timegrp:D2", "site2:D2", "site3:D2",
  #                                 "site4:D2", "site5:D2", "site6:D2", "site7:D2", "site8:D2", "site9:D2", "site10:D2"),
  #                               c("D1", "D1:grp", "D1:time", "D1:cigs0.10", "D1:pack0.20",
  #                                 "D1:female", "D1:age.10", "D1:bmi0.5", "D1:bmi.change.5",
  #                                 "D1:timefemale", "D1:timeage", "D1:timegrp", "D1:site2", "D1:site3",
  #                                 "D1:site4", "D1:site5", "D1:site6", "D1:site7", "D1:site8", "D1:site9", "D1:site10",
  #                                 "D2", "grp:D2", "time:D2", "cigs0.10:D2", "pack0.20:D2",
  #                                 "female:D2", "age.10:D2", "bmi0.5:D2", "bmi.change.5:D2",
  #                                 "timefemale:D2", "timeage:D2", "timegrp:D2", "site2:D2", "site3:D2",
  #                                 "site4:D2", "site5:D2", "site6:D2", "site7:D2", "site8:D2", "site9:D2", "site10:D2")])
  # 
  # LME.cov <- c(vcov(lme.ML)["D1:time", "D1:timegrp"], vcov(lme.ML)["D1:time", "D1:timefemale"],
  #              vcov(lme.ML)["D1:time", "D1:timeage"], vcov(lme.ML)["time:D2", "timegrp:D2"],
  #              vcov(lme.ML)["time:D2", "timefemale:D2"],
  #              vcov(lme.ML)["time:D2", "timeage:D2"])
  
  # save results in list
  #CompleteCase.MLE  <- list(coefficients = c(LME.b.y1, LME.b.y2), variance = LME.var, covariance = LME.cov)
  
  #####################################################################################################################
  ############################### BDS and ODS with ACL ################################################################
  #####################################################################################################################
  
  #cat("Fit ACML with ODS and BDS sampling")
  
  #params <- CompleteCase.MLE$coefficients
  
  # random effects covariance matrix
  #r.eff.cov   <- getVarCov(lme.ML)
  #trueval.biv <- c(params, log(diag(r.eff.cov)), fisherZ(rep(0.1, 6)),
  #                 log(summary(lme.ML)$sigma), log(summary(lme.ML)$sigma))
  #names(trueval.biv) <- NULL
  
  
  ############################# ODS Sampling ##########################################################################
  
  # dat.longODS <- fulldat %>% filter(SampledSlpODS == 1) %>%
  #        gather(., key = "y_name", value = "Y", Y1, Y2) %>%
  #        mutate(D1 = as.integer(y_name == "Y1"), D2 = as.integer(y_name == "Y2")) %>%
  #        arrange(id)
  #   
  # dat.longODS$NoWeighting <- 1
  # NoWeighting             <- dat.longODS$NoWeighting
  #  
  # # get the variables in the format needed for ACML
  # dat.longODS$g1             <- dat.longODS$D1*dat.longODS$grp
  # dat.longODS$t1             <- dat.longODS$D1*dat.longODS$time
  # dat.longODS$cigs0.101      <- dat.longODS$D1*dat.longODS$cigs0.10
  # dat.longODS$pack0.201      <- dat.longODS$D1*dat.longODS$pack0.20
  # dat.longODS$female1        <- dat.longODS$D1*dat.longODS$female
  # dat.longODS$age.101        <- dat.longODS$D1*dat.longODS$age.10
  # dat.longODS$bmi0.51        <- dat.longODS$D1*dat.longODS$bmi0.5
  # dat.longODS$bmi.change.51  <- dat.longODS$D1*dat.longODS$bmi.change.5
  # dat.longODS$tf1            <- dat.longODS$D1*dat.longODS$timefemale
  # dat.longODS$ta1            <- dat.longODS$D1*dat.longODS$timeage
  # dat.longODS$tg1            <- dat.longODS$D1*dat.longODS$timegrp
  # dat.longODS$site21         <- dat.longODS$D1*dat.longODS$site2
  # dat.longODS$site31         <- dat.longODS$D1*dat.longODS$site3
  # dat.longODS$site41         <- dat.longODS$D1*dat.longODS$site4
  # dat.longODS$site51         <- dat.longODS$D1*dat.longODS$site5
  # dat.longODS$site61         <- dat.longODS$D1*dat.longODS$site6
  # dat.longODS$site71         <- dat.longODS$D1*dat.longODS$site7
  # dat.longODS$site81         <- dat.longODS$D1*dat.longODS$site8
  # dat.longODS$site91         <- dat.longODS$D1*dat.longODS$site9
  # dat.longODS$site101        <- dat.longODS$D1*dat.longODS$site10
  # 
  # dat.longODS$g2             <- dat.longODS$D2*dat.longODS$grp
  # dat.longODS$t2             <- dat.longODS$D2*dat.longODS$time
  # dat.longODS$cigs0.102      <- dat.longODS$D2*dat.longODS$cigs0.10
  # dat.longODS$pack0.202      <- dat.longODS$D2*dat.longODS$pack0.20
  # dat.longODS$female2        <- dat.longODS$D2*dat.longODS$female
  # dat.longODS$age.102        <- dat.longODS$D2*dat.longODS$age.10
  # dat.longODS$bmi0.52        <- dat.longODS$D2*dat.longODS$bmi0.5
  # dat.longODS$bmi.change.52  <- dat.longODS$D2*dat.longODS$bmi.change.5
  # dat.longODS$tf2            <- dat.longODS$D2*dat.longODS$timefemale
  # dat.longODS$ta2            <- dat.longODS$D2*dat.longODS$timeage
  # dat.longODS$tg2            <- dat.longODS$D2*dat.longODS$timegrp
  # dat.longODS$site22         <- dat.longODS$D2*dat.longODS$site2
  # dat.longODS$site32         <- dat.longODS$D2*dat.longODS$site3
  # dat.longODS$site42         <- dat.longODS$D2*dat.longODS$site4
  # dat.longODS$site52         <- dat.longODS$D2*dat.longODS$site5
  # dat.longODS$site62         <- dat.longODS$D2*dat.longODS$site6
  # dat.longODS$site72         <- dat.longODS$D2*dat.longODS$site7
  # dat.longODS$site82         <- dat.longODS$D2*dat.longODS$site8
  # dat.longODS$site92         <- dat.longODS$D2*dat.longODS$site9
  # dat.longODS$site102        <- dat.longODS$D2*dat.longODS$site10
  # 
  # # Sampling Probabilities for ODS only
  # dat.longODS$ProbLow       <- dat.samp$SampProbsODS[[1]][1]*(dat.longODS$var.sampled == "Y1") +
  #     dat.samp$SampProbsODS[[2]][1]*(dat.longODS$var.sampled == "Y2")
  # dat.longODS$ProbMid       <- dat.samp$SampProbsODS[[1]][2]*(dat.longODS$var.sampled == "Y1") +
  #     dat.samp$SampProbsODS[[2]][2]*(dat.longODS$var.sampled == "Y2")
  # dat.longODS$ProbHigh      <- dat.samp$SampProbsODS[[1]][3]*(dat.longODS$var.sampled == "Y1") +
  #     dat.samp$SampProbsODS[[2]][3]*(dat.longODS$var.sampled == "Y2")
  # SampProbMix           <- cbind(dat.longODS$ProbLow, dat.longODS$ProbMid, dat.longODS$ProbHigh)
  #  
  # # Cutoff for ODS only
  # dat.longODS$CutoffLow    <- dat.samp$CutoffsODS[[1]][1]*(dat.longODS$var.sampled == "Y1") +
  #     dat.samp$CutoffsODS[[2]][1]*(dat.longODS$var.sampled == "Y2")
  # dat.longODS$CutoffHig    <- dat.samp$CutoffsODS[[1]][2]*(dat.longODS$var.sampled == "Y1") +
  #     dat.samp$CutoffsODS[[2]][2]*(dat.longODS$var.sampled == "Y2")
  # cutpointsMix          <- cbind(dat.longODS$CutoffLow, dat.longODS$CutoffHig)
  #  
  # dat.longODS$w.function <- ifelse(dat.longODS$var.sampled == "Y1", "slope1", "slope2")
  # w.functionMixSlp    <-  dat.longODS$w.function
  #  
  # # fit model
  # acml.ods <- acml.lmem2(formula.fixed = Y ~ 0 + D1 + g1 + t1 + cigs0.101 + pack0.201 +
  #                            female1 + age.101 +
  #                            bmi0.51 + bmi.change.51 + tf1 + ta1 + tg1 + site21 + site31 + site41 + site51 +
  #                            site61 + site71 + site81 + site91 + site101 +
  #                            D2 + g2 + t2 + cigs0.102 + pack0.202 + female2 + age.102 +
  #                            bmi0.52 + bmi.change.52 + tf2 + ta2 + tg2 + site22 + site32 + site42 + site52 +
  #                            site62 + site72 + site82 + site92 + site102,
  #                           formula.random = ~ 0 + D1 + t1 + D2 + t2,
  #                           data = dat.longODS, id = id, InitVals = trueval.biv,
  #                           SampProb = SampProbMix, cutpoints = cutpointsMix,
  #                           Weights = NoWeighting, w.function = w.functionMixSlp)
  # # save results
  # acml.ods  <- list(coefficients = acml.ods$coefficients[1:42], 
  #                   variance = diag(acml.ods$covariance)[1:42],
  #                    covariance = c(acml.ods$covariance[2, 11], acml.ods$covariance[2, 9],
  #                                   acml.ods$covariance[2, 10], acml.ods$covariance[23, 32],
  #                                   acml.ods$covariance[23, 30], acml.ods$covariance[23, 31]))
  
  #####################################################################################################################
  ############################### BDS sampling ########################################################################
  #####################################################################################################################
  
  # fit FC with phase one only variables to get the coefficients of the first phase
  # trueval.phase1 <- trueval.biv[-c(2,12,23,33)]
  # dat.longFC <- fulldat %>% 
  #   gather(., key = "y_name", value = "Y", Y1, Y2) %>%
  #   mutate(D1 = as.integer(y_name == "Y1"), D2 = as.integer(y_name == "Y2")) %>%
  #   arrange(id)
  # 
  # dat.longFC$NoWeighting <- 1
  # NoWeighting          <- dat.longFC$NoWeighting
  # 
  # # # get the variables in the format needed for ACML
  # dat.longFC$g1             <- dat.longFC$D1*dat.longFC$grp
  # dat.longFC$t1             <- dat.longFC$D1*dat.longFC$time
  # dat.longFC$cigs0.101      <- dat.longFC$D1*dat.longFC$cigs0.10
  # dat.longFC$pack0.201      <- dat.longFC$D1*dat.longFC$pack0.20
  # dat.longFC$female1        <- dat.longFC$D1*dat.longFC$female
  # dat.longFC$age.101        <- dat.longFC$D1*dat.longFC$age.10
  # dat.longFC$bmi0.51        <- dat.longFC$D1*dat.longFC$bmi0.5
  # dat.longFC$bmi.change.51  <- dat.longFC$D1*dat.longFC$bmi.change.5
  # dat.longFC$tf1            <- dat.longFC$D1*dat.longFC$timefemale
  # dat.longFC$ta1            <- dat.longFC$D1*dat.longFC$timeage
  # dat.longFC$tg1            <- dat.longFC$D1*dat.longFC$timegrp
  # dat.longFC$site21         <- dat.longFC$D1*dat.longFC$site2
  # dat.longFC$site31         <- dat.longFC$D1*dat.longFC$site3
  # dat.longFC$site41         <- dat.longFC$D1*dat.longFC$site4
  # dat.longFC$site51         <- dat.longFC$D1*dat.longFC$site5
  # dat.longFC$site61         <- dat.longFC$D1*dat.longFC$site6
  # dat.longFC$site71         <- dat.longFC$D1*dat.longFC$site7
  # dat.longFC$site81         <- dat.longFC$D1*dat.longFC$site8
  # dat.longFC$site91         <- dat.longFC$D1*dat.longFC$site9
  # dat.longFC$site101        <- dat.longFC$D1*dat.longFC$site10
  # 
  # dat.longFC$g2             <- dat.longFC$D2*dat.longFC$grp
  # dat.longFC$t2             <- dat.longFC$D2*dat.longFC$time
  # dat.longFC$cigs0.102      <- dat.longFC$D2*dat.longFC$cigs0.10
  # dat.longFC$pack0.202      <- dat.longFC$D2*dat.longFC$pack0.20
  # dat.longFC$female2        <- dat.longFC$D2*dat.longFC$female
  # dat.longFC$age.102        <- dat.longFC$D2*dat.longFC$age.10
  # dat.longFC$bmi0.52        <- dat.longFC$D2*dat.longFC$bmi0.5
  # dat.longFC$bmi.change.52  <- dat.longFC$D2*dat.longFC$bmi.change.5
  # dat.longFC$tf2            <- dat.longFC$D2*dat.longFC$timefemale
  # dat.longFC$ta2            <- dat.longFC$D2*dat.longFC$timeage
  # dat.longFC$tg2            <- dat.longFC$D2*dat.longFC$timegrp
  # dat.longFC$site22         <- dat.longFC$D2*dat.longFC$site2
  # dat.longFC$site32         <- dat.longFC$D2*dat.longFC$site3
  # dat.longFC$site42         <- dat.longFC$D2*dat.longFC$site4
  # dat.longFC$site52         <- dat.longFC$D2*dat.longFC$site5
  # dat.longFC$site62         <- dat.longFC$D2*dat.longFC$site6
  # dat.longFC$site72         <- dat.longFC$D2*dat.longFC$site7
  # dat.longFC$site82         <- dat.longFC$D2*dat.longFC$site8
  # dat.longFC$site92         <- dat.longFC$D2*dat.longFC$site9
  # dat.longFC$site102        <- dat.longFC$D2*dat.longFC$site10
  # 
  # cutpointsFC  <- matrix(c(-70, 80), ncol = 2, byrow = TRUE, nrow = nrow(dat.longFC))
  # SampProbFC   <- matrix(1, ncol=3, nrow = nrow(dat.longFC))
  # w.functionFC <- rep("intercept1",  nrow(dat.longFC))   
  # 
  # mod.FC.phase1 <- acml.lmem2(formula.fixed = Y ~ 0 + D1 + t1 + cigs0.101 + pack0.201 + female1 + age.101 +
  #                               bmi0.51 + bmi.change.51 + tf1 + ta1 + site21 + site31 + site41 + site51 +
  #                               site61 + site71 + site81 + site91 + site101 +
  #                               D2 + t2 + cigs0.102 + pack0.202 + female2 + age.102 +
  #                               bmi0.52 + bmi.change.52 + tf2 + ta2 + site22 + site32 + site42 + site52 +
  #                               site62 + site72 + site82 + site92 + site102,
  #                             formula.random = ~ 0 + D1 + t1 + D2 + t2,
  #                             data = dat.longFC, id = id, InitVals = trueval.phase1,
  #                             SampProb = SampProbFC, cutpoints = cutpointsFC,
  #                             Weights = NoWeighting, w.function = w.functionFC)
  # CoefPhase1 <- coef(mod.FC.phase1)
  # 
  # dat.longBDS <- fulldat %>% filter(SampledSlpBLUP == 1) %>%
  #   gather(., key = "y_name", value = "Y", Y1, Y2) %>%
  #   mutate(D1 = as.integer(y_name == "Y1"), D2 = as.integer(y_name == "Y2")) %>%
  #   arrange(id)
  # 
  # dat.longBDS$NoWeighting <- 1
  # NoWeighting          <- dat.longFC$NoWeighting
  # 
  # # # get the variables in the format needed for ACML
  # dat.longBDS$g1             <- dat.longBDS$D1*dat.longBDS$grp
  # dat.longBDS$t1             <- dat.longBDS$D1*dat.longBDS$time
  # dat.longBDS$cigs0.101      <- dat.longBDS$D1*dat.longBDS$cigs0.10
  # dat.longBDS$pack0.201      <- dat.longBDS$D1*dat.longBDS$pack0.20
  # dat.longBDS$female1        <- dat.longBDS$D1*dat.longBDS$female
  # dat.longBDS$age.101        <- dat.longBDS$D1*dat.longBDS$age.10
  # dat.longBDS$bmi0.51        <- dat.longBDS$D1*dat.longBDS$bmi0.5
  # dat.longBDS$bmi.change.51  <- dat.longBDS$D1*dat.longBDS$bmi.change.5
  # dat.longBDS$tf1            <- dat.longBDS$D1*dat.longBDS$timefemale
  # dat.longBDS$ta1            <- dat.longBDS$D1*dat.longBDS$timeage
  # dat.longBDS$tg1            <- dat.longBDS$D1*dat.longBDS$timegrp
  # dat.longBDS$site21         <- dat.longBDS$D1*dat.longBDS$site2
  # dat.longBDS$site31         <- dat.longBDS$D1*dat.longBDS$site3
  # dat.longBDS$site41         <- dat.longBDS$D1*dat.longBDS$site4
  # dat.longBDS$site51         <- dat.longBDS$D1*dat.longBDS$site5
  # dat.longBDS$site61         <- dat.longBDS$D1*dat.longBDS$site6
  # dat.longBDS$site71         <- dat.longBDS$D1*dat.longBDS$site7
  # dat.longBDS$site81         <- dat.longBDS$D1*dat.longBDS$site8
  # dat.longBDS$site91         <- dat.longBDS$D1*dat.longBDS$site9
  # dat.longBDS$site101        <- dat.longBDS$D1*dat.longBDS$site10
  # 
  # dat.longBDS$g2             <- dat.longBDS$D2*dat.longBDS$grp
  # dat.longBDS$t2             <- dat.longBDS$D2*dat.longBDS$time
  # dat.longBDS$cigs0.102      <- dat.longBDS$D2*dat.longBDS$cigs0.10
  # dat.longBDS$pack0.202      <- dat.longBDS$D2*dat.longBDS$pack0.20
  # dat.longBDS$female2        <- dat.longBDS$D2*dat.longBDS$female
  # dat.longBDS$age.102        <- dat.longBDS$D2*dat.longBDS$age.10
  # dat.longBDS$bmi0.52        <- dat.longBDS$D2*dat.longBDS$bmi0.5
  # dat.longBDS$bmi.change.52  <- dat.longBDS$D2*dat.longBDS$bmi.change.5
  # dat.longBDS$tf2            <- dat.longBDS$D2*dat.longBDS$timefemale
  # dat.longBDS$ta2            <- dat.longBDS$D2*dat.longBDS$timeage
  # dat.longBDS$tg2            <- dat.longBDS$D2*dat.longBDS$timegrp
  # dat.longBDS$site22         <- dat.longBDS$D2*dat.longBDS$site2
  # dat.longBDS$site32         <- dat.longBDS$D2*dat.longBDS$site3
  # dat.longBDS$site42         <- dat.longBDS$D2*dat.longBDS$site4
  # dat.longBDS$site52         <- dat.longBDS$D2*dat.longBDS$site5
  # dat.longBDS$site62         <- dat.longBDS$D2*dat.longBDS$site6
  # dat.longBDS$site72         <- dat.longBDS$D2*dat.longBDS$site7
  # dat.longBDS$site82         <- dat.longBDS$D2*dat.longBDS$site8
  # dat.longBDS$site92         <- dat.longBDS$D2*dat.longBDS$site9
  # dat.longBDS$site102        <- dat.longBDS$D2*dat.longBDS$site10
  # 
  # # BLUP Slope Sampling
  # dat.longBDS$MixProbLow  <- dat.samp$SampProbsBDS[[1]][1]*(dat.longBDS$var.sampled == "Y1") +
  #   dat.samp$SampProbsBDS[[2]][1]*(dat.longBDS$var.sampled == "Y2")
  # dat.longBDS$MixProbMid  <- dat.samp$SampProbsBDS[[1]][2]*(dat.longBDS$var.sampled == "Y1") +
  #   dat.samp$SampProbsBDS[[2]][2]*(dat.longBDS$var.sampled == "Y2")
  # dat.longBDS$MixProbHigh <- dat.samp$SampProbsBDS[[1]][3]*(dat.longBDS$var.sampled == "Y1") +
  #   dat.samp$SampProbsBDS[[2]][3]*(dat.longBDS$var.sampled == "Y2")
  # SampProbMixBDS <- cbind(dat.longBDS$MixProbLow, dat.longBDS$MixProbMid, dat.longBDS$MixProbHigh)
  # 
  # dat.longBDS$MixCutoff1  <- dat.samp$CutoffsBDS[[1]][1]*(dat.longBDS$var.sampled == "Y1") +
  #   dat.samp$CutoffsBDS[[2]][1]*(dat.longBDS$var.sampled == "Y2")
  # dat.longBDS$MixCutoff2 <- dat.samp$CutoffsBDS[[1]][2]*(dat.longBDS$var.sampled == "Y1") +
  #   dat.samp$CutoffsBDS[[2]][2]*(dat.longBDS$var.sampled == "Y2")
  # cutpointsMixBDS <- cbind(dat.longBDS$MixCutoff1, dat.longBDS$MixCutoff2)
  # 
  # dat.longBDS$MixSlp  <- ifelse(dat.longBDS$var.sampled == "Y1", "blup.slope1", "blup.slope2")
  # w.functionMixSlpBDS <- dat.longBDS$MixSlp
  # 
  # # fit model
  # acml.bds <- acml.lmem2(formula.fixed = Y ~ 0 + D1 + g1 + t1 + cigs0.101 + pack0.201 +
  #                          female1 + age.101 +
  #                          bmi0.51 + bmi.change.51 + tf1 + ta1 + tg1 + site21 + site31 + site41 + site51 +
  #                          site61 + site71 + site81 + site91 + site101 +
  #                          D2 + g2 + t2 + cigs0.102 + pack0.202 + female2 + age.102 +
  #                          bmi0.52 + bmi.change.52 + tf2 + ta2 + tg2 + site22 + site32 + site42 + site52 +
  #                          site62 + site72 + site82 + site92 + site102,
  #                        formula.random = ~ 0 + D1 + t1 + D2 + t2,
  #                        data = dat.longBDS, id = id, InitVals = trueval.biv,
  #                        SampProb = SampProbMixBDS, cutpoints = cutpointsMixBDS,
  #                        Weights = NoWeighting, w.function = w.functionMixSlpBDS,
  #                        xcol.phase1 = c(1, 3:11, 13:22, 24:32, 34:42), 
  #                        ests.phase1 = CoefPhase1)
  # # save results
  # acml.bds  <- list(coefficients = acml.bds$coefficients[1:42], 
  #                    variance = diag(acml.bds$covariance)[1:42],
  #                   covariance = c(acml.bds$covariance[2, 11], acml.bds$covariance[2, 9],
  #                                  acml.bds$covariance[2, 10], acml.bds$covariance[23, 32],
  #                                  acml.bds$covariance[23, 30], acml.bds$covariance[23, 31]))
  # 
  
  
  #####################################################################################################################
  ############################### Imputation ##########################################################################
  #####################################################################################################################
  
  cat("Fit Imputation iteratively")
  
  dat.long <- fulldat %>%
    gather(., key = "y_name", value = "Y", Y1, Y2) %>%
    mutate(D1 = as.integer(y_name == "Y1"), D2 = as.integer(y_name == "Y2")) %>%
    arrange(id)
  
  dat.long$g1             <- dat.long$D1*dat.long$grp
  dat.long$t1             <- dat.long$D1*dat.long$time
  dat.long$cigs0.101      <- dat.long$D1*dat.long$cigs0.10
  dat.long$pack0.201      <- dat.long$D1*dat.long$pack0.20
  dat.long$female1        <- dat.long$D1*dat.long$female
  dat.long$age.101        <- dat.long$D1*dat.long$age.10
  dat.long$bmi0.51        <- dat.long$D1*dat.long$bmi0.5
  dat.long$bmi.change.51  <- dat.long$D1*dat.long$bmi.change.5
  dat.long$tf1            <- dat.long$D1*dat.long$timefemale
  dat.long$ta1            <- dat.long$D1*dat.long$timeage
  dat.long$tg1            <- dat.long$D1*dat.long$timegrp
  dat.long$site21         <- dat.long$D1*dat.long$site2
  dat.long$site31         <- dat.long$D1*dat.long$site3
  dat.long$site41         <- dat.long$D1*dat.long$site4
  dat.long$site51         <- dat.long$D1*dat.long$site5
  dat.long$site61         <- dat.long$D1*dat.long$site6
  dat.long$site71         <- dat.long$D1*dat.long$site7
  dat.long$site81         <- dat.long$D1*dat.long$site8
  dat.long$site91         <- dat.long$D1*dat.long$site9
  dat.long$site101        <- dat.long$D1*dat.long$site10
  dat.long$g2             <- dat.long$D2*dat.long$grp
  dat.long$t2             <- dat.long$D2*dat.long$time
  dat.long$cigs0.102      <- dat.long$D2*dat.long$cigs0.10
  dat.long$pack0.202      <- dat.long$D2*dat.long$pack0.20
  dat.long$female2        <- dat.long$D2*dat.long$female
  dat.long$age.102        <- dat.long$D2*dat.long$age.10
  dat.long$bmi0.52        <- dat.long$D2*dat.long$bmi0.5
  dat.long$bmi.change.52  <- dat.long$D2*dat.long$bmi.change.5
  dat.long$tf2            <- dat.long$D2*dat.long$timefemale
  dat.long$ta2            <- dat.long$D2*dat.long$timeage
  dat.long$tg2            <- dat.long$D2*dat.long$timegrp
  dat.long$site22         <- dat.long$D2*dat.long$site2
  dat.long$site32         <- dat.long$D2*dat.long$site3
  dat.long$site42         <- dat.long$D2*dat.long$site4
  dat.long$site52         <- dat.long$D2*dat.long$site5
  dat.long$site62         <- dat.long$D2*dat.long$site6
  dat.long$site72         <- dat.long$D2*dat.long$site7
  dat.long$site82         <- dat.long$D2*dat.long$site8
  dat.long$site92         <- dat.long$D2*dat.long$site9
  dat.long$site102        <- dat.long$D2*dat.long$site10
  
  dat.rs         <- dat.long
  dat.rs$grp     <- ifelse(dat.rs$SampledRan == 1, dat.rs$grp, NA)
  dat.rs$timegrp <- ifelse(dat.rs$SampledRan == 1, dat.rs$timegrp, NA)
  dat.rs$g1      <- ifelse(dat.rs$SampledRan == 1, dat.rs$g1, NA)
  dat.rs$g2      <- ifelse(dat.rs$SampledRan == 1, dat.rs$g2, NA)
  dat.rs$gt1     <- ifelse(dat.rs$SampledRan == 1, dat.rs$tg1, NA)
  dat.rs$gt2     <- ifelse(dat.rs$SampledRan == 1, dat.rs$tg2, NA)
  
  dat.ods         <- dat.long
  dat.ods$grp     <- ifelse(dat.ods$SampledSlpODS == 1, dat.ods$grp, NA)
  dat.ods$timegrp <- ifelse(dat.ods$SampledSlpODS == 1, dat.ods$timegrp, NA)
  dat.ods$g1      <- ifelse(dat.ods$SampledSlpODS == 1, dat.ods$g1, NA)
  dat.ods$g2      <- ifelse(dat.ods$SampledSlpODS == 1, dat.ods$g2, NA)
  dat.ods$gt1     <- ifelse(dat.ods$SampledSlpODS == 1, dat.ods$tg1, NA)
  dat.ods$gt2     <- ifelse(dat.ods$SampledSlpODS == 1, dat.ods$tg2, NA)
  
  dat.bds         <- dat.long
  dat.bds$grp     <- ifelse(dat.bds$SampledSlpBLUP == 1, dat.bds$grp, NA)
  dat.bds$timegrp <- ifelse(dat.bds$SampledSlpBLUP == 1, dat.bds$timegrp, NA)
  dat.bds$g1      <- ifelse(dat.bds$SampledSlpBLUP == 1, dat.bds$g1, NA)
  dat.bds$g2      <- ifelse(dat.bds$SampledSlpBLUP == 1, dat.bds$g2, NA)
  dat.bds$gt1     <- ifelse(dat.bds$SampledSlpBLUP == 1, dat.bds$tg1, NA)
  dat.bds$gt2     <- ifelse(dat.bds$SampledSlpBLUP == 1, dat.bds$tg2, NA)
  
  
  model.rs   <- IIM(formula.fixed = Y ~ 0 + D1 + g1 + t1 + cigs0.101 + pack0.201 +
                      female1 + age.101 +
                      bmi0.51 + bmi.change.51 + tf1 + ta1 + gt1 + site21 + 
                      site31 + site41 + site51 +
                      site61 + site71 + site81 + site91 + site101 +
                      D2 + g2 + t2 + cigs0.102 + pack0.202 + female2 + age.102 +
                      bmi0.52 + bmi.change.52 + tf2 + ta2 + gt2 + site22 +
                      site32 + site42 + site52 +
                      site62 + site72 + site82 + site92 + site102,
                    formula.random = ~ 0 + D1 + t1 + D2 + t2 | id, 
                    formula.imp = grp ~ cigs0.10 + pack0.20 + female + 
                      age.10 + bmi0.5 + site2 + site3 + site4 + site5 + 
                      site6 + site7 + site8 + 
                      site9 + site10, 
                    data = dat.rs, n.imp = 5, iter = 5, id = "id", grp = c("g1", "g2"), 
                    timegrp = c("gt1", "gt2"), time = c("t1", "t2"), yname = "y_name")
  
  # model.ods   <- IIM(formula.fixed = Y ~ 0 + D1 + g1 + t1 + cigs0.101 + pack0.201 +
  #                     female1 + age.101 +
  #                     bmi0.51 + bmi.change.51 + tf1 + ta1 + tg1 + site21 + 
  #                     site31 + site41 + site51 +
  #                     site61 + site71 + site81 + site91 + site101 +
  #                     D2 + g2 + t2 + cigs0.102 + pack0.202 + female2 + age.102 +
  #                     bmi0.52 + bmi.change.52 + tf2 + ta2 + tg2 + site22 +
  #                     site32 + site42 + site52 +
  #                     site62 + site72 + site82 + site92 + site102,
  #                   formula.random = ~ 0 + D1 + t1 + D2 + t2 | id, 
  #                   formula.imp = grp ~ cigs0.10 + pack0.20 + female + 
  #                     age.10 + bmi0.5 + site2 + site3 + site4 + site5 + 
  #                     site6 + site7 + site8 + 
  #                     site9 + site10, 
  #                   data = dat.ods, n.imp = 5, iter = 5, id = "id", grp = c("g1", "g2"), 
  #                   timegrp = c("gt1", "gt2"), time = c("t1", "t2"), yname = "y_name")
  
  # model.bds   <- IIM(formula.fixed = Y ~ 0 + D1 + g1 + t1 + cigs0.101 + pack0.201 +
  #                     female1 + age.101 +
  #                     bmi0.51 + bmi.change.51 + tf1 + ta1 + tg1 + site21 + 
  #                     site31 + site41 + site51 +
  #                     site61 + site71 + site81 + site91 + site101 +
  #                     D2 + g2 + t2 + cigs0.102 + pack0.202 + female2 + age.102 +
  #                     bmi0.52 + bmi.change.52 + tf2 + ta2 + tg2 + site22 +
  #                     site32 + site42 + site52 +
  #                     site62 + site72 + site82 + site92 + site102,
  #                   formula.random = ~ 0 + D1 + t1 + D2 + t2 | id, 
  #                   formula.imp = grp ~ cigs0.10 + pack0.20 + female + 
  #                     age.10 + bmi0.5 + site2 + site3 + site4 + site5 + 
  #                     site6 + site7 + site8 + 
  #                     site9 + site10, 
  #                   data = dat.bds, n.imp = 5, iter = 5, id = "id", grp = c("g1", "g2"), 
  #                   timegrp = c("gt1", "gt2"), time = c("t1", "t2"), yname = "y_name")
  
  # get the results
  #cat("Saving Results")
  fits <- list(model.rs)
  results <- cbind(Sampling = c("IIM RS"), ests(fits), vars(fits))
  # save the results
  save(results, file=paste0("outputLHS/run-", run, "-", count, ".RData"))
}
