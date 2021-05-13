#########################################################################################################################
####### Generate the Data and Sample in a Bivariate Setting #############################################################
#########################################################################################################################

## Function to generate the data (int + time + time*grp + grp + conf + conf*time)
## INPUT:  N (# of subjects), ni (# observations for subjects), prob (prevalence of the expensive covariate),
##         c.parm (vector of parameters linking the expensive covariate with the confounder), beta_y1, beta_y2 
##         (true values of the fixed effects coefficients), sigma2.a1, sigma2.a2, sigma2.b1, sigma2.b2, sigma.a1a2, 
##         sigma.a1b1, sigma.a1b2, sigma.a2b1, sigma.a2b2, sigma.b1b2 (values needed to build the variance and covariance
##         matrix of the random effects), sigma2.e1, sigma2.e2, rho.e (values needed to build the variance and covariance
##         matrix of the error), type (can be balanced.complete, balanced.incomplete or unbalanced)
## OUTPUT: Dataset to be used in the analysis
datGenBiv <- function(N, ni, prob, c.parm,
                      beta_y1, beta_y2,
                      sigma2.a1, sigma2.a2, sigma2.b1, sigma2.b2, sigma.a1a2, 
                      sigma.a1b1,
                      sigma.a1b2, sigma.a2b1, sigma.a2b2, 
                      sigma.b1b2, sigma2.e1, sigma2.e2, rho.e, type){
  
  ###################################################################################################
  # maximum number of observations
  # generate maximum number pe ID, and drop the information we do not need for the
  # balanced and incomplete and unbalanced case
  m <- max(ni)
  # generate id
  id   <- rep(1:N, each = m)
  # balanced data
  if(type == "balanced.complete" | type == "balanced.incomplete") {
    time <- rep(c(0:(m - 1)), N)
  }
  # unbalanced data
  if(type == "unbalanced") {
    time <- runif(N*m, 0, m)
  }
  # generate group and confounder
  grp.tmp  <- rbinom(N, 1, prob)
  conf.tmp <- rnorm(N, c.parm[1]+grp.tmp*c.parm[2], 1) 
  grp      <- rep(grp.tmp, each = m)
  conf     <- rep(conf.tmp, each = m)
  
  # create X matrix
  X         <- data.frame(id = id, grp = grp,  conf = conf,  time = time)
  # add the time*group interaction
  X$timegrp  <- X$grp*X$time
  # add the time*conf interaction
  X$timeconf <- X$conf*X$time
  # sort by time and ID
  X <- X[order(X$id, X$time),]
  # remove ID from X
  X$id <- NULL

  ###################################################################################################
  # generate joint distribution of the RE
  # first list all element of G
  G.el <- c(sigma2.a1, sigma.a1b1, sigma2.b1, sigma.a1a2, sigma.a2b1, sigma2.a2, 
            sigma.a1b2, sigma.b1b2,
            sigma.a2b2, sigma2.b2) 
  # create a symmetric matrix
  G     <- matrix(0, ncol = 4, nrow = 4)
  G[upper.tri(G, diag = TRUE)] <- G.el
  G.low <- c(sigma.a1b1, sigma.a1a2, sigma.a1b2, sigma.a2b1, sigma.b1b2, sigma.a2b2)
  G[lower.tri(G, diag = FALSE)] <- G.low
  re    <- mvrnorm(n = N,  mu = rep(0, 4), Sigma = G)
  # re for the first outcome
  a1 <- rep(re[, 1], each = m)
  b1 <- rep(re[, 2], each = m)
  # re for the secod outcome
  a2 <- rep(re[, 3], each = m)
  b2 <- rep(re[, 4], each = m)
  ###################################################################################################
  # generate error term
  # Error normally distributed
  R   <-  matrix(c(sigma2.e1, rho.e*sqrt(sigma2.e1*sigma2.e2), 
                   rho.e*sqrt(sigma2.e1*sigma2.e2), sigma2.e2), nrow = 2, byrow = TRUE)
  err <-  mvrnorm(n = N*m,  mu = rep(0, 2), Sigma = R)
  # generate the two outcomes
  Y1 <- cbind(1, as.matrix(X))%*%beta_y1 + a1 + b1*X[, "time"] + err[, 1]
  Y2 <- cbind(1, as.matrix(X))%*%beta_y2 + a2 + b2*X[, "time"] + err[, 2]
  # put results in a dataset to use for the analysis
  dat <- data.frame(id = rep(1:N, each = m), grp = X$grp, conf = X$conf, time = X$time, 
                    timegrp = X$timegrp, timeconf = X$timeconf, Y1 = Y1, Y2 = Y2)
  # sort by time and ID
  dat <- dat[order(dat$id, dat$time),]
  
  # if we need to generate balanced and incomplete
  if(type == "balanced.incomplete"){
    droptime <- rep(sample(seq(ni[1], ni[2]), N, replace=TRUE), each=ni[2]) - 1
    dat <- dat[dat$time<=droptime,]
  }
  # if we need to generate unbalanced
  if(type == "unbalanced"){
    dat <- create.unbalanced(dat, ni)
    dat <- dat[order(dat$id, dat$time),]
    dat$n <- NULL
  }
  
  # output 
  out <- dat
  return(out)
}



## Function to generate unbalanced data 
## INPUT:  d (dataset with balanced and complete data), ni (vector of interger within whom time points are generated)
## OUTPUT: dataset to be used in the analysis 
## Note:   input data needs to have a variable named id
create.unbalanced <- function(d, ni){
  obseach <- sample(min(ni):max(ni), length(unique(d$id)), replace = TRUE)
  d <- d %>% group_by(id) %>% nest() %>% ungroup() %>%  mutate(n = obseach) %>% 
    mutate(samp = map2(data, n, sample_n)) %>% dplyr::select(-data) %>% unnest(samp) %>% ungroup()
  return(d)
}



## Function for simple random sampling
## INPUT:  (id.long) id in long format, (n) number of subject to sample using srs
## OUTPUT: id of the sampled subjects
random.sampling <- function(id.long, n){
  s <- sample(unique(id.long), n)
  Sampled <- as.integer(id.long %in% s)
  return(Sampled)
}



## Function to get subject-specific intercepts and slopes for ODS sampling for subject i
## INPUT:  dataset
## OUTPUT: estimated intercep and slope
## Note:   data must contain variables named Y and time
LinRegFn <- function(data){  
  X  <- cbind(1, data$time)
  Xt <- t(X)
  solve(Xt %*% X) %*% Xt %*% data$Y
}



## Calculate subject-specific intercepts and slopes for every subjects
## INPUT: Y (vector of outcome variable in the dataset), time (time varianble in the dataset) 
##        and id (ID variable in the dataset)
## OUTPUT: subject-specific intercept and slopes for all subjects in the dataset
CalcSSIntSlp <- function(Y, time, id){
  data.tmp  <- data.frame(id = id, Y = Y, time = time)
  data.list <- split(data.tmp, id)
  L.id      <- c(unlist(tapply(id,id,length)))
  mtx       <- matrix(unlist(lapply(data.list, LinRegFn)), byrow=TRUE, ncol=2)
  out       <- list(Int = rep(mtx[,1], L.id), Slp = rep(mtx[,2], L.id))
  return(out)
}



## Function to generate the cutoffs for ODS sampling
## INPUT: Y (vector of outcome variable in the dataset), time variable in the dataset (vector of time variable in the 
##        dataset), id (vector of subject ID in the dataset), PropInCentralRegion (the proportion of subjects we want 
##        in the central region)
## OUTPUT: Cutoff points for ODS intercept, and ODS slope 
est.cutoffsODS <- function(Y, time, id, PropInCentralRegion){
  p         <- PropInCentralRegion
  data.tmp  <- data.frame(id=id, Y=Y, time=time)
  data.list <- split(data.tmp, id)
  print("Running individual regressions")
  out      <- matrix(unlist(lapply(data.list, LinRegFn)), byrow=TRUE, ncol=2)
  
  Ints <- quantile(out[,1], c((1-p)/2,(1+p)/2))
  Slps <- quantile(out[,2], c((1-p)/2,(1+p)/2))
  
  q1 <- .99
  Del <- 1
  while (Del>0.003){
    q1 <- q1-.001
    Del <- abs(mean(out[,1] > quantile(out[,1], probs=1-q1) & out[,1] < quantile(out[,1], probs=q1) &
                      out[,2] > quantile(out[,2], probs=1-q1) & out[,2] < quantile(out[,2], probs=q1)) - PropInCentralRegion)
  }
  
  out <- list(IntCutUniv = Ints,
              SlpCutUniv = Slps,
              IntCutBiv   = c(quantile(out[,1], probs=1-q1), quantile(out[,1], probs=q1)),
              SlpCutBiv   = c(quantile(out[,2], probs=1-q1), quantile(out[,2], probs=q1)))
  out
}



## Function to generate the cutoffs for BDS slope sampling
## INPUT:  Y (vector of outcome variable in the dataset), time variable in the dataset (vector of time variable in the 
##         dataset), id (vector of subject ID in the dataset), conf (vector of confounders in the dataset) 
##         PropInCentralRegion (the proportion of subjects we want in the central region)  
## OUTPUT: Cutoff points for BDS slope sampling
est.cutoffsBDS <- function(Y, time, id, conf, PropInCentralRegion){
   p         <- PropInCentralRegion
   
   # lme model with time and confounder
   print("Running LME model")
   lme.mod  <- lme(Y ~ time*conf, random = ~ 1 + time | id, control = lmeControl(maxIter = 1e9, opt = c("optim")))
   out      <- ranef(lme.mod)
   Slps <- quantile(out[,2], c((1-p)/2, (1+p)/2))
   out <- list(SlpCutUniv = Slps)
   out
}



## Function to identify stratum for each subject i (indentify three possible strata)
## INPUT:  w.function (whehter we want to sample based on intercept, slope or bivariate sampling), cutpoints (cutoff for define
##         extreme), Int (value of subject i estimated intercept), Slp (value of subject i estimated slope)
## OUTPUT: wheter a subject is in startum 1, 2 or 3
identify.stratum <- function(w.function, cutpoints, Int, Slp){
  if(w.function == "intercept"){
    stratum <- ifelse( Int <  cutpoints[1], 1, ifelse( Int >= cutpoints[2], 3, 2))}
  if(w.function == "slope"){
    stratum <- ifelse( Slp <  cutpoints[1], 1, ifelse( Slp >= cutpoints[2], 3, 2))}
  if(w.function=="bivar"){
    stratum <- ifelse(Int >= cutpoints[1] & Int < cutpoints[2] & Slp >=cutpoints[3] & 
                        Slp <cutpoints[4], 1, 2)}
  return(stratum)
}



## Function for outcome dependent sampling
## INPUT:  id.long (vector of ID in long format), stratum.long (stratum identifier in long format), 
##         SamplingStrategy (either dependent or independent ODS), NsPerStratum (vector indicating how many subject
##         we want to sample in each of the three strata)
## OUTPUT: List of three element (indicator vector of whether a subject has been sampled, Sampling probability
##         for each subject, sampling probability in each strata)
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


## Function for outcome dependent sampling of a single outcome (used to show how to re-use data)
## INPUT:  d (dataset), NsPerStratum (vector indicating how many subject
##         we want to sample in each of the three strata), Y1.cutoffODS (lower cutoff for ODS sampling),
##         Y1.cutoffBDS (lower cutoff for BDS sampling)
## OUTPUT: d (original dataset with additional columns of indicator variables showing whether each subject has been sampled
##         SampProbsODS (sampling probabilities for ODS), SampProbsBDS (sampling probabilities for BDS)
## Note:   This function only does slope sampling for Y1
SampledFunY1 <- function(d, NsPerStratumUniv, Y1.cutoffsODS, Y1.cutoffsBDS){

  d$SampledSlpODS  <- 0
  d$SampledSlpBDS  <- 0
  SampProbsODS     <- SampProbsBDS  <- list()
  # calculate subject specific intercept and slope for the first outcome
  
  #### ODS SAMPLING ####
  IntSlps    <- CalcSSIntSlp(Y = d$Y1, time = d$time, id = d$id)
  Slp        <- IntSlps[[2]]
  # identify stratum membership
  StratSlp   <- identify.stratum(w.function="slope", cutpoints = Y1.cutoffsODS$SlpCutUniv, Slp = Slp)
  # identify those sampled along with individual sampling probs and stratum sampling probs
  SampledSlp <- ods.sampling(id.long = d$id, stratum.long = StratSlp, SamplingStrategy = "IndepODS",
                             NsPerStratum = NsPerStratumUniv)
  d$SampledSlpODS         <- SampledSlp[[1]]
  SampProbsODS[[1]]       <- SampledSlp$SampProbs

  #### BDS SAMPLING ####
  mod          <- lme(Y1 ~ time*conf, random = ~ 1 + time | id, data = d,
                      control = lmeControl(maxIter = 1e9, opt = c("optim")))
  random.slope <- ranef(mod)[,2]
  StratSlpBDS  <- identify.stratum(w.function="slope", cutpoints = Y1.cutoffsBDS$SlpCutUniv, 
                                   Slp = random.slope)
    # each number of StratSlpBLUP represent wheter a subject was in stratum 1, 2, and 3
    # replicate the number per number of time a subject is observed
    # compute number of observations per subject
  nobs         <- d %>% group_by(id) %>% summarise(nobs = n()) %>% pull(nobs)
  StratSlpBDS  <- rep(StratSlpBDS, times = nobs)
  SampledSlp   <- ods.sampling(id.long = d$id, stratum.long = StratSlpBDS, SamplingStrategy = "IndepODS",
                               NsPerStratum = NsPerStratumUniv)
  d$SampledSlpBDS   <- SampledSlp$Sampled
  SampProbsBDS[[1]] <- SampledSlp$SampProbs

  #### Output ####
  out <- list(dat = d, SampProbsODS = SampProbsODS, SampProbsBDS = SampProbsBDS)
  return(out)
}

## Function for outcome dependent sampling when we sample each outcome independently
## INPUT:  d (dataset), Nsample (number of subjects we want to sample), NsPerStratum (vector indicating how many subject
##         we want to sample in each of the three strata), Y1.cutoffODS (Y_1 cutoff for ODS sampling), Y2.cutoffODS
##         (Y_2 cutoff for ODS sampling), Y1.cutoffBDS (Y_1 cutoff for BDS sampling), Y2.cutoffBDS 
##         (Y_2 cutoff for BDS sampling)
## OUTPUT: d (original dataset with additional columns of indicator variables showing whether each subject has been sampled
##         SampProbsODS (sampling probabilities for ODS), SampProbsBDS (sampling probabilities for BDS)
## Note:   This function only does slope sampling
SampledFun <- function(d, Nsample, NsPerStratumUniv, Y1.cutoffsODS, Y2.cutoffsODS, 
                       Y1.cutoffsBDS, Y2.cutoffsBDS){

  Nsample.sim           <- 2*sum(NsPerStratumUniv)

  SampledRan            <- random.sampling(id.long = d$id, n = Nsample.sim)
  d$SampledRan          <- SampledRan

   # divide N/2 to be sampled based on Y_1 and N/2 to be sampled based on Y_2
   nobs             <- d %>% group_by(id) %>% summarise(nobs = n()) %>% pull(nobs)
   d$var.sampled    <- rep(sample(c("Y1", "Y2"), size = length(unique(d$id)), replace = TRUE,
                                  prob = c(0.5, 0.5)), times = nobs)

   #### slope sampling ####
   SampProbsODS     <- SampProbsBDS <- list()

   # calculate subject specific intercept and slope for the first outcome
   # get all of Y_1 that do not have missing data?
   dY1              <- subset(d, d$var.sampled == "Y1")

   #### ODS SAMPLING ####
   IntSlps <- CalcSSIntSlp(Y = dY1$Y1, time = dY1$time, id = dY1$id)
   Slp <- IntSlps[[2]]
   # identify stratum membership
   StratSlp <- identify.stratum(w.function="slope", cutpoints = Y1.cutoffsODS$SlpCutUniv, Slp = Slp)
   # identify those sampled along with individual sampling probs and stratum sampling probs
   SampledSlp <- ods.sampling(id.long = dY1$id, stratum.long = StratSlp, SamplingStrategy = "IndepODS",
                              NsPerStratum = NsPerStratumUniv)
   dY1$SampledSlpODS       <- SampledSlp[[1]]
   SampProbsODS[[1]]       <- SampledSlp$SampProbs

   #### BLUP SAMPLING ####
   mod <- lme(Y1 ~ time*conf, random = ~ 1 + time | id, data = dY1,
              control = lmeControl(maxIter = 1e9, opt = c("optim")))
   random.slope <- ranef(mod)[,2]
   StratSlpBDS <- identify.stratum(w.function="slope", cutpoints = Y1.cutoffsBDS$SlpCutUniv, 
                                   Slp = random.slope)
   # each number of StratSlpBDS represent wheter a subject was in stratum 1, 2, and 3
   # replicate the number per number of time a subject is observed
   # compute number of observations per subject
   nobs         <- dY1 %>% group_by(id) %>% summarise(nobs = n()) %>% pull(nobs)
   StratSlpBDS  <- rep(StratSlpBDS, times = nobs)
   SampledSlp   <- ods.sampling(id.long = dY1$id, stratum.long = StratSlpBDS, 
                                SamplingStrategy = "IndepODS",
                                NsPerStratum = NsPerStratumUniv)
   dY1$SampledSlpBDS  <- SampledSlp$Sampled
   SampProbsBDS[[1]]  <- SampledSlp$SampProbs

   # repeat the same for Y2
   #### ODS SAMPLING ####
   dY2            <- subset(d, d$var.sampled == "Y2")
   IntSlps        <- CalcSSIntSlp(Y = dY2$Y2, time = dY2$time, id = dY2$id)
   Slp            <- IntSlps[[2]]
   # identify stratum membership
   StratSlp <- identify.stratum(w.function="slope", cutpoints = Y2.cutoffsODS$SlpCutUniv, Slp = Slp)
   # identify those sampled along with individual sampling probs and stratum sampling probs
   SampledSlp <- ods.sampling(id.long = dY2$id, stratum.long = StratSlp, 
                              SamplingStrategy = "IndepODS",
                              NsPerStratum = NsPerStratumUniv)
   dY2$SampledSlpODS     <- SampledSlp[[1]]
   SampProbsODS[[2]]        <- SampledSlp$SampProbs

   #### BLUP SAMPLING ####
   mod <- lme(Y2 ~ time*conf, random = ~ 1 + time | id, data = dY2,
              control = lmeControl(maxIter = 1e9, opt = c("optim")))
   random.slope <- ranef(mod)[,2]
   StratSlpBDS <- identify.stratum(w.function="slope", cutpoints = Y2.cutoffsBDS$SlpCutUniv, 
                                   Slp = random.slope)
   # each number of StratSlpBDS represent wheter a subject was in stratum 1, 2, and 3
   # replicate the number per number of time a subject is observed
   # compute number of observations per subject
   nobs         <- dY2 %>% group_by(id) %>% summarise(nobs = n()) %>% pull(nobs)
   StratSlpBDS  <- rep(StratSlpBDS, times = nobs)
   SampledSlp   <- ods.sampling(id.long = dY2$id, stratum.long = StratSlpBDS, 
                              SamplingStrategy = "IndepODS",
                              NsPerStratum = NsPerStratumUniv)
   dY2$SampledSlpBDS  <- SampledSlp[[1]]
   SampProbsBDS[[2]]  <- SampledSlp$SampProbs

   # put the data back together and sort by ID
   d            <- bind_rows(dY1, dY2) %>% arrange(id)

   #### output ####
   out <- list(dat = d, SampProbsODS = SampProbsODS, SampProbsBDS = SampProbsBDS)
   return(out)
}



## Inverse logit function
## INPUT:  vector
## OUTPUT: inverse logit transformation of the input
expit <- function(x) {
  exp(x)/(1+exp(x))
}
