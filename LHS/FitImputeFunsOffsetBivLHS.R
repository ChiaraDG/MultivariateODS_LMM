#########################################################################################################################
####### Functions Offset Imputation with Bivariate Data #################################################################
#########################################################################################################################


## Function to get individual intercept and slope 
# Important to note that data must contain variables Y and time
# INPUT: dataset
# OUTPUT: linear regression function
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


## Function for simple random sampling
# INPUT: id in long format, number of subject to sample
# OUTPUT: id of sampled subjects
random.sampling <- function(id.long, n){
  s <- sample(unique(id.long), n)
  Sampled <- as.integer(id.long %in% s)
  return(Sampled)
}


## find the inverse of variance matrix for each subject:
# INPUT: time variable, the matrix of variance/covariance for the random effects and the variance for the error term
InvVarMat <- function(time, refmat, se2, Y1, Y2){
  # if none of the values of Y are missing
  if(suppressWarnings(!is.na(all(Y1)) & !is.na(all(Y2)))){
    Zi      <- cbind(1, time)
    Z       <- as.matrix(bdiag(Zi, Zi))
    se1     <- diag(se2[1], nrow(Zi))
    se2     <- diag(se2[2], nrow(Zi))
    out <- solve(Z%*%refmat%*%t(Z) + as.matrix(bdiag(se1, se2)))
  }
  # if Y1 is missing completely
  if(suppressWarnings(is.na(all(Y1)))){
    Z     <- cbind(0, 0, 1, time)
    se2   <- diag(se2[2], nrow(Z))
    out   <- solve(Z%*%refmat%*%t(Z) + se2)
  }
  # if only Y2 is completely missing
  if(suppressWarnings(is.na(all(Y2)))){
    Z     <- cbind(1, time, 0, 0)
    se2   <- diag(se2[1], nrow(Z))
    out   <- solve(Z%*%refmat%*%t(Z) + se2)
  }
  return(out)
}

## find the offset
#INPUT: dataframe, subject specific weights, model from previous iteration
#OUTPUT: dataframe with ID and offset to be used during the imputation
offsetSearchBiv <- function(d, weights, betas, var.names = NULL){

  # d need to be a dataframe. Convert to dataframe
  if(!is.data.frame(d)) {d <- as.data.frame(d)}
  # weights need to be a list. If not stop and return an error
  if(!is.list(weights)) {stop("Error: weights need to be a list")}

  # estimate the offset
  d.y1  <- subset(d, D1 == 1)
  d.y2  <- subset(d, D2 == 1)
  
  X0.y1 <- cbind(1, 0, d.y1$time, d.y1$cigs0.10, d.y1$pack0.20, d.y1$female, d.y1$age.10,
                 d.y1$bmi0.5, d.y1$bmi.change.5, d.y1$time*d.y1$female, d.y1$time*d.y1$age.10,  
                 0*d.y1$time, d.y1$site2, d.y1$site3, d.y1$site4, d.y1$site5, d.y1$site6, d.y1$site7,
                 d.y1$site8, d.y1$site9, d.y1$site10)
  X0.y2 <- cbind(1, 0, d.y2$time, d.y2$cigs0.10, d.y2$pack0.20, d.y2$female, d.y2$age.10,
                 d.y2$bmi0.5, d.y2$bmi.change.5, d.y2$time*d.y2$female, d.y2$time*d.y2$age.10,  
                 0*d.y2$time, d.y2$site2, d.y2$site3, d.y2$site4, d.y2$site5, d.y2$site6, d.y2$site7,
                 d.y2$site8, d.y2$site9, d.y2$site10)
  d$m0  <- as.matrix(bdiag(X0.y1, X0.y2))%*%betas
  
  X1.y1 <- cbind(1, 1, d.y1$time, d.y1$cigs0.10, d.y1$pack0.20, d.y1$female, d.y1$age.10,
                 d.y1$bmi0.5, d.y1$bmi.change.5, d.y1$time*d.y1$female, d.y1$time*d.y1$age.10,  
                 1*d.y1$time, d.y1$site2, d.y1$site3, d.y1$site4, d.y1$site5, d.y1$site6, d.y1$site7,
                 d.y1$site8, d.y1$site9, d.y1$site10)
  X1.y2 <- cbind(1, 1, d.y2$time, d.y2$cigs0.10, d.y2$pack0.20, d.y2$female, d.y2$age.10,
                 d.y2$bmi0.5, d.y2$bmi.change.5, d.y2$time*d.y2$female, d.y2$time*d.y2$age.10,  
                 1*d.y2$time, d.y2$site2, d.y2$site3, d.y2$site4, d.y2$site5, d.y2$site6, d.y2$site7,
                 d.y2$site8, d.y2$site9, d.y2$site10)
  d$m1  <- as.matrix(bdiag(X1.y1, X1.y2))%*%betas
  
  # compute (a) and (b) terms in the offset
  split_dat <- split(d, d$id)
  # get the y and the mu for each person
  ys <- lapply(split_dat, '[[', "value")
  m0 <- lapply(split_dat, '[[', "m0")
  m1 <- lapply(split_dat, '[[', "m1")
  # set up the terms for the offset
  a.term <- Map(function(y, w, m1, m0){t(y)%*%w%*%(m1 - m0)}, ys, weights, m1, m0)
  b.term <- Map(function(w, m1, m0){(-0.5)*(t(m1)%*%w%*%m1 -  t(m0)%*%w%*%m0)}, weights, m1, m0)
  offset <- Map("+", a.term, b.term)
  out <- data.frame(id = unique(d$id), offset = unlist(offset))
  return(out)
}

DirectImputationOffsetBiv <- function(d, weights, betas, sampling){
  
  attributes(d$Y1) <- NULL
  attributes(d$Y2) <- NULL
  
  if(sampling == "random"){
    # take those selected with random sampling
    dat.rs <- d
    dat.rs$grp <- ifelse(d$SampledRan == 1, d$grp, NA)  
    dat.rs$timegrp <- ifelse(d$SampledRan == 1, d$timegrp, NA)  
    # transform to long format to help with the lme fitting
    dat.long <- dat.rs %>% 
      gather(., key = "y_name", value = "value", Y1, Y2) %>%
      mutate(D1 = as.integer(y_name == "Y1"), D2 = as.integer(y_name == "Y2"))
    # compute the offset
    dat.summ <- offsetSearchBiv(dat.long, weights = weights, betas = betas)
    # add other variables of interest
    dat.wide <- dat.rs %>% 
      dplyr::select(id, grp, cigs0.10, pack0.20, female, age.10, bmi0.5, site2, 
                    site3, site4, site5, site6, site7, site8, site9, site10) %>%
      distinct %>% 
      inner_join(., dat.summ, by = "id")
  }
  
  if(sampling == "SlpODS"){
    # take those selected with random sampling
    dat.slpods <- d
    dat.slpods$grp <- ifelse(d$SampledSlpODS == 1, d$grp, NA)  
    dat.slpods$timegrp <- ifelse(d$SampledSlpODS == 1, d$timegrp, NA)  
    # transform to long format to help with the lme fitting
    dat.long <- dat.slpods %>% 
      gather(., key = "y_name", value = "value", Y1, Y2) %>%
      mutate(D1 = as.integer(y_name == "Y1"), D2 = as.integer(y_name == "Y2"))
    # compute the offset
    dat.summ <- offsetSearchBiv(dat.long, weights = weights, betas = betas)
    # add other variables of interest
    dat.wide <- dat.slpods %>% 
      dplyr::select(id, grp, cigs0.10, pack0.20, female, age.10, bmi0.5, site2, 
                    site3, site4, site5, site6, site7, site8, site9, site10) %>% 
      distinct %>% 
      inner_join(., dat.summ, by = "id")
  }
  
  if(sampling == "SlpBLUP"){
    # take those selected with random sampling
    dat.slpblup <- d
    dat.slpblup$grp <- ifelse(d$SampledSlpBLUP == 1, d$grp, NA)  
    dat.slpblup$timegrp <- ifelse(d$SampledSlpBLUP == 1, d$timegrp, NA)  
    # transform to long format to help with the lme fitting
    dat.long <- dat.slpblup %>% 
      gather(., key = "y_name", value = "value", Y1, Y2) %>%
      mutate(D1 = as.integer(y_name == "Y1"), D2 = as.integer(y_name == "Y2"))
    # compute the offset
    dat.summ <- offsetSearchBiv(dat.long, weights = weights, betas = betas)
    # add other variables of interest
    dat.wide <- dat.slpblup %>% 
      dplyr::select(id, grp, cigs0.10, pack0.20, female, age.10, bmi0.5, site2, 
                    site3, site4, site5, site6, site7, site8, site9, site10) %>% 
      distinct %>% 
      inner_join(., dat.summ, by = "id")
  }
  
  # fit imputation model
  mod.imp <- glm(grp ~  cigs0.10 + pack0.20 + female + age.10 + bmi0.5 + site2 + site3 + site4 + site5 + 
                   site6 + site7 + site8 + site9 + site10 + offset(offset), 
                 data = dat.wide, family = "binomial")
  
  # imputation
  SampleParams <- rmvnorm(1, mod.imp$coef, summary(mod.imp)$cov.unscaled)
  pred.imp <- expit(cbind(1, dat.wide$cigs0.10, dat.wide$pack0.20, dat.wide$female, 
                          dat.wide$age.10, dat.wide$bmi0.5, dat.wide$site2, 
                          dat.wide$site3, dat.wide$site4, dat.wide$site5, dat.wide$site6, 
                          dat.wide$site7, dat.wide$site8, dat.wide$site9, dat.wide$site10)%*%t(SampleParams) + 
                      dat.wide$offset)
  imp.grp  <- rbinom(nrow(dat.wide), 1, pred.imp)
  # dataset with imputed group
  Completed        <- dat.wide
  Completed$grpall <- ifelse(is.na(Completed$grp), imp.grp, Completed$grp)
  Completed        <- Completed %>% dplyr::select(id, grpall)
  # merge with information on y and x: this reshape back to long format
  long.dat <- merge(dat.long, Completed, by = "id")
  # add interaction time*grp and sort by id
  long.dat$timegrpall <- long.dat$time*long.dat$grpall 
  long.dat$timefemale <- long.dat$time*long.dat$female
  long.dat$timeage    <- long.dat$time*long.dat$age.10
  long.dat            <- long.dat[order(long.dat$id, long.dat$y_name),]
  
  ## do an lme on the fully imputed data
  mod.imp <- lme(value ~ -1 + D1 + D1:grpall + D1:time + D1:cigs0.10 + D1:pack0.20 + 
                   D1:female + D1:age.10 + D1:bmi0.5 +  D1:bmi.change.5 + 
                   D1:timefemale + D1:timeage + D1:timegrpall + D1:site2 + D1:site3 + D1:site4 + D1:site5 +
                   D1:site6 + D1:site7 + D1:site8 + D1:site9 + D1:site10 +
                   D2 + D2:grpall + D2:time + D2:cigs0.10 + D2:pack0.20 + 
                   D2:female + D2:age.10 + D2:bmi0.5 + D2:bmi.change.5 + 
                   D2:timefemale + D2:timeage + D2:timegrpall + D2:site2 + D2:site3 + D2:site4 + D2:site5 +
                   D2:site6 + D2:site7 + D2:site8 + D2:site9 + D2:site10, data = long.dat, 
                 random = ~ 0 + D1 + D1:time + D2 + D2:time | id, 
                 weights = varIdent(form = ~1 | y_name), 
                 control = lmeControl(maxIter = 1e9, opt = c("optim")))
  
  # get the variance and covariance matrix of random effects
  VarCov.Mat.Random <- getVarCov(mod.imp)[c("D1", "D1:time", "D2", "time:D2"),
                                          c("D1", "D1:time", "D2", "time:D2")]
  # Take the sigma for the residuals, one for each of the subjects.
  SE   <- (summary(mod.imp)$sigma)
  w.Y2 <- coef(mod.imp$modelStruct$varStruct, uncons = FALSE)
  SE2  <- c(SE, SE*w.Y2)^2
  
  
  # get a set of new betas to re-compute offset
  b.y1   <- fixef(mod.imp)[c("D1", "D1:grpall", "D1:time", "D1:cigs0.10", "D1:pack0.20", 
                             "D1:female", "D1:age.10", "D1:bmi0.5", "D1:bmi.change.5", 
                             "D1:timefemale", "D1:timeage", "D1:timegrpall", "D1:site2", "D1:site3", 
                             "D1:site4", "D1:site5", "D1:site6", "D1:site7", "D1:site8", "D1:site9", "D1:site10")]
  b.y2   <- fixef(mod.imp)[c("D2", "grpall:D2", "time:D2", "cigs0.10:D2", "pack0.20:D2", 
                             "female:D2", "age.10:D2", "bmi0.5:D2", "bmi.change.5:D2", 
                             "timefemale:D2", "timeage:D2", "timegrpall:D2", "site2:D2", "site3:D2", "site4:D2",
                             "site5:D2", "site6:D2", "site7:D2", "site8:D2", "site9:D2", "site10:D2")]
  betas  <- c(b.y1, b.y2)
  
  sigmas <- vcov(mod.imp)[c("D1", "D1:grpall", "D1:time", "D1:cigs0.10", "D1:pack0.20", 
                            "D1:female", "D1:age.10", "D1:bmi0.5", "D1:bmi.change.5", 
                            "D1:timefemale", "D1:timeage", "D1:timegrpall", "D1:site2", "D1:site3", 
                            "D1:site4", "D1:site5", "D1:site6", "D1:site7", "D1:site8", "D1:site9", "D1:site10", 
                            "D2", "grpall:D2", "time:D2", "cigs0.10:D2", "pack0.20:D2", 
                            "female:D2", "age.10:D2", "bmi0.5:D2", "bmi.change.5:D2", 
                            "timefemale:D2", "timeage:D2", "timegrpall:D2", "site2:D2", "site3:D2", "site4:D2",
                            "site5:D2", "site6:D2", "site7:D2", "site8:D2", "site9:D2", "site10:D2"), 
                          c("D1", "D1:grpall", "D1:time", "D1:cigs0.10", "D1:pack0.20", 
                            "D1:female", "D1:age.10", "D1:bmi0.5", "D1:bmi.change.5", 
                            "D1:timefemale", "D1:timeage", "D1:timegrpall", "D1:site2", "D1:site3", 
                            "D1:site4", "D1:site5", "D1:site6", "D1:site7", "D1:site8", "D1:site9", "D1:site10",
                            "D2", "grpall:D2", "time:D2", "cigs0.10:D2", "pack0.20:D2", 
                            "female:D2", "age.10:D2", "bmi0.5:D2", "bmi.change.5:D2", 
                            "timefemale:D2", "timeage:D2", "timegrpall:D2", "site2:D2", "site3:D2", "site4:D2",
                            "site5:D2", "site6:D2", "site7:D2", "site8:D2", "site9:D2", "site10:D2")]
  # get a set of new betas to re-compute offset
  b.new <- as.vector(rmvnorm(1, as.vector(betas), as.matrix(sigmas)))
  
  # output
  out <- list(new.VarCov = VarCov.Mat.Random, new.se2 = SE2, est.new = b.new, lme.imp = mod.imp)
  return(out)
}





ImputationBiv <- function(d, imp, burn, sampling){
  if(burn >= imp){stop("Burn needs to be less than the number of replicates")}
  
  if(sampling == "random"){
    # take those selected with random sampling
    dat <- d
    # set to NA both the group and the group*time interaction
    dat$grp     <- ifelse(d$SampledRan == 1, d$grp, NA) 
    dat$timegrp <- ifelse(d$SampledRan == 1, d$timegrp, NA) 
  }
  
  if(sampling == "SlpODS"){
    # take those selected with ODS sampling
    dat <- d
    # set to NA both the group and the group*time interaction
    dat$grp     <- ifelse(d$SampledSlpODS == 1, d$grp, NA) 
    dat$timegrp <- ifelse(d$SampledSlpODS == 1, d$timegrp, NA) 
  }
  
  if(sampling == "SlpBLUP"){
    # take those selected with BLUP sampling
    dat <- d
    # set to NA both the group and the group*time interaction
    dat$grp     <- ifelse(d$SampledSlpBLUP == 1, d$grp, NA) 
    dat$timegrp <- ifelse(d$SampledSlpBLUP == 1, d$timegrp, NA) 
  }
  
  # transform data in long format to be able to run the algorithm
  dat.long <- dat %>% 
    gather(., key = "y_name", value = "value", Y1, Y2) %>%
    mutate(D1 = as.integer(y_name == "Y1"), D2 = as.integer(y_name == "Y2")) %>% 
    arrange(id)
  lme.mod <- lme(value ~ -1 + D1 + D1:grp + D1:time + D1:cigs0.10 + D1:pack0.20 + 
                   D1:female + D1:age.10 + D1:bmi0.5 +  D1:bmi.change.5 + 
                   D1:timefemale + D1:timeage + D1:timegrp + D1:site2 + D1:site3 + D1:site4 + D1:site5 +
                   D1:site6 + D1:site7 + D1:site8 + D1:site9 + D1:site10 +
                   D2 + D2:grp + D2:time + D2:cigs0.10 + D2:pack0.20 + 
                   D2:female + D2:age.10 + D2:bmi0.5 + D2:bmi.change.5 + 
                   D2:timefemale + D2:timeage + D2:timegrp + D2:site2 + D2:site3 + D2:site4 + D2:site5 +
                   D2:site6 + D2:site7 + D2:site8 + D2:site9 + D2:site10,
                 data = dat.long, random = ~ 0 + D1 + D1:time + D2 + D2:time| id, 
                 weights = varIdent(form = ~1 | y_name), na.action = na.omit,
                 control = lmeControl(maxIter = 1e9, opt = c("optim")))
  # get the variance and covariance matrix of random effects
  VarCov.Mat.Random <- getVarCov(lme.mod)[c("D1", "D1:time", "D2", "time:D2"),
                                          c("D1", "D1:time", "D2", "time:D2")]
  # Take the sigma for the residuals, one for each of the subjects.
  SE   <- (summary(lme.mod)$sigma)
  w.Y2 <- coef(lme.mod$modelStruct$varStruct, uncons = FALSE)
  SE2  <- c(SE, SE*w.Y2)^2
  # get the weight matrix
  # If there are not missing Y we could just use this code
  w.new <- tapply(d$time, d$id, InvVarMat, refmat = VarCov.Mat.Random, se2 = SE2, Y1 = d$Y1, Y2 = d$Y2)
  # get a set of new betas to re-compute offset
  b.new.y1 <- fixef(lme.mod)[c("D1", "D1:grp", "D1:time", "D1:cigs0.10", "D1:pack0.20", 
                               "D1:female", "D1:age.10", "D1:bmi0.5", "D1:bmi.change.5", 
                               "D1:timefemale", "D1:timeage", "D1:timegrp", 
                               "D1:site2", "D1:site3", 
                               "D1:site4", "D1:site5", "D1:site6", "D1:site7", "D1:site8", "D1:site9", "D1:site10")]
  b.new.y2 <- fixef(lme.mod)[c("D2", "grp:D2", "time:D2", "cigs0.10:D2", "pack0.20:D2", 
                               "female:D2", "age.10:D2", "bmi0.5:D2", "bmi.change.5:D2", 
                               "timefemale:D2", "timeage:D2", "timegrp:D2", "site2:D2", "site3:D2", "site4:D2",
                               "site5:D2", "site6:D2", "site7:D2", "site8:D2", "site9:D2", "site10:D2")]
  b.new    <- c(b.new.y1, b.new.y2)
  
  # start imputing
  i <- 1
  
  Ests.Imp <- Covs.Imp <- list()
  while(i <= imp){
    print(i)
    imputation <- DirectImputationOffsetBiv(d, weights = w.new, betas = b.new, sampling = sampling)
    # If there are not missing Y we could just use this code
    w.new  <- tapply(d$time, d$id, InvVarMat, refmat = imputation$new.VarCov, se2 = imputation$new.se2, 
                     Y1 = d$Y1, Y2 = d$Y2)
    # save the new betas
    b.new <- imputation$est.new
    # store the results from imputation
    Ests.Imp[[i]] <- fixef(imputation$lme.imp)
    Covs.Imp[[i]] <- vcov(imputation$lme.imp)
    i <- i + 1
  }
  # burn the first few models
  Ests.Imp <- Ests.Imp[-c(1:burn)]
  Covs.Imp <- Covs.Imp[-c(1:burn)]
  # put results together
  coeff.res <- MIcombine(Ests.Imp, Covs.Imp)
  out <- list(coefficients = coeff.res$coefficients[c("D1", "D1:grpall", "D1:time", "D1:cigs0.10", "D1:pack0.20", 
                                                      "D1:female", "D1:age.10", "D1:bmi0.5", "D1:bmi.change.5", 
                                                      "D1:timefemale", "D1:timeage", "D1:timegrpall", "D1:site2", "D1:site3", 
                                                      "D1:site4", "D1:site5", "D1:site6", "D1:site7", "D1:site8", "D1:site9", "D1:site10",
                                                      "D2", "grpall:D2", "time:D2", "cigs0.10:D2", "pack0.20:D2", 
                                                      "female:D2", "age.10:D2", "bmi0.5:D2", "bmi.change.5:D2", 
                                                      "timefemale:D2", "timeage:D2", "timegrpall:D2", "site2:D2", "site3:D2", "site4:D2",
                                                      "site5:D2", "site6:D2", "site7:D2", "site8:D2", "site9:D2", "site10:D2")],
              
              variance = diag(coeff.res$variance[c("D1", "D1:grpall", "D1:time", "D1:cigs0.10", "D1:pack0.20", 
                                                     "D1:female", "D1:age.10", "D1:bmi0.5", "D1:bmi.change.5", 
                                                     "D1:timefemale", "D1:timeage", "D1:timegrpall", "D1:site2", "D1:site3", 
                                                   "D1:site4", "D1:site5", "D1:site6", "D1:site7", "D1:site8", "D1:site9", "D1:site10",
                                                     "D2", "grpall:D2", "time:D2", "cigs0.10:D2", "pack0.20:D2", 
                                                     "female:D2", "age.10:D2", "bmi0.5:D2", "bmi.change.5:D2", 
                                                     "timefemale:D2", "timeage:D2", "timegrpall:D2", "site2:D2", "site3:D2", "site4:D2",
                                                   "site5:D2", "site6:D2", "site7:D2", "site8:D2", "site9:D2", "site10:D2"), 
                                                   c("D1", "D1:grpall", "D1:time", "D1:cigs0.10", "D1:pack0.20", 
                                                     "D1:female", "D1:age.10", "D1:bmi0.5", "D1:bmi.change.5", 
                                                     "D1:timefemale", "D1:timeage", "D1:timegrpall", "D1:site2", "D1:site3", 
                                                     "D1:site4", "D1:site5", "D1:site6", "D1:site7", "D1:site8", "D1:site9", "D1:site10",
                                                     "D2", "grpall:D2", "time:D2", "cigs0.10:D2", "pack0.20:D2", 
                                                     "female:D2", "age.10:D2", "bmi0.5:D2", "bmi.change.5:D2", 
                                                     "timefemale:D2", "timeage:D2", "timegrpall:D2", "site2:D2", "site3:D2", "site4:D2",
                                                     "site5:D2", "site6:D2", "site7:D2", "site8:D2", "site9:D2", "site10:D2")]),
              covariance = c(coeff.res$variance["D1:time", "D1:timegrpall"], coeff.res$variance["D1:time", "D1:timefemale"],
                             coeff.res$variance["D1:time", "D1:timeage"], coeff.res$variance["time:D2", "timegrpall:D2"], 
                             coeff.res$variance["time:D2", "timefemale:D2"],
                             coeff.res$variance["time:D2", "timeage:D2"]))
  #Ests.Imp <- do.call(rbind, Ests.Imp)
  #out <- list(coefficients = Ests.Imp)
  return(out)
}