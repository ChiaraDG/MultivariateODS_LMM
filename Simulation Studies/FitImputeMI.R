#########################################################################################################################
####### Functions Offset Imputation with Bivariate Data #################################################################
#########################################################################################################################

## Find the inverse of variance matrix for each subject:
## Formula: ZDZ' + \sigma^2I
## INPUT:   the matrix of random effects Z
##          the matrix of variance/covariance for the random effects 
##          and the standard deviation for the error term (assume that this is a diagonal matrix)
## OUTPUT:  inverse of the variance matrix for one subject
## Notes:   For this function to work data need one line per id per time per outcome.

vinv.i <- function(Zi, D, sigma.e){
  # stopping rule to check inputs are given correctly
  if(!is.matrix(D)){stop("D needs to be a matrix")}
  
  # from the vector sigma.e create the matrix Sigma
  nZi      <- nrow(Zi)
  nsigma.e <- length(sigma.e)
  Sigma    <- diag(rep(sigma.e^2, each = nZi/nsigma.e))
  Zi       <- as.matrix(Zi)
  
  # the Vi matrix create 
  vi <- Zi%*%D%*%t(Zi) + Sigma
  
  # return the inverse of Vi
  out <- solve(vi)
  return(out)
}

## Find the offset for each person i to be used for the imputation part
## Formula: Y'V^{-1}*(mu1 - mu0) - (0.5)*(mu1'V^{-1}mu1 - mu0'V^{-1}mu0)
## INPUT:   the vector yi of outcomes for subject i, a dataframe with column names of fixed effects X.i matrix, 
##          the inverse of the variance matrix resulting from vinv.i,
##          the coefficients beta computed in a previous iteration (or an initial value),
##          grp indicates how the expensive variable is called in the dataset we are using for the analysis
##          timegrp indicates how the time by grp is called in the dataset we are using for the analysis
##          time    indicates how the time by grp is called in the dataset we are using for the analysis
## OUPUT:   offset for subject i
## Notes:   X.i is the design matrix so the first column needs to be
##          a column of 1s he coefficients' betas order here needs to be the same of the order X.i in the model

offset.i <- function(subjectData, betas, D, sigma.e, grp, timegrp, time){
  
  # stopping rules to check the outputs are given correctly
  if(!is.character(grp)){stop("grp needs to be a character vector")}
  if(!is.character(timegrp)){stop("timegrp needs to be a character vector")}
  if(!is.character(time)){stop("time needs to be a character vector")}
  
  Xi       <- subjectData[["Xi"]]
  Zi       <- subjectData[["Zi"]]
  yi       <- subjectData[["yi"]]
  
  vinv     <- vinv.i(Zi = Zi, D, sigma.e)
  
  # first we need to compute mu0 and mu1
  X0.i                                <- Xi
  X0.i[, colnames(X0.i) %in% grp]     <- 0
  X0.i[, colnames(X0.i) %in% timegrp] <- 0
  mu0.i                               <- as.matrix(X0.i)%*%betas
  
  X1.i                                <- Xi
  X1.i[, colnames(X1.i) == grp[1]]    <- c(rep(1, nrow(Xi)/2), rep(0, nrow(Xi)/2))
  X1.i[, colnames(X1.i) == grp[2]]    <- c(rep(0, nrow(Xi)/2), rep(1, nrow(Xi)/2))
  X1.i[, colnames(X1.i) %in% timegrp] <- X1.i[, colnames(X1.i) %in% time]
  mu1.i                               <- as.matrix(X1.i)%*%betas
  
  # now we can compute the offset
  offset <- t(yi)%*%vinv%*%(mu1.i - mu0.i) - (0.5)*(t(mu1.i)%*%vinv%*%mu1.i - t(mu0.i)%*%vinv%*%mu0.i)
  
  as.numeric(offset)
}


## Code to run one iteration of the imputation
## INPUT:   formula.fixed (formula for the fixed effects in our model of interest), formula.random
##          (formula for the random effects in our model of interest). formula.imp 
##          (formula for the logistic model used for impuation), data (dataset), D (matrix of 
##          variance and covariance matrix of random effects), sigma.e (std of error term),
##          betas (vector of parameters beta computed in a previous iteration or an initial value),
##          id, grp, timegrp, time indicates how the id, expensive exposure, time by expensive exposure,
##          and time variables are called in the dataset we are using for analysis. yname indicates
##          the name of a variable indicating whether the outcome in each line is Y_1 or Y_2
## OUPUT:   new.VarCov (matrix of estimated variance and covariance of random effects) 
##          new.se (estimated std of the error terms), est.new (coefficients to be used in the next step), 
##          lme.imp (the output of nlme model used to compute the final estimates with mitool)
## Notes:   data needs to be in the format used for acml. Each row is a single observation 
##          (i.e.if subject i has 3 Y_1 and 3 Y_2, we will have 6 rows)

DirectImputationOffsetBiv <- function(formula.fixed, formula.random, formula.imp, 
                                      data, D, sigma.e, betas, id, grp, timegrp, time, yname){
  
  ## get the terms of the model: outcome, matrix of random effects, matrix of fixed effects,
  fixed.f  <- model.frame(formula.fixed, data, na.action = NULL)
  y        <- model.response(fixed.f, 'numeric')
  x        <- model.matrix(formula.fixed, fixed.f)
  randef   <- data[, all.vars(formula.random)]
  randef.f <- randef[, 1:(ncol(randef) - 1)]
  all.id   <- data$id
  z        <- as.matrix(randef.f)
  
  ## Compute the offset for each subject i
  SubjectDat  <- CreateSubjectData2(id = all.id, y = y, x = x, z = z)
  offset      <- lapply(SubjectDat, offset.i, D = D, sigma.e = sigma.e, betas = betas, 
                        grp = grp, timegrp = timegrp, time = time)
  offset.dat  <- data.frame(id = unique(all.id), offset = unlist(offset))
  # get the data in the wide format and merge with offset data
  dat.wide    <- data[,c("id",  all.vars(formula.imp))] %>% distinct %>% 
    inner_join(.,  offset.dat, by = "id")
  # get the design matrix for the imputation model
  imp.f  <- model.frame(formula.imp, dat.wide, na.action = NULL)
  x.imp  <- model.matrix(formula.imp, imp.f)
  
  ## Fit imputation model
  mod.imp <- glm(formula.imp, offset = offset, data = dat.wide, family = "binomial")
  # imputation
  SampleParams <- rmvnorm(1, mod.imp$coef, summary(mod.imp)$cov.unscaled)
  pred.imp     <- expit(cbind(x.imp)%*%t(SampleParams) + dat.wide$offset)
  imp.grp      <- rbinom(nrow(dat.wide), 1, pred.imp)
  # dataset with imputed group
  Completed        <- dat.wide
  Completed$grpall <- ifelse(is.na(Completed$grp), imp.grp, Completed$grp)
  Completed        <- Completed %>% dplyr::select(id, grpall)
  # merge with information on y and x: this reshape back to long format
  long.dat           <- merge(data, Completed, by = "id")
  long.dat[,timegrp] <- NULL
  long.dat[,grp]     <- NULL
  # make sure the format is correct for nlme
  n.var              <- ncol(x)/2
  groups             <- long.dat$grpall*cbind(x[,1], x[,(n.var + 1)])
  colnames(groups)   <- grp
  # add interaction time*grp and sort by id
  timegrpall          <- long.dat[, names(long.dat) %in% time]*long.dat$grpall
  colnames(timegrpall)<- timegrp
  long.dat            <- cbind(long.dat, groups, timegrpall)
  
  ## do an lme on the fully imputed data 
  mod.imp <- lme(formula.fixed, data = long.dat, 
                 random = formula.random, weights = varIdent(form = ~1 | y_name), 
                 control = lmeControl(maxIter = 1e9, opt = c("nlminb"), returnObject = TRUE))
  
  # get the variance and covariance matrix of random effects
  VarCov.Mat.Random <- getVarCov(mod.imp)
  # Take the sigma for the residuals, one for each of the subjects.
  SE   <- (summary(mod.imp)$sigma)
  w.Y2 <- coef(mod.imp$modelStruct$varStruct, uncons = FALSE)
  SE  <- c(SE, SE*w.Y2)
  # get a set of new betas to re-compute offset
  betas   <- fixef(mod.imp)
  sigmas  <- vcov(mod.imp)
  b.new   <- as.vector(rmvnorm(1, as.vector(betas), as.matrix(sigmas)))
  
  # output
  out <- list(new.VarCov = VarCov.Mat.Random, new.se = SE, est.new = b.new,
              data = long.dat, b.sim = betas)
  return(out)
}



## Code to run the MI
## INPUT:   formula.fixed (formula for the fixed effects in our model of interest), formula.random
##          (formula for the random effects in our model of interest). formula.imp 
##          (formula for the logistic model used for imputation), data (dataset), n.imp (number of imputations),
##          n.burs (number of imputation to discard), id, grp, timegrp, time indicates how the id, 
##          expensive exposure, time by expensive exposure,
##          and time variables are called in the dataset we are using for analysis. yname indicates
##          the name of a variable indicating whether the outcome in each line is Y_1 or Y_2
## OUPUT:   coefficients (estimated coefficients), covariance (estimated variance of the fixed effects)
## Notes:   data needs to be in the format used for acml. Each row is a single observation 
##          (i.e.if subject i has 3 Y_1 and 3 Y_2, we will have 6 rows)

MI.Imputation <- function(formula.fixed, formula.random, formula.imp, data, iter, id,
                grp, timegrp, time, yname){
  
  names(data)[names(data) == id]     <- "id"
  data$y_name  <- data[,yname]
  
  ## fit the model with available data in order to get estimated parameter and start IIM
  mod.init   <- lme(formula.fixed, data = data, 
                    random = formula.random, 
                    weights = varIdent(form = ~1 | y_name), na.action = na.omit,
                    control = lmeControl(maxIter = 1e9, opt = c("optim"), returnObject = TRUE))
  # get the variance and covariance matrix of random effects
  VarCov.Mat <- getVarCov(mod.init)
  # Take the sigma for the residuals, one for each of the subjects.
  SE   <- (summary(mod.init)$sigma)
  w.Y2 <- coef(mod.init$modelStruct$varStruct, uncons = FALSE)
  SE   <- c(SE, SE*w.Y2)
  # get the betas to compute the offset
  b.new    <- fixef(mod.init)
  
  ###############################################
  
  # start imputing
  i <- 1
  
  while(i <= iter){
    print(i)
    imputation <- DirectImputationOffsetBiv(formula.fixed, formula.random, formula.imp, 
                                            data = data, D = VarCov.Mat, 
                                            sigma.e = SE, betas = b.new, 
                                            id = id, grp = grp, timegrp = timegrp, time = time, 
                                            yname = yname)
    
    # save the updated parameters needed to compute the offset
    b.new      <- imputation$est.new
    VarCov.Mat <- imputation$new.VarCov
    SE         <- imputation$new.se
    i <- i + 1
  }
  
  
  # get the dataset from the last impuation
  new.dat  <- imputation$data
  # fit the model and save the results
  mod.final   <- lme(formula.fixed, data = new.dat, 
                    random = formula.random, 
                    weights = varIdent(form = ~1 | y_name), na.action = na.omit,
                    control = lmeControl(maxIter = 1e9, opt = c("optim"), returnObject = TRUE))
  # get the betas and variance anc covariance matrix
  b.new    <- fixef(mod.final)
  sigma    <- vcov(mod.final)
  varcov   <- getVarCov(mod.final)
  # Take the sigma for the residuals, one for each of the subjects.
  se   <- (summary(mod.final)$sigma)
  w.Y2 <- coef(mod.final$modelStruct$varStruct, uncons = FALSE)
  se   <- c(se, se*w.Y2)
  
  # varcomp_vcov(mod.init) this gets the information matrix of the variance components
  
  out <- list(coefficients = b.new, covariance = sigma, D = as.vector(t(varcov)),
              Sigma = se)
  
  return(out)

}


MI <- function(formula.fixed, formula.random, formula.imp, data, n.imp, id,
                           grp, timegrp, time, yname, iter){
  
  Ests.Imp  <- Covs.Imp <- list()
  D         <- matrix(NA, ncol = 16, nrow = n.imp)
  Sigma     <- matrix(NA, ncol = 2, nrow = n.imp)
  
  for(i in 1:n.imp){
    print(i)
    mod  <- MI.Imputation(formula.fixed, formula.random, formula.imp, data, iter, id,
                           grp, timegrp, time, yname)
    # save results
    Ests.Imp[[i]]   <- mod$coefficients
    Covs.Imp[[i]]   <- mod$covariance
    D[i,]           <- mod$D
    Sigma[i,]       <- mod$Sigma
  }
  
  
  # put results together using Rubin's rule
  coeff.res <- MIcombine(Ests.Imp, Covs.Imp)
  # save D and Sigma
  D         <- apply(D, 2, mean)
  Sigma     <- apply(Sigma, 2, mean)
  out       <- list(coefficients = coeff.res$coefficients, 
                     covariance = diag(coeff.res$variance),
                    D = D, Sigma = Sigma)
  
  return(out)
  
}

# create subject data
CreateSubjectData2 <- function(id, y, x, z){
  
  # split data to have yi, xi and zi
  split_y <- split(y, id)
  split_x <- split(x, id)
  split_z <- split(z, id)
  
  # create list with subject data
  n.id        <- length(unique(id))
  subjectData <- list()
  
  for(i in 1:n.id){
    yi             <- split_y[[i]]
    Zi             <- matrix(split_z[[i]], nrow = length(yi))
    colnames(Zi)   <- colnames(z)
    Xi             <- matrix(split_x[[i]], nrow = length(yi))
    colnames(Xi)   <- colnames(x)
    
    subjectData[[i]] <- list(yi = yi, Xi = Xi, Zi = Zi)
    
  }
  
  # return a list of the outcome, fixed effects and random effects
  subjectData
}
