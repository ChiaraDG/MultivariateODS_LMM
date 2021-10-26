# Two-Phase Studies with Multivariate Longitudinal Data

The R scripts in this repository contain the code necessary to reproduce simulations and data analysis of the manuscript "Two-Phase Studies with Multivariate Longitudinal Data" (Di Gravio, Tao, Schildcrout). The repository contains two folder:

* [Simulation Studies](https://github.com/ChiaraDG/MultivariateODS_LMM/tree/main/Simulation%20Studies): contains the code necessary to replicate the simulations in the paper

* [LHS](https://github.com/ChiaraDG/MultivariateODS_LMM/tree/main/LHS): contains the code necessary to replicate the analysis of the Lung Health Study

More detailed instruction on how to reproduce the results are provided in each folder

## Example

For each simulation scenario, there is a separate R file that specify the parameters used in the simulation. For instance:

```{r, eval = FALSE}
N.sim               <- 2400
prev.grp.sim        <- 0.3
conf.param.sim      <- c(-0.15, -0.05) 
ni.sim              <- c(4, 6)
p.central           <- 0.70
# design
type                 <- "balanced.incomplete"
NsPerStratumUniv.sim <- c(150, 100, 150)

# fixed effects
beta_y1.sim <- c(80, 0.5, -2.5, -1.5, -0.25, -0.10)
beta_y2.sim <- c(65, -0.6, -2, -1, -0.15, -0.15)

# variance component of random effects
#sigma.a1.sim <- 81/4; sigma.a2.sim <- 36/4; sigma.b1.sim <- 4/4; sigma.b2.sim <- 1/4
#sigma.a1a2.sim <- 30/4; sigma.a1b1.sim <- 1/4; sigma.a1b2.sim <- 0.5/4; sigma.a2b1.sim <- 3/4
#sigma.a2b2.sim <- 1.5/4; sigma.b1b2.sim <- 1/4                                           

# error term
rho.e.sim <- 0
sigma.e1.sim <- 9/4; sigma.e2.sim <- 4/4
```

The data, the cut-offs to be used under ODS and/or BDS scenario and the selection of informative individuals can be done using the following code:

```{r, eval = FALSE}
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
  
dat.samp <- SampledFun(d = dat, Nsample = Nsample.sim, NsPerStratumUniv = NsPerStratumUniv.sim,
                         Y1.cutoffsODS =  cutoffs.Y1.ODS, Y2.cutoffsODS = cutoffs.Y2.ODS,
                         Y1.cutoffsBDS =  cutoffs.Y1.BDS, Y2.cutoffsBDS = cutoffs.Y2.BDS)
fulldat <- dat.samp$dat
```

Next, the user can decide which design and analysis procedure they can run. The code in this repository allows for the following:

1. **Random Sampling with Maximum Likelihood**

```{r, eval = FALSE}
dat.longRan <- fulldat %>% filter(SampledRan == 1) %>% gather(., key = "y_name", value = "value", Y1, Y2) %>%
               mutate(D1 = as.integer(y_name == "Y1"), D2 = as.integer(y_name == "Y2")) %>%
               arrange(id)
     
lme.ML      <- lme(value ~ 0 + D1 + D1:grp + D1:time + D1:conf + D1:timegrp + D1:timeconf +
                           D2 + D2:grp + D2:time + D2:conf + D2:timegrp + D2:timeconf, data = dat.longRan,
                           random = ~ 0 + D1 + D1:time + D2 + D2:time | id,
                           na.action = na.omit, weights = varIdent(form = ~1 | y_name),
                           control = lmeControl(maxIter = 1e9, opt = c("optim")))
```

2. **ODS Sampling with Ascertainment Corrected Maximum Likelihood**

```{r, eval = F}
dat.longODS <- fulldat %>% filter(SampledSlpODS == 1) %>% gather(., key = "y_name", value = "Y", Y1, Y2) %>%
               mutate(D1 = as.integer(y_name == "Y1"), D2 = as.integer(y_name == "Y2")) %>%
               arrange(id)
     
acml.ods    <- acml.lmem2(formula.fixed = Y ~ 0 + D1 + g1 + c1 + t1 + gt1 + ct1 + D2 + g2 + c2 + t2 + gt2 + ct2, 
                          formula.random = ~ 0 + D1 + t1 + D2 + t2,
                          data = dat.longODS, id = id, InitVals = truevals.biv.ods,
                          SampProb = SampProbMixODS, cutpoints = cutpointsMixODS,
                          Weights = NoWeighting, w.function = w.functionMixSlpODS)
```

3. **BDS Sampling with Ascertainment Corrected Maximum Likelihood**

```{r, eval = FALSE}
dat.longBDS <- fulldat %>% filter(SampledSlpBDS == 1) %>% gather(., key = "y_name", value = "Y", Y1, Y2) %>%
               mutate(D1 = as.integer(y_name == "Y1"), D2 = as.integer(y_name == "Y2")) %>%
               arrange(id)
               
acml.bds    <- acml.lmem2(formula.fixed = Y ~ 0 + D1 + g1 + c1 + t1 + gt1 + ct1 + D2 + g2 + c2 + t2 + gt2 + ct2,
                          formula.random = ~ 0 + D1 + t1 + D2 + t2,
                          data = dat.longBDS, id = id, InitVals = truevals.biv.ods,
                          SampProb = SampProbMixBDS, cutpoints = cutpointsMixBDS,
                          Weights = NoWeighting, w.function = w.functionMixSlpBDS,
                          xcol.phase1 = c(1, 3, 4, 6, 7, 9, 10, 12), ests.phase1 = CoefPhase1)
```

4. **Random Sampling with Multiple Imputation**

```{r, eval = FALSE}
# prepare the data for MI
dat.long       <- fulldat %>% gather(., key = "y_name", value = "Y", Y1, Y2) %>%
                  mutate(D1 = as.integer(y_name == "Y1"), D2 = as.integer(y_name == "Y2")) %>% arrange(id)
dat.long$t1    <- dat.long$D1*dat.long$time
dat.long$g1    <- dat.long$D1*dat.long$grp
dat.long$c1    <- dat.long$D1*dat.long$conf
dat.long$gt1   <- dat.long$D1*dat.long$timegrp
dat.long$ct1   <- dat.long$D1*dat.long$timeconf
dat.long$t2    <- dat.long$D2*dat.long$time
dat.long$g2    <- dat.long$D2*dat.long$grp
dat.long$c2    <- dat.long$D2*dat.long$conf
dat.long$gt2   <- dat.long$D2*dat.long$timegrp
dat.long$ct2   <- dat.long$D2*dat.long$timeconf

dat.rs         <- dat.long
dat.rs$grp     <- ifelse(dat.rs$SampledRan == 1, dat.rs$grp, NA)
dat.rs$timegrp <- ifelse(dat.rs$SampledRan == 1, dat.rs$timegrp, NA)
dat.rs$g1      <- ifelse(dat.rs$SampledRan == 1, dat.rs$g1, NA)
dat.rs$g2      <- ifelse(dat.rs$SampledRan == 1, dat.rs$g2, NA)
dat.rs$gt1     <- ifelse(dat.rs$SampledRan == 1, dat.rs$gt1, NA)
dat.rs$gt2     <- ifelse(dat.rs$SampledRan == 1, dat.rs$gt2, NA)


model.rs.mi   <- MI(formula.fixed = Y ~ 0 + D1 + g1 + c1 + t1 + gt1 + ct1 + D2 + g2 + c2 + t2 + gt2 + ct2, 
                   formula.random = ~0 + D1 + t1 + D2 + t2 | id, 
                   formula.imp = grp ~ conf, data = dat.rs, 
                   n.imp = 70, iter = 10, id = "id", grp = c("g1", "g2"), 
                   timegrp = c("gt1", "gt2"), time = c("t1", "t2"), yname = "y_name")
```

4. **ODS Sampling with Multiple Imputation**

```{r, eval = FALSE}
# prepare the data for MI
dat.long       <- fulldat %>% gather(., key = "y_name", value = "Y", Y1, Y2) %>%
                  mutate(D1 = as.integer(y_name == "Y1"), D2 = as.integer(y_name == "Y2")) %>% arrange(id)
dat.long$t1    <- dat.long$D1*dat.long$time
dat.long$g1    <- dat.long$D1*dat.long$grp
dat.long$c1    <- dat.long$D1*dat.long$conf
dat.long$gt1   <- dat.long$D1*dat.long$timegrp
dat.long$ct1   <- dat.long$D1*dat.long$timeconf
dat.long$t2    <- dat.long$D2*dat.long$time
dat.long$g2    <- dat.long$D2*dat.long$grp
dat.long$c2    <- dat.long$D2*dat.long$conf
dat.long$gt2   <- dat.long$D2*dat.long$timegrp
dat.long$ct2   <- dat.long$D2*dat.long$timeconf

dat.ods         <- dat.long
dat.ods$grp     <- ifelse(dat.ods$SampledSlpODS == 1, dat.ods$grp, NA)
dat.ods$timegrp <- ifelse(dat.ods$SampledSlpODS == 1, dat.ods$timegrp, NA)
dat.ods$g1      <- ifelse(dat.ods$SampledSlpODS == 1, dat.ods$g1, NA)
dat.ods$g2      <- ifelse(dat.ods$SampledSlpODS == 1, dat.ods$g2, NA)
dat.ods$gt1     <- ifelse(dat.ods$SampledSlpODS == 1, dat.ods$gt1, NA)
dat.ods$gt2     <- ifelse(dat.ods$SampledSlpODS == 1, dat.ods$gt2, NA)


model.ods.mi   <- MI(formula.fixed = Y ~ 0 + D1 + g1 + c1 + t1 + gt1 + ct1 + D2 + g2 + c2 + t2 + gt2 + ct2, 
                   formula.random = ~0 + D1 + t1 + D2 + t2 | id, 
                   formula.imp = grp ~ conf, data = dat.ods, 
                   n.imp = 70, iter = 10, id = "id", grp = c("g1", "g2"), 
                   timegrp = c("gt1", "gt2"), time = c("t1", "t2"), yname = "y_name")
```

5. **BDS Sampling with Multiple Imputation**

```{r, eval = FALSE}
# prepare the data for MI
dat.long       <- fulldat %>% gather(., key = "y_name", value = "Y", Y1, Y2) %>%
                  mutate(D1 = as.integer(y_name == "Y1"), D2 = as.integer(y_name == "Y2")) %>% arrange(id)
dat.long$t1    <- dat.long$D1*dat.long$time
dat.long$g1    <- dat.long$D1*dat.long$grp
dat.long$c1    <- dat.long$D1*dat.long$conf
dat.long$gt1   <- dat.long$D1*dat.long$timegrp
dat.long$ct1   <- dat.long$D1*dat.long$timeconf
dat.long$t2    <- dat.long$D2*dat.long$time
dat.long$g2    <- dat.long$D2*dat.long$grp
dat.long$c2    <- dat.long$D2*dat.long$conf
dat.long$gt2   <- dat.long$D2*dat.long$timegrp
dat.long$ct2   <- dat.long$D2*dat.long$timeconf

dat.bds         <- dat.long
dat.bds$grp     <- ifelse(dat.bds$SampledSlpBDS == 1, dat.bds$grp, NA)
dat.bds$timegrp <- ifelse(dat.bds$SampledSlpBDS == 1, dat.bds$timegrp, NA)
dat.bds$g1      <- ifelse(dat.bds$SampledSlpBDS == 1, dat.bds$g1, NA)
dat.bds$g2      <- ifelse(dat.bds$SampledSlpBDS == 1, dat.bds$g2, NA)
dat.bds$gt1     <- ifelse(dat.bds$SampledSlpBDS == 1, dat.bds$gt1, NA)
dat.bds$gt2     <- ifelse(dat.bds$SampledSlpBDS == 1, dat.bds$gt2, NA)

model.bds.mi   <- MI(formula.fixed = Y ~ 0 + D1 + g1 + c1 + t1 + gt1 + ct1 + D2 + g2 + c2 + t2 + gt2 + ct2, 
                   formula.random = ~0 + D1 + t1 + D2 + t2 | id, 
                   formula.imp = grp ~ conf, data = dat.bds, 
                   n.imp = 70, iter = 10, id = "id", grp = c("g1", "g2"), 
                   timegrp = c("gt1", "gt2"), time = c("t1", "t2"), yname = "y_name")
```
