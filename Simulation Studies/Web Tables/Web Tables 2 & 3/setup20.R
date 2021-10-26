#########################################################################################################################
#### Parameter simulation ###############################################################################################
#########################################################################################################################

N.sim               <- 2400
prev.grp.sim        <- 0.3
conf.param.sim      <- c(-0.15, -0.05) 
ni.sim              <- c(4, 6)
p.central           <- 0.70
# design
type                <- "balanced.incomplete"
# number of people sampled per stratum
NsPerStratumUniv.sim <- c(150, 100, 150)

# fixed effects
beta_y1.sim <- c(80, 0.5, -2.5, -1.5, -0.25, -0.10)
beta_y2.sim <- c(65, -0.6, -2, -1, -0.15, -0.15)

# variance component of random effects
sigma.a1.sim <- 81/4; sigma.a2.sim <- 36/4; sigma.b1.sim <- 4/4; sigma.b2.sim <- 1/4
sigma.a1a2.sim <- 30/4; sigma.a1b1.sim <- 1/4; sigma.a1b2.sim <- 0.5/4; sigma.a2b1.sim <- 3/4
sigma.a2b2.sim <- 1.5/4; sigma.b1b2.sim <- 1/4

# error term
rho.e.sim <- 0
sigma.e1.sim <- 9/4; sigma.e2.sim <- 4/4

# numebr of people sampled 
Nsample.sim <- 400

# number of iteration for MI
m <- 20

# for ACCRE
n                   <- 1
