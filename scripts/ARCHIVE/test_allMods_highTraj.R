# load data
scenarios <- readRDS(here("data", "scenarios.RDS"))
high.lam.combos <- readRDS(here("data", "high.lam.combos.RDS"))

# functions
source(here("scripts", "current version",
            "1 - simulating data", "IPM_sim_2.0function.R"))
source(here("scripts", "current version",
            "2 - models", "IPM_marray.R"))
source(here("scripts", "current version",
            "4 - run models", "run_scenarios_helperFns.R"))

high.comb <- high.lam.combos[sample(1:5000, 1), 1:3]

highpopTraj <- simPopTrajectory(n.years=15,
                               n.data.types=c(0.25,0.25,0.25),
                               age.init=c(150,150),
                               phi.1=as.numeric(high.comb[2]),
                               phi.ad=as.numeric(high.comb[3]),
                               f=as.numeric(high.comb[1]))

# simulate data

detect.l <- 0.3
detect.m <- 0.5
detect.h <- 0.8

highpopDat <- simData (indfates = highpopTraj$indfates, #highpopTraj$indfates
                      n.years = 15,
                      n.data.types = c(0.25,0.25,0.25),
                      ADonly = T,
                      p.1 = detect.h,
                      p.ad = detect.h,
                      p.count = detect.h,
                      p.prod = detect.h,
                      BinMod = T,
                      n.sam = 3,
                      sig = 0,
                      productivity = T)

# RUN MODELS ###########


######### IPM ######

#########

#### DATA ####
dat1 <- list(y = highpopDat$SUR,
             ch.y = highpopDat$ch,
             OBS_nestlings = highpopDat$OBS_nestlings,
             R_obs = highpopDat$R_obs
)


#### CONSTANTS ####

const1 <- list(nyears = ncol(highpopDat$ch),
               n.sam = nrow(highpopDat$SUR),
               n.ind = nrow(highpopDat$ch),
               first = highpopDat$firstobs)

#### INITIAL VALUES ####
z.state <- state.data(highpopDat$ch)

inits1 <- list(
  mean.phi = c(high.comb$phi1, high.comb$phiad),#c(detect.h, detect.h),
  mean.p = detect.h,
  p.surv = detect.h,
  fec = high.comb$fec,#detect.h,
  #mean.phi = runif(2,0,1),#c(detect.h, detect.h),
  #mean.p = runif(1,0,1),#detect.h,
  #fec = runif(1,0,5),#detect.h,
  z=z.state,
  n1.start=highpopTraj$Nouts[1,1], #HAS changed this to just pull from popTraj
  nad.start=highpopTraj$Nouts[2,1]
)

#### PARAMETERS TO MONITOR ####
params1 <- c("p.surv", "mean.phi","mean.p", "fec", "lambda") #,"Ntot","N1","Nad","f","rho")#0.3764911
#dont need all of these later, but wanted to look at them now

#### MCMC SETTINGS ####
nb <- 100000#0 #burn-in
ni <- nb + nb #total iterations
nt <- 10  #thin
nc <- 3  #chains

#### COMPILE CONFIGURE AND BUILD ####
Rmodel1 <- nimbleModel(code = IPMmod, constants = const1, data = dat1,
                       check = FALSE, calculate = FALSE, inits = inits1)
conf1 <- configureMCMC(Rmodel1, monitors = params1)#, thin = nt,
#control = list(maxContractions = 1000))
# lots of initial model checking you can do by exploring conf1
# if you wanted to change samplers this is where you would do that
Rmcmc1 <- buildMCMC(conf1)
Cmodel1 <- compileNimble(Rmodel1, showCompilerOutput = FALSE)
Cmcmc1 <- compileNimble(Rmcmc1, project = Rmodel1)
library(beepr)
beep(sound = 2)

#### RUN MCMC ####
t.start <- Sys.time()
#sink("sad_output.txt")
outIPM <- runMCMC(Cmcmc1, niter = ni , nburnin = nb , nchains = nc, inits = inits1,
                  setSeed = FALSE, progressBar = TRUE, samplesAsCodaMCMC = TRUE)
#sink()
t.end <- Sys.time()
(runTime <- t.end - t.start)
beep(sound = 3)

########HAS: if it doesnt run and there are slice sampler or other issues, use:
Cmcmc1$mvSamples[["Ntot"]][[1]] # or any other parameters to help diagnose
#I believe it is the starting values for the chain, so you can (hopefully) find the problem

####

# run each model

library(coda)

# check that results are unbiased
summary(outIPM)
high.comb

# check model convergence

gelman.diag(outIPM)
#plot(out2)

######### NO NESTS ######

#########

#### DATA ####
dat1 <- list(y = highpopDat$SUR,
             ch.y = highpopDat$ch
)


#### CONSTANTS ####

const1 <- list(nyears = ncol(highpopDat$ch),
               n.sam = nrow(highpopDat$SUR),
               n.ind = nrow(highpopDat$ch),
               first = highpopDat$firstobs)

#### INITIAL VALUES ####
z.state <- state.data(highpopDat$ch)

inits1 <- list(mean.phi = c(high.comb$phi1, high.comb$phiad),#c(detect.h, detect.h),
               mean.p = detect.h,
               p.surv = detect.h,
               fec = high.comb$fec,#detect.h,
               #mean.phi = runif(2,0,1),#c(detect.h, detect.h),
               #mean.p = runif(1,0,1),#detect.h,
               #fec = runif(1,0,5),#detect.h,
               z=z.state,
               n1.start=highpopTraj$Nouts[1,1], #HAS changed this to just pull from popTraj
               nad.start=highpopTraj$Nouts[2,1]
)

#### PARAMETERS TO MONITOR ####
params1 <- c("p.surv", "mean.phi","mean.p", "fec", "lambda") #,"Ntot","N1","Nad","f","rho")#0.3764911
#dont need all of these later, but wanted to look at them now

#### MCMC SETTINGS ####
nb <- 100000#0 #burn-in
ni <- nb + nb #total iterations
nt <- 10  #thin
nc <- 3  #chains

#### COMPILE CONFIGURE AND BUILD ####
Rmodel1 <- nimbleModel(code = nonests, constants = const1, data = dat1,
                       check = FALSE, calculate = FALSE, inits = inits1)
conf1 <- configureMCMC(Rmodel1, monitors = params1)#, thin = nt,
#control = list(maxContractions = 1000))
# lots of initial model checking you can do by exploring conf1
# if you wanted to change samplers this is where you would do that
Rmcmc1 <- buildMCMC(conf1)
Cmodel1 <- compileNimble(Rmodel1, showCompilerOutput = FALSE)
Cmcmc1 <- compileNimble(Rmcmc1, project = Rmodel1)
library(beepr)
beep(sound = 2)

#### RUN MCMC ####
t.start <- Sys.time()
#sink("sad_output.txt")
#changed to checking to just see a matrix, since it is working!
outnonests <- runMCMC(Cmcmc1, niter = ni , nburnin = nb , nchains = nc, inits = inits1,
                      setSeed = FALSE, progressBar = TRUE, samplesAsCodaMCMC = TRUE)
#sink()
t.end <- Sys.time()
(runTime <- t.end - t.start)
beep(sound = 3)

########HAS: if it doesnt run and there are slice sampler or other issues, use:
Cmcmc1$mvSamples[["Ntot"]][[1]] # or any other parameters to help diagnose
#I believe it is the starting values for the chain, so you can (hopefully) find the problem

####

# run each model

library(coda)

# check that results are unbiased
summary(outnonests)
high.comb

# check model convergence

gelman.diag(outnonests)
#plot(out2)

######### NO MR ######

#########

#### DATA ####
dat1 <- list(y = highpopDat$SUR,
             OBS_nestlings = highpopDat$OBS_nestlings,
             R_obs = highpopDat$R_obs
)


#### CONSTANTS ####

const1 <- list(nyears = ncol(highpopDat$ch),
               n.sam = nrow(highpopDat$SUR))

#### INITIAL VALUES ####
z.state <- state.data(highpopDat$ch)

inits1 <- list(mean.phi = c(high.comb$phi1, high.comb$phiad),#c(detect.h, detect.h),
               #mean.p = detect.h,
               p.surv = detect.h,
               fec = high.comb$fec,#detect.h,
               #mean.phi = runif(2,0,1),#c(detect.h, detect.h),
               #mean.p = runif(1,0,1),#detect.h,
               #fec = runif(1,0,5),#detect.h,
               #z=z.state,
               n1.start=highpopTraj$Nouts[1,1], #HAS changed this to just pull from popTraj
               nad.start=highpopTraj$Nouts[2,1]
)

#### PARAMETERS TO MONITOR ####
params1 <- c("p.surv", "mean.phi", "fec", "lambda") #,"Ntot","N1","Nad","f","rho")#0.3764911
#dont need all of these later, but wanted to look at them now

#### MCMC SETTINGS ####
nb <- 100000#0 #burn-in
ni <- nb + nb #total iterations
nt <- 10  #thin
nc <- 3  #chains

#### COMPILE CONFIGURE AND BUILD ####
Rmodel1 <- nimbleModel(code = nomr, constants = const1, data = dat1,
                       check = FALSE, calculate = FALSE, inits = inits1)
conf1 <- configureMCMC(Rmodel1, monitors = params1)#, thin = nt,
#control = list(maxContractions = 1000))
# lots of initial model checking you can do by exploring conf1
# if you wanted to change samplers this is where you would do that
Rmcmc1 <- buildMCMC(conf1)
Cmodel1 <- compileNimble(Rmodel1, showCompilerOutput = FALSE)
Cmcmc1 <- compileNimble(Rmcmc1, project = Rmodel1)
library(beepr)
beep(sound = 2)

#### RUN MCMC ####
t.start <- Sys.time()
#sink("sad_output.txt")
outnomr <- runMCMC(Cmcmc1, niter = ni , nburnin = nb , nchains = nc, inits = inits1,
                   setSeed = FALSE, progressBar = TRUE, samplesAsCodaMCMC = TRUE)
#sink()
t.end <- Sys.time()
(runTime <- t.end - t.start)
beep(sound = 3)

# AEB note - this model is the fastest to run by a LOT

########HAS: if it doesnt run and there are slice sampler or other issues, use:
Cmcmc1$mvSamples[["Ntot"]][[1]] # or any other parameters to help diagnose
#I believe it is the starting values for the chain, so you can (hopefully) find the problem

####

# run each model

library(coda)

# check that results are unbiased
summary(outnomr)
high.comb

# check model convergence

gelman.diag(outnomr)
#plot(out2)

######### ABUND ONLY ######

#########

#### DATA ####
dat1 <- list(y = highpopDat$SUR
)


#### CONSTANTS ####

const1 <- list(nyears = ncol(highpopDat$ch),
               n.sam = nrow(highpopDat$SUR))

#### INITIAL VALUES ####
z.state <- state.data(highpopDat$ch)

inits1 <- list(mean.phi = c(high.comb$phi1, high.comb$phiad),#c(detect.h, detect.h),
               #mean.p = detect.h,
               p.surv = detect.h,
               fec = high.comb$fec,#detect.h,
               #mean.phi = runif(2,0,1),#c(detect.h, detect.h),
               #mean.p = runif(1,0,1),#detect.h,
               #fec = runif(1,0,5),#detect.h,
               #z=z.state,
               n1.start=highpopTraj$Nouts[1,1], #HAS changed this to just pull from popTraj
               nad.start=highpopTraj$Nouts[2,1]
)

#### PARAMETERS TO MONITOR ####
params1 <- c("p.surv", "mean.phi", "fec", "lambda") #,"Ntot","N1","Nad","f","rho")#0.3764911
#dont need all of these later, but wanted to look at them now

#### MCMC SETTINGS ####
nb <- 100000#0 #burn-in
ni <- nb + nb #total iterations
nt <- 10  #thin
nc <- 3  #chains

#### COMPILE CONFIGURE AND BUILD ####
Rmodel1 <- nimbleModel(code = abundonly, constants = const1, data = dat1,
                       check = FALSE, calculate = FALSE, inits = inits1)
conf1 <- configureMCMC(Rmodel1, monitors = params1)#, thin = nt,
#control = list(maxContractions = 1000))
# lots of initial model checking you can do by exploring conf1
# if you wanted to change samplers this is where you would do that
Rmcmc1 <- buildMCMC(conf1)
Cmodel1 <- compileNimble(Rmodel1, showCompilerOutput = FALSE)
Cmcmc1 <- compileNimble(Rmcmc1, project = Rmodel1)
library(beepr)
beep(sound = 2)

#### RUN MCMC ####
t.start <- Sys.time()
#sink("sad_output.txt")
outabund <- runMCMC(Cmcmc1, niter = ni , nburnin = nb , nchains = nc, inits = inits1,
                    setSeed = FALSE, progressBar = TRUE, samplesAsCodaMCMC = TRUE)
#sink()
t.end <- Sys.time()
(runTime <- t.end - t.start)
beep(sound = 3)

# AEB note
# this model is also very very fast

########HAS: if it doesnt run and there are slice sampler or other issues, use:
Cmcmc1$mvSamples[["Ntot"]][[1]] # or any other parameters to help diagnose
#I believe it is the starting values for the chain, so you can (hopefully) find the problem

####

# run each model

library(coda)

# check that results are unbiased
summary(outabund)
high.comb

# check model convergence

gelman.diag(outabund)
#plot(out2)
