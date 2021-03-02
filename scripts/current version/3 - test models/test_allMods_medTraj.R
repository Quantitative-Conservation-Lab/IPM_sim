# load data
scenarios <- readRDS(here("data", "scenarios.RDS"))
med.lam.combos <- readRDS(here("data", "med.lam.combos.RDS"))

# functions
#source(here("scripts", "in progress scripts", "simHelperFns.R"))
source(here("scripts", "in progress scripts", "IPM_sim_2.0function.R"))
source(here("scripts", "in progress scripts", "compute_time_calc.R"))
source(here("scripts", "in progress scripts", "IPMNimble_v2.0.R"))

med.comb <- med.lam.combos[sample(1:5000, 1), 1:3]

medpopTraj <- simPopTrajectory(n.years=15, 
                               n.data.types=c(0.25,0.25,0.25),
                               age.init=c(150,150), 
                               phi.1=as.numeric(med.comb[2]), 
                               phi.ad=as.numeric(med.comb[3]), 
                               f=as.numeric(med.comb[1]))

# simulate data

detect.m <- 0.5

medpopDat <- simData (indfates = medpopTraj$indfates, #medpopTraj$indfates
                      n.years = 15, 
                      n.data.types = c(0.25,0.25,0.25), 
                      ADonly = T, 
                      p.1 = detect.m, 
                      p.ad = detect.m, 
                      p.count = detect.m,
                      p.prod = detect.m,
                      BinMod = T, 
                      n.sam = 3,  
                      sig = 0, 
                      productivity = T)

# RUN MODELS ###########


######### IPM ######

#########

#### DATA ####
dat1 <- list(y = medpopDat$SUR, 
             ch.y = medpopDat$ch, 
             OBS_nestlings = medpopDat$OBS_nestlings, 
             R_obs = medpopDat$R_obs
)


#### CONSTANTS ####

# AEB TODO
# this does not work well...
#HAS: we dont need age if they are all adults at their first resight?
#age = ageunknown(medpopDat$age_ch)
const1 <- list(nyears = ncol(medpopDat$ch), 
               n.sam = nrow(medpopDat$SUR),
               n.ind = nrow(medpopDat$ch),
               first = medpopDat$firstobs)
#age = age) age doesnt matter, so removing it

#### INITIAL VALUES ####
z.state <- state.data(medpopDat$ch)

inits1 <- list(#p.surv = runif(1,0,1),
  # AEB - ok I am an idiot - was putting the detection parameters as initial values
  # on the demographic parameters!!!!!!!!! ugh
  
  mean.phi = c(med.comb$phi1, med.comb$phiad),#c(detect.m, detect.m),
  mean.p = detect.m,
  p.surv = detect.m,
  fec = med.comb$fec,#detect.m, 
  #mean.phi = runif(2,0,1),#c(detect.m, detect.m),
  #mean.p = runif(1,0,1),#detect.m,
  #fec = runif(1,0,5),#detect.m, 
  z=z.state,
  # TODO
  # why are these 150 - should start ad stable age distrib
  # need to double check this
  #n1.start= sum(medpopTraj$indfates[1, 1, ], na.rm = TRUE),
  #nad.start= sum(medpopTraj$indfates[2, 1, ], na.rm = TRUE)
  n1.start=medpopTraj$Nouts[1,1], #HAS changed this to just pull from popTraj
  nad.start=medpopTraj$Nouts[2,1] 
  #n1.start=medpopTraj$Nouts[1,1], #HAS changed this to just pull from popTraj
  #nad.start=medpopTraj$Nouts[2,1] 
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
#changed to checking to just see a matrix, since it is working!
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
med.comb

# check model convergence

# AEB TODO
# what the heck why is it mixing so terribly?
gelman.diag(outIPM)
#plot(out2)

######### NO NESTS ######

#########

#### DATA ####
dat1 <- list(y = medpopDat$SUR, 
             ch.y = medpopDat$ch
)


#### CONSTANTS ####

# AEB TODO
# this does not work well...
#HAS: we dont need age if they are all adults at their first resight?
#age = ageunknown(medpopDat$age_ch)
const1 <- list(nyears = ncol(medpopDat$ch), 
               n.sam = nrow(medpopDat$SUR),
               n.ind = nrow(medpopDat$ch),
               first = medpopDat$firstobs)
#age = age) age doesnt matter, so removing it

#### INITIAL VALUES ####
z.state <- state.data(medpopDat$ch)

inits1 <- list(#p.surv = runif(1,0,1),
  # AEB - ok I am an idiot - was putting the detection parameters as initial values
  # on the demographic parameters!!!!!!!!! ugh
  
  mean.phi = c(med.comb$phi1, med.comb$phiad),#c(detect.m, detect.m),
  mean.p = detect.m,
  p.surv = detect.m,
  fec = med.comb$fec,#detect.m, 
  #mean.phi = runif(2,0,1),#c(detect.m, detect.m),
  #mean.p = runif(1,0,1),#detect.m,
  #fec = runif(1,0,5),#detect.m, 
  z=z.state,
  # TODO
  # why are these 150 - should start ad stable age distrib
  # need to double check this
  #n1.start= sum(medpopTraj$indfates[1, 1, ], na.rm = TRUE),
  #nad.start= sum(medpopTraj$indfates[2, 1, ], na.rm = TRUE)
  n1.start=medpopTraj$Nouts[1,1], #HAS changed this to just pull from popTraj
  nad.start=medpopTraj$Nouts[2,1] 
  #n1.start=medpopTraj$Nouts[1,1], #HAS changed this to just pull from popTraj
  #nad.start=medpopTraj$Nouts[2,1] 
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
med.comb

# check model convergence

# AEB TODO
# what the heck why is it mixing so terribly?
gelman.diag(outnonests)
#plot(out2)

######### NO MR ######

#########

#### DATA ####
dat1 <- list(y = medpopDat$SUR, 
             OBS_nestlings = medpopDat$OBS_nestlings, 
             R_obs = medpopDat$R_obs
)


#### CONSTANTS ####

# AEB TODO
# this does not work well...
#HAS: we dont need age if they are all adults at their first resight?
#age = ageunknown(medpopDat$age_ch)
const1 <- list(nyears = ncol(medpopDat$ch), 
               n.sam = nrow(medpopDat$SUR))
#age = age) age doesnt matter, so removing it

#### INITIAL VALUES ####
z.state <- state.data(medpopDat$ch)

inits1 <- list(#p.surv = runif(1,0,1),
  # AEB - ok I am an idiot - was putting the detection parameters as initial values
  # on the demographic parameters!!!!!!!!! ugh
  
  mean.phi = c(med.comb$phi1, med.comb$phiad),#c(detect.m, detect.m),
  #mean.p = detect.m,
  p.surv = detect.m,
  fec = med.comb$fec,#detect.m, 
  #mean.phi = runif(2,0,1),#c(detect.m, detect.m),
  #mean.p = runif(1,0,1),#detect.m,
  #fec = runif(1,0,5),#detect.m, 
  #z=z.state,
  # TODO
  # why are these 150 - should start ad stable age distrib
  # need to double check this
  #n1.start= sum(medpopTraj$indfates[1, 1, ], na.rm = TRUE),
  #nad.start= sum(medpopTraj$indfates[2, 1, ], na.rm = TRUE)
  n1.start=medpopTraj$Nouts[1,1], #HAS changed this to just pull from popTraj
  nad.start=medpopTraj$Nouts[2,1] 
  #n1.start=medpopTraj$Nouts[1,1], #HAS changed this to just pull from popTraj
  #nad.start=medpopTraj$Nouts[2,1] 
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
#changed to checking to just see a matrix, since it is working!
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
med.comb

# check model convergence

# AEB TODO
# what the heck why is it mixing so terribly?
gelman.diag(outnomr)
#plot(out2)

######### ABUND ONLY ######

#########

#### DATA ####
dat1 <- list(y = medpopDat$SUR
)


#### CONSTANTS ####

# AEB TODO
# this does not work well...
#HAS: we dont need age if they are all adults at their first resight?
#age = ageunknown(medpopDat$age_ch)
const1 <- list(nyears = ncol(medpopDat$ch), 
               n.sam = nrow(medpopDat$SUR))
#age = age) age doesnt matter, so removing it

#### INITIAL VALUES ####
z.state <- state.data(medpopDat$ch)

inits1 <- list(#p.surv = runif(1,0,1),
  # AEB - ok I am an idiot - was putting the detection parameters as initial values
  # on the demographic parameters!!!!!!!!! ugh
  
  mean.phi = c(med.comb$phi1, med.comb$phiad),#c(detect.m, detect.m),
  #mean.p = detect.m,
  p.surv = detect.m,
  fec = med.comb$fec,#detect.m, 
  #mean.phi = runif(2,0,1),#c(detect.m, detect.m),
  #mean.p = runif(1,0,1),#detect.m,
  #fec = runif(1,0,5),#detect.m, 
  #z=z.state,
  # TODO
  # why are these 150 - should start ad stable age distrib
  # need to double check this
  #n1.start= sum(medpopTraj$indfates[1, 1, ], na.rm = TRUE),
  #nad.start= sum(medpopTraj$indfates[2, 1, ], na.rm = TRUE)
  n1.start=medpopTraj$Nouts[1,1], #HAS changed this to just pull from popTraj
  nad.start=medpopTraj$Nouts[2,1] 
  #n1.start=medpopTraj$Nouts[1,1], #HAS changed this to just pull from popTraj
  #nad.start=medpopTraj$Nouts[2,1] 
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
#changed to checking to just see a matrix, since it is working!
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
med.comb

# check model convergence

# AEB TODO
# what the heck why is it mixing so terribly?
gelman.diag(outabund)
#plot(out2)