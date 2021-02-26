# todo

# simulate trajectory

# load data
scenarios <- readRDS(here("data", "scenarios.RDS"))
low.lam.combos <- readRDS(here("data", "low.lam.combos.RDS"))
med.lam.combos <- readRDS(here("data", "med.lam.combos.RDS"))
high.lam.combos <- readRDS(here("data", "high.lam.combos.RDS"))

# functions
source(here("scripts", "in progress scripts", "simHelperFns.R"))
source(here("scripts", "in progress scripts", "IPM_sim_2.0function.R"))
source(here("scripts", "in progress scripts", "compute_time_calc.R"))

low.comb <- low.lam.combos[sample(1:5000, 1), 1:3]
med.comb <- med.lam.combos[sample(1:5000, 1), 1:3]
high.comb <- high.lam.combos[sample(1:5000, 1), 1:3]

# thought
# make adult survival be higher?? force it

lowpopTraj <- simPopTrajectory(n.years=15, 
                               n.data.types=c(0.25,0.25,0.25),
                               age.init=c(150,150), 
                               phi.1=as.numeric(low.comb[2]), 
                               phi.ad=as.numeric(low.comb[3]), 
                               f=as.numeric(low.comb[1]))
# returns indfates and n

medpopTraj <- simPopTrajectory(n.years=15, 
                               n.data.types=c(0.25,0.25,0.25),
                               age.init=c(150,150), 
                               phi.1=as.numeric(med.comb[2]), 
                               phi.ad=as.numeric(med.comb[3]), 
                               f=as.numeric(med.comb[1]))

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

# next steps
# abby needs to make the repository more readable

lowpopDat <- simData (indfates = lowpopTraj$indfates, 
                      n.years = 15, 
                      n.data.types = c(0.25,0.25,0.25), 
                      ADonly = T, 
                      p.1 = detect.l, 
                      p.ad = detect.l, 
                      p.count = detect.l,
                      p.prod = detect.l,
                      BinMod = T, 
                      n.sam = 3,  
                      sig = 0, 
                      productivity = T)


medpopDat <- simData (indfates = medpopTraj$indfates, #lowpopTraj$indfates
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
  
highpopDat <- simData (indfates = lowpopTraj, 
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

#####

#### DATA ####
dat1 <- list(y = lowpopDat$SUR, 
             ch.y = lowpopDat$ch, 
             OBS_nestlings = lowpopDat$OBS_nestlings, 
             R_obs = lowpopDat$R_obs
)


#### CONSTANTS ####

# AEB TODO
# this does not work well...
#HAS: we dont need age if they are all adults at their first resight?
#age = ageunknown(lowpopDat$age_ch)
const1 <- list(nyears = ncol(lowpopDat$ch), 
               n.sam = nrow(lowpopDat$SUR),
               n.ind = nrow(lowpopDat$ch),
               first = lowpopDat$firstobs)
               #age = age) age doesnt matter, so removing it

#### INITIAL VALUES ####
z.state <- state.data(lowpopDat$ch)
inits1 <- list(p.surv = detect.l,
               mean.phi = runif(2,0,1),#c(detect.l, detect.l),
               mean.p = runif(1,0,1),#detect.l,
               fec = runif(1,0,1),#detect.l, 
               z=z.state,
               # TODO
               # why are these 150 - should start ad stable age distrib
               # need to double check this
               #n1.start= sum(lowpopTraj$indfates[1, 1, ], na.rm = TRUE),
               #nad.start= sum(lowpopTraj$indfates[2, 1, ], na.rm = TRUE)
               n1.start=lowpopTraj$Nouts[1,1], #HAS changed this to just pull from popTraj
               nad.start=lowpopTraj$Nouts[2,1] 
               )

#### PARAMETERS TO MONITOR ####
params1 <- c("p.surv", "mean.phi", "fec", "lambda","Ntot","N1","Nad","f","rho")#0.3764911

#### MCMC SETTINGS ####
nb <- 10000#000 #burn-in
ni <- nb + nb #total iterations
nt <- 10  #thin
nc <- 1  #chains

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
beep(sound = 2)

#### RUN MCMC ####
t.start <- Sys.time()
sink("sad_output.txt")
out2 <- runMCMC(Cmcmc1, niter = ni , nburnin = nb , nchains = nc, inits = inits1,
                setSeed = FALSE, progressBar = TRUE, samplesAsCodaMCMC = TRUE)  
sink()
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
summary(out2)
low.comb

# check model convergence

# AEB TODO
# what the heck why is it mixing so terribly?
gelman.diag(out2)
plot(out2)



