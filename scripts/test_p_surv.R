# todo

# simulate trajectory

# load data
scenarios <- readRDS(here("data", "scenarios.RDS"))
low.lam.combos <- readRDS(here("data", "low.lam.combos.RDS"))
med.lam.combos <- readRDS(here("data", "med.lam.combos.RDS"))
high.lam.combos <- readRDS(here("data", "high.lam.combos.RDS"))

# functions
source(here("scripts", "current version",
            "1 - simulating data", "IPM_sim_2.0function.R"))
source(here("scripts", "current version",
            "2 - models", "IPM_marray.R"))
source(here("scripts", "current version",
            "4 - run models", "run_scenarios_helperFns.R"))


 low.comb <- low.lam.combos[sample(1:5000, 1), 1:3]
 med.comb <- med.lam.combos[sample(1:5000, 1), 1:3]
 high.comb <- high.lam.combos[sample(1:5000, 1), 1:3]
# 
# # thought
# # make adult survival be higher?? force it
# 
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


#just load data from files
#is the true param values in scenarios.RDS? Where to pull those from
# low.lam.combos[1,] #for the traj below?
# lowpopTraj<-readRDS("~/Desktop/IPMSim/data/lowTrajectories/lowpopTraj-1-1.RDS")
# med.lam.combos[1,]
# medpopTraj<-readRDS("~/Desktop/IPMSim/data/medTrajectories/medpopTraj-1-1.RDS")
# high.lam.combos[1,]
# highpopTraj<-readRDS("~/Desktop/IPMSim/data/highTrajectories/highpopTraj-1-1.RDS")

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
                      p.1 = detect.h, 
                      p.ad = detect.h, 
                      p.count = detect.m,
                      p.prod = detect.h,
                      BinMod = T, 
                      n.sam = 3,  
                      sig = 0, 
                      productivity = T)

highpopDat <- simData (indfates = highpopTraj$indfates, 
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
#testing each scenario, make data for the models
# dat1l <- list(y = lowpopDat$SUR, 
#              marr = marray(lowpopDat$ch), 
#              R=rowSums(marray(lowpopDat$ch)),
#              OBS_nestlings = lowpopDat$OBS_nestlings, 
#              R_obs = lowpopDat$R_obs
# )
# dat1m <- list(y = medpopDat$SUR, 
#               marr = marray(medpopDat$ch), 
#               R=rowSums(marray(medpopDat$ch)),
#               OBS_nestlings = medpopDat$OBS_nestlings, 
#               R_obs = medpopDat$R_obs
# )
# dat1h <- list(y = highpopDat$SUR, 
#               marr = marray(highpopDat$ch), 
#               R=rowSums(marray(highpopDat$ch)),
#               OBS_nestlings = highpopDat$OBS_nestlings, 
#               R_obs = highpopDat$R_obs
# )

datfn<-function(popDat){
  datout<-list(
               y = popDat$SUR, 
               marr = marray(popDat$ch), 
               R=rowSums(marray(popDat$ch)),
               OBS_nestlings = popDat$OBS_nestlings, 
               R_obs = popDat$R_obs)
}

dat1l<-datfn(lowpopDat)
dat1m<-datfn(medpopDat)
dat1h<-datfn(highpopDat)
#### CONSTANTS ####

# AEB TODO
# this does not work well...
#HAS: we dont need age if they are all adults at their first resight?
#age = ageunknown(lowpopDat$age_ch)
#constants for each model for each scenario - will be the same, but whatev
const1l <- list(nyears = ncol(lowpopDat$SUR), 
               n.sam = nrow(lowpopDat$SUR))
const1m <- list(nyears = ncol(medpopDat$SUR), 
                n.sam = nrow(medpopDat$SUR))
const1h <- list(nyears = ncol(highpopDat$SUR), 
                n.sam = nrow(highpopDat$SUR))


#### INITIAL VALUES ####

initsfn<-function(poptraj_Noutsmat, d, combo){ #put in the population traj Nouts for initial values
  initsout <- list(
    p.surv = runif(1,0,1),
    mean.phi = c(combo$phi1, combo$phiad),#runif(2,0,1),
    mean.p = d, #runif(1,0,1),#detect.h,
    p.surv = d,#runif(1,0,1),#detect.h,
    fec = combo$fec,#runif(1,0,5),#high.comb$fec,#detect.l,
    n1.start=poptraj_Noutsmat[1,1],#highpopTraj$Nouts[1,1], #HAS changed this to just pull from popTraj
    nad.start=poptraj_Noutsmat[2,1])#highpopTraj$Nouts[2,1]
    #              # p.surv = runif(1,0,1),
                 # mean.phi = runif(2,0,1),
                 # mean.p = runif(1,0,1),#detect.h,
                 # p.surv = runif(1,0,1),#detect.h,
                 # fec =runif(1,0,5),#high.comb$fec,#detect.l,
                 # n1.start=poptraj_Noutsmat[1,1],#highpopTraj$Nouts[1,1], #HAS changed this to just pull from popTraj
                 # nad.start=poptraj_Noutsmat[2,1])#highpopTraj$Nouts[2,1]
  return(initsout)
}
initslow<-initsfn(lowpopTraj$Nouts, d=detect.l, combo=low.lam.combos[1,])
initsmed<-initsfn(medpopTraj$Nouts, d=detect.m, combo=med.lam.combos[1,])
initshigh<-initsfn(highpopTraj$Nouts, d=detect.h, combo=high.lam.combos[1,])
# inits1 <- list(p.surv = runif(1,0,1),
#                # AEB - ok I am an idiot - was putting the detection parameters as initial values
#                # on the demographic parameters!!!!!!!!! ugh
#                mean.phi =runif(2,0,1),
#                #mean.phi = c(high.comb$phi1, high.comb$phiad),#c(detect.l, detect.l),
#                mean.p = runif(1,0,1),#detect.h,
#                p.surv = runif(1,0,1),#detect.h,
#                fec = high.comb$fec,#detect.l, 
#                #mean.phi = runif(2,0,1),#c(detect.l, detect.l),
#                #mean.p = runif(1,0,1),#detect.l,
#                #fec = runif(1,0,5),#detect.l, 
#                # z=z.state,
#                # TODO
#                # why are these 150 - should start ad stable age distrib
#                # need to double check this
#                #n1.start= sum(lowpopTraj$indfates[1, 1, ], na.rm = TRUE),
#                #nad.start= sum(lowpopTraj$indfates[2, 1, ], na.rm = TRUE)
#                n1.start=highpopTraj$Nouts[1,1], #HAS changed this to just pull from popTraj
#                nad.start=highpopTraj$Nouts[2,1] 
#                #n1.start=medpopTraj$Nouts[1,1], #HAS changed this to just pull from popTraj
#                #nad.start=medpopTraj$Nouts[2,1] 
# )

#### PARAMETERS TO MONITOR ####
params1 <- c("p.surv", "mean.phi","mean.p", "fec", "lambda") #,"Ntot","N1","Nad","f","rho")#0.3764911
#dont need all of these later, but wanted to look at them now

#### MCMC SETTINGS ####
nb <- 100000 #burn-in
ni <- nb + nb #total iterations
nt <- 10  #thin
nc <- 1  #chains

#### lazy testing ####
#low
Rmodellow <- nimbleModel(code = IPMmod, constants = const1l, data = dat1l, 
                       check = FALSE, calculate = FALSE, inits = initslow)
conflow <- configureMCMC(Rmodellow, monitors = params1)
Rmcmclow <- buildMCMC(conflow)  
Cmodellow <- compileNimble(Rmodellow, showCompilerOutput = FALSE)
Cmcmclow <- compileNimble(Rmcmclow, project = Rmodellow)
Cmcmclow$run(niter=ni, nburnin=nb, thin=nt, reset=T)
out_low<-as.data.frame(as.matrix(Cmcmclow$mvSamples))
mean(out_low$p.surv)
summary(out_low$p.surv)

runlow<-runIPMmod(nb=nb, ni=ni, nt=nt, nc=5, popDat = lowpopDat, 
                   popTraj = lowpopTraj, comb=low.comb, detect=0.3)
summary(runlow[,19])
plot(runlow[,19])

#medium
Rmodelmed <- nimbleModel(code = IPMmod, constants = const1m, data = dat1m, 
                         check = FALSE, calculate = FALSE, inits = initsmed)# inits1)
confmed <- configureMCMC(Rmodelmed, monitors = params1)
Rmcmcmed <- buildMCMC(confmed)  
Cmodelmed <- compileNimble(Rmodelmed, showCompilerOutput = FALSE)
Cmcmcmed <- compileNimble(Rmcmcmed, project = Rmodelmed)
Cmcmcmed$run(niter=ni, nburnin=nb, thin=1, reset=T)
out_med<-as.data.frame(as.matrix(Cmcmcmed$mvSamples))
mean(out_med$p.surv)
summary(out_med$p.surv)

runmed<-runIPMmod(nb=nb, ni=ni, nt=nt, nc=5, popDat = medpopDat, 
                  popTraj = medpopTraj, comb=med.comb, detect=0.5)
plot(runmed[,19])
summary(runmed[,19])

#high
Rmodelhigh <- nimbleModel(code = IPMmod, constants = const1h, data = dat1h, 
                         check = FALSE, calculate = FALSE, inits = initshigh)# inits1)
confhigh <- configureMCMC(Rmodelhigh, monitors = params1)
Rmcmchigh <- buildMCMC(confhigh)  
Cmodelhigh <- compileNimble(Rmodelhigh, showCompilerOutput = FALSE)
Cmcmchigh <- compileNimble(Rmcmchigh, project = Rmodelhigh)
Cmcmchigh$run(niter=ni, nburnin=nb)
out_high<-as.data.frame(as.matrix(Cmcmchigh$mvSamples))
mean(out_high$p.surv)
summary(out_high$p.surv)


runhigh<-runIPMmod(nb=nb, ni=ni, nt=nt, nc=5, popDat = highpopDat, 
                   popTraj = highpopTraj, comb=high.comb, detect=0.8)
summary(runhigh[,19])
plot(runhigh[,19])
# 
# Cmcmc1$run(niter=ni, nburnin=nb)
# out22<-as.data.frame(as.matrix(Cmcmc1$mvSamples))
# out2 <- runMCMC(Cmcmc1, niter = ni , nburnin = nb , nchains = nc, inits = inits1,
#                 setSeed = FALSE, progressBar = TRUE, samplesAsCodaMCMC = TRUE)
# #sink()
# t.end <- Sys.time()
# (runTime <- t.end - t.start)
# beep(sound = 3)
# 
# ########HAS: if it doesnt run and there are slice sampler or other issues, use:
# Cmcmc1$mvSamples[["Ntot"]][[1]] # or any other parameters to help diagnose
# #I believe it is the starting values for the chain, so you can (hopefully) find the problem
# 
# ####
# 
# # run each model
# 
# library(coda)
# 
# # check that results are unbiased
# summary(out2)
# low.comb
# 
# # check model convergence
# 
# # AEB TODO
# # what the heck why is it mixing so terribly?
# gelman.diag(out2)
# plot(out2)
# 
# 
# # LOOKS GOOD
# # still need to check
# # H M L
# # with all combinations of
