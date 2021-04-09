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


 #low.comb <- low.lam.combos[sample(1:5000, 1), 1:3]
 #med.comb <- med.lam.combos[sample(1:5000, 1), 1:3]
 #high.comb <- high.lam.combos[sample(1:5000, 1), 1:3]
# 
# # thought
# # make adult survival be higher?? force it
# 
# lowpopTraj <- simPopTrajectory(n.years=15,
#                                n.data.types=c(0.25,0.25,0.25),
#                                age.init=c(150,150),
#                                phi.1=as.numeric(low.comb[2]),
#                                phi.ad=as.numeric(low.comb[3]),
#                                f=as.numeric(low.comb[1]))
# # returns indfates and n
# 
# medpopTraj <- simPopTrajectory(n.years=15,
#                                n.data.types=c(0.25,0.25,0.25),
#                                age.init=c(150,150),
#                                phi.1=as.numeric(med.comb[2]),
#                                phi.ad=as.numeric(med.comb[3]),
#                                f=as.numeric(med.comb[1]))
# 
# highpopTraj <- simPopTrajectory(n.years=15,
#                                 n.data.types=c(0.25,0.25,0.25),
#                                 age.init=c(150,150),
#                                 phi.1=as.numeric(high.comb[2]),
#                                 phi.ad=as.numeric(high.comb[3]),
#                                 f=as.numeric(high.comb[1]))


#just load data from files
#is the true param values in scenarios.RDS? Where to pull those from
low.lam.params <- readRDS("~/Desktop/IPMSim/low.lam.params.RDS")
low.lam.params[1,] #for the traj below?
lowpopTraj<-readRDS("~/Desktop/IPMSim/data/lowTrajectories/lowpopTraj-1-1.RDS")
med.lam.params <- readRDS("~/Desktop/IPMSim/med.lam.params.RDS")
med.lam.params[1,]
medpopTraj<-readRDS("~/Desktop/IPMSim/data/medTrajectories/medpopTraj-1-1.RDS")
high.lam.params <- readRDS("~/Desktop/IPMSim/high.lam.params.RDS")
high.lam.params[1,]
highpopTraj<-readRDS("~/Desktop/IPMSim/data/highTrajectories/highpopTraj-1-1.RDS")

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
    p.surv = d,
    mean.phi = c(combo$phi1, combo$phiad),
    mean.p = d, 
    fec = combo$fec,
    n1.start=poptraj_Noutsmat[1,1],
    nad.start=poptraj_Noutsmat[2,1])
  return(initsout)
}
initslow<-initsfn(lowpopTraj$Nouts, d=detect.l, combo=low.lam.params[1,])
initsmed<-initsfn(medpopTraj$Nouts, d=detect.m, combo=med.lam.params[1,])
initshigh<-initsfn(highpopTraj$Nouts, d=detect.h, combo=high.lam.params[1,])


#### PARAMETERS TO MONITOR ####
params1 <- c("p.surv", "mean.phi","mean.p", "fec", "lambda") 

#### MCMC SETTINGS ####
nb <- 10000 #burn-in
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
summary(out_low)
low.lam.params[1,]

runlow<-runIPMmod(nb=nb, ni=ni, nt=nt, nc=5, popDat = lowpopDat,
                   popTraj = lowpopTraj, comb=low.lam.params[1,], detect=0.3)
summary(runlow[,19])
plot(runlow[,19])
hist(runlow[,19])
low.lam.params[1,]

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
summary(out_med)
med.lam.params[1,]

# runmed<-runIPMmod(nb=nb, ni=ni, nt=nt, nc=5, popDat = medpopDat, 
#                   popTraj = medpopTraj, comb=med.comb, detect=0.5)
# plot(runmed[,19])
# summary(runmed[,19])

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
summary(out_high)
high.lam.params[1,]

# runhigh<-runIPMmod(nb=nb, ni=ni, nt=nt, nc=5, popDat = highpopDat, 
#                    popTraj = highpopTraj, comb=high.comb, detect=0.8)
# summary(runhigh[,19])
# plot(runhigh[,19])
